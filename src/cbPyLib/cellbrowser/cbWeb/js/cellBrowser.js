// A viewer for (x,y) scatter plots of small circles
// shows associated meta data (usually key->val attributes, one mapping per circle)
// and expression data (string -> float, one mapping per circle)

/* jshint -W097 */
/* jshint -W117 */  // don't complain about unknown classes, like Set()
/* jshint -W104 */  // allow 'const'
/* jshint -W069 */  // object access with ["xx"]

/* TODO:
 * - status bug - reset last expression array when coloring by meta */

"use strict";

var cellbrowser = function() {
    var db = null; // the cbData object from cbData.js. Loads coords,
                   // annotations and gene expression vectors

    var gVersion = "$VERSION$"; // cellbrowser.py:copyStatic will replace this with the pip version or git release
    var gCurrentCoordName = null; // currently shown coordinates

    // object with all information needed to map to the legend colors
    var gLegend = null;

    // optional second legend, for split screen mode
    var gOtherLegend = null;

    // all info about the current legend. gLegend.rows is:
    // [ colorHexStr, defaultColorHexStr, label, count, internalKey, uniqueKey ]
    // internalKey can be int or str, depending on current coloring mode.
    // E.g. in meta coloring, it's the metadata string value.
    // When doing expression coloring, it's the expression bin index.
    // uniqueKey is used to save manually defined colors to localStorage

    var renderer = null;

    var background = null;

    // last 10 genes
    var gRecentGenes = [];

    // -- CONSTANTS
    var gTitle = "UCSC Cell Browser";
    var COL_PREFIX = "col_";

    var gOpenDataset = null; // while navigating the open dataset dialog, this contains the current name
        // it's a global variable as the dialog is not a class (yet?) and it's the only piece of data
        // it is a subset of dataset.json , e.g. name, description, cell count, etc.

    // depending on the type of data, single cell or bulk RNA-seq, we call a circle a
    // "sample" or a "cell". This will adapt help menus, menus, etc.
    var gSampleDesc = "cell";

    var gLabelCoordCache = {};

    // width of left meta bar in pixels
    var metaBarWidth = 250;
    // margin between left meta bar and drawing canvas
    var metaBarMargin = 0;
    // width of legend, pixels
    var legendBarWidth = 200;
    var legendBarMargin = 0;
    // width of the metaBar tooltip (histogram)
    var metaTipWidth = 400;
    // height of pull-down menu bar at the top, in pixels
    var menuBarHeight = null; // defined at runtime by div.height
    // height of the toolbar, in pixels
    var toolBarHeight = 28;
    // position of first combobox in toolbar from left, in pixels
    var toolBarComboLeft = metaBarWidth;
    var toolBarComboTop   = 2;
    // width of the collection combobox
    var collectionComboWidth = 200;
    var layoutComboWidth = 200;
    // width of a single gene cell in the meta gene bar tables
    //var gGeneCellWidth = 66;

    // height of bottom gene bar
    var geneBarHeight = 100;
    var geneBarMargin = 5;
    // color for missing value when coloring by expression value
    //var cNullColor = "CCCCCC";
    //const cNullColor = "DDDDDD";
    //const cNullColor = "95DFFF"; //= light blue
    const cNullColor = "e1f6ff"; //= light blue

    const cDefGradPalette = "tol-sq-blue";  // default legend gradient palette for gene expression
    // this is a special palette, tol-sq with the first entry being a light blue, so 0 stands out a bit more
    const cDefGradPaletteHeat = "tol-sq";  // default legend gradient palette for the heatmap
    const cDefQualPalette  = "rainbow"; // default legend palette for categorical values

    var datasetGradPalette = cDefGradPalette;
    var datasetQualPalette = cDefQualPalette;


    const exprBinCount = 10; //number of expression bins for genes
    // has to match cbData.js.exprBinCount - TODO - share the constant between these two files

    var HIDELABELSNAME = "Hide labels";
    var SHOWLABELSNAME = "Show labels";
    var METABOXTITLE   = "By Annotation";

    // maximum number of distinct values that one can color on
    const MAXCOLORCOUNT = 500;

    // histograms show only the top X values and summarize the rest into "other"
    var HISTOCOUNT = 12;
    // the sparkline is a bit shorter
    var SPARKHISTOCOUNT = 12;

    // links to various external databases
    var dbLinks = {
        "HPO" : "https://hpo.jax.org/app/browse/gene/", // entrez ID
        "OMIM" : "https://omim.org/entry/", // OMIM ID
        "COSMIC" : "http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=", // gene symbol
        "SFARI" : "https://gene.sfari.org/database/human-gene/", // gene symbol
        "BrainSpLMD" : "http://www.brainspan.org/lcm/search?exact_match=true&search_type=gene&search_term=", // entrez
        "BrainSpMouseDev" : "http://developingmouse.brain-map.org/gene/show/", // internal Brainspan ID
        "Eurexp" : "http://www.eurexpress.org/ee/databases/assay.jsp?assayID=", // internal ID
        "LMD" : "http://www.brainspan.org/lcm/search?exact_match=true&search_type=gene&search_term=" // entrez
    };

    var DEBUG = true;

    function _dump(o) {
    /* for debugging */
        console.log(JSON.stringify(o));
    }

    function formatString (str) {
        /* Stackoverflow code https://stackoverflow.com/a/18234317/233871 */
        /* "a{0}bcd{1}ef".formatUnicorn("foo", "bar"); // yields "aFOObcdBARef" */
            if (arguments.length) {
                var t = typeof arguments[0];
                var key;
                var args = ("string" === t || "number" === t) ?
                    Array.prototype.slice.call(arguments)
                    : arguments[0];

                for (key in args) {
                    str = str.replace(new RegExp("\\{" + key + "\\}", "gi"), args[key]);
                }
            }
            return str;
        }

    // Median of medians: https://en.wikipedia.org/wiki/Median_of_medians
    // find median in an unsorted array, worst-case complexity O(n).
    // from https://gist.github.com/wlchn/ee15de1da59b8d6981a400eee4376ea4
    const selectMedian = (arr, compare) => {
        return _selectK(arr, Math.floor(arr.length / 2), compare);
    };

    const _selectK = (arr, k, compare) => {
        if (!Array.isArray(arr) || arr.length === 0 || arr.length - 1 < k) {
            return;
        }
        if (arr.length === 1) {
            return arr[0];
        }
        let idx = _selectIdx(arr, 0, arr.length - 1, k, compare || _defaultCompare);
        return arr[idx];
    };

    const _partition = (arr, left, right, pivot, compare) => {
        let temp = arr[pivot];
        arr[pivot] = arr[right];
        arr[right] = temp;
        let track = left;
        for (let i = left; i < right; i++) {
            // if (arr[i] < arr[right]) {
            if (compare(arr[i], arr[right]) === -1) {
                let t = arr[i];
                arr[i] = arr[track];
                arr[track] = t;
                track++;
            }
        }
        temp = arr[track];
        arr[track] = arr[right];
        arr[right] = temp;
        return track;
    };

    const _selectIdx = (arr, left, right, k, compare) => {
        if (left === right) {
            return left;
        }
        let dest = left + k;
        while (true) {
            let pivotIndex =
                right - left + 1 <= 5
                ? Math.floor(Math.random() * (right - left + 1)) + left
                : _medianOfMedians(arr, left, right, compare);
            pivotIndex = _partition(arr, left, right, pivotIndex, compare);
            if (pivotIndex === dest) {
                return pivotIndex;
            } else if (pivotIndex < dest) {
                left = pivotIndex + 1;
            } else {
                right = pivotIndex - 1;
            }
        }
    };

    const _medianOfMedians = (arr, left, right, compare) => {
        let numMedians = Math.ceil((right - left) / 5);
        for (let i = 0; i < numMedians; i++) {
        let subLeft = left + i * 5;
        let subRight = subLeft + 4;
        if (subRight > right) {
            subRight = right;
        }
        let medianIdx = _selectIdx(arr, subLeft, subRight, Math.floor((subRight - subLeft) / 2), compare);
        let temp = arr[medianIdx];
        arr[medianIdx] = arr[left + i];
        arr[left + i] = temp;
        }
        return _selectIdx(arr, left, left + numMedians - 1, Math.floor(numMedians / 2), compare);
    };

    const _defaultCompare = (a, b) => {
        return a < b ? -1 : a > b ? 1 : 0;
    };
    // End median of medians

    function debug(msg, args) {
        if (DEBUG) {
            console.log(formatString(msg, args));
        }
    }

    function warn(msg) {
        alert(msg);
    }

    function getById(query) {
        return document.getElementById(query);
    }

    function cloneObj(d) {
    /* returns a copy of an object, wasteful */
        // see http://stackoverflow.com/questions/122102/what-is-the-most-efficient-way-to-deep-clone-an-object-in-javascript
        return JSON.parse(JSON.stringify(d));
    }

    function cloneArray(a) {
    /* returns a copy of an array */
        return a.slice();
    }

    function copyNonNull(srcArr, trgArr) {
    /* copy non-null values to trgArr */
        if (srcArr.length!==trgArr.length)
            alert("warning - copyNonNull - target and source array have different sizes.");

        for (var i = 0; i < srcArr.length; i++) {
            if (srcArr[i]!==null)
                trgArr[i] = srcArr[i];
        }
        return trgArr;
    }

    function isEmpty(obj) {
        for(var key in obj) {
            if(obj.hasOwnProperty(key))
                return false;
        }
        return true;
    }

    function allEmpty(arr) {
        /* return true if all members of array are white space only strings */
        var newArr = arr.filter(function(str) { return /\S/.test(str); });
        return (newArr.length===0);
    }

    function copyNonEmpty(srcArr, trgArr) {
    /* copy from src to target array if value is not "". Just return trgArr is srcArr is null or lengths don't match.  */
        if (!srcArr || (srcArr.length!==trgArr.length))
            return trgArr;

        for (var i = 0; i < srcArr.length; i++) {
            if (srcArr[i]!=="")
                trgArr[i] = srcArr[i];
        }
        return trgArr;
    }

    function keys(o) {
    /* return all keys of object as an array */
        var allKeys = [];
        for(var k in o) allKeys.push(k);
        return allKeys;
    }

    function trackEvent(eventName, eventLabel) {
    /* send an event to google analytics */
        if (typeof gtag !== 'function')
            return;
        gtag('event', eventName, eventLabel);
    }

    function trackEventObj(eventName, obj) {
    /* send an event obj to google analytics */
        if (typeof gtag !== 'function')
            return;
        gtag('event', obj);
    }

    function classAddListener(className, type, listener) {
        /* add an event listener for all elements of a class */
        var els = document.getElementsByClassName(className);
        for (let el of els) {
            el.addEventListener(type, listener);
        }
    }

    function capitalize(s) {
        return s[0].toUpperCase() + s.slice(1);
    }

    function cleanString(s) {
        /* make sure that string only contains normal characters. Good when printing something that may contain
         * dangerous */
        if (s===undefined)
            return undefined;
        return s.replace(/[^0-9a-zA-Z _-]/g, '');
    }

    function findMetaValIndex(metaInfo, value) {
        /* return the index of the value of an enum meta field */
        var valCounts = metaInfo.valCounts;
        for (var valIdx = 0; valIdx < valCounts.length; valIdx++) {
            if (valCounts[valIdx][0]===value)
                return valIdx;
        }
    }

    function intersectArrays(arrList) {
        /* return the intersection of all arrays as an array. Non-IE11? */
        var smallSet = new Set(arrList[0]);
        for (var i=1; i < arrList.length; i++) {
            var otherSet = new Set(arrList[i]);
            smallSet = new Set([...smallSet].filter(x => otherSet.has(x)));
        }
        var newArr = Array.from(smallSet);
        // alternative without spread:
        //function intersection(setA, setB) {
          //  var _intersection = new Set();
          //  for (var elem of setB) {
          //      if (setA.has(elem)) {
          //          _intersection.add(elem);
          //      }
          //  }
          //  return _intersection;
        //}
        return newArr;
    }

    function saveToUrl(key, value, defaultValue) {
    /* save a value in both localStorage and the URL. If the value is defaultValue or null, remove it */
        if (value===defaultValue || value===null) {
            localStorage.removeItem(key);
            delState(key);
        }
        else {
            localStorage.setItem(key, value);
            addStateVar(key, value);
        }
    }

    function getFromUrl(key, defaultValue) {
    /* get a value from localStorage or the current URL or return the default if not defined in either place.
     * The URL overrides localStorage. */
        var val = getVar(key);
        if (val!==undefined)
            return val;

        val = localStorage.getItem(key);
        if (val===null)
            return defaultValue
        else
            return val;
    }

    function getBaseUrl() {
       /* return URL of current page, without args or query part */
       var myUrl = window.location.href;
       myUrl = myUrl.replace("#", "");
       var urlParts = myUrl.split("?");
       var baseUrl = urlParts[0];
       return baseUrl;
    }

    function copyToClipboard(element) {
    /* https://stackoverflow.com/questions/22581345/click-button-copy-to-clipboard-using-jquery */
        var $temp = $("<input>");
        $("body").append($temp);
        $temp.val($(element).text()).select();
        document.execCommand("copy");
        $temp.remove();
    }

    function iWantHue(n) {
        /* a palette as downloaded from iwanthue.com - not sure if this is better. Ellen likes it */
        var colList = ["7e4401", "244acd", "afc300", "a144cb", "00a13e",
            "f064e5", "478700", "727eff", "9ed671", "b6006c", "5fdd90", "f8384b",
            "00b199", "bb000f", "0052a3", "fcba56", "005989", "c57000", "7a3a78",
            "ccca76", "ff6591", "265e1c", "ff726c", "7b8550", "923223", "9a7e00",
            "ffa9ad", "5f5300", "ff9d76", "b3885f"];
        var colList2 = ["cd6a00", "843dc3", "c9cd31", "eda3ff", "854350"];
        if (n<=5)
            colList = colList2;
        return colList.slice(0, n);
    }

    function activateTooltip(selector) {
        // noconflict in html, I had to rename BS's tooltip to avoid overwrite by jquery
        var ttOpt = {
            "html": true,
            "animation": false,
            "delay": {"show":350, "hide":100},
            "trigger" : "hover",
            container:"body"
        };
        $(selector).bsTooltip(ttOpt);
    }

    function menuBarHide(idStr) {
    /* hide a menu bar selector */
         $(idStr).parent().addClass("disabled").css("pointer-events", "none");
    }

    function menuBarShow(idStr) {
    /* show a menu bar entry given its selector */
         $(idStr).parent().removeClass("disabled").css("pointer-events", '');
    }

    function updateMenu() {
    /* deactivate menu options based on current variables */
     // the "hide selected" etc menu options are only shown if some cells are selected
     if (renderer.selCells.length===0) {
         menuBarHide("#tpOnlySelectedButton");
         menuBarHide("#tpFilterButton");
     }
     else {
         menuBarShow("#tpOnlySelectedButton");
         menuBarShow("#tpFilterButton");
     }

     // the "show all" menu entry is only shown if some dots are actually hidden
     //if ((pixelCoords!==null && shownCoords!==null) && pixelCoords.length===shownCoords.length)
         //menuBarHide("#tpShowAllButton");
     //else
         //menuBarShow("#tpShowAllButton");

     //$("#tpTrans"+(transparency*100)).addClass("active");
     //$("#tpSize"+circleSize).addClass("active");

     // the "hide labels" menu entry is only shown if there are labels
     //if (gCurrentDataset.labelField === null)
         //menuBarHide("#tpHideLabels");
     //if (gCurrentDataset.showLabels===true)
         //$("#tpHideLabels").text(HIDELABELSNAME);
     //else
         //$("#tpHideLabels").text(SHOWLABELSNAME);
    }

    function prettySeqDist(count, addSign) {
        /* create human-readable string from chrom distance */
        var f = count;
        var sign = "";
        if (addSign && count > 0)
            sign = "+";

        if (Math.abs(count)>=1000000) {
            f = (count / 1000000);
            return sign+f.toFixed(3)+"Mbp";
        }
        if (Math.abs(count)>=10000) {
            f = (count / 1000);
            return sign+f.toFixed(2)+"kbp";
        }
        if (Math.abs(count)>=1000) {
            f = (count / 1000);
            return sign+f.toFixed(2)+"kbp";
        }
        return sign+f+"bp";
    }

    function prettyNumber(count, isBp) /*str*/ {
        /* convert a number to a shorter string, e.g. 1200 -> 1.2k, 1200000 -> 1.2M, etc */
        var f = count;
        if (count>1000000) {
            f = (count / 1000000);
            return f.toFixed(1)+"M";
        }
        if (count>10000) {
            f = (count / 1000);
            return f.toFixed(0)+"k";
        }
        if (count>1000) {
            f = (count / 1000);
            return f.toFixed(1)+"k";
        }
        return f;
    }

    function addMd5(url, md5s, key) {
        /* lookup key in md5s and add value to url separate by ? */
        if (md5s && md5s[key])
            url += "?"+md5s[key];
        return url;
    }

    function preloadImage(url) {
        let img= new Image();
        img.src = url;
    }

    function openDatasetLoadPane(datasetInfo) {
        /* open dataset dialog: load html into the three panes  */
        //var datasetName = datasetInfo.name;
        //var md5 = datasetInfo.md5;
        // the UCSC apache serves latin1, so we force it back to utf8
        gOpenDataset = datasetInfo; // for click handlers in the right panel
        var thumbUrl = cbUtil.joinPaths([datasetInfo.name, "thumb.png"]);
        preloadImage(thumbUrl); // many datasets have thumb.png, so preload it now

        $.ajaxSetup({
            'beforeSend' : function(xhr) {
                if (xhr && xhr.overrideMimeType)
                    xhr.overrideMimeType('text/html; charset=utf8');
            },
        });

        let datasetName = datasetInfo.name;
        let md5 = datasetInfo.md5;
        if (datasetInfo.hasFiles && datasetInfo.hasFiles.indexOf("datasetDesc")!==-1) {
            // description is not through html files but a json file
            var jsonUrl = cbUtil.joinPaths([datasetName, "desc.json"]) +"?"+md5;
            fetch(jsonUrl)
              .then(function(response) {
                if(!response.ok) {
                    throw new Error('Could not find desc.json file');
                }
                return response.json();
              })
              .catch(function(err) {
                  var msg = "File "+jsonUrl+" was not found but datasetDesc.json has 'datasetDesc' in hasFiles. Internal error. Please contact the site admin or cells@ucsc.edu";
                  $( "#pane1" ).html(msg);
                  $( "#pane2" ).html(msg);
                  $( "#pane3" ).html(msg);
                  $( "#pane3" ).show();
              })
              .then(function(desc) {
                  datasetDescToHtml(datasetInfo, desc);
              });
        }
        else {
          var message = "This dataset does not seem to have a desc.conf file. Please "+
              "read https://cellbrowser.readthedocs.io/en/master/dataDesc.html or run 'cbBuild --init' to create one";
          if (datasetInfo.abstract)
              // the top-level non-hierarchy dataset.conf has a message in it. Use it here, as a fallback.
              message = datasetInfo.abstract;

          $( "#pane1" ).html(message);
          $( "#pane2" ).hide();
          $( "#tabLink2" ).hide();
          $( "#pane3" ).hide();
          $( "#tabLink3" ).hide();
        }
        $("#tpOpenDialogTabs").tabs("refresh").tabs("option", "active", 0);
    }

    let descLabels = {
        "paper_url":"Publication",
        "other_url" : "Website",
        "geo_series" : "NCBI GEO Series", // = CIRM tagsV5
        "sra" : "NCBI Short Read Archive",
        "pmid" : "PubMed Abstract",
        "pmcid" : "PubMed Fulltext",
        "sra_study" : "NCBI Short-Read Archive",
        "ega_study" : "European Genotype-Phenot. Archive",
        "bioproject" : "NCBI Bioproject",
        "dbgap" : "NCBI DbGaP",
        "biorxiv_url" : "BioRxiv preprint",
        "doi" : "Publication Fulltext",
        "arrayexpress" : "ArrayExpress",
        "ena_project" : "European Nucleotide Archive",
        "cirm_dataset" : "California Institute of Regenerative Medicine Dataset",
    };

    let descUrls = {
        "geo_series" : "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
        "sra_study" : "https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=",
        "bioproject" : "https://www.ncbi.nlm.nih.gov/bioproject/",
        "ega_study" : "https://ega-archive.org/studies/",
        "pmid" : "https://www.ncbi.nlm.nih.gov/pubmed/",
        "pmcid" : "https://www.ncbi.nlm.nih.gov/pmc/articles/",
        "dbgap" : "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=",
        "doi" : "http://dx.doi.org/",
        "ena_project" : "https://www.ebi.ac.uk/ena/data/view/",
        "cirm_dataset" : "https://cirm.ucsc.edu/d/",
        "arrayexpress" : "https://www.ebi.ac.uk/arrayexpress/experiments/",
    }

    function htmlAddLink(htmls, desc, key, linkLabel) {
        /* add a link to html on a new line. if desc[key] includes a space, the part after it is the link label. */
        if (!desc[key])
            return;

        let label = "Link";
        if (linkLabel)
            label = linkLabel;
        else
            label = descLabels[key];

        htmls.push("<b>");
        htmls.push(label);
        htmls.push(": </b>");

        // for cases where more than one ID is needed, this function also accepts a list in the object
        // for 99% of the cases, it'll be a string though
        let urls = desc[key];
        if (!(urls instanceof Array))
            urls = [urls];

        let frags = []; // html fragments, one per identifier
        for (let url of urls) {
            url = url.toString(); // in case it's an integer or float
            let urlLabel = url;
            let spcPos = url.indexOf(" ");
            if (spcPos!==-1) {
                urlLabel = url.slice(spcPos+1);
                url = url.slice(0,spcPos);
            }

            if (!url.startsWith("http"))
                url = descUrls[key]+url;

            let parts = []
            parts.push("<a target=_blank href='");
            parts.push(url);
            parts.push("'>");
            parts.push(urlLabel);
            parts.push("</a>");
            let htmlLine = parts.join("");
            frags.push(htmlLine);
        }
        htmls.push(frags.join(", "));
        htmls.push("<br>");
    }

    function buildDownloadsPane(datasetInfo, desc) {
        var htmls = [];
        if (datasetInfo.name==="") { // the top-level desc page has no methods/downloads, it's just a dataset list
            $( "#pane3" ).hide();
            $( "#tabLink3" ).hide();
        } else {
            if (desc.coordFiles===undefined) {
                htmls.push("To download the data for datasets in this collection: open the collection, ");
                htmls.push("select a dataset in the list to the left, and navigate to the 'Data Download' tab. ");
                htmls.push("This information can also be accessed while viewing a dataset by clicking the 'Info &amp; Downloads' button.");
            } else if (desc.hideDownload===true || desc.hideDownload=="True" || desc.hideDownload=="true") {
                htmls.push("The downloads section has been deactivated by the authors.");
                htmls.push("Please contact the dataset authors to get access.");
            } else {
                if (desc.matrixFile!==undefined && desc.matrixFile.endsWith(".mtx.gz")) {
                    htmls.push("<p><b>Expression in MTX format:</b> <a href='"+datasetInfo.name);
                    htmls.push("/matrix.mtx.gz'>matrix.mtx.gz</a>");
                    htmls.push(", <a href='"+datasetInfo.name);
                    htmls.push("/features.tsv.gz'>features.tsv.gz</a>");
                    htmls.push(", <a href='"+datasetInfo.name);
                    htmls.push("/barcodes.tsv.gz'>barcodes.tsv.gz</a>");
                } else {
                htmls.push("<p><b>Expression matrix:</b> <a href='"+datasetInfo.name);
                htmls.push("/exprMatrix.tsv.gz'>exprMatrix.tsv.gz</a>");
                }
                if (desc.unitDesc)
                    htmls.push("<br>Values are: "+desc.unitDesc);
                htmls.push("</p>");

                if (desc.rawMatrixFile) {
                    htmls.push("<p><b>Raw expression matrix:</b> <a href='"+datasetInfo.name);
                    htmls.push("/"+desc.rawMatrixFile+"'>"+desc.rawMatrixFile+"</a>");
                    if (desc.rawMatrixNote)
                        htmls.push("<br>"+desc.rawMatrixNote);
                    htmls.push("</p>");
                }

                htmls.push("<p><i><a style='float:right; padding-left: 100px'; target=_blank href='https://cellbrowser.readthedocs.io/en/master/load.html'>Help: Load matrix/meta into Seurat or Scanpy</a></i></p>");

                htmls.push("<p><b>Cell meta annotations:</b> <a target=_blank href='"+datasetInfo.name);
                htmls.push("/meta.tsv'>meta.tsv</a>");
                if (desc.metaNote)
                    htmls.push("<br>"+desc.metaNote);
                htmls.push("</p>");

                htmls.push("<p><b>Dimensionality reduction coordinates:</b><br>");
                for (let fname of desc.coordFiles)
                    htmls.push("<a target=_blank href='"+datasetInfo.name+"/"+fname+"'>"+fname+"</a><br>");
                htmls.push("</p>");

                if (desc.supplFiles) {
                    let supplFiles = desc.supplFiles;
                    for (let suppFile of supplFiles) {
                        let label = suppFile.label;
                        let fname = suppFile.file;
                        htmls.push("<p><b>"+label+":</b> <a href='"+datasetInfo.name);
                        htmls.push("/"+fname+"'>"+fname+"</a>");
                        htmls.push("</p>");
                    }
                }

                htmls.push("<p><b>Dataset description</b>: ");
                htmls.push("<a target=_blank href='"+datasetInfo.name+"/desc.json'>desc.json</a></p>");

                htmls.push("<p><b>Cell Browser configuration</b>: ");
                htmls.push("<a target=_blank href='"+datasetInfo.name+"/dataset.json'>dataset.json</a></p>");

                $( "#pane3" ).html(htmls.join(""));
                $( "#pane3" ).show();
                $( "#tabLink3" ).show();
            }

        }
    }

    function buildImagesPane(datasetInfo, desc) {
        if (!desc.imageSets) {
            $( "#tabLinkImg" ).hide();
            $( "#paneImg" ).hide();
            return;
        }

        let htmls = [];
        htmls.push("<h4>Microscopy images</h4>");
        // TOC
        let catIdx = 0;
        htmls.push("<div style='padding-bottom:8px'>Jump to: ");
        for (let catInfo of desc.imageSets) {
            htmls.push("<a style='padding-left:12px' href='#imgCat"+catIdx+"'>"+catInfo.categoryLabel+"</a>");
            catIdx++;
        }
        htmls.push("</div>");

        if (desc.imageSetNote)
            htmls.push("<p>"+desc.imageSetNote+"<p>");

        // actual HTML
        catIdx = 0;
        for (let catInfo of desc.imageSets) {
            htmls.push("<div style='padding-top:6px; padding-bottom:4px' class='tpImgCategory'>");
            htmls.push("<a name='imgCat"+catIdx+"'></a>");
            catIdx++;
            htmls.push("<b>"+catInfo.categoryLabel+":</b><br>");
            let imgDir = datasetInfo.name+"/images/";
            htmls.push("<div style='padding-left:1em; padding-top:4px' class='tpImgSets'>");

            for (let imgSet of catInfo.categoryImageSets) {
                let imgLinks = [];
                if (imgSet.setLabel)
                    htmls.push("<b>"+imgSet.setLabel+"</b><br>");
                htmls.push("<div style='padding-left:1em;' class='tpImgSetLinks'>");
                for (let img of imgSet.images) {
                    imgLinks.push("Show: <a target=_blank href='"+imgDir+img.file+"'>"+img.label+"</a>");
                }
                htmls.push(imgLinks.join(", "));

                if (imgSet.downloads) {
                    let dlLinks = [];
                    for (let dl of imgSet.downloads) {
                        //dlLinks.push("<a href='"+imgDir+dl.file+"' download><span style='font-size:12px' class='material-icons-round'>get_app</span>"+dl.label+"</a>");
                        dlLinks.push("<a href='"+imgDir+dl.file+"' download>"+dl.label+"</a>");
                    }
                    //htmls.push("<br><div style='padding-left: 1em'>Download: ");
                    htmls.push("<br><div>Download: ");
                    htmls.push(dlLinks.join(", "));
                    htmls.push("</div>");
                }
                //htmls.push("<br>");
                htmls.push("</div>"); //  tpImgSetLinks
            }
            htmls.push("</div>"); //  tpImgSets
            htmls.push("</div>"); //  tpImgCategory
        }
        //htmls.push("</ul>");
        $( "#paneImg" ).html(htmls.join(""));
        $( "#paneImg" ).show();
        $( "#tabLinkImg" ).show();
    }

    function buildMethodsPane(datasetInfo, desc) {
        // methods panel
        //
        var htmls = [];
        if (desc.methods) {
            htmls.push("<p>");
            htmls.push(desc.methods);
            htmls.push("</p>");
        }
        if (desc.algParams) {
            htmls.push("<p><b>Algorithm parameters: </b>");
            let algParams = desc.algParams;
            if (algParams instanceof Object)
                algParams = Object.entries(algParams);

            for (let i=0; i<algParams.length; i++) {
                let key = algParams[i][0];
                let val = algParams[i][1];
                htmls.push(key+"="+val+", ");
            }
            htmls.push("</p>");
        }
        if (htmls.length!==0) {
            $( "#pane2" ).html(htmls.join(""));
            $( "#pane2" ).show();
            $( "#tabLink2" ).show();
        } else {
            $( "#pane2" ).hide();
            $( "#tabLink2" ).hide();
        }
    }

    function pageAtUcsc() {
        // return true if current page is at ucsc.edu
        return (window.location.hostname.endsWith("ucsc.edu"));
    }

    function datasetDescToHtml(datasetInfo, desc) {
        /* given an object with keys title, abstract, pmid, etc, fill the dataset description tabs with html */
        if (!desc) // http errors call this with undefined
            return;

        let htmls = [];

        if (datasetInfo.name==="") // the root dataset
            $('#tabLink1').text("Overview");
        else
            $('#tabLink1').text("Abstract");

        if (desc.title) {
            htmls.push("<h4>");
            htmls.push(desc.title);
            htmls.push("</h4>");
        }
        if (desc.image) {
            htmls.push("<img style='float:right; padding-left:5px' src='");
            htmls.push(datasetInfo.name+"/"+desc.image[0]+"'");
            if (desc.imageMap)
                htmls.push(" usemap='#clickmap'");
            htmls.push(" width='"+desc.image[1]+"' height='"+desc.image[2]+"'>");
        }
        if (desc.imageMap) {
            htmls.push('<map name="clickmap">');
            htmls.push(desc.imageMap);
            htmls.push('</map>');
        }

        if (desc.abstract) {
            htmls.push("<p>");
            htmls.push(desc.abstract);
            htmls.push("</p>");
        }
        else {
            // the top-level hardcoded dataset for non-hierarchy mode has the abstract in the
            // dataset config. It's a lot easier this way, so just pull it in here.
            htmls.push("<p>");
            htmls.push(datasetInfo.abstract);
            htmls.push("</p>");
        }

        if (desc.author) {
            htmls.push("<b>Author: </b> "+desc.author);
            htmls.push("<br>");
        }

        if (desc.authors) {
            htmls.push("<b>Authors: </b> "+desc.authors);
            htmls.push("<br>");
        }

        if (desc.lab) {
            htmls.push("<b>Lab: </b> "+desc.lab);
            htmls.push("<br>");
        }
        if (desc.institution) {
            htmls.push("<b>Institution: </b> "+desc.institution);
            htmls.push("<br>");
        }


        htmlAddLink(htmls, desc, "biorxiv_url");
        htmlAddLink(htmls, desc, "paper_url");
        htmlAddLink(htmls, desc, "other_url");
        htmlAddLink(htmls, desc, "geo_series");
        htmlAddLink(htmls, desc, "pmid");
        htmlAddLink(htmls, desc, "dbgap");
        htmlAddLink(htmls, desc, "sra_study");
        htmlAddLink(htmls, desc, "bioproject");
        htmlAddLink(htmls, desc, "sra");
        htmlAddLink(htmls, desc, "doi");
        htmlAddLink(htmls, desc, "arrayexpress");
        htmlAddLink(htmls, desc, "cirm_dataset");
        htmlAddLink(htmls, desc, "ega_study");
        htmlAddLink(htmls, desc, "ena_project");

        if (desc.urls) {
            for (let key in desc.urls)
                htmlAddLink(htmls, desc.urls, key, key);
        }

        if (desc.custom) {
            for (let key in desc.custom) {
                htmls.push("<b>"+key+": </b> "+desc.custom[key]);
                htmls.push("<br>");
            }
        }

        if (desc.submitter) {
            htmls.push("<b>Submitted by: </b> "+desc.submitter);
            if (desc.submission_date) {
                htmls.push(" ("+desc.submission_date);
                htmls.push(")");
            }
            if (desc.version)
                htmls.push(", Version "+desc.version);
            htmls.push("<br>");
        }

        if (desc.shepherd) {
            htmls.push("<b>UCSC Data Shepherd: </b> "+desc.shepherd);
            htmls.push("<br>");
        }
        if (desc.wrangler) {
            htmls.push("<b>UCSC Data Wrangler: </b> "+desc.wrangler);
            htmls.push("<br>");
        }

        let topName = datasetInfo.name.split("/")[0];
        if (pageAtUcsc()) {
            if (datasetInfo.name!=="") {
                if ((datasetInfo.parents) && (datasetInfo.parents.length > 1)) {
                    // if the dataset is a collection
                    htmls.push("<b>Direct link to this collection for manuscripts: </b> https://"+topName+".cells.ucsc.edu");
                    htmls.push("<br>");
                }
                else {
                    htmls.push("<b>Direct link to this plot for manuscripts: </b> https://"+topName+".cells.ucsc.edu");
                    htmls.push("<br>");
                }
            }
        }
        htmls.push("<p style='padding-top: 15px'><small>Cell Browser dataset ID: "+datasetInfo.name+"</small></p>");

        $( "#pane1" ).html(htmls.join(""));

        buildMethodsPane(datasetInfo, desc);
        buildDownloadsPane(datasetInfo, desc);
        buildImagesPane(datasetInfo, desc);

        $("#tpOpenDialogTabs").tabs("refresh");
        //.tabs("option", "active", 0) does not do the color change of the tab so doing this instead
        $("#tabLink1").click();
        $("area").click( function(ev) {
            var dsName = ev.target.href.split("/").pop();
            loadDataset(gOpenDataset.name+"/"+dsName, true);
            $(".ui-dialog-content").dialog("close");
            ev.preventDefault();
        });

    }

    function buildListPanel(datasetList, noteSpace, listGroupHeight, leftPaneWidth, htmls, selName) {
        /* make a dataset list and append its html lines to htmls */
        htmls.push("<div id='tpDatasetList' class='list-group' style='width:400px; position:absolute; top:"+noteSpace+"; height:"+listGroupHeight+"px; overflow-y:scroll; width:"+leftPaneWidth+"px'>");
        if (!datasetList || datasetList.length===0) {
            alert("No datasets are available. Please make sure that at least one dataset does not set visibility=hide "+
                " or that at least one collection is defined. Problems? -> cells@ucsc.edu");
            return;
        }

        var selIdx = 0;
        for (var i = 0; i < datasetList.length; i++) {
            var dataset = datasetList[i];

            var clickClass = "tpDatasetButton";
            if (dataset.isCollection)
                clickClass = "tpCollectionButton";
            if (dataset.name===selName || (selName===undefined && i===0)) {
                clickClass += " active";
                selIdx = i;
            }

            var bodyPartStr = "";
            if (dataset.body_parts) {
                bodyPartStr = dataset.body_parts.join("|");
            }

            var line = "<a id='tpDatasetButton_"+i+"' data-body-parts='"+bodyPartStr+"' role='button' class='tpListItem list-group-item "+clickClass+"' data-datasetid='"+i+"'>"; // bootstrap seems to remove the id
            htmls.push(line);

            if (!dataset.isSummary)
                htmls.push('<button type="button" class="btn btn-primary btn-xs load-dataset" data-placement="bottom">Open</button>');

            if (dataset.sampleCount!==undefined) {
                var countDesc = prettyNumber(dataset.sampleCount);
                htmls.push("<span class='badge' style='background-color: #888'>"+countDesc+"</span>");
            }

            if (dataset.datasetCount!==undefined) {
                htmls.push("<span class='badge' style='background-color: #28a745'>"+dataset.datasetCount+" datasets</span>");
            }

            if (dataset.collectionCount!==undefined) {
                htmls.push("<span class='badge' style='background-color: #188725'>"+dataset.collectionCount+" collections</span>");
            }

            if (dataset.tags!==undefined) {
                for (var tagI = 0; tagI < dataset.tags.length; tagI++) {
                var tag = dataset.tags[tagI];
                htmls.push("<span class='badge'>"+tag+"</span>");
                }
            }
            htmls.push(dataset.shortLabel+"</a>");
            //if (db!==null && db.name===dataset.name)
                //activeIdx = i;
        }
        htmls.push("</div>"); // list-group
        return selIdx;
    }

    function getBodyParts(datasets) {
        /* return list of body_parts given a dataset array */
        var bpObj = {};
        for (let i=0; i < datasets.length; i++) {
            let ds = datasets[i];
            if (ds.body_parts===undefined)
                continue
            for (let bp of ds.body_parts)
                bpObj[bp] = true;
        }

        let bodyParts = keys(bpObj);
        bodyParts.sort();
        return bodyParts;
    }

    function filterDatasetsDom(onlyBps) {
        /* keep only datasets with a a body_tag in filtNames */

        let elList = $(".tpListItem");
        for (let el of elList) {
            let bpStr = el.getAttribute("data-body-parts");
            let found = false;
            if (!onlyBps || onlyBps.length==0)
                // if no filtering is specified just show everything
                found = true;
            else
                if (bpStr && bpStr!=="") {
                    let bps = bpStr.split("|");
                    if (onlyBps.length===0)
                        found = true;
                    else {
                        if (bps.indexOf("summary")!==-1) { // never filter the summary
                            found = true;
                        }
                        //if (bps.indexOf("all")!==-1) { // always match datasets with "all"
                            //found = true;
                        //}
                        else
                            for (let bp of onlyBps) {
                                if (bps.indexOf(bp)!==-1) {
                                    found = true;
                                    break;
                                }
                            }
                    }
                }

            if (found)
                el.style.display="";
            else
                el.style.display="none";
        }
    }

    function openDatasetDialog(openDsInfo, selName) {
    /* build dataset open dialog,
     * - openDsInfo is the currently open object or a collection.
     * - selName is the currently selected dataset in this list
     */

        var datasetList = [];
        var listGroupHeight = 0;
        var leftPaneWidth = 400;
        var title = "Choose Cell Browser Dataset";
        var noteSpace = "3em"; // space from top of dialog to info pane and tabs

        // inline functions
        function openCollOrDataset(selDatasetIdx) {
            /* click handler, opens either a collection or a dataset */
            var dsInfo = datasetList[selDatasetIdx];
            var datasetName = dsInfo.name;
            if (dsInfo.isCollection)
                showCollectionDialog(datasetName);
            else
                loadDataset(datasetName, true, dsInfo.md5);
            $(".ui-dialog-content").dialog("close");
            //changeUrl({"bp":null});
        }

        function connectOpenPane(selDatasetIdx, datasetList) {
            /* set all the click handlers for the left open dataset pane */
            $("button.list-group-item").eq(selDatasetIdx).css("z-index", "1000"); // fix up first overlap
            $("button.list-group-item").keypress(function(e) {
                // load the current dataset when the user presses Return
                if (e.which === '13') {
                    openCollOrDataset(selDatasetIdx);
                }
            });
            $(".list-group-item").click( function (ev) {
                selDatasetIdx = parseInt($(ev.target).data('datasetid')); // index of clicked dataset
                $(".list-group-item").removeClass("active");
                $('#tpDatasetButton_'+selDatasetIdx).bsButton("toggle"); // had to rename .button() in index.html
                var datasetInfo = datasetList[selDatasetIdx];
                openDatasetLoadPane(datasetInfo);
            });
            $(".list-group-item").dblclick( function(ev) {
                selDatasetIdx = parseInt($(this).data('datasetid'));
                openCollOrDataset(selDatasetIdx);
            });
            $(".load-dataset").click( function (ev) {
                ev.preventDefault();
                ev.stopPropagation();
                selDatasetIdx = parseInt($(this).parents('.list-group-item').data('datasetid'));
                openCollOrDataset(selDatasetIdx);
                return false;
            });
            $(".list-group-item").focus( function (event) {
                selDatasetIdx = parseInt($(event.target).data('datasetid')); // index of clicked dataset
                // bootstrap has a bug where the blue selection frame is hidden by neighboring buttons
                // Working around this here by bumping up the current z-index.
                $("button.list-group-item").css("z-index", "0");
                $("button.list-group-item").eq(selDatasetIdx).css("z-index", "1000");
            });
        }

        function onBodyPartChange(ev) {
            /* called when user changes body part list */
            let filtNames = $("#tpBodyCombo").val();
            // change the URL
            let filtArg = filtNames.join("_");
            changeUrl({"bp":filtArg});
            filterDatasetsDom(filtNames);
        }

        // -- end inline functions

        gOpenDataset = openDsInfo;
        var activeIdx = 0;
        var onlyInfo = false;

        datasetList = openDsInfo.datasets;

        if (datasetList===undefined)
            onlyInfo = true;

        var noteLines = [];

        // if this is a collection, not a dataset, change descriptive text in dialog
        if (datasetList && gOpenDataset.name!=="") {
            let dsCount = datasetList.length;
            title = 'Select one dataset from the collection "'+openDsInfo.shortLabel+'"';
            title = title.replace(/'/g, "&apos;");
            noteLines.push( "<p>The collection '"+openDsInfo.shortLabel+"' contains "+dsCount+" datasets. " +
                "Double-click or click 'Open' below.<br>To move between datasets later in the cell browser, " +
                "use the 'Collection' dropdown. </p>");

            changeUrl({"ds":openDsInfo.name.replace(/\//g, " ")}); // + is easier to type
        }

        let doFaceting = false;
        let filtList = [];
        if (openDsInfo.parents === undefined && openDsInfo.datasets !== undefined) {
            //noteLines.push("<span>Filter:</span>");
            let bodyParts = getBodyParts(openDsInfo.datasets);
            if (bodyParts.length!==0) {
                noteLines.push("<span style='margin-right:5px'>Filter datasets by organ:</span>");
                doFaceting = true;
                // some mirrors don't use the "body_parts" statement and don't need the faceting
                let selPar = getVarSafe("bp");
                if (selPar && selPar!=="")
                    filtList = selPar.split("_");

                buildComboBox(noteLines, "tpBodyCombo", bodyParts, filtList, "select organs...", 200, {multi:true});
            }
        }

        // create links to the parents of the dataset
        if (openDsInfo && openDsInfo.parents && !onlyInfo) {

            noteLines.push("Go back to: " );
            // make the back links
            let backLinks = [];
            let parents = openDsInfo.parents;
            for (let i=0; i<parents.length; i++) {
                let parentInfo = parents[i];
                let parName = parentInfo[0];
                let parLabel = parentInfo[1];
                let childName = null;
                if (i === parents.length-1)
                    childName = openDsInfo.name;
                else
                    childName = parents[i+1][0];
                backLinks.push("<span class='tpBackLink link' data-open-dataset='"+parName+"' data-sel-dataset='"+childName+"'>"+parLabel+"</span>");
            }
            noteLines.push(backLinks.join("&nbsp;&gt;&nbsp;"));
            noteSpace = "5em"; // TODO: redesign dialog to not have hard-coded spacing
        }

        if (onlyInfo)
            title = "Dataset Information";
        else {
            datasetList.unshift( {
                shortLabel:"Overview",
                name:openDsInfo.name,
                hasFiles:openDsInfo.hasFiles,
                body_parts:["summary"],
                isSummary:true,
                abstract:openDsInfo.abstract
            });
        }

        var winWidth = window.innerWidth - 0.05*window.innerWidth;
        var winHeight = window.innerHeight - 0.05*window.innerHeight;
        var tabsWidth = winWidth - leftPaneWidth - 70;
        listGroupHeight = winHeight - 100;

        var htmls = ["<div style='line-height: 1.1em'>"];
        htmls.push(noteLines.join(""));
        htmls.push("</div>");

        if (onlyInfo)
            leftPaneWidth = 0;
        else
            activeIdx = buildListPanel(datasetList, noteSpace, listGroupHeight, leftPaneWidth, htmls, selName);

        htmls.push("<div id='tpOpenDialogDatasetDesc' style='width:"+tabsWidth+"px; position:absolute; left: " + (leftPaneWidth + 20) + "px; top: "+noteSpace+"; border: 0'>");
        htmls.push("<div id='tpOpenDialogTabs' style='border: 0'>");
        htmls.push("<ul class='nav nav-tabs'>");
        htmls.push("<li class='active'><a class='tpDatasetTab' id='tabLink1' data-toggle='tab' href='#pane1'>Abstract</a></li>");
        htmls.push("<li><a class='tpDatasetTab' id='tabLink2' data-toggle='tab' href='#pane2'>Methods</a></li>");
        htmls.push("<li><a class='tpDatasetTab' id='tabLink3' data-toggle='tab' style='display:none' href='#pane3'>Data Download</a></li>");
        htmls.push("<li><a class='tpDatasetTab' id='tabLinkImg' data-toggle='tab' href='#paneImg'>Images</a></li>");
        htmls.push("</ul>");

        htmls.push("<div id='pane1' class='tpDatasetPane tab-pane'>");
        htmls.push("<p>Loading abstract...</p>");
        htmls.push("</div>");

        htmls.push("<div id='pane2' class='tpDatasetPane tab-pane'>");
        htmls.push("<p>Loading methods...</p>");
        htmls.push("</div>");

        htmls.push("<div id='pane3' class='tpDatasetPane tab-pane'>");
        htmls.push("<p>Loading download instructions...</p>");
        htmls.push("</div>");

        htmls.push("<div id='paneImg' class='tpDatasetPane tab-pane'>");
        htmls.push("<p>Loading image data...</p>");
        htmls.push("</div>");

        htmls.push("</div>"); // tpOpenDialogTabs

        htmls.push("</div>"); // tpOpenDialogDatasetDesc

        //htmls.push("<div id='tpSelectedId' data-selectedid='0'>"); // store the currently selected datasetId in the DOM
        var selDatasetIdx = 0;

        var buttons = [];
        if (db!==null) {
            var cancelLabel = "Cancel";
            if (onlyInfo)
                cancelLabel = "Close";
            buttons.push( {
                text: cancelLabel,
                click: function() {
                    $( this ).dialog( "close" );
                    if (openDsInfo.isCollection)
                        openDatasetDialog(openDsInfo, null); // show top-level dialog
                }
            });
        }

        $(".ui-dialog-content").dialog("close"); // close the last dialog box

        showDialogBox(htmls, title, {width: winWidth, height:winHeight, buttons: buttons});

        $("#tpOpenDialogTabs").tabs();

        if (doFaceting) {
            activateCombobox("tpBodyCombo", 200);
            $("#tpBodyCombo").change( onBodyPartChange );
        }

        $('.tpBackLink').click( function(ev) {
            let openDatasetName = $(ev.target).attr('data-open-dataset');
            let selDatasetName = $(ev.target).attr('data-sel-dataset');
            loadCollectionInfo(openDatasetName, function(newCollInfo) {
                openDatasetDialog(newCollInfo, selDatasetName);
            });
            changeUrl({"ds":openDatasetName.replace(/\//g, " ")});
        });

        var focused = document.activeElement;
        var scroller = $("#tpDatasetList").overlayScrollbars({ });
        $(focused).focus();


        $("#tabLink1").tab("show");

        if (activeIdx!==null && !onlyInfo) {
            if (activeIdx!==0)
                scroller.scroll($("#tpDatasetButton_"+activeIdx)); // scroll left pane to current button
            $("tpDatasetButton_"+activeIdx).addClass("active");
        }

        if (getVarSafe("ds")===undefined) // only filter on the top level
            filterDatasetsDom(filtList);
        connectOpenPane(selDatasetIdx, datasetList);
        // finally, activate the default pane and load its html
        openDatasetLoadPane(openDsInfo);
    }


    function onSelChange(selection) {
    /* called each time when the selection has been changed */
        var cellIds = [];
        selection.forEach(function(x) {cellIds.push(x)});
        $("#tpSetBackground").parent("li").removeClass("disabled");

        if (cellIds.length===0 || cellIds===null) {
            clearMetaAndGene();
            clearSelectionState();
            $("#tpSetBackground").parent("li").addClass("disabled");
        } else if (cellIds.length===1) {
            $("#tpHoverHint").hide();
            $("#tpSelectHint").show();
            var cellId = cellIds[0];
            var cellCountBelow = cellIds.length-1;
            updateMetaBarCustomFields(cellId);
            db.loadMetaForCell(cellId, function(ci) {
                updateMetaBarOneCell(ci, cellCountBelow);
            }, onProgress);
        } else {
            $("#tpHoverHint").hide();
            $("#tpSelectHint").show();
            updateMetaBarManyCells(cellIds);
        }

        updateGeneTableColors(cellIds);
        if ("geneSym" in gLegend)
            buildViolinPlot();

        var cols = renderer.col.arr;
        var selectedLegends = {};
        for (var i = 0; i < gLegend.rows.length; i++) {
            selectedLegends[i] = 0;
        }
        selection.forEach(function(cellId) {
            selectedLegends[cols[cellId]]++;
        });
        for (var i = 0; i < gLegend.rows.length; i++) {
            if (selectedLegends[i] == gLegend.rows[i].count) {
                $("#tpLegendCheckbox_" + i).prop("checked", true);
            } else {
                $("#tpLegendCheckbox_" + i).prop("checked", false);
            }
        }
        updateLegendGrandCheckbox();
    }

    function onSaveAsClick() {
    /* File - Save Image as ... */
        var canvas = $("canvas")[0];
        canvas.toBlob(function(blob) { saveAs( blob , "cellBrowser.png"); } , "image/png");
    }

    function onSelectAllClick() {
    /* Edit - select all visible*/
        clearSelectionState();
        renderer.selectClear();
        renderer.selectVisible();
        renderer.drawDots();
    }

    function onSelectNoneClick() {
    /* Edit - Select None */
        clearSelectionState();
        renderer.selectClear();
        renderer.drawDots();
    }

    function onSelectInvertClick() {
    /* Edit - Invert selection */
        clearSelectionState();
        renderer.selectInvert();
        renderer.drawDots();
    }

    function buildOneComboboxRow(htmls, comboWidth, rowIdx, queryExpr) {
        /* create one row of combobox elements in the select dialog */
        htmls.push('<div class="tpSelectRow" id="tpSelectRow_'+rowIdx+'">');
        // or &#x274e; ?
        htmls.push('<span style="cursor: default; font-size:1.2em" id="tpSelectRemove_'+rowIdx+'">&#xd7;</span>&nbsp;'); // unicode: mult sign
        htmls.push('<select name="type" id="tpSelectType_'+rowIdx+'">');
        htmls.push('<option value="meta">Cell annotation field</option><option value="expr">Expression of gene </option>');
        htmls.push('</select>');

        buildMetaFieldCombo(htmls, "tpSelectMetaComboBox_"+rowIdx, "tpSelectMetaCombo_"+rowIdx, 0);

        var id = "tpSelectGeneCombo_"+rowIdx;
        htmls.push('<select style="width:'+comboWidth+'px" id="'+id+'" placeholder="gene search..." class="tpCombo">');
        htmls.push('</select>');

        htmls.push('<select name="operator" id="tpSelectOperator_'+rowIdx+'">');
        htmls.push('<option value="eq" selected>is equal to</option>');
        htmls.push('<option value="neq" selected>is not equal to</option>');
        htmls.push('<option value="gt">is greater than</option>');
        htmls.push('<option value="lt">is less than</option>');
        htmls.push('</select>');

        htmls.push('<input id="tpSelectValue_'+rowIdx+'" type="text" name="exprValue">');
        htmls.push('</input>');

        htmls.push('<select style="max-width: 300px" id="tpSelectMetaValueEnum_'+rowIdx+'" name="metaValue">');
        // if ("m" in queryExpr) {
          //   var op = getQueryOp(queryExpr);
            // var fieldName = query["m"];
            //var fieldValue = queryExpr[op];
            //var metaInfo = findMetaInfo(fieldName);
            //for (var valIdx = 0; valIdx < metaInfo.valCounts.length; valIdx++) {
            //    htmls.push('<option value="'+valIdx+'">'+metaInfo.valCounts[i][0]+'</option>');
           // }
       // }
        htmls.push('</select>');
        htmls.push('</div>'); // tpSelectRow_<rowIdx>

        htmls.push('<p>');
    }

    function findCellsUpdateMetaCombo(rowIdx, fieldIdx) {
        /* given the row and the ID name of the field, setup the combo box row */
        var metaInfo = db.getMetaFields()[fieldIdx];
        var valCounts = metaInfo.valCounts;
        var shortLabels = metaInfo.ui.shortLabels;
        $('#tpSelectMetaCombo_'+rowIdx).val("tpMetaVal_"+fieldIdx).trigger('chosen:updated'); // update the meta dropdown

        if (valCounts===undefined) {
            // this is a numeric field
            $('#tpSelectValue_'+rowIdx).val("");
            $('#tpSelectValue_'+rowIdx).show();
            $('#tpSelectMetaValueEnum_'+rowIdx).hide();
        } else {
            // it's an enum field
            $('#tpSelectValue_'+rowIdx).hide();
            $('#tpSelectMetaValueEnum_'+rowIdx).empty();
            for (var i = 0; i < valCounts.length; i++) {
                //var valName = valCounts[i][0];
                var valLabel = shortLabels[i];
                $('#tpSelectMetaValueEnum_'+rowIdx).append("<option value='"+i+"'>"+valLabel+"</option>");
            }
            $('#tpSelectMetaValueEnum_'+rowIdx).show();
        }
    }

    function findCellsUpdateRowType(rowIdx, rowType) {
        if (rowType === "meta") {
            $("#tpSelectType_"+rowIdx).val("meta");
            $("#tpSelectGeneCombo_"+rowIdx).next().hide();
            $("#tpSelectValue_"+rowIdx).hide();
            $("#tpSelectMetaComboBox_"+rowIdx).show();
            $("#tpSelectMetaValueEnum_"+rowIdx).show();
        } else {
            $("#tpSelectType_"+rowIdx).val("expr");
            $("#tpSelectGeneCombo_"+rowIdx).next().show();
            $("#tpSelectValue_"+rowIdx).show();
            $("#tpSelectMetaComboBox_"+rowIdx).hide();
            $("#tpSelectMetaValueEnum_"+rowIdx).hide();
        }
    }

    function connectOneComboboxRow(comboWidth, rowIdx, query) {
        /* Filter dialog. Call the jquery inits and setup the change listeners for a combobox row */
        /* Yes, a UI framework, like react or angular, would be very helpful here */

        // first of all: check if the meta name actually exists in this dataset still
        let metaInfo = null;
        var metaName = query["m"];
        if (metaName!==undefined) {
            metaInfo = db.findMetaInfo(metaName);
            if (metaInfo===null)
                return;
        }

        // auto-suggest for gene searches
        $('#tpSelectGeneCombo_'+rowIdx).selectize({
                "labelField" : 'text',
                "valueField" : 'id',
                "searchField" : 'text',
                "load" : comboLoadGene,
       });
        activateCombobox("tpSelectMetaCombo_"+rowIdx, comboWidth);

        $('#tpSelectRemove_'+rowIdx).click( function(ev) {
            //console.log(ev);
            var rowToDel = (this.id).split("_")[1];
            $("#tpSelectRow_"+rowToDel).remove();
        });

        //$("#tpSelectGeneCombo_"+rowIdx).next().hide();
        //$("#tpSelectValue_"+rowIdx).hide();


        var rowType = "gene";
        var op = getQueryOp(query);
        if (metaName===undefined) {
            // this is a gene query
            findCellsUpdateRowType(rowIdx, rowType);
            selectizeSetValue("#tpSelectGeneCombo_"+rowIdx, query["g"]);
            $("#tpSelectValue_"+rowIdx).val(query[op]);
        } else {
            // it's a meta query
            rowType = "meta";
            findCellsUpdateRowType(rowIdx, rowType);
            findCellsUpdateMetaCombo(rowIdx, metaInfo.index);
            var enumIdx = findMetaValIndex(metaInfo, query[op]);
            $("#tpSelectMetaValueEnum_"+rowIdx).val(enumIdx);
        }
        $("#tpSelectOperator_"+rowIdx).val(op);

        $('#tpSelectMetaCombo_'+rowIdx).change(function(ev) {
            // when the user changes the meta field, update the list of meta field values in the dropdown
            var selVal = this.value;
            var fieldIdx = parseInt(selVal.split("_")[1]);
            findCellsUpdateMetaCombo(rowIdx, fieldIdx);
        });

        $('#tpSelectType_'+rowIdx).change(function(ev) {
            // when the user changes the gene expression / meta dropdown, hide/show the
            // respective other dropdowns
            var rowType = this.value;
            findCellsUpdateRowType(rowIdx, rowType);
        });
    }

    function readSelectForm() {
        /* convert the current state of the dialog box to a short string and return it */
        // example: [{"g":"PITX2", "gt":0.05}, {"m":"Cluster", "eq":"cluster 2"}]
        // XX TODO: non-enum meta data fields ???
        var queries = [];

        var rowCount = $(".tpSelectRow").length;
        for (var rowIdx=0;  rowIdx<rowCount; rowIdx++) {
            var query = {};
            var op = $('#tpSelectOperator_'+rowIdx).val();

            var queryType = $('#tpSelectType_'+rowIdx).val();
            if (queryType==="expr") {
                var gene = $('#tpSelectGeneCombo_'+rowIdx).val();
                var val = $('#tpSelectValue_'+rowIdx).val();
                query["g"] = gene;
                query[op] = val;
            }
            else {
                var metaValTag = $('#tpSelectMetaCombo_'+rowIdx).val();
                var metaIdx = parseInt(metaValTag.split("_")[1]);
                var metaInfo = db.conf.metaFields[metaIdx];
                var metaName = metaInfo.name;
                query["m"] = metaName;

                var opVal = null;
                var selOp = null;
                if (metaInfo.type==="enum") {
                    let selVal = $('#tpSelectMetaValueEnum_'+rowIdx).val();
                    var valIdx = parseInt(selVal);
                    opVal = db.conf.metaFields[metaIdx].valCounts[valIdx][0];
                } else {
                    let selVal = $('#tpSelectValue_'+rowIdx).val();
                    opVal = parseFloat(selVal);
                }

                query[op] = opVal;
            }
            queries.push(query);
        }
        return queries;
    }

    function greaterThan(x, y) { return (x>y); }
    function lessThan(x, y) { return (x<y); }
    function equals(x, y) { return (x===y); }
    function notequals(x, y) { return (x!==y); }

    function getQueryOp(query) {
        if ("eq" in query) return "eq";
        if ("neq" in query) return "neq";
        if ("gt" in query) return "gt";
        if ("lt" in query) return "lt";
    }

    function makeFuncAndVal(query) {
        /* return a comparator function and the value given a query object */
        var compFunc = equals;
        var val = query["eq"];

        if ("lt" in query) {
            compFunc = lessThan;
            val = query["lt"];
        }
        if ("gt" in query) {
            compFunc = greaterThan;
            val = query["gt"];
        }
        if ("neq" in query) {
            compFunc = notequals;
            val = query["neq"];
        }

        return [compFunc, val];
    }

    function searchArrayForFuncAndVal(arr, funcAndVal) {
        /* given an array and function and a value, return an array with the indices the matching array elements */
        var compFunc = funcAndVal[0];
        var compVal  = funcAndVal[1];
        var selCells = [];
        for (var i=0; i<arr.length; i++) {
            if (compFunc(arr[i], compVal))
                selCells.push(i);
        }
        return selCells;
    }

    function findCellsMatchingQueryList(queries, onDone) {
        /* given a list of dicts, return the identifiers of the matching cells */
        /* example: [{"g":"PITX2", "gt":0.05}, {"m":"Cluster", "eq":"cluster 2"}] */

        var doneQueries = 0;
        var queryResults = [];

        function allQueriesDone() {
            // XX this could use a promise framework
            if (queryResults.length===1)
                onDone(queryResults[0]);
            else
                onDone(intersectArrays(queryResults));
        }

        function gotGeneVec(exprArr, sym, desc, funcVal) {
            funcVal[1] = Number(funcVal[1]); // for gene queries, must be a number, not string
            var selCells = searchArrayForFuncAndVal(exprArr, funcVal);
            queryResults.push(selCells);
            doneQueries++;
            if (doneQueries===queries.length)
                allQueriesDone();
        }

        function gotMetaArr(metaArr, metaInfo, funcVal) {
            /* filter meta data array with funcVal = function + value  */
            if (metaInfo.origVals)
                metaArr = metaInfo.origVals; // numerical meta fields have the original numbers stored
            var selCells = searchArrayForFuncAndVal(metaArr, funcVal);
            queryResults.push(selCells);
            doneQueries++;
            if (doneQueries===queries.length)
                allQueriesDone();
        }

        var selCells = [];
        for (var i=0; i < queries.length; i++) {
            var query = queries[i];
            var funcVal = makeFuncAndVal(query); // [0] = function to compare, [1] = value for comparison
            if ("g" in query) {
                db.loadExprVec(query.g, gotGeneVec, null, funcVal);
            }
            else {
                var fieldName = query.m;
                var fieldIdx = cbUtil.findIdxWhereEq(db.conf.metaFields, "name", fieldName);
                var findVal = funcVal[1];

                var metaInfo = db.getMetaFields()[fieldIdx];
                if (metaInfo.type==="enum")
                    findVal = findMetaValIndex(metaInfo, findVal);

                let searchDesc = [funcVal[0], findVal];

                if (metaInfo.origVals)
                    // for numeric fields, the raw data is already in memory
                    gotMetaArr(metaInfo.origVals, metaInfo, searchDesc)
                else
                    // other fields may not be loaded yet
                    db.loadMetaVec(metaInfo, gotMetaArr, null, searchDesc);
            }
        }
    }

    function makeSampleQuery() {
        /* find first enum meta field and return a query for its name and its first value */
        var metaFieldInfo = db.conf.metaFields;
        for (var i = 0; i < metaFieldInfo.length; i++) {
            var field = metaFieldInfo[i];
            if (field.type!=="enum")
                continue;
            var fieldName = field.name;
            var val1 = field.valCounts[0][0];
            return {"m":fieldName, "eq":val1};
        }
    }

    function addNewAnnotation(fieldLabel, newMetaValue, cellIds) {
        var metaInfo;
        let cellCount = db.conf.sampleCount;
        if (!db.getMetaFields()[0].isCustom) {
            // add a new enum meta field
            metaInfo = {
                name:"custom",
                label: fieldLabel,
                type: "enum",
                arr : Array.from(new Uint8Array(cellCount)), // cannot JSON-serialize Typed Arrays
                ui : {
                    shortLabels : ["No annotation"]
                },
                valCounts : [ ["No annotation", cellCount]]
            }
            db.addCustomMetaField(metaInfo);
            rebuildMetaPanel();
            activateTab("meta");
        } else
            metaInfo = db.getMetaFields()[0];

        metaInfo.ui.shortLabels.push( newMetaValue );
        let newValIdx = metaInfo.valCounts.length;
        metaInfo.valCounts.push( [newMetaValue, cellIds.length]);
        // update the "No annotation" count
        let noAnnotCount = metaInfo.valCounts[0][1];
        metaInfo.valCounts[0][1] = noAnnotCount - cellIds.length;

        let arr = metaInfo.arr;
        for (let i=0; i<cellIds.length; i++)
            arr[cellIds[i]] = newValIdx;
        // need to update the value histogram
        db.metaHist[metaInfo.name] = makeFieldHistogram(metaInfo, cellIds, arr);
        updateMetaBarManyCells(cellIds); // redo, as we rebuilt the meta panel

        var jsonStr = JSON.stringify(metaInfo);
        var comprStr = LZString.compress(jsonStr);
        localStorage.setItem(db.name+"|custom", comprStr);
        }

    function onSelectNameClick() {
        /* Edit > Name selection */

        let title = "Annotate selected "+gSampleDesc+"s";

        let htmls = [];
        htmls.push('<p>There are '+renderer.getSelection().length+' '+gSampleDesc+' in the current selection.</p>');

        htmls.push('<p><b>Name of annotation field</b>:<br>');
        htmls.push('<input class="tpDialogInput" id="tpFieldLabel" type="text" value="My custom annotations"></p>');
        htmls.push('<p><b>Annotate selected cells as:</b><br>');
        htmls.push('<input class="tpDialogInput" id="tpMetaVal" type="text"></p>');
        htmls.push('<p>Remove annotations later by clicking <b>Tools > Remove all annotations</b>.</p>');

        var dlgHeight = 400;
        var dlgWidth = 800;
        var buttons = [
            //"Close and remove all annotations" : function() {
                //db.getMetaFields().shift();
                //rebuildMetaPanel();
                //localStorage.removeItem(db.name+"|custom");
                //resetCustomAnnotations();
            //},
        {text:"OK", click: function() {
                let fieldLabel = $('#tpFieldLabel').val();
                if (fieldLabel==="")
                    fieldLabel = "My custom annotations";
                let newMetaValue = $("#tpMetaVal").val();
                if (newMetaValue==="")
                    return;

                addNewAnnotation(fieldLabel, newMetaValue, renderer.getSelection());
                $( this ).dialog( "close" );
                colorByMetaField("custom");
            }
        }];
        showDialogBox(htmls, title, {showClose:true, height:dlgHeight, width:dlgWidth, buttons:buttons});
        $("#tpMetaVal").focus();
        return true;
    }

    function onBackgroudSetClick() {
// Tools -> Set cells as background
        if ($("#tpSetBackground").parent("li").hasClass("disabled")) {
            return;
        }

        background = renderer.getSelection();
        $("#tpResetBackground").parent("li").removeClass("disabled");
        if ("geneSym" in gLegend)
            buildViolinPlot();
    }

    function onBackgroudResetClick() {
// Tools -> Reset background cells
        background = null;
        $("#tpResetBackground").parent("li").addClass("disabled");
        if ("geneSym" in gLegend)
            buildViolinPlot();
    }

    function saveQueryList(queryList) {
        var queryStr = JSURL.stringify(queryList);
        changeUrl({'select':queryStr});
        localStorage.setItem(db.name+"|select", queryStr);
    }

    function selectByQueryList(queryList) {
        /* select cells defined by query list, save to local storage and URL and redraw */
            findCellsMatchingQueryList(queryList, function(cellIds) {
                if (cellIds.length===0) {
                    alert("No matching "+gSampleDesc+"s.");
                } else {
                    renderer.selectSet(cellIds);
                    saveQueryList(queryList);
                }
            });
    }

    function onFindCellsClick() {
    /* Edit - Find cells */

        var dlgHeight = 400;
        var dlgWidth = 800;
        var buttons =
        [
            //{
                //text:"Cancel",
                //click:function() {
                    //// save state even if user presses cancel  - good idea?
                    //var queryStr = JSURL.stringify(queryList);
                    //var queryList = readSelectForm();
                    //localStorage.setItem(db.name+"|select", queryStr);
                    //$(this).dialog("close");
                //}
            //},
            {
                text:"OK",
                click: function() {
                    var queryList = readSelectForm();
                    selectByQueryList(queryList);
                    $("#tpDialog").dialog("close");
                    renderer.drawDots();
                }
            }
        ];

        var htmls = [];

        // build from current query or create a sample query
        var queries = [];
        var queryStr = getVar("select");

        if (queryStr===undefined)
            queryStr = localStorage.getItem(db.name+"|select");

        if (queryStr===undefined || queryStr===null)
            queries = [makeSampleQuery()];
        else {
            queries = JSURL.parse(queryStr);
        }

        var comboWidth = 250;

        var query;
        for (var i=0; i < queries.length; i++) {
            query = queries[i];
            buildOneComboboxRow(htmls, comboWidth, i, query);
        }

        htmls.push("<div id='tpSelectAddRowLink' class='link'>Add another search criterion</div>");

        showDialogBox(htmls, "Find cells based on annotation or gene expression", {showClose:true, height:dlgHeight, width:dlgWidth, buttons:buttons});


        for (i=0; i < queries.length; i++) {
            query = queries[i];
            connectOneComboboxRow(comboWidth, i, query);
        }

        var rowIdx = queries.length+1;
        $('#tpSelectAddRowLink').click( function(ev) {
            var htmls = [];
            var rowIdx = $(".tpSelectRow").length;
            var newRowQuery = makeSampleQuery();
            buildOneComboboxRow(htmls, comboWidth, rowIdx, newRowQuery);
            $(htmls.join("")).insertBefore("#tpSelectAddRowLink");
            connectOneComboboxRow(comboWidth, rowIdx, newRowQuery);
        });
    }

    function onMarkClick() {
        /* Edit - mark selection (with arrows) */
        if (gCurrentDataset.markedCellIds===undefined)
            gCurrentDataset.markedCellIds = {};

        var markedIds = gCurrentDataset.markedCellIds;

        var selIds = keys(gSelCellIds);
        if (selIds.length>100) {
            warn("You cannot mark more than 100 "+gSampleDesc+"s");
            return;
        }

        for (var i = 0; i < selIds.length; i++) {
            var selId = selIds[i];
            markedIds[selId] = true;
        }

        plotDots();
        renderer.render(stage);
    }

    function onMarkClearClick() {
        gCurrentDataset.markedCellIds = {};
        plotDots();
        renderer.render(stage);
    }

    function cartSave(db) {
        /* save db.cart dataset long-term changes that may be shared with others, like colors, labels, annotations, etc
         * For now, save to localStorage and also to the URL. */

        var datasetName = db.name;
        var data = db.cart;
        var key = "cart";

        var fullKey = datasetName+"###"+key;
        var jsonStr = JSON.stringify(data);
        var comprStr = LZString.compress(jsonStr);
        var uriStr = LZString.compressToEncodedURIComponent(jsonStr);
        localStorage.setItem(fullKey, comprStr);
        var urlData = {};
        if (isEmpty(data))
            uriStr = null;
        urlData[key] = uriStr;
        changeUrl(urlData);
        console.log("Saving state: ", data);
    }

    function createMetaUiFields(db) {
        /* This function changes db.metaFields[fieldName],
         * it adds: .ui.shortLabels, .ui.longLabels, ui.palette and .ui.colors;
         * ui info fields hold the final data as shown in the ui, they're calculated when the cart is loaded.
         * apply changes like labels/color/etc stored the userMeta object to db.conf.metaFields.
         * Potentially clean up the changes and recreate the cart object.
         * Currently, this only does something for enum fields.
         * */

        if (db.cart===undefined)
            db.cart = {};
        var userMeta = db.cart;
        if (userMeta===null)
            alert("the 'cart' argument in the URL is invalid. Please remove cart=xxxx from the URL and reload");
        var metaFields = db.conf.metaFields;

        for (var metaIdx = 0; metaIdx < metaFields.length; metaIdx++) {
            var metaInfo = metaFields[metaIdx];
            var fieldChanges = userMeta[metaInfo.name] || {};

            if (metaInfo.type!=="enum") {
                metaInfo.ui = {};
                continue;
            }

            // create shortLabels
            var shortLabels = null;
            var oldCounts = metaInfo.valCounts;
            if (oldCounts) {
                shortLabels = [];
                for (var i = 0; i < oldCounts.length; i++)
                    shortLabels.push(oldCounts[i][0]);
                var newLabels = fieldChanges.shortLabels;
                shortLabels = copyNonEmpty(newLabels, shortLabels);
            }

            // create the long labels
            var longLabels = [];
            if ("longLabels" in metaInfo)
                longLabels = cloneArray(metaInfo.longLabels);
            else
                longLabels = cloneArray(shortLabels);
            longLabels = copyNonEmpty(fieldChanges.longLabels, longLabels);

            // create the colors: configured colors override default colors and cart overrides those
            var colors = makeColorPalette(cDefQualPalette, metaInfo.valCounts.length);
            if ("colors" in metaInfo)
                copyNonNull(metaInfo.colors, colors);
            var newColors = fieldChanges.colors;
            colors = copyNonEmpty(newColors, colors);

            var ui = {};
            ui.colors = newColors;
            ui.longLabels = longLabels;
            ui.shortLabels = shortLabels;
            metaInfo.ui = ui;
        }

        //var delFields = [];
        //var delAttrs = [];
        //if (metaInfo===null) { // field does not exist anymore
            //delFields.push(fieldName);
            //continue;
        //}

        // remove all fields and attributes that were found to be invalid in the current state
        // so we don't accumulate crap
        //var cleanedFields = [];
        //for (var i = 0; i < delAttrs.length; i++) {
            //var fieldName = delAttrs[i][0];
            //var attrName = delAttrs[i][1];
            //delete userMeta[fieldName][attrName];
            //cleanedFields.push(fieldName);
        //}

        //for (var i = 0; i < delFields.length; i++) {
            //var fieldName = delFields[i];
            //delete userMeta[fieldName];
            //cleanedFields.push(fieldName);
        //}
        //if (delAttrs.length!==0)
            //warn("You had previously changed labels or colors or annotations but the dataset has been updated since then. "+
                //"As a result, your annotations had to be removed. This concerned the following annotation fields: "+
                //cleanedFields.join(", "));
    }

    function cartFieldArrayUpdate(db, metaInfo, key, pos, val) {
        /* write val into db.cart[fieldName][key][pos], save the dataset cart and apply it.
         * db.cart[fieldName][key] is an array of arrLen
         * */
        var cart = db.cart;
        var fieldName = metaInfo.name;
        var arrLen = metaInfo.valCounts.length;

        if (!(fieldName in cart))
            cart[fieldName] = {};

        // init array
        if (!(key in cart[fieldName])) {
            var emptyArr = [];
            for (var i=0; i<arrLen; i++)
                emptyArr.push("");
            cart[fieldName][key] = emptyArr;
        }

        cart[fieldName][key][pos] = val;
        cartSave(db);
        createMetaUiFields(db);
    }

    function cartOverwrite(db, key, val) {
        /* write val into db.cart[key][attr], save the dataset cart and apply it.
         * attrName can be null in which case db.cart[key] will be replaced with val  */
        var cart = db.cart;
        if (!(key in db.cart))
            cart[key] = {};

        cart[key] = val;

        cartSave(db);
        createMetaUiFields(db);
    }

    function cartLoad(db) {
        /* load the entire dataset cart from the URL or - if there is no cart on the URL - from localStorage */
        var datasetName = db.name;
        var cart = {};
        var key = "cart";
        var uriStr = getVar(key, null);
        var jsonStr;
        if (uriStr!==null) {
            jsonStr = LZString.decompressFromEncodedURIComponent(uriStr);
            cart = JSON.parse(jsonStr);
            console.log("Loading cart from URL: ", cart);
        }
        else {
            var fullKey = datasetName+"###"+key;
            var comprStr = localStorage.getItem(fullKey);
            if (comprStr) {
                jsonStr = LZString.decompress(comprStr);
                cart = JSON.parse(jsonStr);
                console.log("Loading cart from local storage: ", cart);
            }
        }
        db.cart = cart;
        createMetaUiFields(db);

    }

    function onRunClusteringClick() {
        /* not used yet */
        var myArray = new ArrayBuffer(512);
        var longInt8View = new Uint8Array(myArray);

        // generate some data
        for (var i=0; i< longInt8View.length; i++) {
        longInt8View[i] = i % 256;
        }

        //let url = "http://localhost:5050/upload";
        let url = "http://localhost:5050/bin";
        var xhr = new XMLHttpRequest();
        //xhr.open("POST", url, false);
        xhr.open("GET", url, true);
        xhr.responseType = "arraybuffer";
        xhr.onload = function() { var buf = xhr.response; console.log(buf)};
        xhr.send(null);
        //xhr.send(myArray);
    }

    function resetCustomAnnotations() {
        db.removeAllCustomAnnots();
        rebuildMetaPanel();
        changeUrl({"meta":null});
        colorByDefaultField();
        localStorage.removeItem(db.name+"|custom");
    }

    function onCustomAnnotationsClick() {
        /* */
        if (!db.getMetaFields()[0].isCustom) {
            alert("You have currently no custom annotations. Select a few cells and use Edit > Name selection "+
                    "to create a custom annotation field.");
        } else {
            resetCustomAnnotations();
        }
    }

    function onRenameClustersClick() {
    /* Tools - Rename Clusters */
        var htmls = [];
        htmls.push("<p>Change labels below. To keep the old name, leave the 'New Name' cell empty. You cannot modify the 'Orig. Name' column.</p>");
        //htmls.push('<p>To rename a single cluster without this dialog: click onto it in the legend, then click its label.</p>');
        htmls.push('<div id="tpGrid" style="width:100%;height:500px;"></div>');
        var title = "Rename Clusters";
        var dlgHeight = window.innerHeight - 100;
        var dlgWidth = 650;

        var clusterField = db.conf.labelField;

	var data = [];

        var buttons = [

        { text: "Empty All",
          click : function() {
                for (var i = 0; i < data.length; i++) {
                    var row = data[i];
                    row["newName"] = "";
                    row["mouseOver"] = "";
                    row["color"] = "";
                }
                grid.invalidate();
            }
        } ,
        {
            text: "OK",
            click: function() {
                Slick.GlobalEditorLock.commitCurrentEdit(); // save currently edited cell to data

                var shortLabels = [];
                var longLabels = [];
                var colors = [];

                for (var i = 0; i < data.length; i++) {
                    var row = data[i];
                    if (row["newName"]===row["origName"])
                        shortLabels.push("");
                    else
                        shortLabels.push(row["newName"]);

                    longLabels.push(row["mouseOver"]);
                    colors.push(row["color"]);
                }

                var fieldMeta = {};

                if (!allEmpty(shortLabels))
                    fieldMeta["shortLabels"] = shortLabels;
                if (!allEmpty(longLabels))
                    fieldMeta["longLabels"] = longLabels;
                if (!allEmpty(colors))
                    fieldMeta["colors"] = colors;

                cartOverwrite(db, clusterField, fieldMeta);
                var metaInfo = db.findMetaInfo(clusterField);

                renderer.setLabels(metaInfo.ui.shortLabels);

                // only need to update the legend if the current field is shown
                if (gLegend.type==="meta" && gLegend.metaInfo.name===clusterField) {
                    //var shortLabels = findMetaInfo(clusterField).ui.shortLabels;
                    legendUpdateLabels(clusterField);
                    buildLegendBar();
                }

                $( this ).dialog( "close" );
                renderer.drawDots();
            }
        },

        ];

        showDialogBox(htmls, title, {showClose:true, height:dlgHeight, width:dlgWidth, buttons:buttons});

        var columns = [
        	{id: "origName", width:100, maxWidth: 200, name: "Orig. Name", field: "origName"},
        	{id: "newName", width:150, maxWidth: 200, name: "New Name", field: "newName", editor: Slick.Editors.Text},
        	{id: "mouseOver", width:200, maxWidth: 300, name: "Mouseover Label", field: "mouseOver", editor: Slick.Editors.Text},
        	//{id: "color", width: 80, maxWidth: 120, name: "Color Code", field: "color", editor: Slick.Editors.Text}
        ];

	var options = {
            editable: true,
	    enableCellNavigation: true,
	    enableColumnReorder: false,
            enableAddRow: true,
            //asyncEditorLoading: false,
            autoEdit: true
	};

        var metaInfo = db.findMetaInfo(clusterField);

        var fieldChanges = db.cart[clusterField] || {};

        var shortLabels = fieldChanges["shortLabels"];
        var longLabels = fieldChanges["longLabels"];
        var colors = fieldChanges["colors"];

        for (var i = 0; i < metaInfo.valCounts.length; i++) {
            var shortLabel = "";
            if (shortLabels)
                shortLabel = shortLabels[i];

            var longLabel = "";
            if (longLabels)
                longLabel = longLabels[i];

            var color = "";
            if (colors)
                color = colors[i];

            var valCount = metaInfo.valCounts[i];
            var valRow = {
                "origName":valCount[0],
                "newName": shortLabel,
                "mouseOver": longLabel,
                "color": color
                };
            data.push(valRow);
        }

	var grid = new Slick.Grid("#tpGrid", data, columns, options);
        grid.setSelectionModel(new Slick.CellSelectionModel());
        grid.onAddNewRow.subscribe(function (e, args) {
          var item = args.item;
          grid.invalidateRow(data.length);
          data.push(item);
          grid.updateRowCount();
          grid.render();
        });
    }

    function onSelectByIdClick() {
        /* Edit - select cells by ID */

        function onSearchDone(searchRes) {
            var idxArr = searchRes[0];
            var notFoundIds = searchRes[1];

            if (notFoundIds.length!==0) {
                $('#tpNotFoundIds').text("Could not find these IDs: "+notFoundIds.join(", "));
                $('#tpNotFoundHint').text("Please fix them and click the OK button to try again.");
            } else
                $( "#tpDialog" ).dialog( "close" );

            renderer.selectSet(idxArr);
            renderer.drawDots();
        }

        var dlgHeight = 500;
        var dlgWidth = 500;
        var htmls = [];
        var buttons =
        [
            {
            text:"OK",
            click : function() {
                var idListStr = $("#tpIdList").val();
                idListStr = idListStr.trim().replace(/\r\n/g,"\n");
                var idList = idListStr.split("\n");
                var re = new RegExp("\\*");

                var hasWildcards = $("#tpHasWildcard")[0].checked;
                db.loadFindCellIds(idList, onSearchDone, onProgressConsole, hasWildcards);
                }
            }
        ];

        htmls.push("<textarea id='tpIdList' style='height:320px;width:400px;display:block'>");
        htmls.push("</textarea><div id='tpNotFoundIds'></div><div id='tpNotFoundHint'></div>");
        htmls.push("<input id='tpHasWildcard' type='checkbox' style='margin-right: 10px' /> Allow RegEx search, e.g. enter '^TH' to find all IDs that <br>start with 'TH' or '-1$' to find all IDs that end with '-1'");
        var title = "Paste a list of IDs (one per line) to select "+gSampleDesc+"s";
        showDialogBox(htmls, title, {showClose:true, height:dlgHeight, width:dlgWidth, buttons:buttons});
    }

    function onExportIdsClick() {
        /* Edit - Export cell IDs */
        var selCells = renderer.getSelection();

        function buildExportDialog(idList) {
            /* callback when cellIds have arrived */
            var dlgHeight = 500;

            var htmls = [];
            if (selCells.length===0)
                htmls.push("Shown below are the identifiers of all "+idList.length+" cells in the dataset.<p><p>");

            var idListEnc = encodeURIComponent(idList.join("\n"));
            htmls.push("<textarea style='height:320px;width:350px;display:block'>");
            htmls.push(idList.join("\n"));
            htmls.push("</textarea>");

            var buttons =
            [
                {
                    text:"Download as file",
                    click : function() {
                        var blob = new Blob([idList.join("\n")], {type: "text/plain;charset=utf-8"});
                        saveAs(blob, "identifiers.txt");
                    },
                },
                {
                    text:"Copy to clipboard",
                    click :function() {
                        $("textarea").select();
                        document.execCommand('copy');
                        $( this ).dialog( "close" );
                    }
                }
            ];

            let title = "List of "+idList.length+" selected IDs";
            if (selCells.length===0)
                title = "No cells selected";

            showDialogBox(htmls, title,
                {showClose:true, height:dlgHeight, width:500, buttons:buttons}
                );
        }

        db.loadCellIds(selCells, buildExportDialog);
    }


    function onAboutClick() {
        /* user clicked on help > about */
        var dlgHeight = 500;

        var htmls = [];
        var title = "UCSC Cell Browser";

        htmls.push("<p><b>Version:</b> "+gVersion+"</p>");
        htmls.push("<p><b>Written by:</b> Maximilian Haeussler, Nikolay Markov (U Northwestern), Brian Raney, Lucas Seninge</p>");
        htmls.push("<p><b>Testing / User interface / Documentation / Data import / User support:</b> Matt Speir</p>");
        htmls.push("<p><b>Code contributions by:</b> Pablo Moreno (EBI, UK)</p>");
        htmls.push("<p><b>Documentation:</b> <a target=_blank href='https://cellbrowser.readthedocs.io/'>Readthedocs</a></p>");
        htmls.push("<p><b>Github Repo: </b><a target=_blank href='https://github.com/maximilianh/cellBrowser/'>cellBrowser</a></p>");

        showDialogBox(htmls, title, {showClose:true, height:dlgHeight, width:500});
    }

    function buildMenuBar() {
        /* draw the menubar at the top */
       var htmls = [];
       htmls.push("<div style='width:"+menuBarHeight+"px' id='tpMenuBar'>");
       htmls.push('<nav class="navbar navbar-default navbar-xs">');

       htmls.push('<div class="container-fluid">');

         htmls.push('<div class="navbar-header">');
           htmls.push('<a class="navbar-brand" href="#">'+gTitle+'</a>');
         htmls.push('</div>');

         htmls.push('<ul class="nav navbar-nav">');

         htmls.push('<li class="dropdown">');
           htmls.push('<a href="#" class="dropdown-toggle" data-toggle="dropdown" data-submenu role="button" aria-haspopup="true" aria-expanded="false">File</a>');
           htmls.push('<ul class="dropdown-menu">');
             htmls.push('<li><a href="#" id="tpOpenDatasetLink"><span class="dropmenu-item-label">Open dataset...</span><span class="dropmenu-item-content">o</span></a></li>');
             //htmls.push('<li class="dropdown-submenu"><a tabindex="0" href="#">Download Data</a>');
               //htmls.push('<ul class="dropdown-menu" id="tpDownloadMenu">');
                 //htmls.push('<li><a href="#" id="tpDownload_matrix">Gene Expression Matrix</a></li>');
                 //htmls.push('<li><a href="#" id="tpDownload_meta">Cell Metadata</a></li>');
                 //htmls.push('<li><a href="#" id="tpDownload_coords">Visible coordinates</a></li>');
               //htmls.push('</ul>'); // Download sub-menu
             htmls.push('<li><a href="#" id="tpSaveImage">Download current image</a></li>');
             htmls.push('</li>');   // sub-menu container

           htmls.push('</ul>'); // File menu
         htmls.push('</li>'); // File dropdown

         htmls.push('<li class="dropdown">');
         htmls.push('<a href="#" class="dropdown-toggle" data-toggle="dropdown" data-submenu role="button" aria-haspopup="true" aria-expanded="false">Edit</a>');
         htmls.push('<ul class="dropdown-menu">');
         htmls.push('<li><a id="tpSelectAll" href="#"><span class="dropmenu-item-label">Select all visible</span><span class="dropmenu-item-content">s a</span></a></li>');
         htmls.push('<li><a id="tpSelectNone" href="#"><span class="dropmenu-item-label">Select none</span><span class="dropmenu-item-content">s n</span></a></li>');
         htmls.push('<li><a id="tpSelectInvert" href="#"><span class="dropmenu-item-label">Invert selection</span><span class="dropmenu-item-content">s i</span></a></li>');
         htmls.push('<li><a id="tpSelectName" href="#"><span class="dropmenu-item-label">Name selection...</span><span class="dropmenu-item-content">s s</span></a></li>');
         htmls.push('<li><a id="tpExportIds" href="#">Export selected...</a></li>');
         htmls.push('<li><a id="tpSelectComplex" href="#"><span class="dropmenu-item-label">Find cells...</span><span class="dropmenu-item-content">f c</span></a></li>');
         //htmls.push('<li><a id="tpMark" href="#"><span class="dropmenu-item-label">Mark selected</span><span class="dropmenu-item-content">h m</span></a></li>');
         //htmls.push('<li><a id="tpMarkClear" href="#"><span class="dropmenu-item-label">Clear marks</span><span class="dropmenu-item-content">c m</span></a></li>');
         htmls.push('<li><a id="tpSelectById" href="#">Find by ID...<span class="dropmenu-item-content">f i</span></a></li>');
         htmls.push('</ul>'); // View dropdown
         htmls.push('</li>'); // View dropdown

         htmls.push('<li class="dropdown">');
         htmls.push('<a href="#" class="dropdown-toggle" data-toggle="dropdown" data-submenu role="button" aria-haspopup="true" aria-expanded="false">View</a>');
         htmls.push('<ul class="dropdown-menu">');

         htmls.push('<li><a href="#" id="tpZoomPlus"><span class="dropmenu-item-label">Zoom in</span><span class="dropmenu-item-content">+</span></a></li>');
         htmls.push('<li><a href="#" id="tpZoomMinus"><span class="dropmenu-item-label">Zoom out</span><span class="dropmenu-item-content">-</span></a></li>');
         htmls.push('<li><a href="#" id="tpZoom100Menu"><span class="dropmenu-item-label">Zoom 100%</span><span class="dropmenu-item-content">space</span></a></li>');
         htmls.push('<li><a href="#" id="tpSplitMenu"><span id="tpSplitMenuEntry" class="dropmenu-item-label">Split screen</span><span class="dropmenu-item-content">t</span></a></li>');
         htmls.push('<li><a href="#" id="tpHeatMenu"><span id="tpHeatMenuEntry" class="dropmenu-item-label">Toggle Heatmap</span><span class="dropmenu-item-content">h</span></a></li>');

         htmls.push('<li><hr class="half-rule"></li>');

         //htmls.push('<li><a href="#" id="tpOnlySelectedButton">Show only selected</a></li>');
         //htmls.push('<li><a href="#" id="tpFilterButton">Hide selected '+gSampleDesc+'s</a></li>');
         //htmls.push('<li><a href="#" id="tpShowAllButton">Show all '+gSampleDesc+'</a></li>');
         htmls.push('<li><a href="#" id="tpHideShowLabels"><span id="tpHideMenuEntry">Hide labels</span><span class="dropmenu-item-content">c l</span></a></li>');
         //htmls.push('<li><hr class="half-rule"></li>');

         //htmls.push('<li class="dropdown-submenu"><a tabindex="0" href="#">Transparency</a>');
           //htmls.push('<ul class="dropdown-menu" id="tpTransMenu">');
             //htmls.push('<li id="tpTrans0"><a href="#">0%</a></li>');
             //htmls.push('<li id="tpTrans40"><a href="#">40%</a></li>');
             //htmls.push('<li id="tpTrans60"><a href="#">60%</a></li>');
             //htmls.push('<li id="tpTrans80"><a href="#">80%</a></li>');
           //htmls.push('</ul>'); // Transparency sub-menu
         //htmls.push('</li>');   // sub-menu container

         //htmls.push('<li class="dropdown-submenu"><a tabindex="0" href="#">Circle size</a>');
           //htmls.push('<ul class="dropdown-menu" id="tpSizeMenu">');
             //htmls.push('<li id="tpSize1"><a href="#">1 px</a></li>');
             //htmls.push('<li id="tpSize2"><a href="#">2 px</a></li>');
             //htmls.push('<li id="tpSize3"><a href="#">3 px</a></li>');
             //htmls.push('<li id="tpSize4"><a href="#">4 px</a></li>');
             //htmls.push('<li id="tpSize5"><a href="#">5 px</a></li>');
             //htmls.push('<li id="tpSize6"><a href="#">6 px</a></li>');
             //htmls.push('<li id="tpSize7"><a href="#">7 px</a></li>');
             //htmls.push('<li id="tpSize8"><a href="#">8 px</a></li>');
           //htmls.push('</ul>'); // Circle size sub-menu
         //htmls.push('</li>');   // sub-menu container

         htmls.push('</ul>'); // View dropdown-menu
         htmls.push('</li>'); // View dropdown container

         htmls.push('<li class="dropdown">');
         htmls.push('<a href="#" class="dropdown-toggle" data-toggle="dropdown" data-submenu role="button" aria-haspopup="true" aria-expanded="false">Tools</a>');
         htmls.push('<ul class="dropdown-menu">');
         //htmls.push('<li><a href="#" id="tpRenameClusters">Rename clusters...<span class="dropmenu-item-content"></span></a></li>');
         htmls.push('<li><a href="#" id="tpCustomAnnots">Remove all custom annotations<span class="dropmenu-item-content"></span></a></li>');
         //htmls.push('<li><a href="#" id="tpCluster">Run clustering...<span class="dropmenu-item-content"></span></a></li>');
         htmls.push('<li class="disabled"><a href="#" id="tpSetBackground">Set as background cells<span class="dropmenu-item-content">b s</span></a></li>');
         htmls.push('<li class="disabled"><a href="#" id="tpResetBackground">Reset background cells<span class="dropmenu-item-content">b r</span></a></li>');
         htmls.push('</ul>'); // Tools dropdown-menu
         htmls.push('</li>'); // Tools dropdown container


         htmls.push('<li class="dropdown">');
         htmls.push('<a href="#" class="dropdown-toggle" data-toggle="dropdown" data-submenu role="button" aria-haspopup="true" aria-expanded="false">Help</a>');
         htmls.push('<ul class="dropdown-menu">');
         htmls.push('<li><a href="#" id="tpAboutButton">About</a></li>');
         htmls.push('<li><a href="#" id="tpTutorialButton">Tutorial</a></li>');
         htmls.push('<li><a target=_blank href="https://github.com/maximilianh/cellBrowser#readme" id="tpGithubButton">Setup your own cell browser</a></li>');
         htmls.push('</ul>'); // Help dropdown-menu
         htmls.push('</li>'); // Help dropdown container

       htmls.push('</ul>'); // navbar-nav

       htmls.push('</div>'); // container
       htmls.push('</nav>'); // navbar
       htmls.push('</div>'); // tpMenuBar

       $(document.body).append(htmls.join(""));

       $('#tpTransMenu li a').click( onTransClick );
       $('#tpSizeMenu li a').click( onSizeClick );
       //$('#tpFilterButton').click( onHideSelectedClick );
       //$('#tpOnlySelectedButton').click( onShowOnlySelectedClick );
       $('#tpZoom100Menu').click( onZoom100Click );
       $('#tpSplitMenu').click( onSplitClick );
       $('#tpHeatMenu').click( onHeatClick );
       $('#tpZoomPlus').click( onZoomInClick );
       $('#tpZoomMinus').click( onZoomOutClick );
       //$('#tpShowAllButton').click( onShowAllClick );
       $('#tpHideShowLabels').click( onHideShowLabelsClick );
       $('#tpExportIds').click( onExportIdsClick );
       $('#tpSelectById').click( onSelectByIdClick );
       $('#tpMark').click( onMarkClick );
       $('#tpMarkClear').click( onMarkClearClick );
       $('#tpTutorialButton').click( function()  { showIntro(false); } );
       $('#tpAboutButton').click( onAboutClick );
       $('#tpOpenDatasetLink').click( openCurrentDataset );
       $('#tpSaveImage').click( onSaveAsClick );
       $('#tpSelectAll').click( onSelectAllClick );
       $('#tpSelectNone').click( onSelectNoneClick );
       $('#tpSelectInvert').click( onSelectInvertClick );
       $('#tpSelectName').click( onSelectNameClick );
       $('#tpSelectComplex').click( onFindCellsClick );


       $('#tpRenameClusters').click( onRenameClustersClick );
       $('#tpCustomAnnots').click( onCustomAnnotationsClick );
       $('#tpSetBackground').click( onBackgroudSetClick );
       $('#tpResetBackground').click( onBackgroudResetClick );
       //$('#tpCluster').click( onRunClusteringClick );

       // This version is more like OSX/Windows:
       // - menus only open when you click on them
       // - once you have clicked, they start to open on hover
       // - a click anywhere else will stop the hovering
       var doHover = false;
       $(".nav > .dropdown").click( function(){ doHover = !doHover; return true;} );
       $(".nav > .dropdown").hover(
           function(event) {
               if (doHover) {
                   $(".dropdown-submenu").removeClass("open"); $(".dropdown").removeClass("open"); $(this).addClass('open');
               }
           });

       $(document).click ( function() { doHover= false; });

       // when user releases the mouse outside the canvas, remove the zooming marquee
       $(document).mouseup ( function(ev) { if (ev.target.nodeName!=="canvas") { renderer.resetMarquee(); }} );

       $('[data-submenu]').submenupicker();

    }

    function resizeDivs(skipRenderer) {
       /* resize all divs and the renderer to current window size */
       var rendererLeft = metaBarWidth+metaBarMargin;
       var rendererHeight  = window.innerHeight - menuBarHeight - toolBarHeight;

       var rendererWidth = window.innerWidth - legendBarWidth - rendererLeft;
       var legendBarLeft = rendererWidth+metaBarMargin+metaBarWidth;

       var heatWidth, heatHeight;
       if (db && db.heatmap) {
            heatWidth = rendererWidth;
            heatHeight = db.heatmap.height;
            rendererHeight = rendererHeight - heatHeight;
            db.heatmap.setSize(heatWidth, heatHeight);
            let heatTop = window.innerHeight - heatHeight;
            db.heatmap.div.style.top = heatTop+"px";
            db.heatmap.draw();
       }

       $("#tpToolBar").css("width", rendererWidth+"px");

       $("#tpToolBar").css("height", toolBarHeight+"px");
       $("#tpLeftSidebar").css("height", (window.innerHeight - menuBarHeight)+"px");

       // when this is run the first time, these elements don't exist yet.
       // Note that the whole concept of forcing these DIVs to go up to the screen size is strange, but
       // I have not found a way in CSS to make them go to the end of the screen. They need to have a fixed size,
       // as otherwise the scroll bars of tpLegendBar and tpMetaPanel won't appear
       if ($('#tpMetaPanel').length!==0)
           $("#tpMetaPanel").css("height", (window.innerHeight - $('#tpMetaPanel').offset().top)+"px");
       if ($('#tpLegendRows').length!==0)
           $("#tpLegendBar").css("height", (window.innerHeight - $('#tpLegendBar').offset().top)+"px");
       $('#tpLegendBar').css('left', legendBarLeft+"px");

       if (skipRenderer!==true)
           renderer.setSize(rendererWidth, rendererHeight, true);

    }

    var progressUrls = {};

    function onProgressConsole(ev) {
        //console.log(ev);
    }

    function onProgress(ev) {
        /* update progress bars. The DOM elements of these were added in maxPlot (not optimal?)  */
        var url = ev.currentTarget.responseURL;
        url = url.split("?")[0]; // strip off the md5 checksum

        if (url.search("exprMatrix.bin")!==-1) // never show progress bar for single gene vector requests
            return;

        var progressRowIdx = progressUrls[url]; // there can be multiple progress bars
        if (progressRowIdx===undefined) {
            // if there is none yet, find the first free index
            progressRowIdx = 0;
            for (var oldUrl in progressUrls) {
                progressRowIdx = Math.max(progressRowIdx, progressUrls[oldUrl]);

            }
            progressRowIdx++;
            progressUrls[url] = progressRowIdx;
        }

        var label = url;
        if (url.endsWith("coords.bin"))
            label = "Loading Coordinates";
        else if (url.endsWith(".bin"))
            label = "Loading cell annotations";

        var labelId = "#mpProgressLabel"+progressRowIdx;
        $(labelId).html(label);

        var percent = Math.round(100 * (ev.loaded / ev.total));

        if (percent>=99) {
            $("#mpProgress"+progressRowIdx).css("width", percent+"%");
            $("#mpProgress"+progressRowIdx).show(0);
            //progressUrls.splice(index, 1);
            delete progressUrls[url];
            $("#mpProgressDiv"+progressRowIdx).css("display", "none");
        }
        else {
            $("#mpProgress"+progressRowIdx).css("width", percent+"%");
            $("#mpProgressDiv"+progressRowIdx).css("display", "inherit");
        }
    }


    function colorByMetaField(fieldName, doneLoad) {
        /* load the meta data for a field, setup the colors, send it all to the renderer and call doneLoad */

        function onMetaArrLoaded(metaArr, metaInfo) {
            gLegend = buildLegendForMeta(metaInfo);
            buildLegendBar();
            var renderColors = legendGetColors(gLegend.rows);
            renderer.setColors(renderColors);
            renderer.setColorArr(metaArr);
            metaInfo.arr = metaArr;
            doneLoad();
        }

       if (doneLoad===undefined)
           doneLoad = function() { renderer.drawDots(); };

       if (fieldName===null || fieldName===undefined) {
           // obscure hacky option: you can set the default color field to "None"
           // so there is no coloring at all on startup
           renderer.setColors(["black"]);
           var cellCount = db.conf.sampleCount;
           renderer.setColorArr(new Uint8Array(cellCount));
           gLegend.rows = [];
           gLegend.title = "Nothing selected";
           gLegend.rows.push( {
                   color:"000000", defColor:null, label:"No Value",
                   count:cellCount, intKey:0, strKey:null
           } );
           doneLoad();
           return;
       }

       var metaInfo  = db.findMetaInfo(fieldName);
       console.log("Color by meta field "+fieldName);

       // cbData always keeps the most recent expression array. Reset it now.
       if (db.lastExprArr)
           delete db.lastExprArr;

       var defaultMetaField = db.getDefaultColorField()[1];

       // internal field names cannot contain non-alpha chars, so tolerate user errors here
       // otherwise throw an error
       if (metaInfo === null && fieldName!==undefined) {
           metaInfo = db.findMetaInfo(fieldName.replace(/[^0-9a-z]/gi, ''));
           if (metaInfo === null) {
               alert("The field "+fieldName+" does not exist in the sample/cell annotations. Cannot color on it.");
               metaInfo = db.findMetaInfo(defaultMetaField);
           }
       }

       if (metaInfo.type==="uniqueString") {
           warn("This field contains a unique identifier. You cannot color on such a field. However, you can search for values in this field using 'Edit > Find by ID'.");
           return null;
       }

       if (metaInfo.diffValCount > MAXCOLORCOUNT && metaInfo.type==="enum") {
           warn("This field has "+metaInfo.diffValCount+" different values. Coloring on a field that has more than "+MAXCOLORCOUNT+" different values is not supported.");
           return null;
       }


        if (fieldName===defaultMetaField)
           changeUrl({"meta":null, "gene":null});
        else
           changeUrl({"meta":fieldName, "gene":null});

        if (metaInfo.arr) // eg custom fields
            onMetaArrLoaded(metaInfo.arr, metaInfo);
        else
            db.loadMetaVec(metaInfo, onMetaArrLoaded, onProgress);


       changeUrl({"pal":null});
       // clear the gene search box
       var select = $('#tpGeneCombo')[0].selectize.clear();
    }

    function activateTab(name) {
        /* activate a tab on the left side */
        var idx = 0;
        if (name==="gene")
            idx = 1;

        $( "#tpLeftTabs" ).tabs( "option", "active", idx );
    }

    function doLog2(arr) {
        /* take log2(x+1) for all values in array and return the result */
        var arr2 = new Float64Array(arr.length);
        for (var i = 0; i < arr.length; i++) {
            arr2.push(Math.log2(arr[i]+1));
        }
        return arr2;
    }

    function splitExprByMeta(exprVec, splitArr, selCells) {
        /* split the expression vector into two vectors. splitArr is an array with 0/1, indicates where values go.
         * if selCells is not null, restrict the splitting to just indices in selCells.
         * Returns array of the two arrays.
         * */
        console.time("findCellsWithMeta");
        if (exprVec.length!==splitArr.length) {
            warn("internal error - splitExprByMeta: exprVec has diff length from splitArr");
        }

        var arr1 = [];
        var arr2 = [];

        // code duplication, not very elegant, but avoids creating an array just for the indices
        if (selCells.length===0)
            // the version if no cells are selected
            for (var cellIdx = 0; cellIdx < exprVec.length; cellIdx++) {
                var val = exprVec[cellIdx];
                if (splitArr[cellIdx]===0)
                    arr1.push(val);
                else
                    arr2.push(val);
            }
        else
            // the version with a cell selection
            for (var i = 0; i < selCells.length; i++) {
                let cellIdx = selCells[i];
                let val = exprVec[cellIdx];
                if (splitArr[cellIdx]===0)
                    arr1.push(val);
                else
                    arr2.push(val);
            }

        if (db.conf.violinDoLog2) {
            console.time("log2");
            arr1 = doLog2(arr1);
            arr2 = doLog2(arr2);
            console.timeEnd("log2");
        }

        console.timeEnd("findCellsWithMeta");
        return [arr1, arr2];
    }

    function splitExpr(exprVec, selCells) {
        /* split the expression vector into two vectors, one for selected and one for unselected cells */
        console.time("splitExpr");
        var selMap = {};
        for (var i = 0; i < selCells.length; i++) {
            selMap[selCells[i]] = null;
        }

        var sel = [];
        var unsel = [];
        for (i = 0; i < exprVec.length; i++) {
            if (i in selMap)
                sel.push(exprVec[i]);
            else
                unsel.push(exprVec[i]);
        }

        console.timeEnd("splitExpr");
        return [sel, unsel];
    }

    function buildViolinFromValues(labelList, dataList) {
        /* make a violin plot given the labels and the values for them */
        if ("violinChart" in window)
            window.violinChart.destroy();

        var labelLines = [];
        labelLines[0] = labelList[0].split("\n");
        labelLines[0].push(dataList[0].length);
        if (dataList.length > 1) {
            labelLines[1] = labelList[1].split("\n");
            labelLines[1].push(dataList[1].length);
        }

        const ctx = getById("tpViolinCanvas").getContext("2d");

 	var violinData = {
          labels : labelLines,
	  datasets: [{
            data : dataList,
	    label: 'Mean',
	    backgroundColor: 'rgba(255,0,0,0.5)',
	    borderColor: 'red',
	    borderWidth: 1,
	    outlierColor: '#999999',
	    padding: 7,
	    itemRadius: 0
	}]
	};

        var optDict = {
            maintainAspectRatio: false,
            legend: { display: false },
	    title: { display: false }
        };

        var yLabel = null;
        if (db.conf.unit===undefined && db.conf.matrixArrType==="Uint32")
            yLabel = "read/UMI count";
        if (db.conf.unit!==undefined)
            yLabel = db.conf.unit;

        if (yLabel!==null)
            optDict.scales = {
                yAxes: [{
                    scaleLabel: {
                        display: true,
                        labelString: yLabel
                    }
                    }]
            };

        window.setTimeout(function() {
            console.time("violinDraw");
            window.violinChart = new Chart(ctx, {
                type: 'violin',
                data: violinData,
                options: optDict
            });
            console.timeEnd("violinDraw");
        }, 10);
    }


    function buildViolinFromMeta(exprVec, metaName, selCells) {
        /* load a binary meta vector, split the exprVector by it and make two violin plots, one meta value vs the other.  */
        var metaInfo = db.findMetaInfo(metaName);
        if (metaInfo.valCounts.length!==2) {
            alert("Config error: meta field in 'violinField', '"+db.conf.violinField+"' does not have two distinct values.");
            return;
        }

        var labelList = [metaInfo.valCounts[0][0], metaInfo.valCounts[1][0]];
        db.loadMetaVec( metaInfo,
                function(metaArr) {
                    var dataList = splitExprByMeta(exprVec, metaArr, selCells);
                    buildViolinFromValues(labelList, dataList);
                },
                null);
    }

    //function removeViolinPlot() {
        /* destroy the violin plot */
        //if ("violinChart" in window)
            //window.violinChart.destroy();
        //$('#tpViolinCanvas').remove();
    //}

    function buildViolinPlot() {
        /* create the violin plot, depending on the current selection and the violinField config */
        var exprVec = gLegend.exprVec;
        if (exprVec===undefined)
            return;

        var dataList = [];
        var labelList = [];
        var selCells = renderer.getSelection();

        // filter exprVec by background
        if (background !== null) {
            var ourSelCells = {};
            for (var i = 0; i < selCells.length; i++) {
                ourSelCells[selCells[i]] = true;
            }
            var ourCells = {};
            for (i = 0; i < background.length; i++) {
                ourCells[background[i]] = true;
            }

            var result = [];
            var renamedSelCells = [];
            for (i = 0; i < exprVec.length; i++) {
                if (i in ourSelCells) {
                    renamedSelCells.push(result.length);
                    result.push(exprVec[i]);
                } else if (i in ourCells) {
                    result.push(exprVec[i]);
                }
            }
            exprVec = result;
            selCells = renamedSelCells;
        }
        // if we have a violin meta field to split on, make two violin plots, one meta vs, the other meta
        // restrict the plot to the selected cells, if any
        if (db.conf.violinField!==undefined) {
            buildViolinFromMeta(exprVec, db.conf.violinField, selCells);
        } else {
            // there is no violin field
            if (selCells.length===0) {
                // no selection, no violinField: default to a single violin plot
                dataList = [Array.prototype.slice.call(exprVec)];
                if (background === null) {
                    labelList = ['All cells'];
                } else {
                    labelList = ['Background\ncells'];
                }
                buildViolinFromValues(labelList, dataList);
            } else {
                // cells are selected and no violin field: make two violin plots, selected against other cells
                dataList = splitExpr(exprVec, selCells);
                if (background === null) {
                    labelList = ['Selected', 'Others'];
                } else {
                    labelList = ['Selected', 'Background'];
                }
                if (dataList[1].length===0) {
                    dataList = [dataList[0]];
                    labelList = ['All Selected'];
                }
                buildViolinFromValues(labelList, dataList);
            }
        }
    }

    function selectizeSetValue(selector, name) {
        /* little convenience method to set a selective dropdown to a given
         * value. does not trigger the change event. */
        if (name===undefined)
            return;
        var sel = $(selector)[0].selectize;
        sel.addOption({id: name, text: name});
        sel.setValue(name, 1); // 1 = do not fire change
    }

    function colorByLocus(locusStr, onDone) {
        /* color by a gene or peak, load the array into the renderer and call onDone or just redraw 
         * peak can be in format: +chr1:1-1000
         * gene can be in format: geneSym or geneSym=geneId
         * */
        if (onDone===undefined)
            onDone = function() { renderer.drawDots(); };

        function gotGeneVec(exprArr, decArr, locusStr, geneDesc, binInfo) {
            /* called when the expression vector has been loaded and binning is done */
            if (decArr===null)
                return;
            console.log("Received expression vector, for "+locusStr+", desc: "+geneDesc);
            // update the URL and possibly the gene combo box
            var fullLocusStr = locusStr;
            if (locusStr.indexOf("|") > -1) {
                fullLocusStr = peakListSerialize(); // the locus str can be +chr1:1-1000 to add a single gene, but we want all
                changeUrl({"locus":locusStr, "meta":null});
            } else {
                changeUrl({"gene":locusStr, "meta":null});
            }

            makeLegendExpr(fullLocusStr, geneDesc, binInfo, exprArr, decArr);
            renderer.setColors(legendGetColors(gLegend.rows));
            renderer.setColorArr(decArr);
            buildLegendBar();
            onDone();

            // update the "recent genes" div
            for (var i = 0; i < gRecentGenes.length; i++) {
                // remove previous gene entry with the same symbol
                if (gRecentGenes[i][0]===locusStr) {
                    gRecentGenes.splice(i, 1);
                    break;
                }
            }
            gRecentGenes.unshift([locusStr, geneDesc]); // insert at position 0
            gRecentGenes = gRecentGenes.slice(0, 9); // keep only nine last
            buildGeneTable(null, "tpRecentGenes", null, null, gRecentGenes);
            $('#tpRecentGenes .tpGeneBarCell').click( onGeneClick );
        }

        // clear the meta combo
        $('#tpMetaCombo').val(0).trigger('chosen:updated');

        console.log("Loading gene expression vector for "+locusStr);
        db.loadExprAndDiscretize(locusStr, gotGeneVec, onProgress);

    }

    function gotCoords(coords, info, clusterInfo, newRadius) {
        /* called when the coordinates have been loaded */
        if (coords.length===0)
            alert("cellBrowser.js/gotCoords: coords.bin seems to be empty");
        var opts = {};
        if (newRadius)
            opts["radius"] = newRadius;

        // label text can be overriden by the user cart
        var labelField = db.conf.labelField;

        if (clusterInfo) {
            var origLabels = [];
            var clusterMids = clusterInfo.labels;
            // old-style files contain just coordinates, no order
            if (clusterMids === undefined) {
                clusterMids = clusterInfo;
            }

            gLabelCoordCache[labelField] = clusterMids;
            for (var i = 0; i < clusterMids.length; i++) {
                origLabels.push(clusterMids[i][2]);
            }
            renderer.origLabels = origLabels;
            clusterMids = [];
         }

        opts["lines"] = clusterInfo.lines;
        opts["lineWidth"] = db.conf.lineWidth;

        renderer.setCoords(coords, clusterMids, info.minX, info.maxX, info.minY, info.maxY, opts);
    }

    function computeAndSetLabels(values, metaInfo) {
        var labelCoords;
        if (gLabelCoordCache[metaInfo.label] !== undefined) {
            labelCoords = gLabelCoordCache[metaInfo.label];
        } else {
            var coords = renderer.coords.orig;
            var names = null;
            if (metaInfo.type !== "float" && metaInfo.type !== "int") {
                var names = metaInfo.ui.shortLabels;
            }

            // console.log(metaInfo);

            console.time("cluster centers");
            var calc = {};
            for (var i = 0, I = values.length; i < I; i++) {
                if (names) {
                    var label = names[values[i]];
                } else {
                    var label = metaInfo.origVals[i].toFixed(2);
                }
                if (calc[label] === undefined) {
                    calc[label] = [[], [], 0]; // all X, all Y, count
                }
                calc[label][0].push(coords[i * 2]);
                calc[label][1].push(coords[i * 2 + 1]);
                calc[label][2] += 1;
            }
            labelCoords = [];
            for (label in calc) {
                var midX = selectMedian(calc[label][0]);
                var midY = selectMedian(calc[label][1]);
                labelCoords.push([midX, midY, label]);
            }
            console.timeEnd("cluster centers");
            gLabelCoordCache[metaInfo.label] = labelCoords;
        }
        var shouldRedraw = renderer.setLabelCoords(labelCoords);
        if (shouldRedraw) {
            renderer.drawDots();
        } else {
            renderer.redrawLabels();
        }
    }

    function setLabelField(labelField) {
        var metaInfo = db.findMetaInfo(labelField);
        db.conf.activeLabelField = labelField;

        if (metaInfo.arr) // preloaded
            computeAndSetLabels(metaInfo.arr, metaInfo);
        else
            db.loadMetaVec(metaInfo, computeAndSetLabels);
    }

   function colorByDefaultField(onDone) {
       /* get the default coloring field from the config or the URL and start coloring by it.
        * Call onDone() when done. */
       var colorByInfo = db.getDefaultColorField();
       var colorType = colorByInfo[0];
       var colorBy = colorByInfo[1];

       // allow to override coloring by URL args
       if (getVar("gene")!==undefined) {
           colorType = "gene";
           colorBy = getVar("gene");
           activateTab("gene");
       }
       else if (getVar("meta")!==undefined) {
           colorType = "meta";
           colorBy = getVar("meta");
           activateTab("meta");
       } else if (getVar("locus")!==undefined) {
           colorType = "locus";
           colorBy = getVar("locus");
           activateTab("gene");
       }

       gLegend = {};
       if (colorType==="meta") {
           colorByMetaField(colorBy, onDone);
           // update the meta field combo box
           var fieldIdx  = db.fieldNameToIndex(colorBy);
           if (fieldIdx===null) {
               alert("Default coloring is configured to be on field "+colorBy+
                    " but cannot find a field with this name, using field 1 instead.");
               fieldIdx = 1;
           }

           $('#tpMetaCombo').val(fieldIdx).trigger('chosen:updated');
           $('#tpMetaBox_'+fieldIdx).addClass('tpMetaSelect');
       }
       else if (colorType==="locus") {
           colorByLocus(colorBy, onDone);
           peakListSetStatus(colorBy);
       } else
           // must be gene then
           colorByLocus(colorBy, onDone);
    }

   function makeFullLabel(db) {
       /* return full name of current dataset, including parent names */
       var nameParts = [];
       var parents = db.conf.parents;
       if (parents)
           for (var i=0; i < parents.length; i++)
               if (parents[i][0]!="") // "" is the root dataset = no need to add
                   nameParts.push( parents[i][1] );

       nameParts.push( db.conf.shortLabel );
       var datasetLabel = nameParts.join(" - ");
       return datasetLabel;
   }

    function renderData() {
    /* init the basic UI parts, the main renderer in the center, start loading and draw data when ready
     */
       var forcePalName = getVar("pal", null);

       var loadsDone = 0;

       var selList = null; // search expression to select, in format accepted by findCellsMatchingQueryList()

       function doneOnePart() {
       /* make sure renderer only draws when both coords and other data have loaded */
           loadsDone +=1;
           if (loadsDone===2) {
               buildLegendBar();
               setLabelField(db.conf.labelField);

               if (forcePalName!==null) {
                   legendChangePaletteAndRebuild(forcePalName);
               }
               else
                   renderer.setColors(legendGetColors(gLegend.rows));


               renderer.setTitle("Dataset: "+makeFullLabel(db));

               if (selList)
                   findCellsMatchingQueryList(selList, function (cellIds) {
                        renderer.selectSet(cellIds);
                        renderer.drawDots();
                   });
               else
                   renderer.drawDots();
           }
       }

       function guessRadiusAlpha(dotCount) {
           /* return reasonable radius and alpha values for a number of dots */
           if (dotCount<3000)
               return [4, 0.7];
           if (dotCount<6000)
               return [4, 0.6];
           if (dotCount<10000)
               return [3, 0.5];
           if (dotCount<35000)
               return [2, 0.3];
           if (dotCount<60000)
               return [1, 0.5];
           // everything else
           return [0, 0.3];
       }

       function makeRendConf(dbConf, dotCount) {
           /* return the 'args' object for the renderer, based on dbConf and count of dots */
           var rendConf = {};

           var radius = dbConf.radius;
           var alpha = dbConf.alpha;

           if (radius===undefined || alpha===undefined) {
               var radiusAlpha = guessRadiusAlpha(dotCount);

               if (radius===undefined)
                   radius = radiusAlpha[0];
               if (alpha===undefined)
                   alpha = radiusAlpha[1];
            }

           rendConf["radius"] = radius;
           rendConf["alpha"] = alpha;

           rendConf["mode"] = "move"; // default mode, one of "move", "select" or "zoom"

           return rendConf;
        }

       function gotFirstCoords(coords, info, clusterMids) {
           /* XX very ugly way to implement promises. Need a better approach one day. */
           gotCoords(coords, info, clusterMids);
           doneOnePart();
       }

       var rendConf = makeRendConf(db.conf, db.conf.sampleCount);
       renderer.initPlot(rendConf);

       if (db.conf.showLabels===false)
           renderer.setShowLabels(false);

       buildLeftSidebar();
       buildToolBar(db.conf.coords, db.conf, metaBarWidth+metaBarMargin, menuBarHeight);

       db.loadCoords(0, gotFirstCoords, onProgress);

       if (getVar("select")!==undefined) {
           selList = JSURL.parse(getVar("select"));
       }

       colorByDefaultField(doneOnePart);

       if (db.conf.atacSearch)
           db.loadGeneLocs(db.conf.atacSearch, db.conf.fileVersions.geneLocs);
       // pre-load the dataset description file, as the users will often go directly to the info dialog
       // and the following pre-loads risk blocking this load.
       var jsonUrl = cbUtil.joinPaths([db.conf.name, "desc.json"]) +"?"+db.conf.md5;
       fetch(jsonUrl);

       //if (db.conf.sampleCount < 50000) {
           if (db.conf.quickGenes)
               db.preloadGenes(db.conf.quickGenes, function() {
                   updateGeneTableColors(null);
                   if (getVar("heat")==="1")
                       onHeatClick();
                }, onProgressConsole);
           db.preloadAllMeta();
        //}
    }

    function onTransClick(ev) {
    /* user has clicked transparency menu entry */
        var transText = ev.target.innerText;
        var transStr = transText.slice(0, -1); // remove last char
        var transFloat = 1.0 - (parseFloat(transStr)/100.0);
        transparency = transFloat;
        plotDots();
        $("#tpTransMenu").children().removeClass("active");
        $("#tpTrans"+transStr).addClass("active");
        renderer.render(stage);
    }

    function legendSort(sortBy) {
        /* sort the legend by "name" or "count" */
        var rows = gLegend.rows;

        if (sortBy==="name") {
            // index 2 is the label
            rows.sort(function(a, b) { return naturalSort(a.label, b.label); });
        }
        else {
            // sort this list by count = index 3
            rows.sort(function(a, b) { return b.count - a.count; }); // reverse-sort by count
        }
        buildLegendBar();
    }

    //function filterCoordsAndUpdate(cellIds, mode) {
    /* hide/show currently selected cell IDs or "show all". Rebuild the legend and the coloring. */
        //if (mode=="hide")
            //shownCoords = removeCellIds(shownCoords, cellIds);
        //else if (mode=="showOnly")
            //shownCoords = showOnlyCellIds(shownCoords, cellIds);
        //else
            //shownCoords = allCoords.slice();

        //pixelCoords = scaleData(shownCoords);

        //makeLegendObject();
        //buildLegendBar();
        //gSelCellIds = {};
        //plotDots();
        //renderer.render(stage);
        //menuBarShow("#tpShowAllButton");
    //}

    /* function onHideSelectedClick(ev) {
     user clicked the hide selected button
        filterCoordsAndUpdate(gSelCellIds, "hide");
        menuBarHide("#tpFilterButton");
        menuBarHide("#tpOnlySelectedButton");
        menuBarShow("#tpShowAllButton");
        ev.preventDefault();
    }


    function onShowOnlySelectedClick(ev) {
     // user clicked the only selected button
        filterCoordsAndUpdate(gSelCellIds, "showOnly");
        menuBarHide("#tpFilterButton");
        menuBarHide("#tpOnlySelectedButton");
        menuBarShow("#tpShowAllButton");
        ev.preventDefault();
    }

    function onShowAllClick(ev) {
    // user clicked the show all menu entry
        //gSelCellIds = {};
        filterCoordsAndUpdate(gSelCellIds, "showAll");
        shownCoords = allCoords.slice(); // complete copy of list, fastest in Blink
        pixelCoords = scaleData(shownCoords);
        makeLegendObject();
        buildLegendBar();
        gClasses = assignCellClasses();
        plotDots();
        menuBarHide("#tpFilterButton");
        menuBarHide("#tpOnlySelectedButton");
        menuBarHide("#tpShowAllButton");
        gLegend.lastClicked = null;
        renderer.render(stage);
    } */

    function onHideShowLabelsClick(ev) {
    /* user clicked the hide labels / show labels menu entry */
        if ($("#tpHideMenuEntry").text()===SHOWLABELSNAME) {
            renderer.setShowLabels(true);
            $("#tpHideMenuEntry").text(HIDELABELSNAME);
        }
        else {
            renderer.setShowLabels(false);
            $("#tpHideMenuEntry").text(SHOWLABELSNAME);
        }

        renderer.drawDots();
    }

    function onSizeClick(ev) {
    /* user clicked circle size menu entry */
        var sizeText = ev.target.innerText;
        var sizeStr = sizeText.slice(0, 1); // keep only first char
        circleSize = parseInt(sizeStr);
        $("#tpSizeMenu").children().removeClass("active");
        $("#tpSize"+circleSize).addClass("active");
        plotDots();
        renderer.render(stage);
    }

    function onZoom100Click(ev) {
    /* in addition to zooming (done by maxPlot already), reset the URL */
        changeUrl({'zoom':null});
        renderer.zoom100();
        renderer.drawDots();
        $("#tpZoom100Button").blur(); // remove focus
        ev.preventDefault();
        return false;
    }

    function activateMode(modeName) {
        renderer.activateMode(modeName);
    }

    function onZoomOutClick(ev) {
        var zoomRange = renderer.zoomBy(0.8);
        pushZoomState(zoomRange);
        renderer.drawDots();
        ev.preventDefault();
    }

    function onZoomInClick(ev) {
        var zoomRange = renderer.zoomBy(1.2);
        pushZoomState(zoomRange);
        renderer.drawDots();
        ev.preventDefault();
    }

    function onWindowResize(ev) {
        /* called when window is resized by user */
        if (ev.target.id==="tpHeat") // hack: do not do anything if jquery resizable() started this.
            return;
        resizeDivs();
    }

    function onColorPaletteClick(ev) {
        /* called when users clicks a color palette */
        var palName = ev.target.getAttribute("data-palette");
        if (palName==="default")
            legendSetColors(gLegend, null) // reset the colors
        legendChangePaletteAndRebuild(palName);
        renderer.drawDots();
    }

    function buildEmptyLegendBar(fromLeft, fromTop) {
        // create an empty right side legend bar
        var htmls = [];
        htmls.push("<div id='tpLegendBar' style='position:absolute;top:"+fromTop+"px;left:"+fromLeft+"px; width:"+legendBarWidth+"px'>");
        htmls.push("<div class='tpSidebarHeader'><div style='float:left'>Legend</div>");

        //htmls.push("<div id='tpToolbarButtons' style='padding-bottom: 2px'>");
        htmls.push("<div style='float:right' class='btn-group btn-group-xs'>");
            htmls.push("<button type='button' class='btn btn-default dropdown-toggle' data-toggle='dropdown' aria-haspopup='true' aria-expanded='false' id='tpChangeColorScheme'>Colors&nbsp;<span class='caret'> </span></button>");
            htmls.push('<ul class="dropdown-menu pull-right">');
            htmls.push('<li><a class="tpColorLink" data-palette="default" href="#">Reset to Default</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="rainbow" href="#">Qualitative: Rainbow</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="tol-dv" href="#">Qualitative: Paul Tol&#39;s</a></li>');
            //htmls.push('<li><a class="tpColorLink" data-palette="cb-Paired" href="#">Qualitative: paired</a></li>');
            //htmls.push('<li><a class="tpColorLink" data-palette="cb-Set3" href="#">Qualitative: pastel</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="blues" href="#">Gradient: shades of blue</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="reds" href="#">Gradient: shades of red</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="tol-sq-blue" href="#">Gradient: beige to red</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="tol-rainbow" href="#">Gradient: blue to red</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="viridis" href="#">Gradient: Viridis</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="magma" href="#">Gradient: Magma</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="inferno" href="#">Gradient: Inferno</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="plasma" href="#">Gradient: Plasma</a></li>');
            htmls.push('</ul>');
        htmls.push("</div>"); // btn-group
        //htmls.push("</div>"); // tpToolbarButtons

        htmls.push("</div>"); // tpSidebarHeader

        //htmls.push("<div id='tpLegendTitleBox' style='position:relative; width:100%; height:1.5em'>");
        htmls.push("<div id='tpLegendContent'>");
        htmls.push("</div>"); // content
        htmls.push("</div>"); // bar
        $(document.body).append(htmls.join(""));

        $(".tpColorLink").click( onColorPaletteClick );
    }

    function getTextWidth(text, font) {
        // re-use canvas object for better performance
        // http://stackoverflow.com/questions/118241/calculate-text-width-with-javascript
        var canvas = getTextWidth.canvas || (getTextWidth.canvas = document.createElement("canvas"));
        var context = canvas.getContext("2d");
        context.font = font;
        var metrics = context.measureText(text);
        return metrics.width;
    }

    function populateTable(table, rows, cells, content) {
    /* build table from DOM objects */
        var is_func = (typeof content === 'function');
        if (!table) table = document.createElement('table');
        for (var i = 0; i < rows; i++) {
            var row = document.createElement('tr');
            for (var j = 0; j < cells; j++) {
                row.appendChild(document.createElement('td'));
                var text = !is_func ? (content + '') : content(table, i, j);
                row.cells[j].appendChild(document.createTextNode(text));
            }
            table.appendChild(row);
        }
        return table;
    }

    function legendGetColors(rows) {
    /* go over the legend lines: create an array of colors in the order of their meta value indexes.
     * (the values in the legend may be sorted not in the order of their internal indices) */
        if (rows===undefined)
            return [];

        var colArr = [];
        for (var i = 0; i < rows.length; i++) {
            var row = rows[i];
            var col = row.color;
            if (col===null)
                col = row.defColor; // only use default color if nothing else set

            var idx = row.intKey;
            colArr[idx] = col; // 0 = color
        }

        return colArr;
    }

    function legendUpdateLabels(fieldName) {
        /* re-copy the labels into the legend rows */
        // rows have attributes like these: defColor, currColor, label, count, valueIndex, uniqueKey
        var shortLabels = db.findMetaInfo(fieldName).ui.shortLabels;
        var rows = gLegend.rows;
        for (var i = 0; i < rows.length; i++) {
            var row = rows[i];
            var valIdx = row.intKey;
            var shortLabel = shortLabels[valIdx];
            row.label = shortLabel;
        }
    }

    function legendRemoveManualColors(gLegend) {
        /* remove all manually defined colors from the legend */
        // reset the legend object
        legendSetColors(gLegend, null, "color");

        // reset the URL and local storage settings
        var rows = gLegend.rows;
        var urlChanges = {};
        for (var i = 0; i < rows.length; i++) {
            var row = rows[i];
            var urlVar = COL_PREFIX+row.strKey;
            localStorage.removeItem(urlVar);
            urlChanges[urlVar] = null;

        }
        changeUrl(urlChanges);
    }

    function legendChangePaletteAndRebuild(palName, resetManual) {
        /* change the legend color palette and put it into the URL */
        var success = legendSetPalette(gLegend, palName);
        if (success) {
            if (palName==="default") {
                legendRemoveManualColors(gLegend);
                changeUrl({"pal":null});
            }
            else
                changeUrl({"pal":palName});
            buildLegendBar();
            var colors = legendGetColors(gLegend.rows);
            renderer.setColors(colors);
        }
    }

    function legendSetColors(legend, colors, keyName) {
        /* set the colors for all legend rows, keyName can be "color" or "defColor", depending on
         * whether the current row color or the row default color should be changed.
         * colors can also be null to reset all values to null. */
        if (!keyName)
            keyName = "color";
        var rows = legend.rows;
        for (let i = 0; i < rows.length; i++) {
            var colorVal = null;
            if (colors)
                colorVal = colors[i];

            var legendRow = rows[i];
            if (legendRow.label == "0" && legend.type=="expr")
                colorVal = cNullColor;
            legendRow[keyName] = colorVal;
        }
    }

    function legendSetPalette(legend, origPalName) {
    /* update the defColor [1] attribute of all legend rows. pal is an array of hex colors.
     * Will use the predefined colors that are
     * in the legend.metaInfo.colors configuration, if present.
     * */
        var palName = origPalName;
        if (origPalName==="default") {
            if (legend.rowType==="category")
                palName = datasetQualPalette;
            else
                palName = datasetGradPalette;
        }

        var rows = legend.rows;
        var n = rows.length;
        var pal = null;
        var usePredefined = false;

        pal = makeColorPalette(palName, n);
        // if this is a field for which colors were defined manually, use them
        if (legend.metaInfo!==undefined && legend.metaInfo.colors!==undefined && origPalName==="default") {
            copyNonNull(legend.metaInfo.colors, pal);
            usePredefined = true;
        } else
            pal = makeColorPalette(palName, n);

        if (pal===null) {
            alert("Sorry, palette '"+palName+"' does not have "+rows.length+" different colors");
            return false;
        }

        legendSetColors(legend, pal, "defColor");
        legend.palName = palName;

        // update the dropdown menu
        $('.tpColorLink').parent().removeClass("active");
        // force the menu to the "defaults" entry if we're using predefined colors
        if (usePredefined)
            palName = "default";
        $('.tpColorLink[data-palette="'+palName+'"]').parent().addClass("active");
        return true;
    }

     function labelForBinMinMax(binMin, binMax, isAllInt) {
        /* given the min/max of a numeric value bin, return a good legend for it */
        // pretty print the numbers
        var minDig = 2;
        //if (binMin % 1 === 0) // % 1 = fractional part
            //minDig = 0

        var maxDig = 2;
        //if (binMin % 1 === 0)
         //   maxDig = 0


        if (isAllInt) {
            minDig = 0;
            maxDig = 0
        }

        var legLabel = "";
        if (binMax===0 && binMax===0)
            legLabel = "0";
        else if (binMin==="Unknown")
            legLabel = "Unknown";
        else if (binMin!==binMax) {
            if (Math.abs(binMin) > 1000000)
                binMin = binMin.toPrecision(4);
            if (Math.abs(binMax) > 1000000)
                binMax = binMax.toPrecision(4);
            if (typeof(binMin)=== 'number')
                binMin = binMin.toFixed(minDig);
            if (typeof(binMax)=== 'number')
                binMax = binMax.toFixed(minDig);

            legLabel = binMin+'–'+binMax;
        }
        else
            legLabel = binMin.toFixed(minDig);
        return legLabel;
    }

    function makeLegendRowsNumeric(binInfo) {
        /* return an array of legend lines given bin info from gene expression or a numeric meta field  */
        var legendRows = [];

        // figure out if all our ranges are integers
        var isAllInt = true;
        for (var binIdx = 0; binIdx < binInfo.length; binIdx++) {
            let oneBin = binInfo[binIdx];
            var binMin = oneBin[0];
            var binMax = oneBin[1];

            var restMin = binMin - Math.trunc(binMin);
            var restMax = binMax - Math.trunc(binMax);
            if (restMin!==0 || restMax!==0)
                isAllInt = false;
        }

        var colIdx = 0;
        for (var binIdx = 0; binIdx < binInfo.length; binIdx++) {
            let oneBin = binInfo[binIdx];

            var binMin = oneBin[0];
            var binMax = oneBin[1];
            var count  = oneBin[2];

            var legendId = binIdx;

            var legLabel = labelForBinMinMax(binMin, binMax, isAllInt);

            var uniqueKey = legLabel;

            // override any color with the color specified in the current URL
            var savKey = COL_PREFIX+legLabel;
            var legColor = getVar(savKey, null);

            if (binMin===0 && binMax===0) {
                uniqueKey = "noExpr";
                legColor = cNullColor;
            }
            else if (binMin==="Unknown" && binMax==="Unknown") {
                uniqueKey = "noExpr";
                legColor = cNullColor;
            }
            else
                colIdx++;

            legendRows.push( {
                "color": legColor,
                "defColor":null,
                "label":legLabel,
                "count":count,
                "intKey":binIdx,
                "strKey":uniqueKey
            });
        }
        return legendRows;
    }

    function makeLegendExpr(geneSym, mouseOver, binInfo, exprVec, decExprVec) {
        /* build gLegend object for coloring by expression
         * return the colors as an array of hex codes */

        activateTooltip("#tpGeneSym");

        var legendRows = makeLegendRowsNumeric(binInfo);

        gLegend = {};
        gLegend.type = "expr";
        gLegend.rows = legendRows;
        var subTitle = null;
        if (db.conf.atacSearch)
            gLegend.title = ("sum of "+geneSym.split("+").length) + " peaks";
        else {
            //  make a best effort to find the gene sym and gene ID
            var geneInfo = db.getGeneInfo(geneSym);
            geneSym = geneInfo.sym;
            subTitle = geneInfo.id;
            gLegend.title = getGeneLabel()+": "+geneSym;
        }

        gLegend.titleHover = mouseOver;
        gLegend.geneSym = geneSym;
        gLegend.subTitle = mouseOver;
        gLegend.rowType = "range";
        gLegend.exprVec = exprVec; // raw expression values, e.g. floats
        gLegend.decExprVec = decExprVec; // expression values as deciles, array of bytes
        gLegend.selectionDirection = "all";
        var oldPal = getVar("pal", "default")
        legendSetPalette(gLegend, oldPal);

        var colors = legendGetColors(legendRows);
        return colors;
    }

    function alphaNumSearch(genes, saneSym) {
        /* use alpha-num-only search in gene list for saneSym */
        for (var i=0; i<genes.length; i++) {
            var geneSym = genes[i][0];
            if (saneSym===onlyAlphaNum(geneSym))
                return geneSym;
        }
        return null;
    }

    function phoneHome() {
      /* add empty javascript to document so we can count usage at UCSC for grant reports */
      /* -> does this pose a data protection issue? Need to document. Does it require opt-in ? */
      var s = document.createElement( 'script' );
      s.setAttribute( 'src', "https://cells.ucsc.edu/js/cbTrackUsage.js" );
      s.async = true;
      document.body.appendChild( s );
    }

    function onGeneClick (event) {
    /* user clicked on a gene in the gene table */
        var saneSym = event.target.id.split("_")[1]; // the symbol of the gene, as only-alphaNum chars
        $('.tpMetaBox').removeClass('tpMetaSelect');
        $('.tpGeneBarCell').removeClass("tpGeneBarCellSelected");
        $('#tpGeneBarCell_'+saneSym).addClass("tpGeneBarCellSelected");

        // search through both quick and recent gene lists to find the real gene symbol
        var geneSym = null;
        var quickGenes = db.conf.quickGenes;
        if (quickGenes)
            geneSym = alphaNumSearch(quickGenes, saneSym);
        if (geneSym===null)
            geneSym = alphaNumSearch(gRecentGenes, saneSym);

        colorByLocus(geneSym);
        event.stopPropagation();
    }

    function showDialogBox(htmlLines, title, options) {
        /* show a dialog box with html in it */
        $('#tpDialog').remove();

        var addStr = "";
        if (options.width!==undefined)
            addStr = "max-width:"+options.width+"px;";
        var maxHeight = $(window).height()-200;
        // unshift = insert at pos 0
        //htmlLines.unshift("<div style='display:none;"+addStr+"max-height:"+maxHeight+"px' id='tpDialog' title='"+title+"'>");
        htmlLines.unshift("<div style='display:none;"+addStr+"' id='tpDialog' title='"+title+"'>");
        htmlLines.push("</div>");
        $(document.body).append(htmlLines.join(""));

        var dialogOpts = {modal:true, closeOnEscape:true};
        if (options.width!==undefined)
            dialogOpts["width"] = options.width;
        if (options.height!==undefined)
            dialogOpts["height"] = options.height;
        //dialogOpts["maxHeight"] = maxHeight;
        if (options.buttons!==undefined)
            dialogOpts["buttons"] =  options.buttons;
        else
            dialogOpts["buttons"] =  [];

        if (options.showClose!==undefined)
            dialogOpts["buttons"].unshift({ text:"Cancel", click:function() { $( this ).dialog( "close" ); }});
        if (options.showOk!==undefined)
            dialogOpts["buttons"].push({ text: "OK", click:function() { $( this ).dialog( "close" ); }});
        //dialogOpts["position"] = "center";
        //dialogOpts["height"] = "auto";
        //dialogOpts["width"] = "auto";

        $( "#tpDialog" ).dialog(dialogOpts);
    }

    function onChangeGenesClick(ev) {
    /* called when user clicks the "change" button in the gene list */
        var htmls = [];
        htmls.push("<p style='padding-bottom: 5px'>Enter a list of gene symbols, one per line:</p>");
        htmls.push("<textarea id='tpGeneListBox' class='form-control' style='height:320px'>");

        var geneFields = gCurrentDataset.preloadGenes.genes;
        for (var i = 0; i < geneFields.length; i++) {
            htmls.push(geneFields[i][1]+"\r\n");
        }
        htmls.push("</textarea>");
        //htmls.push("<p>");
        //htmls.push("<button style='float:center;' id='tpGeneDialogOk' class='ui-button ui-widget ui-corner-all'>OK</button>");
        //htmls.push("</div>");
        //htmls.push("</div>");
        var buttons = [{text:"OK", click:onGeneDialogOkClick}];
        showDialogBox(htmls, "Genes of interest", {height: 500, width:400, buttons:buttons});

        $('#tpGeneDialogOk').click ( onGeneDialogOkClick );
    }

    function onGeneLoadComplete() {
        /* called when all gene expression vectors have been loaded */
        console.log("All genes complete");
        // Close the dialog box only if all genes were OK. The user needs to see the list of skipped genes
        if ( $( "#tpNotFoundGenes" ).length===0 ) {
            $("#tpDialog").dialog("close");
        }

        var cellIds = getSortedCellIds();

        // transform the data from tpLoadGeneExpr as gene -> list of values (one per cellId) to
        // newExprData cellId -> list of values (one per gene)

        // initialize an empty dict cellId = floats, one per gene
        var newExprData = {};
        var geneCount = gLoad_geneList.length;
        var cellIdCount = cellIds.length;
        for (var i = 0; i < cellIdCount; i++) {
            let cellId = cellIds[i];
            newExprData[cellId] = new Float32Array(gLoad_geneList); // XX or a = [] + a.length = x ?
        }

        // do the data transform
        var newGeneFields = [];
        var newDeciles = {};
        for (var geneI = 0; geneI < gLoad_geneList.length; geneI++) {
            var geneSym = gLoad_geneList[geneI];
            var geneInfo = gLoad_geneExpr[geneSym];

            var geneId = geneInfo[0];
            var geneVec = geneInfo[1];
            var deciles = geneInfo[2];
            newGeneFields.push( [geneSym, geneId] );
            newDeciles[geneSym] = deciles;

            for (var cellI = 0; cellI < cellIdCount; cellI++) {
                let cellId = cellIds[cellI];
                newExprData[cellId][geneI] = geneVec[cellI];
            }
        }

        gLoad_geneList = null;
        gLoad_geneExpr = null;

        gCurrentDataset.preloadExpr = {};
        gCurrentDataset.preloadExpr.genes = newGeneFields;
        gCurrentDataset.preloadExpr.cellExpr = newExprData;
        gCurrentDataset.preloadExpr.deciles = newDeciles;

        //tpGeneTable();
    }

    //function onReceiveExprLineProgress(line) {
        /* called when a line of the expression matrix has been loaded: parse line and upd. progressbar */
        //var symbol = this.geneSymbol;
        //console.log("Got gene "+symbol);
        //var exprTuple = parseMatrixLine(line);
        //var exprVec = exprTuple[1];
        //exprTuple.push( getDeciles(exprVec) );
        //gLoad_geneExpr[symbol] = exprTuple;
//
        //var progressbar = $( "#tpGeneProgress" );
        //var val = progressbar.progressbar( "value" ) || 0 ;
        //val++;
        //progressbar.progressbar( "value", val );
        //$( "#tpProgressLabel" ).text(symbol);
//
        //var progrMax = progressbar.progressbar("option", "max");
        //if (val >= progrMax)
            //onGeneLoadComplete();
    //}

    //function singleExprDone() {
    ///* called when both cellIds(header) and single line of the expr matrix have been read */
    //    var sExpr = gCurrentDataset.singleExpr;
    //    if (sExpr.exprVec===undefined || gCurrentDataset.matrixCellIds==undefined)
    //        // don't proceed if one part of the data is not here yet
    //        return;

    //    // convert the expression vector to a dict cellId -> array of float, with a single value
    //    var cellToExpr = {};
    //    var exprVec = sExpr.exprVec;
    //    var cellIds = gCurrentDataset.matrixCellIds;
    //    for (var i = 0; i < cellIds.length; i++) {
    //        var cellId = cellIds[i];
    //        var val = [exprVec[i]];
    //        cellToExpr[cellId] = val;
    //    }

    //    var geneSym = sExpr.symbol;
    //    console.log("XX here with symbol"+geneSym);
    //    if (geneSym !== null) {
    //        console.log("XX adding to combo"+geneSym);
    //        var select = $("#tpGeneCombo").selectize();
    //        var selectize = select[0].selectize;
    //        //var optId = selectize.search(geneSym);
    //        selectize.addOption({text: geneSym, value: geneSym});
    //        selectize.refreshOptions();
    //        selectize.setTextboxValue(geneSym);
    //    }


    //}

    //function onReceiveMatrixSingleExprLine(line) {
    //    /* got a single expression line from the matrix, no progressbar */
    //    var sExpr = gCurrentDataset.singleExpr;
    //    sExpr.symbol = this.geneSymbol;
    //    console.log("Got gene vector for "+sExpr.symbol);
    //    var exprTuple = parseMatrixLine(line);
    //    sExpr.geneId  = exprTuple[0];
    //    sExpr.exprVec = exprTuple[1];
    //    sExpr.deciles = getDeciles(sExpr.exprVec);

    //    singleExprDone();
    //}

    //function onReceiveMatrixHeader(line) {
    //    /* received the */
    //    line = line.trim(); // old scripts always kept the newline
    //    var fields = line.split("\t");
    //    var cellIds = fields.slice(1); // copy all but first
    //    console.log("Got header, cellCount="+cellIds.length);
    //    gCurrentDataset.matrixCellIds = cellIds;
    //    singleExprDone();
    //}

    function onGeneDialogOkClick(ev) {
    /* called the user clicks the OK button on the 'paste your genes' dialog box */
        var genes = $('#tpGeneListBox').val().replace(/\r\n/g, "\n").split("\n");
        $("#tpDialog").remove();

        gLoad_geneList = [];
        gLoad_geneExpr = {};

        var notFoundGenes = [];
        var baseUrl = gCurrentDataset.baseUrl;
        var url = cbUtil.joinPaths([baseUrl, "geneMatrix.tsv"]);
        var validCount = 0; // needed for progressbar later
        var matrixOffsets = gCurrentDataset.matrixOffsets;
        for (var i = 0; i < genes.length; i++) {
            var gene = genes[i];
            if (gene.length===0) // skip empty lines
                continue;
            if (!(gene in matrixOffsets))
                {
                notFoundGenes.push(gene);
                continue;
                }

            gLoad_geneList.push(gene);
            var offsetInfo = matrixOffsets[gene];
            var start = offsetInfo[0];
            var end = start+offsetInfo[1];
            jQuery.ajax( { url: url,
                headers: { Range: "bytes="+start+"-"+end } ,
                geneSymbol : gene,
                success: onReceiveExprLineProgress
            });
        }

        var htmls = [];
        htmls.push("<div id='tpGeneProgress'><div class='tpProgressLabel' id='tpProgressLabel'>Loading...</div></div>");

        if (notFoundGenes.length!==0) {
            htmls.push("<div id='tpNotFoundGenes'>Could not find the following gene symbols: ");
            htmls.push(notFoundGenes.join(", "));
            htmls.push("</div>");
            //htmls.push("<button style='float:center;' id='tpGeneDialogOk' class='ui-button ui-widget ui-corner-all'>OK</button>");
            //for (var i = 0; i < notFoundGenes.length; i++) {
                //htmls.push(notFoundGenes[i]);
            //}
            //htmls.push("<p>Could not find the following gene symbols:</p>");
            //showDialogBox(htmls, "Warning", 400);
        }

        var showOk = (notFoundGenes.length!==0);
        showDialogBox(htmls, "Downloading expression data", {width:350, showOk:true});

        var progressLabel = $( "#tpProgressLabel" );
        $("#tpGeneProgress").progressbar( {
              value: false,
              max  : gLoad_geneList.length
              });

    }

    function buildGeneTable(htmls, divId, title, subtitle, geneInfos, noteStr) {
    /* create gene expression info table. if htmls is null, update DIV with divId in-place. */
        var doUpdate = false;
        if (htmls===null) {
            htmls = [];
            doUpdate = true;
        }

        var tableWidth = metaBarWidth;
        //var cellWidth = gGeneCellWidth;

        if (title) {
            htmls.push("<div style='margin-top:8px' id='"+divId+"_title'>");
            htmls.push("<div style='padding-left:3px; font-weight:bold'>"+title+"</div>");
            if (subtitle) {
                htmls.push('<div style="margin-top:6px" class="tpHint">');
                htmls.push(subtitle);
                htmls.push('</div>');
            }
            htmls.push("</div>"); // divId_title
        }

        if (doUpdate) {
            $('#'+divId).empty();
        }

        htmls.push("<div id='"+divId+"'>");

        if (geneInfos===undefined || geneInfos===null || geneInfos.length===0) {
            if (noteStr!==undefined)
                htmls.push("<div style='font-style:80%'>"+noteStr+"</div>");
            htmls.push("</div>");
            return;
        }

        htmls.push('<table style="margin-top:10px" id="tpGeneTable"><tr>');
        //htmls.push('<td><button id="tpChangeGenes" title="Change the list of genes that are displayed in this table" class = "ui-button ui-widget ui-corner-all" style="width:95%">Change</button></td>');

        // need max length of gene names to make number of columns
        var maxColLen = 0;
        for (var i=0; i < geneInfos.length; i++) {
            var geneId = geneInfos[i][0];
            maxColLen = Math.max(maxColLen, geneId.length);
        }

        var lineLen = 38;
        var colsPerRow = Math.floor(lineLen/maxColLen);
        //var colsPerRow = Math.round(tableWidth / cellWidth);
        var cellWidth = Math.round(tableWidth/colsPerRow);

        var currWidth = 1;
        var i = 0;
        while (i < geneInfos.length) {
            var geneInfo = geneInfos[i];
            var geneId   = geneInfo[0];
            var geneDesc = geneInfo[1];
            if (geneDesc===undefined)
                geneDesc = geneId;

            if (((i % colsPerRow) === 0) && (i!==0)) {
                htmls.push("</tr><tr>");
            }
            if (geneId in db.geneOffsets)
                htmls.push('<td title="'+geneDesc+'" id="tpGeneBarCell_'+onlyAlphaNum(geneId)+'" class="tpGeneBarCell" style="width:'+cellWidth+'px">'+geneId+'</td>');
            i++;
        }
        htmls.push("</tr></table>");
        htmls.push("</div>"); // divId

        if (doUpdate) {
            $('#'+divId).html(htmls.join(""));
        }
    }

    function likeEmptyString(label) {
    /* some special values like "undefined" and empty string get colored in grey  */
        return (label===null || label.trim()==="" || label==="none" || label==="None" || label==="unknown" ||
            label==="nd" || label==="n.d." ||
            label==="Unknown" || label==="NaN" || label==="NA" || label==="undefined" || label==="Na");
    }

    function numMetaToBinInfo(fieldInfo) {
        /* convert a numeric meta field info to look like gene expression info for the legend:
         * an array of [start, end, count] */
        var binInfo = [];
        var binMethod = fieldInfo.binMethod;
        if (binMethod==="uniform") {
            // old method, not used anymore
            let binMin = fieldInfo.minVal;
            var stepSize = fieldInfo.stepSize;
            let binCounts = fieldInfo.binCounts;
            let binCount = fieldInfo.binCounts.length;
            for (var i=0; i<binCount; i++) {
                binInfo.push( [binMin, binMin+stepSize, binCounts[i]] );
                binMin+=stepSize;
            }
        } else if (binMethod==="quantiles") {
            // newer method for pre-quantified fields
            let binMin = fieldInfo.minVal;
            let breaks = fieldInfo.breaks;
            let binCounts = fieldInfo.binCounts;
            let binCountsLen = fieldInfo.binCounts.length;
            let binIdx = 0;
            if (breaks[0]==="Unknown") {
                binInfo.push( ["Unknown", "Unknown", binCounts[0]] );
                binIdx = 1;
            }

            for (; binIdx<binCountsLen; binIdx++) {
                let binMin = breaks[binIdx];
                let binMax = breaks[binIdx+1];
                let binCount = binCounts[binIdx];
                binInfo.push( [binMin, binMax, binCount] );
            }
        }
        else if (binMethod==="raw") {
            // no binning was done at all, this happens when there are not enough distinct values
            let values = fieldInfo.values;
            let binCounts = fieldInfo.binCounts;
            for (let i=0; i<values.length; i++) {
                let value = values[i];
                let valCount = binCounts[i];
                binInfo.push( [value, value, valCount] );
            }
        }
        else
            // these days, we don't pre-quantify anymore. The binning is done on the client now
            // and the meta loading function adds the bin info to the field info object
            //alert("invalid value for meta field binMethod: "+binMethod);
            binInfo = fieldInfo.binInfo;

        return binInfo;
    }

    function getFieldColors(metaInfo) {
        /* return a list with the color codes for a field, taking into account field type and cart settings */

        let palName = metaInfo.ui.palette;
        let colCount;
        if (!palName)
            if (metaInfo.type==="int" || metaInfo.type==="float") {
                palName = datasetGradPalette;
                colCount = exprBinCount+1; // +1 because <unknown> is a special bin
            }
            else {
                palName = cDefQualPalette;
                colCount = metaInfo.valCounts.length;
            }

        var colors = makeColorPalette(palName, colCount);

        let colOverrides = metaInfo.ui.colors;
        if (colOverrides)
            for (let i=0; i<colOverrides.length; i++)
                if (colOverrides[i])
                    colors[i] = colOverrides[i];

        return colors;
    }

    function makeLegendMeta(metaInfo, sortBy) {
    /* Build a new legend object and return it */
        var legend = {};
        legend.type = "meta";
        //legend.metaFieldIdx = metaIndex;
        legend.titleHover = null;

        //var metaInfo = db.getMetaFields()[metaIndex];
        //legend.fieldLabel = metaInfo.label;
        legend.title = metaInfo.label.replace(/_/g, " ");
        legend.metaInfo = metaInfo;

        // numeric meta fields are a special case
        if (metaInfo.type==="int" || metaInfo.type==="float") {
            var binInfo = numMetaToBinInfo(metaInfo);
            legend.rows = makeLegendRowsNumeric(binInfo);
            legend.rowType = "range";
            legendSetPalette(legend, "default");
            return legend;
        }

        if (metaInfo.diffValCount > MAXCOLORCOUNT) {
            warn("This field has "+metaInfo.diffValCount+" different values. Coloring on a "+
                "field that has more than "+MAXCOLORCOUNT+" different values is not supported.");
            return null;
        }

        var metaCounts = metaInfo.valCounts;

        // we are going to sort this list later, so we need to keep track of what the original
        // index in the list was, as every meta data value is stored by its index, not
        // its label. Simply append the index as [2] of the metaCounts array.
        for (var valIdx=0; valIdx < metaCounts.length; valIdx++)
            metaCounts[valIdx].push(valIdx);

        var oldSortBy = getFromUrl("SORT");
        // URL overrides default value
        if (sortBy===undefined && oldSortBy!==undefined)
            sortBy = oldSortBy;

        // default sorting can be specfied with "sortBy" in cellbrowser.conf
        if (sortBy===undefined && metaInfo.sortBy)
            sortBy = metaInfo.sortBy;
        if (sortBy!==undefined && sortBy!=="freq" && sortBy!=="name" && sortBy!=="none") {
            alert("sortBy is '"+cleanString(sortBy)+' but it can only be "freq" or "name"');
            sortBy = undefined;
        }

        var fieldName = metaInfo.label;

        // force field names that look like "cluster" to a rainbow palette
        // even if they are numbers
        if (sortBy===undefined) {
            // should cluster fields be sorted by their name
            if (metaInfo.type==="float" || metaInfo.type==="int" || (metaCounts.length > 60))
                // long lists are easier to grasp if they're sorted by name
                sortBy = "name";
            else if (fieldName.indexOf("luster") || fieldName.indexOf("ouvain") || fieldName.indexOf("res."))
                sortBy = "count";
            else
                sortBy = "count";
        }

        // sort like numbers if the strings are mostly numbers, otherwise sort like strings
        var sortResult = sortPairsBy(metaCounts, sortBy);
        var countListSorted = sortResult.list;

        var useGradient = (metaInfo.type==="float" || metaInfo.type==="int");
        //var defaultColors = makeColorPalette(countListSorted.length, useGradient);

        var rows = [];
        var shortLabels = metaInfo.ui.shortLabels;
        var longLabels = metaInfo.ui.longLabels;
        for (var legRowIdx = 0; legRowIdx < countListSorted.length; legRowIdx++) {
            var legRowInfo = countListSorted[legRowIdx];
            let valIdx = legRowInfo[2]; // index of the original field in fieldInfo
            var label = shortLabels[valIdx];

            var desc  = null;
            if (longLabels)
                desc  = longLabels[valIdx];

            // first use the default palette, then try to get from URL
            var count = legRowInfo[1];
            var uniqueKey = label;
            if (uniqueKey==="")
                uniqueKey = "_EMPTY_";

            if (likeEmptyString(label))
                color = cNullColor;
            // override any color with the color specified in the current URL
            var savKey = COL_PREFIX+uniqueKey;
            var color = getFromUrl(savKey, null);

            rows.push( {
                "color": color,
                "defColor": null,
                "label": label,
                "count": count,
                "intKey":valIdx,
                "strKey":uniqueKey,
                "longLabel" : desc,
            } );
        }

        legend.rows = rows;
        legend.isSortedByName = sortResult.isSortedByName;
        legend.rowType = "category";
        legend.selectionDirection = "all";
        legendSetPalette(legend, "default");
        return legend;
    }

    function legendSetTitle(label) {
        $('#tpLegendTitle').text(label);
    }

    function buildLegendForMeta(metaInfo) {
    /* build the gLegend for a meta field */
        var legend = makeLegendMeta(metaInfo);
        if (legend===null)
            return;

        var metaIdx = db.fieldNameToIndex(metaInfo.name);
        $('.tpMetaBox').removeClass('tpMetaSelect');
        $('.tpMetaValue').removeClass('tpMetaValueSelect');
        $('#tpMetaBox_'+metaIdx).addClass('tpMetaSelect');
        $('#tpMeta_'+metaIdx).addClass('tpMetaValueSelect');
        $('.tpGeneBarCell').removeClass('tpGeneBarCellSelected');
        //$('#tpLegendTitle').text(legend.metaInfo.label.replace(/_/g, " "));
        legendSetTitle(legend.metaInfo.label.replace(/_/g, " "));

        return legend;
    }

    function onMetaClick (event) {
    /* called when user clicks a meta data field or label */
        var fieldName = event.target.dataset.fieldName;
        if (isNaN(fieldName)) {
            // try up one level in the DOM tree
            fieldName = event.target.parentElement.dataset.fieldName;
        }
        colorByMetaField(fieldName);
    }

    function addMetaTipBar(htmls, valFrac, valStr, valFracCategory) {
        /* add another bar to a simple histogram built from divs */
        htmls.push("<div>&nbsp;");
        htmls.push("<div class='tpMetaTipPerc'>"+(100*valFrac).toFixed(1)+"%</div>");
        htmls.push("<div class='tpMetaTipName'>"+valStr);
        if (valFracCategory !== undefined) {
            htmls.push(" <small>(" + (100 * valFracCategory).toFixed(1) + "% of all cells with this value)</small>");
        }
        htmls.push("</div>");
        //htmls.push("<span class='tpMetaTipCount'>"+valCount+"</span>");
        var pxSize = (valFrac * metaTipWidth).toFixed(0);
        htmls.push("<div style='width:"+pxSize+"px' class='tpMetaTipBar'>&nbsp</div>");
        htmls.push("</div>");
    }

    function binInfoToValCounts(binInfo) {
        /* given an array of (start, end, count), return an array of (label, count) */
        var valCounts = [];

        for (var binIdx = 0; binIdx < binInfo.length; binIdx++) {
            var binMin = binInfo[binIdx][0];
            var binMax = binInfo[binIdx][1];
            var count  = binInfo[binIdx][2];
            var label = labelForBinMinMax(binMin, binMax);
            valCounts.push( [label, count] );
        }
        return valCounts;
    }

    function buildMetaTip(metaInfo, valHist, htmls) {
        /* build the content of the tooltip that summarizes the multi selection */
        var valCounts = metaInfo.valCounts;
        var shortLabels = metaInfo.ui.shortLabels;
        if (valCounts===undefined)  // for client-side discretized fields, we have to discretize first
            valCounts = binInfoToValCounts(metaInfo.binInfo);

        var otherCount = 0;
        var totalSum = 0;
        for (var i = 0; i < valHist.length; i++) {
            var valInfo  = valHist[i];
            var valCount = valInfo[0];
            var valFrac  = valInfo[1];
            var valIdx   = valInfo[2];
            var valFracCategory = valInfo[3];
            //var valStr   = valCounts[valIdx][0]; // 0 = label, 1 = count
            var label   = shortLabels[valIdx];

            totalSum += valCount;
            // only show the top values, summarize everything else into "other"
            if (i > HISTOCOUNT) {
                otherCount += valCount;
                continue;
            }

            if (label==="")
                label = "<span style='color:indigo'>(empty)</span>";
            addMetaTipBar(htmls, valFrac, label, valFracCategory);
        }

        if (otherCount!==0) {
            var otherFrac = (otherCount / totalSum);
            addMetaTipBar(htmls, otherFrac, "<span style='color:indigo'>(other)</span>");
        }
        return htmls;
    }

    function metaInfoFromElement(target) {
        /* get the metaInfo object given a DOM element  */
        if (target.dataset.fieldName === undefined)
            target = target.parentNode;
        if (target.dataset.fieldName === undefined)
            target = target.parentNode;
        var fieldName = target.dataset.fieldName;
        var metaInfo = db.findMetaInfo(fieldName);
        return metaInfo;
    }

    function onMetaMouseOver (event) {
    /* called when user hovers over meta element: shows the histogram of selected cells */
        var metaHist = db.metaHist;
        if (metaHist===undefined || metaHist===null)
            return;

        // mouseover over spans or divs will not find the id, so look at their parent, which is the main DIV
        var target = event.target;

        var metaInfo = metaInfoFromElement(target);
        var fieldName = metaInfo.name;

        // change style of this field a little
        var metaSel = "#tpMetaBox_"+metaInfo.index;
        //var backCol = "#666";
        //var foreCol = "#FFF";
        //$(metaSel).css({color: foreCol, backgroundColor: backCol});
        //$(metaSel).children().css({color: foreCol, backgroundColor: backCol});
        //$(metaSel).children().children().css({color: foreCol, backgroundColor: backCol});
        $(metaSel).addClass("tpMetaHover");
        $(metaSel+" .tpMetaValue").addClass("tpMetaHover");

        var htmls = [];

        if (metaInfo.type==="uniqueString")
            htmls.push("<div>Cannot summarize: this is a field with unique values</div>");
        else
            htmls = buildMetaTip(metaInfo, metaHist[fieldName], htmls);

        $('#tpMetaTip').html(htmls.join(""));

        // make sure that tooltip doesn't go outside of screen
        //var tipTop = event.target.offsetTop;
        var tipTop = event.target.getBoundingClientRect().top-8;
        var tipHeight = $('#tpMetaTip').height();
        var screenHeight = $(window).height();
        if (tipTop+tipHeight > screenHeight)
            tipTop = screenHeight - tipHeight - 8;

        $('#tpMetaTip').css({top: tipTop+"px", left: metaBarWidth+"px", width:metaTipWidth+"px"});
        $('#tpMetaTip').show();
    }

    function buildComboBox(htmls, id, entries, selIdx, placeholder, width, opts) {
    /* make html for a combo box and add lines to htmls list.
     * selIdx is an array of values if opt.multi exists, otherwise it's an int or 'undefined' if none. */

        let addStr = "";
        if (opts && opts.multi)
            addStr= " multiple";

        htmls.push('<select style=width:"'+width+'px" id="'+id+'"'+addStr+' data-placeholder="'+placeholder+'" class="tpCombo">');
        for (var i = 0; i < entries.length; i++) {
            var isSelStr = "";
            var entry = entries[i];

            // entries can be either key-val or just values
            var name, label;
            if (Array.isArray(entry)) {
                name = entry[0];
                label = entry[1];
            } else {
                name = label = entry;
            }

            // determine if element is selected
            let isSel = false;
            if (opts && opts.multi) {
                if (selIdx.indexOf(i)!==-1 || selIdx.indexOf(name)!==-1)
                    isSel = true;
            } else if ((selIdx!==undefined && i===selIdx))
                isSel = true;

            if (isSel)
                isSelStr = " selected";

            htmls.push('<option value="'+name+'"'+isSelStr+'>'+label+'</option>');
        }
        htmls.push('</select>');
    }

    function loadCoordSet(coordIdx) {
        /* load coordinates and color by meta data */
        var newRadius = db.conf.coords[coordIdx].radius;
        var colorOnMetaField = db.conf.coords[coordIdx].colorOnMeta;

        db.loadCoords(coordIdx,
                function(coords, info, clusterMids) {
                    gotCoords(coords,info,clusterMids, newRadius);
                    if (colorOnMetaField!==undefined)
                        colorByMetaField(colorOnMetaField);
                    else
                        renderer.drawDots();
                },
                onProgress);
    }

    function onLayoutChange(ev, params) {
        /* user changed the layout in the combobox */
        var coordIdx = parseInt(params.selected);
        loadCoordSet(coordIdx);
        changeUrl({"layout":coordIdx, "zoom":null});
        renderer.coordIdx = coordIdx;
        // remove the focus from the combo box
        removeFocus();
    }

    function onGeneComboChange(ev) {
        /* user changed the gene in the combobox */
        var geneId = ev.target.value;
        if (geneId==="")
            return; // do nothing if user just deleted the current gene
        if (db.conf.atacSearch) {
            updatePeakListWithGene(geneId);
        } else {
            // in the normal, gene-matrix mode.
            var locusStr = null;
            var geneInfo = db.getGeneInfo(geneId);
            //if (geneInfo.id!==geneInfo.sym)
                //locusStr = geneInfo.id;
            //else
                //locusStr = geneInfo.sym
            colorByLocus(geneInfo.id);
        }
    }

    function onMetaComboChange(ev, choice) {
        /* called when user changes the meta field combo box */
        //if (choice.selected==="_none")
        var fieldId = parseInt(choice.selected.split("_")[1]);
        var fieldName = db.getMetaFields()[fieldId].name;
        console.log(choice);
        console.log(ev);

        colorByMetaField(fieldName);
    }

    function onLabelComboChange(ev, choice) {
        /* called when user changes the label field combo box */
        var fieldId = parseInt(choice.selected.split("_")[1]);
        var fieldName = db.getMetaFields()[fieldId].name;

        setLabelField(fieldName);
    }

    function showCollectionDialog(collName) {
        /* load collection with given name and open dialog box for it */
        loadCollectionInfo(collName, function(collData) { openDatasetDialog(collData)});
    }

    function onConfigLoaded(datasetName) {
            // this is a collection if it does not have any field information
            if (db.conf.sampleDesc)
                gSampleDesc = db.conf.sampleDesc;
            else
                gSampleDesc = "cell";

            // allow config to override the default palettes
            datasetGradPalette = cDefGradPalette;
            datasetQualPalette = cDefQualPalette;
            if (db.conf.defQuantPal)
                datasetGradPalette = db.conf.defQuantPal;
            if (db.conf.defCatPal)
                datasetQualPalette = db.conf.defCatPal;


            if (db.conf.metaBarWidth)
                metaBarWidth = db.conf.metaBarWidth;
            else
                metaBarWidth = 250;

            renderer.setPos(null, metaBarWidth+metaBarMargin);

            if (!db.conf.metaFields) {
                // pablo often has single-dataset installations, there is no need to open the
                // dataset selection box then.
                if (db.conf.datasets.length===1 && datasetName==="") // "" is the root dataset
                    loadDataset(db.conf.datasets[0].name, false);
                else
                    showCollectionDialog(datasetName);
                return;
            }

            let binData = localStorage.getItem(db.name+"|custom");
            if (binData) {
                let jsonStr = LZString.decompress(binData);
                let customMeta = JSON.parse(jsonStr);
                db.conf.metaFields.unshift(customMeta);
            }

            cartLoad(db);
            renderData();
            resizeDivs(true);

            if (getVar("openDialog"))
                openDatasetDialog(db.conf, db.name); // open Info dialog
            cartSave(db); // = set the current URL from local storage settings

            // start the tutorial after a while
            var introShownBefore = localStorage.getItem("introShown");
            if (introShownBefore===undefined)
               setTimeout(function(){ showIntro(true); }, 3000); // shown after 5 secs
    }

    function loadDataset(datasetName, resetVars, md5) {
        /* load a dataset and optionally reset all the URL variables.
         * When a dataset is opened through the UI, the variables have to
         * be reset, as their values (gene or meta data) may not exist
         * there. If it's opened via a URL, the variables must stay. */

        // collections are not real datasets, so ask user which one they want
        if (db!==null && db.heatmap)
            removeHeatmap();
        removeSplit();

        db = new CbDbFile(datasetName);
        cellbrowser.db = db; // easier debugging

        var vars;
        if (resetVars)
            vars = {};

        if (datasetName!=="")
            changeUrl({"ds":datasetName.replace(/\//g, " ")}, vars); // + is easier to type than %23

        db.loadConfig(onConfigLoaded, md5);
        trackEvent("open_dataset", datasetName);
        trackEventObj("select_content", {content_type: "dataset", item_id: datasetName});
    }

    function loadCollectionInfo(collName, onDone) {
        /* load collection info and run onDone */
        var jsonUrl = cbUtil.joinPaths([collName, "dataset.json"]);
        cbUtil.loadJson(jsonUrl, onDone);
    }

    function onDatasetChange(ev, params) {
        /* user changed the dataset in the collection dropdown box */
        /* jshint validthis: true */
        $(this).blur();
        removeFocus();

        var parts = params.selected.split("?");
        var datasetName = parts[0];
        var md5 = parts[1];
        loadDataset(datasetName, true, md5);
    }

    function buildLayoutCombo(coordLabel, htmls, files, id, width, left, top) {
        /* files is a list of elements with a shortLabel attribute. Build combobox for them. */
        if (!coordLabel)
            coordLabel = "Layout";

        htmls.push('<div class="tpToolBarItem" style="position:absolute;left:'+left+'px;top:'+top+'px"><label for="'+id+'">');
        htmls.push(coordLabel);
        htmls.push("</label>");

        var entries = [];
        for (var i = 0; i < files.length; i++) {
            var coordFiles = files[i];
            var label = coordFiles.shortLabel;
            if (label===undefined)
                warn("Layout coordinate file "+i+" has no .shortLabel attribute");
            entries.push([i, label]);
        }
        buildComboBox(htmls, id, entries, 0, "Select a layout algorithm...", width);
        htmls.push('</div>');
    }

    function buildCollectionCombo(htmls, id, width, left, top) {
        /* build combobox with shortLabels of all datasets that are part of same collection */
        htmls.push('<div class="tpToolBarItem" style="position:absolute;width:'+width+'px;left:'+left+'px;top:'+top+'px"><label for="'+id+'">Collection</label>');

        var entries = [];
        //var linkedDatasets = parentConf.datasets;
        //for (var i = 0; i < linkedDatasets.length; i++) {
            //var dsInfo = linkedDatasets[i];
            //entries.push( [i, dsInfo.shortLabel] );
        //}

        buildComboBox(htmls, id, entries, 0, "Select a dataset...", width);
        htmls.push('</div>');
    }

    function buildMetaFieldCombo(htmls, idOuter, id, left, selectedField) {
        var metaFieldInfo = db.getMetaFields();
        htmls.push('<div id="'+idOuter+'" style="padding-left:2px; display:inline">');
        var entries = [["_none", ""]];
        var selIdx = 0;
        for (var i = 1; i < metaFieldInfo.length; i++) { // starts at 1, skip ID field
            var field = metaFieldInfo[i];
            var fieldName = field.label;
            //var hasTooManyVals = (field.diffValCount>100);
            //if (hasTooManyVals)
                //continue;
            entries.push( ["tpMetaVal_"+i, fieldName] );
            if (selectedField == fieldName) {
                selIdx = i;
            }
        }

        buildComboBox(htmls, id, entries, selIdx, "select a field...", metaBarWidth+50);
        htmls.push('</div>');
    }

    function getGeneLabel() {
        /* some datasets have data not on genes, but on other things e.g. "lipids". The config can
         * define a label for the rows in the expression matrix */
        var geneLabel = "Gene";
        if (db.conf.atacSearch)
            geneLabel = "Range";
        if (db.conf.geneLabel)
            geneLabel = db.conf.geneLabel;
        return geneLabel;
    }

    function buildGeneCombo(htmls, id, left, width) {
        /* Combobox that allows searching for genes */
        htmls.push('<div class="tpToolBarItem" style="padding-left: 3px">');
        var title = "Color by "+getGeneLabel();
        if (db.conf.atacSearch)
            title = "Find peaks at or close to:"
        htmls.push('<label style="display:block; margin-bottom:8px; padding-top: 8px;" for="'+id+'">'+title+'</label>');
        var geneLabel = getGeneLabel().toLowerCase();
        var boxLabel = 'search for '+geneLabel+'...';
        if (db.conf.atacSearch)
            boxLabel = "enter gene or chrom:start-end";
        htmls.push('<select style="width:'+width+'px" id="'+id+'" placeholder="'+boxLabel+'" class="tpCombo">');
        htmls.push('</select>');
        htmls.push('</div>');
        //htmls.push("<button>Multi-Gene</button>");
    }


    function activateCombobox(id, widthPx) {
        $('#'+id).chosen({
            inherit_select_classes : true,
            disable_search_threshold: 10,
            width : widthPx
        });
    }

    function updateCollectionCombo(id, collData) {
        /* load dataset sibling labels into collection combobox from json */
        var htmls = [];
        var datasets = collData.datasets;
        for (var i = 0; i < datasets.length; i++) {
            var ds = datasets[i];
            var selStr =  "";
            if (ds.name===db.conf.name)
                selStr = "selected";
            var val = ds["name"]+"?"+ds["md5"];
            htmls.push('<option value="'+val+'"'+selStr+'>'+ds.shortLabel+'</option>');
        }
        $('#'+id).html(htmls.join(""));
        $("#"+id).trigger("chosen:updated");
    }

    /* ----- PEAK LIST START ----- */

    function buildPeakList(htmls) {
        /* add a container for the list of peaks to htmls */
        htmls.push("<div id='tpPeakListTitle'>Peaks found</div>");

        htmls.push("<div id='tpPeakList' style='height: 30%'>");
            htmls.push("<span id='noPeaks'>No active genes or ranges</span>");
        htmls.push("</div>");
        htmls.push("<div id='tpPeakListSelector'>");
        //htmls.push("<input id='tpPeakListAuto' style='margin-right: 3px' type='checkbox' checked>");
        htmls.push("<div id='tpPeakListButtons' style='margin-left: 4px'>");
        htmls.push('<button title="Select all peaks in the list above" id="tpPeakListAll">All</button>');
        htmls.push('<button title="Select no peaks in the list above" id="tpPeakListNone">None</button>');
        htmls.push('<button title="Select only peaks within a certain distance upstream from the TSS. Click on the field to edit the distance. Click the button to select the peaks." id="tpPeakListUpstream">');
        htmls.push('<input id="tpPeakListAutoDist" type="text" value="2" style="width:2em;border: 0;height: 0.8em;">');
        htmls.push('kbp upstream</button>');
        htmls.push("<div id='tpPeakListButtons'>");

        //htmls.push("<label for='tpPeakListAuto' style='display:inline; font-weight:normal'>Peaks ");
            //htmls.push("<input id='tpPeakListAutoDist' style='width:2em; height:1.3em; margin-right: 0.3em' type='text' value='2'>");
            //htmls.push("kbp upstream</input>");
            //htmls.push("<button id='tpPeakListUpstream' style='float:right; margin-top: 2px'>Select</button>");
        //htmls.push("</label>");
        htmls.push("</div>");
    }

    function peakListShowTitle(sym, chrom, start, end) {
        /* update the peak list box title */
        var el = getById("tpPeakListTitle");
        if (sym)
            el.innerHTML = "<b>"+sym+"</b>, TSS at <span id='tpTss'>"+chrom+":"+start+"</span>";
        else
            el.innerHTML = "Peaks at <b>"+chrom+":"+prettySeqDist(start)+"-"+prettySeqDist(end)+"</b>";
    }

    function peakListSetStatus(str) {
        /* given a string like "chr1|1000|2000+chr2|3000|4000" set the corresponding checkboxes */
        if (!db.conf.atacSearch)
            return;
        let activePeaks = str.split("+");
        let inEls = document.querySelectorAll(".tpPeak > input");
        for (let el of inEls) {
            let parts = el.id.split(":");
            let peakId = parts[1]+"|"+parts[2]+"|"+parts[3];
            el.checked = (activeRanges.includes(peakId));
        }
    }

    function peakListGetPeaksWith(status) {
        /* return array of objects , e.g. [ {chrom:"chr1", start:1000, end:2000, dist:-70000, el:<domObject>} ]
         * status can be "on" or "off" = will only return ranges that are checked or unchecked.
         * */
        let ranges = [];
        let inEls = document.querySelectorAll(".tpPeak > input");
        for (let el of inEls) {
            if (status==="on" && !el.checked)
                continue;
            else if (status==="off" && el.checked)
                continue;
            let parts = el.id.split(":");
            ranges.push({
                "chrom" : parts[1],
                start : parseInt(parts[2]),
                end : parseInt(parts[3]),
                dist : parseInt(parts[4]),
                el : el,
                locusName : parts[1]+"|"+parts[2]+"|"+parts[3]
            });
        }
        return ranges;
    }

    function peakListSerialize() {
        /* return a summary of all currently selected peaks, e.g. "chr1|1000|2000+chr2|3000|4000" */
        let ranges = [];
        for (let r of peakListGetPeaksWith("on")) {
            let locusId = r.chrom+"|"+r.start+"|"+r.end;
            ranges.push(locusId);
        }
        return ranges.join("+");
    }

    function onPeakChange(ev) {
        /* user checks or unchecks a peak */
        let el = ev.currentTarget.firstChild; // user may have clicked the label
        var isChecked = el.checked;
        var peakInfos = el.id.split(":");
        let chrom = peakInfos[1];
        let start = peakInfos[2];
        let end = peakInfos[3];
        let prefix = "+";
        if (!isChecked)
            prefix = "-";
        let rangeStr = prefix+chrom+"|"+start+"|"+end;
        colorByLocus(rangeStr);
    }

    function peakListShowRanges(chrom, foundRanges, searchStart) {
        /* load a list of ranges (arrays of [start, end] into the peak list box */
        var htmls = [];
        var i = 0;
        for (let rangeInfo of foundRanges) {
            let foundStart = rangeInfo[0];
            let foundEnd = rangeInfo[1];
            //let label = chrom+":"+foundStart+"-"+foundEnd;
            let dist = foundStart-searchStart;
            let label = prettySeqDist(dist, true);
            let regLen = foundEnd-foundStart;
            if (regLen!==0)
                label += ", "+(foundEnd-foundStart)+" bp long";
            let checkBoxId = "range:"+chrom+":"+foundStart+":"+foundEnd+":"+dist;
            htmls.push("<div class='tpPeak'>");
            htmls.push("<input style='margin-right: 4px' id='"+checkBoxId+"' type='checkbox'>");
            htmls.push("<label for='"+checkBoxId+"'>"+label+"</label>");
            htmls.push("</div>");
            i++;
        }
        var divEl = document.getElementById("tpPeakList");
        divEl.innerHTML = htmls.join(""); // set the DIV
        classAddListener("tpPeak", "input", onPeakChange);
    }

    function onPeakAll(ev) {
        /* select all peaks */
        let peaks = peakListGetPeaksWith("off");
        let peakNames = [];
        if (peaks.length>100) {
            alert("More than 100 peaks to select. This will take too long. Please contact us if you need this feature.");
            return;
        }
        for (let p of peaks) {
            peakNames.push(p.locusName);
            p.el.checked = true;
        }
        if (peakNames.length===0)
            return;

        let locusStr = "+"+peakNames.join("+");
        colorByLocus(locusStr);
    }

    function onPeakNone(ev) {
        /* unselect all peaks */
        let peaks = peakListGetPeaksWith("on");
        let peakNames = [];
        for (let p of peaks) {
            peakNames.push(p.locusName);
            p.el.checked = false;
        }
        if (peakNames.length===0)
            return;
        let locusStr = "-"+peakNames.join("-");
        colorByLocus(locusStr);
    }

    function onPeakUpstream(ev) {
        /* select all peaks that are closer than the distance in #tpPeakListAutoDist */
        if (ev.target.id==="tpPeakListAutoDist")
            // user actually clicked the input box = do nothing
            return;
        let maxDistStr = document.getElementById("tpPeakListAutoDist").value;
        let maxDist = parseInt(maxDistStr);
        if (isNaN(maxDist)) {
            alert(maxDistStr+" is not a number");
            return;
        }
        if (maxDist > 2000) {
            alert("The distance filter is too high: "+maxDistStr+". It cannot be larger than 2000, = 2Mbp. If you think this is too restrictive, please contact us.");
            return;
        }
        maxDist = -1*maxDist*1000; // needs to be negative

        let addPeaks = [];
        let peaks = peakListGetPeaksWith("off");
        for (let p of peaks) {
            let dist = p.dist;
            if (dist < 0 && dist > maxDist) {
                addPeaks.push(p.chrom+"|"+p.start+"|"+p.end);
                p.el.checked = true;
            }
        }

        if (addPeaks.length===0) {
            alert("Either there is no active gene search or TSS or no peaks are at "+maxDist+" bp relative to the TSS");
            return;
        }

        let loadStr = addPeaks.join("+")

        colorByLocus(loadStr);
        ev.stopPropagation();
    }

    /* ----- PEAK LIST END ----- */

    function selectizeSendGenes(arr, callback) {
        /* given an array of strings s, return an array of objects with id=s and call callback with it.*/
        let foundArr = [];
        for (let o of arr) {
            var text = o.sym;
            if (o.sym!==o.id)
                text = o.sym + " (" + o.id+")";
            foundArr.push( {"id": o.id, "text": text} );
        }
        callback(foundArr);
    }

    function comboLoadGene(query, callback) {
        /* The load() function for selectize for genes.
         * called when the user types something into the gene box, returns matching gene symbols */
        if (!query.length)
            return callback();
        this.clearOptions();
        var genes = db.findGenes(query);
        selectizeSendGenes(genes, callback);
    }

    function updatePeakListWithGene(geneId) {
        /* update the peak list box with all peaks close to a gene's TSS */
        var peaksInView = db.findRangesByGene(geneId);
        var gene = db.getGeneInfoAtac(geneId);
        peakListShowTitle(gene.sym, gene.chrom, gene.chromStart);
        peakListShowRanges(gene.chrom, peaksInView.ranges, gene.chromStart);
    }

    function comboLoadAtac(query, callback) {
        /* The load() function for selectize for ATAC datasets.
         * called when the user types something into the gene box, calls callback with matching gene symbols or peaks
         * or shows the peaks in the peakList box */
        if (!query.length)
            return callback();

        this.clearOptions();
        this.renderCache = {};

        var range = cbUtil.parseRange(query);
        if (range===null) {
            if (!db.geneToTss || db.geneToTss===undefined)
                db.indexGenesAtac();

            var geneInfos = db.findGenesAtac(query);

            if (geneInfos.length === 0)
                return;

            if (geneInfos.length > 1)
                selectizeSendGenes(geneInfos, callback);
            else {
                var geneId = geneInfos[0].id;
                updatePeakListWithGene(geneId);
            }
        } else {
            // user entered a range e.g. chr1:0-190k or chr1:1m-2m
            let searchStart = range.start;
            let searchEnd = range.end;
            peakListShowTitle(null, range.chrom, range.start, range.end);
            let foundRanges = db.findOffsetsWithinPos(range.chrom, searchStart, searchEnd);
            peakListShowRanges(range.chrom, foundRanges, searchStart);
        }
    }

    function comboRender(item, escape) {
        if (item.dist)
            return '<div style="display:block">'+escape(item.text)+'<div style="text-color: darkgrey; font-size:80%;float:right">+'+prettyNumber(item.dist, "bp")+'</div>';
        else
            return '<div>'+escape(item.text)+'</div>';
    }

    function arrayBufferToString(buf, callback) {
        /* https://stackoverflow.com/questions/8936984/uint8array-to-string-in-javascript */
         var bb = new Blob([new Uint8Array(buf)]);
         var f = new FileReader();
         f.onload = function(e) {
           callback(e.target.result);
         };
         f.readAsText(bb);
    }


    function makeXenaUrl(metaFieldName, geneSyms, geneSym, actMeta) {
        /* return URL to Xena view with this dataset and geneSyms loaded */
        var xenaId = db.conf.xenaId;
        var phenoId = db.conf.xenaPhenoId;
        var browser = 'https://singlecell.xenabrowser.net/';
        var view = [{
                name: phenoId,
                host: 'https://singlecellnew.xenahubs.net',
                fields: metaFieldName,
                columnLabel : 'Cell Annotations',
                fieldLabel: metaFieldName
            }];

        if (geneSyms.length!==0)
            view.push({
                name: xenaId,
                host: 'https://singlecellnew.xenahubs.net',
                fields: geneSyms.join(" "),
                width : 15*geneSyms.length,
                columnLabel : 'Dataset genes',
                fieldLabel: geneSyms.join(" ")
            });

        if (actMeta)
            view.push({
                name: phenoId,
                host: 'https://singlecellnew.xenahubs.net',
                fields: actMeta,
                width : 120,
                columnLabel : 'Current Meta',
                fieldLabel: actMeta
            });

        if (geneSym)
            view.push({
                name: xenaId,
                host: 'https://singlecellnew.xenahubs.net',
                fields: geneSym,
                width : 120,
                columnLabel : 'Current Gene',
                fieldLabel: geneSym
            });

        var heatmap = { "showWelcome" : false };

        var url = browser + 'heatmap/?columns=' +
            encodeURIComponent(JSON.stringify(view)) +
            '&heatmap=' + encodeURIComponent(JSON.stringify(heatmap));
        return url;
    }

    function makeHubUrl(geneSym) {
        /* return URL of the hub.txt file, possibly jumping to a given gene  */
            var hubUrl = db.conf.hubUrl;

            if (hubUrl===undefined)
                return null;

            // we accept full session links in the hubUrl statement and just pass these through
            if (hubUrl && hubUrl.indexOf("genome.ucsc.edu/s/")!==-1)
                return hubUrl;

            var ucscDb = db.conf.ucscDb;
            if (ucscDb===undefined) {
                alert("Internal error: ucscDb is not defined in cellbrowser.conf. Example values: hg19, hg38, mm10, etc. You have to set this variable to make track hubs work.");
                return "";
            }

            var fullUrl = null;
            if (hubUrl.indexOf("/")===-1) {
                // no slash -> it's not a URL at all but just a track name (e.g. "tabulamuris")
                var trackName = hubUrl;
                fullUrl = "https://genome.ucsc.edu/cgi-bin/hgTracks?"+trackName+"=full&genome="+ucscDb;
            } else {
                // it's a url to a hub.txt file: either relative or absolute
                if (!hubUrl.startsWith("http"))
                    // relative URL to a hub.txt file -> make absolute now
                    hubUrl = getBaseUrl()+db.name+"/"+hubUrl
                // URL is an absolute link to a hub.txt URL
                fullUrl = "https://genome.ucsc.edu/cgi-bin/hgTracks?hubUrl="+hubUrl+"&genome="+ucscDb;
            }

            if (geneSym)
                fullUrl += "&position="+geneSym+"&singleSearch=knownCanonical";

            return fullUrl;
    }

    function onGenomeButtonClick(ev) {
        /* run when the user clicks the 'genome browser' button */
        let actSym = null;
        if (gLegend.type==="expr")
            actSym = gLegend.geneSym;
        var fullUrl = makeHubUrl(actSym);
        db.gbWin = window.open(fullUrl, 'gbTab');
        return false;
    }

    function onXenaButtonClick(ev) {
        /* run when the user clicks the xena button */
        var geneInfos = db.conf.quickGenes;
        var syms  = []
        if (geneInfos!=undefined) {
            for (var i = 0; i < geneInfos.length; i++)
                syms.push( geneInfos[i][0] ); // make array of symbols
        }

        var actSym = null;
        var actMeta = null;
        if (gLegend.type==="expr")
            actSym = gLegend.geneSym;
        else
            actMeta = gLegend.metaInfo.label;

        if (syms.length===0 && actSym===null) {
            alert("Sorry, the view is not colored by a gene and there are no 'Dataset genes' "+
                   " defined, so there are no genes active "+
                "that could be shown on the heatmap. Please color by a gene first, "+
                "then click the button again.");
            return;
        }

        var fullUrl = makeXenaUrl(db.conf.labelField, syms, actSym, actMeta);
        //if (!db.xenaWin)
        db.xenaWin = window.open(fullUrl, 'xenaTab');
        //else
            //db.xenaWin.location.href = fullUrl;
        return false;
    }

    function openCurrentDataset() {
            /* open dataset dialog with current dataset highlighted */
            $(this).blur();  // remove focus = tooltip disappears
            var parentNames = db.name.split("/");
            parentNames.pop();
            var newPath = cbUtil.joinPaths([parentNames.join("/"), "dataset.json"]);
            cbUtil.loadJson(newPath, function(parentConf) {
                openDatasetDialog(parentConf, db.name);
            });
    }

    function buildToolBar (coordInfo, dataset, fromLeft, fromTop) {
    /* add the tool bar with icons of tools and add under body to the DOM */
        $("#tpToolBar").remove();

        var htmls = [];

        htmls.push("<div id='tpToolBar' style='position:absolute;left:"+fromLeft+"px;top:"+fromTop+"px'>");
        htmls.push('<button title="More info about this dataset: abstract, methods, data download, etc." id="tpButtonInfo" type="button" class="ui-button tpIconButton" data-placement="bottom">Info &amp; Download</button>');

        if (!getVar("suppressOpenButton", false))
            htmls.push('<button id="tpOpenDatasetButton" class="gradientBackground ui-button ui-widget ui-corner-all" style="margin-top:3px; height: 24px; border-radius:3px; padding-top:3px" title="Open another dataset" data-placement="bottom">Open...</button>');

        var nextLeft = 220;
        if (db.conf.hubUrl!==undefined) {
            htmls.push('<a target=_blank href="#" id="tpOpenGenome" class="gradientBackground ui-button ui-widget ui-corner-all" style="margin-left: 10px; margin-top:3px; height: 24px; border-radius:3px; padding-top:3px" title="Show sequencing read coverage and gene expression on UCSC Genome Browser" data-placement="bottom">Genome Browser</a>');
            nextLeft += 155;
        }

        var xenaId = db.conf.xenaId;
        if (xenaId!==undefined) {
            htmls.push('<a target=_blank href="#" id="tpOpenXena" class="gradientBackground ui-button ui-widget ui-corner-all" style="margin-left: 10px; margin-top:3px; height: 24px; border-radius:3px; padding-top:3px" title="Show gene expression heatmap on UCSC Xena Browser, creates heatmap of current gene (if coloring by gene) and all dataset genes. Click this button also if you have an active Xena window open and want to update the view there." data-placement="bottom">Xena</a>');
            nextLeft += 80;
        }

        if (coordInfo[coordInfo.length-1].shortLabel.length > 20)
            //$('.chosen-drop').css({"width": "300px"});
            layoutComboWidth += 50
        buildLayoutCombo(dataset.coordLabel, htmls, coordInfo, "tpLayoutCombo", layoutComboWidth, nextLeft, 2);
        nextLeft += 65+layoutComboWidth;

        var nameParts = dataset.name.split("/");
        var parentName = null;
        if (nameParts.length > 1) {
            buildCollectionCombo(htmls, "tpCollectionCombo", 330, nextLeft, 0);
            nameParts.pop();
            parentName = nameParts.join("/");
        }

        htmls.push("</div>");

        $(document.body).append(htmls.join(""));

        $('#tpOpenXena').click(onXenaButtonClick);
        $('#tpOpenGenome').click(onGenomeButtonClick);

        activateTooltip('.tpIconButton');
        activateTooltip('#tpOpenUcsc');
        activateTooltip('#tpOpenDatasetButton');

        $('#tpButtonInfo').click( function() { openDatasetDialog(db.conf, db.name) } );

        activateCombobox("tpLayoutCombo", layoutComboWidth);

        if (parentName!==null) {
            activateCombobox("tpCollectionCombo", collectionComboWidth);
            loadCollectionInfo( parentName, function(dataset) {
                updateCollectionCombo("tpCollectionCombo", dataset);
            });
        }

        // selective gene or ATAC Color by search box
        var comboLoad = comboLoadGene;
        if (db.conf.atacSearch) {
            comboLoad = comboLoadAtac;
            getById("tpPeakListUpstream").addEventListener("click", onPeakUpstream);
            getById("tpPeakListAll").addEventListener("click", onPeakAll);
            getById("tpPeakListNone").addEventListener("click", onPeakNone);
            activateTooltip("#tpPeakListButtons > button");

        }

        // This is a hack to deactivate the "sifter" functionality of selectize.
        // It turns out that selectize is a very bad dropdown choice for us,
        // as it makes a strong assumption that matching is done on a keyword basis
        // which is not true for chrom ranges.
        // https://gist.github.com/rhyzx/2281e8d1662b7be21716
        Selectize.prototype.search = function (query) {
            return {
                query: query,
                tokens: [], // disable highlight
                items: $.map(this.options, function (item, key) {
                    return {id: key}
                })
            }
        };

        var select = $('#tpGeneCombo').selectize({
                maxItems: 1,
                valueField : 'id',
                labelField : 'text',
                searchField : 'text',
                closeAfterSelect: true,
                load : comboLoad,
                render : {option : comboRender }
        });

        select.on("change", onGeneComboChange);

        $('#tpCollectionCombo').change(onDatasetChange);
        // update the combobox, select the right dataset
        $('#tpLayoutCombo').change(onLayoutChange);
        $('#tpOpenDatasetButton').click(openCurrentDataset);
    }

    function metaFieldToLabel(fieldName) {
    /* convert the internal meta field string to the label shown in the UI. Fix underscores, _id, etc */
        if (fieldName==="_id")
            fieldName = capitalize(gSampleDesc)+" identifier";
        else
            fieldName = fieldName.replace(/_/g, " ");
        return fieldName;
    }

    function buildMetaPanel(htmls) {
        /* add html strings for the meta panel to the left side bar */
        var metaFields = db.conf.metaFields;
        for (var i = 0; i < metaFields.length; i++) {
            var metaInfo = metaFields[i];
            var fieldLabel = metaInfo.label;
            fieldLabel = fieldLabel.replace(/_/g, " ");
            var fieldMouseOver = metaInfo.desc;

            // fields without binning and with too many unique values are greyed out
            var isGrey = (metaInfo==="enum" && metaInfo.diffValCount>MAXCOLORCOUNT && metaInfo.binMethod===undefined);

            var addClass = "";
            var addTitle="";
            htmls.push("<div class='tpMetaBox' data-field-name='"+metaInfo.name+"' id='tpMetaBox_"+i+"'>");
            if (isGrey) {
                addClass=" tpMetaLabelGrey";
                addTitle=" title='This field contains too many different values. You cannot click it to color on it.'";
            }

            let divId = "tpMetaLabel_"+i;

            htmls.push("<div id='"+divId+"' class='tpMetaLabel"+addClass+"'"+addTitle+">"+fieldLabel);
            if (fieldMouseOver)
                htmls.push('<i title="'+fieldMouseOver+'" '+
                ' class="material-icons tpMetaLabelHelp" style="float:right;font-size:16px">help_outline</i>');
            htmls.push("</div>");

            var styleAdd="";
            if (metaInfo.opt!==undefined) {
                var opt = metaInfo.opt;
                if (opt.fontSize!==undefined)
                    styleAdd = ";font-size:"+metaInfo.opt.fontSize;
            }

            htmls.push("<div class='tpMetaValue' style='width:"+(metaBarWidth-2*metaBarMargin)+"px"+styleAdd+
                "' data-field-name='"+metaInfo.name+"' id='tpMeta_" + i + "'>&nbsp;</div>");
            htmls.push("</div>"); // tpMetaBox
        }
        htmls.push("<div style='background-color:white; float:right' id='tpMetaNote' style='display:none; height:1em'></div>");
    }

    function rebuildMetaPanel() {
        $("#tpMetaPanel").empty();
        let htmls = [];
        buildMetaPanel(htmls);
        $("#tpMetaPanel").html(htmls.join(""));
        connectMetaPanel();
    }

    function connectMetaPanel() {
        activateTooltip(".tpMetaLabelHelp");
        $(".tpMetaLabel").click( onMetaClick );
        $(".tpMetaValue").click( onMetaClick );
        $(".tpMetaValue").mouseover( onMetaMouseOver );
        $(".tpMetaValue").mouseleave ( function() {
            $('#tpMetaTip').hide();
            $('.tpMetaBox').removeClass("tpMetaHover");
            $('.tpMetaBox .tpMetaValue').removeClass("tpMetaHover");
        } );

        // setup the right-click menu
        //var menuItems = [{name: "Use as cluster label"},{name: "Copy field value to clipboard"}];
        var menuItems = [{name: "Copy field value to clipboard"}];
        var menuOpt = {
            selector: ".tpMetaBox",
            items: menuItems,
            className: 'contextmenu-customwidth',
            callback: onMetaRightClick
        };
        $.contextMenu( menuOpt );

        var menuItemsCust = [{name: "Copy field value to clipboard"},
            {name: "Remove custom annotations"}];
        var menuOptCust = {
            selector: "#tpMetaBox_custom",
            items: menuItemsCust,
            className: 'contextmenu-customwidth',
            callback: onMetaRightClick
        };
        $.contextMenu( menuOptCust );
        // setup the tooltips
        //$('[title!=""]').tooltip();
    }

    function buildLeftSidebar () {
    /* add the left sidebar with the meta data fields. db.loadConf
     * must have completed before this can be run, we need the meta field info. */
        $("#tpLeftSidebar").remove();
        // setup the tabs
        var tabsWidth = metaBarWidth;

        var htmls = [];
        htmls.push("<div id='tpMetaTip' style='display:none'></div>");
        htmls.push("<div id='tpLeftSidebar' style='position:absolute;left:0px;top:"+menuBarHeight+"px;width:"+metaBarWidth+"px'>");

        htmls.push("<div class='tpSidebarHeader'>Color By</div>");

        // a bar with the tabs
        htmls.push("<div id='tpLeftTabs'>");
        htmls.push("<ul>");
        htmls.push("<li><a href='#tpAnnotTab'>Annotation</a></li>");
        htmls.push("<li><a href='#tpGeneTab'>"+getGeneLabel()+"</a></li>");
        htmls.push("</ul>");

        htmls.push("<div id='tpAnnotTab'>");
        htmls.push('<label style="padding-left: 2px; margin-bottom:8px; padding-top:8px" for="'+"tpMetaCombo"+'">Color by Annotation</label>');
        buildMetaFieldCombo(htmls, "tpMetaComboBox", "tpMetaCombo", 0);
        htmls.push('<label style="padding-left: 2px; margin-bottom:8px; padding-top:8px" for="tpLabelCombo">Label by Annotation</label>');
        buildMetaFieldCombo(htmls, "tpLabelComboBox", "tpLabelCombo", 0, db.conf.labelField);
        htmls.push('<div style="padding-top:4px; padding-bottom: 4px; padding-left:2px" id="tpHoverHint" class="tpHint">Hover over a '+gSampleDesc+' to update data below</div>');
        htmls.push('<div style="padding-top:4px; padding-bottom: 4px; padding-left:2px; display: none" id="tpSelectHint" class="tpHint">Cells selected. No update on hover.</div>');

        htmls.push("<div id='tpMetaPanel'>");
        buildMetaPanel(htmls);
        htmls.push("</div>"); // tpMetaPanel

        htmls.push("</div>"); // tpAnnotTab

        htmls.push("<div id='tpGeneTab'>");

        buildGeneCombo(htmls, "tpGeneCombo", 0, metaBarWidth-10);

        if (db.conf.atacSearch)
            buildPeakList(htmls);

        var geneLabel = getGeneLabel();
        buildGeneTable(htmls, "tpRecentGenes", "Recent "+geneLabel+"s", "Hover or select cells to update colors", gRecentGenes);

        //var myGenes = loadMyGenes();

        var noteStr = "No genes or peaks defined. Use the setting quickGenesFile in cellbrowser.conf to add a file with gene symbols or peaks that will be shown here";
        buildGeneTable(htmls, "tpGenes", "Dataset "+geneLabel+"s", null, db.conf.quickGenes, noteStr);

        htmls.push("</div>"); // tpGeneTab

        htmls.push("</div>"); // tpLeftSidebar

        $(document.body).append(htmls.join(""));

        $("#tpLeftTabs").tabs();
        $('#tpLeftTabs').tabs("option", "active", 0); // open the first tab

        $('.tpGeneBarCell').click( onGeneClick );
        $('#tpChangeGenes').click( onChangeGenesClick );

        activateCombobox("tpMetaCombo", metaBarWidth-10);
        $("#tpMetaCombo").change( onMetaComboChange );

        activateCombobox("tpLabelCombo", metaBarWidth-10);
        $("#tpLabelCombo").change( onLabelComboChange );
        connectMetaPanel();
    }

    function makeTooltipCont() {
        /* make a div for the tooltips */
        var ttDiv = document.createElement('div');
        ttDiv.id = "tpTooltip";
        ttDiv.style.position = "absolute";
        //ttDiv.style.left = left+"px";
        //ttDiv.style.top = top+"px";
        ttDiv.style["padding"]="2px";
        ttDiv.style["border"]="1px solid black";
        ttDiv.style["border-radius"]="2px";
        ttDiv.style["display"]="none";
        ttDiv.style["cursor"]="pointer";
        ttDiv.style["background-color"]="white";
        ttDiv.style["box-shadow"]="0px 2px 4px rgba(0,0,0,0.3)";
        ttDiv.style["user-select"]="none";
        ttDiv.style["z-index"]="10";
        return ttDiv;
    }

    function showIntro(addFirst) {
        /* add the intro.js data */
        //var htmls = [];
        //htmls.push("<a href='http://example.com/' data-intro='Hello step one!'></a>");
        //$(document.body).append(htmls.join(""));

        localStorage.setItem("introShown", "true");
        activateTab("meta");
        var intro = introJs();
        intro.setOption("hintAnimation", false);
        intro.setOption("exitOnEsc", true);
        intro.setOption("exitOnOverlayClick", true);
        intro.setOption("scrollToElement", false);

        intro.setOption("doneLabel", "Close this window");
        intro.setOption("skipLabel", "Stop the tutorial");

        if (addFirst) {
            intro.setOption("skipLabel", "I know. Close this window.");
            intro.addStep({
                element: document.querySelector('#tpHelpButton'),
                intro: "Are you here for the first time and wondering what this is?<br>The tutorial takes only 1 minute. To skip the tutorial now, click 'I know' below or press Esc.<br>You can always show it again by clicking 'Help > Tutorial'.",
              });
        }

        intro.addSteps(
            [
              {
                intro: "In the middle of the screen, each circle represents a "+gSampleDesc+". You can click the cluster label text to show the marker gene lists of the cluster.",
                element: document.querySelector('#tpCanvas'),
                position: 'auto'
              },
              {
                element: document.querySelector('#tpLeftSidebar'),
                intro: "Info and color control: move the mouse over a circle to show its annotation data.<br>Pick an annotation field or a gene to color on it.<br>",
                position: 'auto'
              },
              //{
                //element: document.querySelector('#tpGeneBar'),
                //intro: "Expression data: when you move the mouse, expression values will be shown here.<br>Click on a gene to color the circles by gene expression level (log'ed).",
                //position: 'top'
              //},
              {
                element: document.querySelector('#tpLegendBar'),
                intro: "Click into the legend to select "+gSampleDesc+"s.<br>Click a color to change it or select a palette from the 'Colors' menu.<br>To setup your own cell browser, see 'Help - Setup your own'",
                position: 'left'
              },
            ]);
        intro.start();
    }

     /**
     https://gist.github.com/mjackson/5311256
     * Converts an HSL color value to RGB. Conversion formula
     * adapted from http://en.wikipedia.org/wiki/HSL_color_space.
     * Assumes h, s, and l are contained in the set [0, 1] and
     * returns r, g, and b in the set [0, 255].
     *
     * @param   Number  h       The hue
     * @param   Number  s       The saturation
     * @param   Number  l       The lightness
     * @return  Array           The RGB representation
     */
    function hue2rgb(p, q, t) {
      if (t < 0) t += 1;
      if (t > 1) t -= 1;
      if (t < 1/6) return p + (q - p) * 6 * t;
      if (t < 1/2) return q;
      if (t < 2/3) return p + (q - p) * (2/3 - t) * 6;
      return p;
    }

    function hslToRgb(h, s, l) {
      var r, g, b;

      if (s === 0) {
	r = g = b = l; // achromatic
      } else {
	var q = l < 0.5 ? l * (1 + s) : l + s - l * s;
	var p = 2 * l - q;

	r = hue2rgb(p, q, h + 1/3);
	g = hue2rgb(p, q, h);
	b = hue2rgb(p, q, h - 1/3);
      }

      return [ r * 255, g * 255, b * 255 ];
    }

    function makeHslPalette(hue, n) {
        /* return a list of n hexcodes from hue to white */
        var pal = [];
        for (var i=1; i<n+1; i++) {
            var c = hslToRgb(hue, 1.0, (0.35+((n-i)/n*0.65)));
            pal.push( palette.rgbColor(c[0]/255, c[1]/255, c[2]/255));
        }
        return pal;
    }

    function isDark(c) {
	/* c is a six-digit hexcode, return true if it's a dark color */
	// from https://stackoverflow.com/questions/12043187/how-to-check-if-hex-color-is-too-black
	var rgb = parseInt(c, 16);   // convert rrggbb to decimal
	var r = (rgb >> 16) & 0xff;  // extract red
	var g = (rgb >>  8) & 0xff;  // extract green
	var b = (rgb >>  0) & 0xff;  // extract blue

	var luma = 0.2126 * r + 0.7152 * g + 0.0722 * b; // per ITU-R BT.709

	return (luma < 40);
    }

    function makePercPalette(palName, n) {
        /* palettes from https://github.com/politiken-journalism/scale-color-perceptual */
        var pal = [];
        var step = 1/n;

        var func = null;
        switch (palName) {
            case 'inferno' : func =  scale.color.perceptual.inferno; break;
            case 'viridis' : func =  scale.color.perceptual.viridis; break;
            case 'magma' : func =  scale.color.perceptual.magma; break;
            case 'plasma' : func =  scale.color.perceptual.plasma; break;
        }

        for (let x=0; x<n; x++) {
            pal.push(func(x*step).substr(1));
        }

        if (pal.length!==n)
            console.log("palette is too small");
        return pal;
    }

    function makeColorPalette(palName, n) {
    /* return an array with n color hex strings */
    /* Use Google's palette functions for now, first Paul Tol's colors, if that fails, use the usual HSV rainbow
     * This code understands our special palette, tol-sq-blue
    */
        var pal = [];
        if (palName==="blues")
            pal = makeHslPalette(0.6, n);
        else if (palName==="magma" || palName==="viridis" || palName==="inferno" || palName=="plasma")
            pal = makePercPalette(palName, n);
        else if (palName==="reds")
            pal = makeHslPalette(0.0, n);
        else {
            if (n===2)
                pal = ["ADD8E6","FF0000"];
            else {
                var realPalName = palName.replace("tol-sq-blue", "tol-sq");
                pal = palette(realPalName, n);
                if (palName==="tol-sq-blue")
                    pal[0]='f4f7ff';
            }
        }

        return pal;
    }

    function colorByCluster() {
    /* called when meta and coordinates have been loaded: scale data and color by meta field  */
        //setZoomRange();
    }

    function loadClusterTsv(fullUrl, func, divName, clusterName) {
    /* load a tsv file relative to baseUrl and call a function when done */
        function conversionDone(data) {
            Papa.parse(data, {
                    complete: function(results, localFile) {
                                func(results, localFile, divName, clusterName);
                            },
                    error: function(err, file) {
                                if (divName!==undefined)
                                    alert("could not load "+fullUrl);
                            }
                    });
        }

        function onTsvLoadDone(res) {
            var data = res.target.response;
            if (res.target.responseURL.endsWith(".gz")) {
                data = pako.ungzip(data);
                //data = String.fromCharCode.apply(null, data); // only good for short strings
                data = arrayBufferToString(data, conversionDone);
            }
            else
                conversionDone(data);
        }

    var req = new XMLHttpRequest();
    req.addEventListener("load", onTsvLoadDone);
    req.open('GET', fullUrl, true);
    req.responseType = "arraybuffer";
    req.send();
    }

    function removeFocus() {
    /* called when the escape key is pressed, removes current focus and puts focus to nothing */
        window.focus();
        if (document.activeElement) {
                document.activeElement.blur();
        }

    }

    function setupKeyboard() {
    /* bind the keyboard shortcut keys */
        phoneHome();
        Mousetrap.bind('o', openCurrentDataset);
        Mousetrap.bind('c m', onMarkClearClick);
        Mousetrap.bind('h m', onMarkClick);

        Mousetrap.bind('space', onZoom100Click);

        Mousetrap.bind('z', function() { activateMode("zoom"); });
        Mousetrap.bind('m', function() { activateMode("move"); });
        Mousetrap.bind('s', function() { activateMode("select"); });

        Mousetrap.bind('-', onZoomOutClick);
        Mousetrap.bind('+', onZoomInClick);
        Mousetrap.bind('s n', onSelectNoneClick);
        Mousetrap.bind('s a', onSelectAllClick);
        Mousetrap.bind('s i', onSelectInvertClick);
        Mousetrap.bind('s s', onSelectNameClick);

        Mousetrap.bind('b s', onBackgroudSetClick);
        Mousetrap.bind('b r', onBackgroudResetClick);

        Mousetrap.bind('m', function() {$('#tpMetaCombo').trigger("chosen:open"); return false;});
        Mousetrap.bind('d', function() {$('#tpDatasetCombo').trigger("chosen:open"); return false;});
        //Mousetrap.bind('l', function() {$('#tpLayoutCombo').trigger("chosen:open"); return false;});
        Mousetrap.bind('g', function() {$("#tpGeneCombo").selectize()[0].selectize.focus(); return false;});
        Mousetrap.bind('c l', onHideShowLabelsClick );
        Mousetrap.bind('f c', onFindCellsClick );
        Mousetrap.bind('f i', function() { onSelectByIdClick(); return false; } );
        Mousetrap.bind('t', onSplitClick );
        Mousetrap.bind('h', onHeatClick );

        Mousetrap.bind('up', function() { renderer.movePerc(0, 0.1); renderer.drawDots(); } );
        Mousetrap.bind('left', function() { renderer.movePerc(-0.1, 0); renderer.drawDots(); } );
        Mousetrap.bind('right', function() { renderer.movePerc(0.1, 0); renderer.drawDots(); } );
        Mousetrap.bind('down', function() { renderer.movePerc(0, -0.1); renderer.drawDots(); } );

        // yay vim
        Mousetrap.bind('i', function() { renderer.movePerc(0, 0.1); renderer.drawDots(); } );
        Mousetrap.bind('j', function() { renderer.movePerc(-0.1, 0); renderer.drawDots(); } );
        Mousetrap.bind('l', function() { renderer.movePerc(0.1, 0); renderer.drawDots(); } );
        Mousetrap.bind('k', function() { renderer.movePerc(0, -0.1); renderer.drawDots(); } );

        //Mousetrap.stopCallback = function(e, element, combo) {
            //var doStop = (element.tagName == 'INPUT' || element.tagName == 'SELECT' || element.tagName == 'TEXTAREA' || (element.contentEditable && element.contentEditable == 'true'));
            //console.log(e, element, combo);
            //return doStop;
        //};
    }

    // https://stackoverflow.com/a/33861088/233871
    function isInt(x) {
        return (typeof x==='number') && x % 1 === 0;
    }

    function naturalSort (a, b) {
    /* copied from https://github.com/Bill4Time/javascript-natural-sort/blob/master/naturalSort.js */
    /* "use strict"; */
    var re = /(^([+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?)?$|^0x[0-9a-f]+$|\d+)/gi,
        sre = /(^[ ]*|[ ]*$)/g,
        dre = /(^([\w ]+,?[\w ]+)?[\w ]+,?[\w ]+\d+:\d+(:\d+)?[\w ]?|^\d{1,4}[\/\-]\d{1,4}[\/\-]\d{1,4}|^\w+, \w+ \d+, \d{4})/,
        hre = /^0x[0-9a-f]+$/i,
        ore = /^0/,
        i = function(s) { return naturalSort.insensitive && ('' + s).toLowerCase() || '' + s; },
        // convert all to strings strip whitespace
        x = i(a).replace(sre, '') || '',
        y = i(b).replace(sre, '') || '',
        // chunk/tokenize
        xN = x.replace(re, '\0$1\0').replace(/\0$/,'').replace(/^\0/,'').split('\0'),
        yN = y.replace(re, '\0$1\0').replace(/\0$/,'').replace(/^\0/,'').split('\0'),
        // numeric, hex or date detection
        xD = parseInt(x.match(hre), 16) || (xN.length !== 1 && x.match(dre) && Date.parse(x)),
        yD = parseInt(y.match(hre), 16) || xD && y.match(dre) && Date.parse(y) || null,
        oFxNcL, oFyNcL;
    // first try and sort Hex codes or Dates
    if (yD) {
        if ( xD < yD ) { return -1; }
        else if ( xD > yD ) { return 1; }
    }
    // natural sorting through split numeric strings and default strings
    for(var cLoc=0, numS=Math.max(xN.length, yN.length); cLoc < numS; cLoc++) {
        // find floats not starting with '0', string or 0 if not defined (Clint Priest)
        oFxNcL = !(xN[cLoc] || '').match(ore) && parseFloat(xN[cLoc]) || xN[cLoc] || 0;
        oFyNcL = !(yN[cLoc] || '').match(ore) && parseFloat(yN[cLoc]) || yN[cLoc] || 0;
        // handle numeric vs string comparison - number < string - (Kyle Adams)
        if (isNaN(oFxNcL) !== isNaN(oFyNcL)) { return (isNaN(oFxNcL)) ? 1 : -1; }
        // rely on string comparison if different types - i.e. '02' < 2 != '02' < '2'
        else if (typeof oFxNcL !== typeof oFyNcL) {
            oFxNcL += '';
            oFyNcL += '';
        }
        if (oFxNcL < oFyNcL) { return -1; }
        if (oFxNcL > oFyNcL) { return 1; }
    }
    return 0;
}

    function sortPairsBy(countList, sortBy) {
    /* sort an array in the format [name, count] by either name (using naturalSort) or count */
        var isSortedByName = null;

        if (sortBy==="name") {
            countList.sort(function(a, b) { return naturalSort(a[0], b[0]); });  // I have a feeling that this is very slow...
            isSortedByName = true;
        }
        else if (sortBy=="count") {
            // sort this list by count
            countList = countList.sort(function(a, b){ return b[1] - a[1]; }); // reverse-sort by count
            isSortedByName = false;
        } else {
            isSortedByName = false;
        }

        var ret = {};
        ret.list = countList;
        ret.isSortedByName = isSortedByName;
        // pallette should be a gradient for data types where this makes sense
        return ret;
    }

    function plotSelection(coords) {
    /* redraw block dots of current selection
       Redrawing makes sure that they are not hidden underneath others.
    */
        for (var i = 0; i < coords.length; i++) {
            var x = coords[i][0];
            var y = coords[i][1];
            var fill = coords[i][2];
            var dot = drawCircle(x, y, fill, 0x000000);
            stage.addChild(dot);
            visibleGlyps.push(dot);
        }
    }

    function removeCellIds(coords, cellIds) {
    /* remove all coords that are in an object of cellIds */
        var newCoords = [];
        for (var i = 0; i < coords.length; i++) {
            var coord = coords[i];
            var cellId = coord[0];
            if (!(cellId in cellIds))
                newCoords.push(coord);
        }
        return newCoords;
    }

    function showOnlyCellIds(coords, cellIds) {
    /* keep only coords that are in an object of cellIds */
        var newCoords = [];
        for (var i = 0; i < coords.length; i++) {
            var coord = coords[i];
            var cellId = coord[0];
            if (cellId in cellIds)
                newCoords.push(coord);
        }
        return newCoords;
    }

    function countValues(arr) {
        var counts = {};
        for (var i = 0; i < arr.length; i++) {
                counts[arr[i]] = 1 + (counts[arr[i]] || 0);
        }
        var countArr = Object.entries(counts);
        return countArr
    }

    function makeLabelRenames(metaInfo) {
        /* return an obj with old cluster name -> new cluster name */
        var valCounts = metaInfo.valCounts;
        if (valCounts===undefined) { // 'int' and 'float' types do not have their values counted yet
            // this doesn't work because the values are not loaded yet, requires moving this call to
            // later
            //metaInfo.valCounts = countValues(metaInfo.arr);
            alert("cannot label on numeric fields, please use the enumFields option in cellbrowser.conf");
        }
        var newLabels = metaInfo.ui.shortLabels;

        var oldToNew = {};
        for (var i = 0; i < valCounts.length; i++) {
            var oldLabel = valCounts[i][0];
            var newLabel = newLabels[i];
            oldToNew[oldLabel] = newLabel;
        }
        return oldToNew;
    }

    function getClusterFieldInfo() {
        var clusterField = db.conf.activeLabelField;
        var clusterMetaInfo = db.findMetaInfo(clusterField);
        return clusterMetaInfo;
    }

    function rendererUpdateLabels(metaInfo) {
        /* update the labels in the renderer from the metaInfo data. */
        var oldToNew = makeLabelRenames(metaInfo);

        var oldRendLabels = null;
        if (renderer.origLabels) {
            oldRendLabels = renderer.origLabels;
        }
        else {
            oldRendLabels = renderer.getLabels();
            renderer.origLabels = oldRendLabels;
        }

        var newRendLabels = [];
        for (var i = 0; i < oldRendLabels.length; i++) {
            var oldLabel = oldRendLabels[i];
            var newLabel = oldToNew[oldLabel];
            newRendLabels.push(newLabel);
        }
        renderer.setLabels(newRendLabels);
    }

    function onLegendLabelClick(ev) {
    /* called when user clicks on legend entry. */

        var legendId = parseInt(ev.target.id.split("_")[1]);
        var colorIndex = gLegend.rows[legendId].intKey;
        $("#tpLegendCheckbox_" + colorIndex).click();
    }

    function onSortByClick (ev) {
    /* flip the current legend sorting */
        var sortBy = null;
        if (ev.target.parentElement.id.endsWith("Col1")) // column 1 is the Name
            sortBy = "name"
        else
            sortBy = "freq";

        saveToUrl("s_"+gLegend.metaInfo.name,sortBy, gLegend.defaultSortBy);
        legendSort(sortBy);
        $(".tooltip").hide();
        buildLegendBar();
    }


    function onMetaRightClick (key, options) {
    /* user right-clicks on a meta field */
        var metaName = options.$trigger[0].id.split("_")[1];

        //if (key==0) {
            //gCurrentDataset.labelField = gCurrentDataset.metaFields[metaIdx];
            //gClusterMids = null; // force recalc
            //plotDots();
            //renderer.render(stage);
            //updateMenu();
        //}
        if (key===0) {
            copyToClipboard("#tpMeta_"+metaName);
        }

    }

    //function onLegendRightClick (key, options) {
    /* user right-clicks on a legend label */
        //var selEls = $(".tpLegendSelect")
        //var legendIds = [];
        //for (var i = 0; i < selEls.length; i++) {
            //var selEl = selEls[i];
            //var legendId = parseInt(selEl.id.split("_")[1]);
            //legendIds.push(legendId);
        //}

        //var cellIds = findCellIdsForLegendIds(gClasses, legendIds);
        //var mode;
        //if (key==0) // hide
            //mode = "hide";
        //else
            //mode = "showOnly";
        //filterCoordsAndUpdate(cellIds, mode);
    //}

    function setLegendHeaders(type) {
    /* set the headers of the right-hand legend */
        if (type==="category") {
            $('#tpLegendCol1').html('<span title="select all checkboxes below" id="tpLegendClear">&#9745;</span><span class="tpLegendHover" title="click to sort by name"> Name<span class="caret"></span></span>');
            $('#tpLegendCol2').html('<span class="tpLegendHover" title="click to sort by frequency"> Frequency<span class="caret"></span></span>');
        }
        else {
            $('#tpLegendCol1').html('<span title="select all checkboxes below" id="tpLegendClear">&#9745;</span> Range<span');
            $('#tpLegendCol2').html('Frequency');
        }
        activateTooltip("#tpLegendClear");
        activateTooltip(".tpLegendHover");
    }

    function updateLegendGrandCheckbox() {
        var checkbox = $("#tpLegendClear");
        var total = renderer.getCount();
        var selected = renderer.selCells.size;
        if (gLegend.selectionDirection == "all" && total == selected) {
            gLegend.selectionDirection = "none";
            checkbox.html("&#9746;");
            // from https://stackoverflow.com/questions/9501921/change-twitter-bootstrap-tooltip-content-on-click
            checkbox.attr('title', "unselect all checkboxes below");
            var tip = checkbox.bsTooltip('fixTitle').data('bs.tooltip').$tip;
            if (tip) {
                tip.find('.tooltip-inner')
                .text("unselect all checkboxes below");
            }
        } else if (gLegend.selectionDirection == "none" && selected === 0) {
            gLegend.selectionDirection = "all";
            checkbox.html("&#9745;");
            checkbox.attr('title', "select all checkboxes below");
            var tip = checkbox.bsTooltip('fixTitle').data('bs.tooltip').$tip;
            if (tip) {
                tip.find('.tooltip-inner')
                .text("select all checkboxes below");
            }
        }
    }

    function onLegendClearClick(ev) {
        /* unselect all checkboxes in the legend and clear the selection */
        if (gLegend.selectionDirection == "all") {
            clearSelectionState();
            renderer.selectAll();
            renderer.drawDots();
        } else {
            onSelectNoneClick();
        }
        ev.stopPropagation();
    }

    function onLegendApplyLimitsClick(ev) {
        /* user clicked the apply button: apply limits to the plot and redraw */
        makeLegendExpr(gLegend.geneSym, gLegend.titleHover, binInfo, exprArr, decArr);
    }

    function onLegendCheckboxClick(ev) {
        /* user clicked the small checkboxes in the legend */
        var valIdx = parseInt(ev.target.getAttribute("data-value-index"));
        if ($(this).is(':checked'))
            renderer.selectByColor(valIdx);
        else
            renderer.unselectByColor(valIdx);
        renderer.drawDots();
        ev.stopPropagation();
        $(this).blur();
    }

    function buildMinMaxPart(htmls) {
        /* create the min/max and apply/reset buttons */
        htmls.push("<div>");
        htmls.push("<table style='margin: 4px; margin-top: 6px'>");
        htmls.push("<tr>");
        htmls.push("<td><label for='exprMin'>Min:</label></td>");
        htmls.push("<td><input name='exprMin' size='8' type='text'></td>");
        htmls.push("<td><button id='tpExprLimitApply' style='border-radius: 4px; margin-left: 4px; padding: 3px 6px 3px 6px' class='ui-button'>Apply</button></td>");
        htmls.push("</tr>");

        htmls.push("<tr>");
        htmls.push("<td><label for='exprMax'>Max:</label></td>");
        htmls.push("<td><input name='exprMax' size='8' type='text'></td>");
        htmls.push("<td><button id='tpExprLimitClear' style='border-radius: 4px; margin-left: 4px; padding: 3px 6px 3px 6px' class='ui-button'>Clear</button></td>");
        htmls.push("</tr>");
        htmls.push("</table>");
        htmls.push("</div>");
    }


    function buildLegendBar() {
    /* draws current legend as specified by gLegend.rows
     * */
        if (gLegend.rows===undefined)
            return;

        $('#tpLegendContent').empty();

        var htmls = [];

        var colors = [];
        var rows = gLegend.rows;

        var legTitle = gLegend.title;
        var subTitle = gLegend.subTitle;

        htmls.push('<span id="tpLegendTitle" title="' +gLegend.titleHover+'">'+legTitle+"</span>");
        if (subTitle)
            htmls.push('<div id="tpLegendSubTitle" >'+subTitle+"</div>");
        htmls.push('<div class="tpHint">Check boxes below to select '+gSampleDesc+'s</small></div>');
        htmls.push("</div>"); // title
        htmls.push('<div id="tpLegendHeader"><span id="tpLegendCol1"></span><span id="tpLegendCol2"></span></div>');
        htmls.push('<div id="tpLegendRows">');

        // get the sum of all, to calculate frequency
        var sum = 0;
        for (var i = 0; i < rows.length; i++) {
            let count = rows[i].count;
            sum += count;
        }

        for (i = 0; i < rows.length; i++) {
            var row = rows[i];
            var colorHex = row.color; // manual color
            if (colorHex===null)
                colorHex = row.defColor; // default color

            var label = row.label;
            var longLabel = row.longLabel;

            let count = row.count;
            var valueIndex = row.intKey;
            var freq  = 100*count/sum;

            if (count===0) // never output categories with 0 count.
                continue;

            var labelClass = "tpLegendLabel";
            label = label.replace(/_/g, " ").replace(/'/g, "&#39;").trim();
            if (longLabel)
                longLabel = longLabel.replace(/_/g, " ").trim();

            if (likeEmptyString(label)) {
                labelClass += " tpGrey";
            }
            if (label==="") {
                label = "(empty)";
            }

            colors.push([i, colorHex]); // save for later

            var mouseOver = "";
            // only show the full value on mouse over if the label is long, "" suppresses mouse over
            if (longLabel && longLabel!=label)
                mouseOver = longLabel;
            else {
                if (label.length > 20)
                    mouseOver = label;
            }

            var classStr = "tpLegend";
            var line = "<div id='tpLegend_" +valueIndex+ "' class='" +classStr+ "'>";
            htmls.push(line);
            htmls.push("<input class='tpLegendCheckbox' data-value-index='"+valueIndex+"' "+
                "id='tpLegendCheckbox_"+valueIndex+"' type='checkbox'>");
            htmls.push("<input class='tpColorPicker' id='tpLegendColorPicker_"+i+"' />");

            htmls.push("<span class='"+labelClass+"' id='tpLegendLabel_"+i+"' data-placement='auto top' title='"+mouseOver+"'>");
            htmls.push(label);
            htmls.push("</span>");
            var prec = 1;
            if (freq<1)
                prec = 2;
            htmls.push("<span class='tpLegendCount' title='"+count+" of "+sum+"'>"+freq.toFixed(prec)+"%</span>");
            htmls.push("</span>");

            htmls.push("</div>");
        }
        htmls.push('</div>'); // tpLegendRows

        // add the div where the violin plot will later be shown
        htmls.push("<div id='tpViolin'>");
        htmls.push("<canvas style='height:200px; padding-top: 10px; padding-bottom:30px' id='tpViolinCanvas'></canvas>");
        htmls.push("</div>"); // violin

        var htmlStr = htmls.join("");
        $('#tpLegendContent').append(htmlStr);
        setLegendHeaders(gLegend.rowType);

        // tpLegendContent has to go only up to the bottom of the screen.
        // This is done again in resizeDivs()
        // I have not found a way to do this in CSS...
        $("#tpLegendBar").css("height", (window.innerHeight - $('#tpLegendBar').offset().top)+"px");

        activateTooltip("#tpResetColors");
        activateTooltip("#tpSortBy");

        $("#tpLegendCol1").click( onSortByClick );
        $("#tpLegendCol2").click( onSortByClick );
        $(".tpLegendCheckbox").click( onLegendCheckboxClick );
        $("#tpLegendClear").click( onLegendClearClick );
        $("#tpExprLimitApply").click( onLegendApplyLimitsClick );

        $('.tpLegend').click( onLegendLabelClick );
        //$('.tpLegendLabel').attr( "title", "Click to select samples with this value. Shift click to select multiple values.");
        activateTooltip(".tpLegendLabel");
        activateTooltip(".tpLegendCount");

        // setup the right-click menu
        //var menuItems = [{name: "Hide "+gSampleDesc+"s with this value"}, {name:"Show only "+gSampleDesc+"s with this value"}];
        //var menuOpt = {
            //selector: ".tpLegend",
            //items: menuItems,
            //className: 'contextmenu-customwidth',
            //callback: onLegendRightClick
        //};
        //$.contextMenu( menuOpt );

        // activate the color pickers
        for (let i = 0; i < colors.length; i++) {
            var colInfo = colors[i];
            var rowIdx = colInfo[0];
            var hexCode = colInfo[1];

            var opt = {
                hideAfterPaletteSelect : true,
                color : hexCode,
                showPalette: true,
                allowEmpty : true,
                showInput: true,
                preferredFormat: "hex",
                change: onColorPickerChange
                }
            $("#tpLegendColorPicker_"+rowIdx).spectrum(opt);
        }

        buildViolinPlot();
    }

    function onColorPickerChange(color, ev) {
        /* called when user manually selects a color in the legend with the color picker */
        console.log(ev);
        /* jshint validthis: true */
        var valueIdx = parseInt(this.id.split("_")[1]);
        var rows = gLegend.rows;
        var clickedRow = rows[valueIdx];
        var oldColorHex = clickedRow.color;
        var defColorHex = clickedRow.defColor;
        var valueKey = clickedRow.strKey;
        if (valueKey==="")
            valueKey="_EMPTY_";

        /* jshint validthis: true */
        var newCol = $(this).spectrum('get');

        var newColHex = "";
        if (newCol===null)
            newColHex = oldColorHex; // if user clicked abort, revert to default color
        else
            newColHex = newCol.toHex();
        clickedRow.color = newColHex;

        // save color to cart if necessary
        saveToUrl(COL_PREFIX+valueKey, newColHex, defColorHex);

        var colors = legendGetColors(gLegend.rows);
        renderer.setColors(colors);
        renderer.drawDots();
    }

    function updateGeneTableColors(cellIds) {
    /* change the colors of the gene table to reflect the expression of cellIds */
        var quickGenes = db.conf.quickGenes;
        if (quickGenes===undefined)
            return;

        if (cellIds===null)
            cellIds = [];

        var pal = makeColorPalette(datasetGradPalette, exprBinCount);

        console.time("avgCalc");
        for (var i=0; i<quickGenes.length; i++) {
            var sym = quickGenes[i][0];
            //console.log("updating colors of "+sym+" for "+cellIds.length+" cells");
            var geneExpr = db.quickExpr[sym];
            if (geneExpr===undefined) { // if any gene is not loaded yet, just quit
                console.log(sym+" is not loaded yet, not updating expr table colors");
                return;
            }
            var vec = geneExpr[0];
            var binInfo = geneExpr[1];
            var sum = 0;

            var avg = 0;
            if (cellIds!==null && cellIds.length!==0) {
                for (var ci=0; ci<cellIds.length; ci++) {
                    sum += vec[cellIds[ci]];
                }
                avg = Math.round(sum / cellIds.length);
            }
            var color = pal[avg];
	    var fontColor = "#333333";
	    if (isDark(color))
		fontColor = "white";
            $("#tpGeneBarCell_"+onlyAlphaNum(sym)).css({"background-color": "#"+color, "color" : fontColor});
        }
        console.timeEnd("avgCalc");
    }

    function makeFieldHistogram(metaInfo, selCellIds, metaVec) {
        /* count the values in metaInfo and return a sorted array of [count, fraction, valIndex] */
        let cellCount = selCellIds.length;
        if (!metaInfo.ui.shortLabels && (metaInfo.type === "float" || metaInfo.type === "int")) {
            // it's a numeric field, so let's create the labels now, they are derived from the bins
            let shortLabels = [];
            for (let bin of metaInfo.binInfo)
                shortLabels.push(labelForBinMinMax(bin[0], bin[1])); // 0,1 is min,max of the bin
            metaInfo.ui.shortLabels = shortLabels;
        }

        var metaCounts = {};
        // make an object with value -> count in the cells
        for (var i = 0; i < cellCount; i++) {
            var cellId = selCellIds[i];
            var metaVal = metaVec[cellId];
            metaCounts[metaVal] = 1 + (metaCounts[metaVal] || 0);
        }
        var categories = metaInfo.binInfo;
        var catCountIdx = 2;
        if (categories === undefined) {
            categories = metaInfo.valCounts;
            catCountIdx = 1;
        }

        // convert the object to an array (count, percent, value) and sort it by count
        var histoList = [];
        for (var key in metaCounts) {
            let intKey = parseInt(key);
            let count = metaCounts[key];
            let frac = (count / cellCount);
            let fracOfCategory;
            if (categories) {
                let catCellCount = categories[intKey][catCountIdx];
                fracOfCategory = count / catCellCount;
            }
            histoList.push([count, frac, intKey, fracOfCategory]);
        }
        if (metaInfo.type !== "float" && metaInfo.type !== "int") {
            histoList = histoList.sort(function (a, b) { return b[0] - a[0]; }); // reverse-sort by count
        }
        return histoList;
    }

    function updateMetaBarManyCells(cellIds) {
    /* update the meta fields on the left to reflect/summarize a list of cellIds */
        var metaFieldInfos = db.getMetaFields();
        //var metaData = gCurrentDataset.metaData;
        var cellCount = cellIds.length;

        if (db.allMeta===undefined) {
            alert("The meta information has not been loaded yet. Please wait and try again in a few seconds.");
            return;
        }

        $('#tpMetaTitle').text("Meta data of "+cellCount+" "+gSampleDesc+"s");

        // for every field...
        var metaHist = {};
        for (var metaIdx = 0; metaIdx < metaFieldInfos.length; metaIdx++) {
            var metaInfo = metaFieldInfos[metaIdx];

            var metaVec = [];
            if (metaInfo.isCustom)
                metaVec = metaInfo.arr;
            else
                metaVec = db.allMeta[metaInfo.name];

            if (metaVec===undefined) {
                var metaMsg = null;
                if (metaInfo.type!=="uniqueString") {
                    console.log("cellBrowser.js:updateMetaBarManyCells - could not find meta info");
                    metaMsg = "(still loading - please wait and retry)";
                }
                else
                    metaMsg = "(unique identifier field)";
                $('#tpMeta_'+metaIdx).html(metaMsg);
                continue;
            }

            var histoList = makeFieldHistogram(metaInfo, cellIds, metaVec); // reverse-sort by count
            metaHist[metaInfo.name] = histoList;

            // make a quick list of the top values for the sparklines, ignore the rest
            var countList = [];
            var otherCount = 0;
            for (let i=0; i < histoList.length; i++) {
                let count = histoList[i][0];
                if (i<SPARKHISTOCOUNT)
                    countList.push(count);
                else
                    otherCount+=count;
            }
            if (otherCount!==0)
                countList.push(otherCount);

            // update the UI

            var topCount = histoList[0][0];
            var topPerc  = histoList[0][1];
            var topVal   = metaInfo.ui.shortLabels[histoList[0][2]];
            var percStr = (100*topPerc).toFixed(1)+"%";

            if (topVal.length > 14)
                topVal = topVal.substring(0, 14)+"...";

            var label = "";
            if (histoList.length === 1) {
                label = topVal;
            } else {
                // numeric type, calculate mean and sd
                if (metaInfo.origVals && (metaInfo.type === "float" || metaInfo.type === "int")) {
                    var sum = 0;
                    for (var i = 0, I = cellIds.length; i < I; i++) {
                        sum += metaInfo.origVals[cellIds[i]];
                    }
                    var mean = sum / cellIds.length;
                    sum = 0;
                    for (i = 0, I = cellIds.length; i < I; i++) {
                        sum += (metaInfo.origVals[cellIds[i]] - mean)**2;
                    }
                    var sd = Math.sqrt(sum / cellIds.length);
                    label = "<span class='tpMetaMultiVal'>mean = " + mean.toFixed(2) + "; sd = " + sd.toFixed(2) + "</span>";
                } else {
                    if (histoList[0][0] === 1) {
                        label = "<span class='tpMetaMultiVal'>" + histoList.length + " unique values</span>";
                    } else {
                        label = "<span class='tpMetaMultiVal'>" + percStr + " " +topVal+"</span>";
                    }
                }
            }

            $('#tpMeta_' + metaIdx).html(label);
        }
        db.metaHist = metaHist;
    }

    function clearMetaAndGene() {
        /* called when user hovers over nothing - clear the meta and gene field field info, hide the tooltip */
        if (db===null) // users moved the mouse while the db is still loading
            return;

        $('#tpMeta_custom').html("");

        var fieldCount = db.getMetaFields();
        for (let metaInfo of db.getMetaFields()) {
            $('#tpMeta_'+metaInfo.name).attr('title', "").html("");
        }
        $('#tpMetaNote').hide();
        updateGeneTableColors(null);
    }

    function updateMetaBarCustomFields(cellId) {
        /* update custom meta fields with custom data */
        if (!db.getMetaFields()[0].isCustom)
            return;

        let metaInfo = db.getMetaFields()[0];
        let intVal = metaInfo.arr[cellId];
        let strVal = metaInfo.ui.shortLabels[intVal];
        let rowDiv = $('#tpMeta_custom').html(strVal);
    }

    function updateMetaBarOneCell(cellInfo, otherCellCount) {
        /* update the meta bar with meta data from a single cellId */
        $('#tpMetaTitle').text(METABOXTITLE);

        let customCount = 0;
        if (db.getMetaFields()[0].isCustom)
            customCount = 1;

        let fieldInfos = db.getMetaFields();

        for (var i = 0; i < cellInfo.length; i++) {
            var fieldValue = cellInfo[i];
            let metaIdx = i + customCount;
            let metaInfo = fieldInfos[metaIdx];

            let rowDiv = $('#tpMeta_'+i);
            if (fieldValue.startsWith("http") && fieldValue.endsWith(".png")) {
                rowDiv.css('height', "40px");
                rowDiv.html("<img src='"+fieldValue+"'></img>");
            } else
                rowDiv.html(fieldValue);
            rowDiv.attr('title', cellInfo[i]);
        }

        if (otherCellCount===0)
            $("#tpMetaNote").hide();
        else {
            $("#tpMetaNote").html("...and "+(otherCellCount)+" other "+gSampleDesc+"s underneath");
            $("#tpMetaNote").show();
        }
    }

    function clearSelectionState() {
        /* clear URL variable with select state, called when user clicks cells or unselects them */
        delState("select");

        $("#tpHoverHint").show();
        $("#tpSelectHint").hide();
    }

    function onCellClickOrHover (cellIds, ev) {
        /* user clicks onto a circle with the mouse or hovers over one.
         * ev is undefined if not a click. cellIds is none if click was on empty background. */

        // do nothing if only hover but something is already selected
        var selCells = renderer.getSelection();
        if (ev===undefined && selCells.length!==0) {
            $("#tpHoverHint").hide();
            $("#tpSelectHint").show();
            return;
            }

        $("#tpHoverHint").show();
        $("#tpSelectHint").hide();

        // if click into background, remove any legend selection
        if (cellIds===null)
            $(".tpLegendLabel").removeClass("tpLegendSelect");

        if (cellIds===null || cellIds.length===0) {
            clearMetaAndGene();
        } else {
            var cellId = cellIds[0];
            var cellCountBelow = cellIds.length-1;
            updateMetaBarCustomFields(cellId);
            db.loadMetaForCell(cellId, function(ci) { updateMetaBarOneCell(ci, cellCountBelow);}, onProgress);
        }

        updateGeneTableColors(cellIds);

        if (ev===undefined) {
            db.metaHist = null;
        } else {
            // it was a click -> we have at least one cell ID
            let cellId = cellIds[0];
            if (!ev.shiftKey && !ev.ctrlKey && !ev.metaKey)
                renderer.selectClear();
            clearSelectionState();
            renderer.selectAdd(cellId);
            renderer.drawDots();
            event.stopPropagation();
        }
    }

    function showTooltip(x, y, labelStr) {
    $("#tpTooltip").css({
        "display":"block",
        "left" : x,
        "top" : y,
       }).html(labelStr);
    }

    function hideTooltip() {
        $("#tpTooltip").hide();
    }

    function onClusterNameHover(clusterName, nameIdx, ev) {
        /* user hovers over cluster label */
        var labelLines = [clusterName];

        var labelField = db.conf.labelField;
        var metaInfo = db.findMetaInfo(labelField);
        var longLabels = metaInfo.ui.longLabels;
        if (longLabels) {
            for (let i=0; i<longLabels.length; i++) {
                let shortLabel = metaInfo.ui.shortLabels[i];
                let longLabel = longLabels[i];
                if (clusterName===shortLabel && longLabel!==shortLabel) {
                    labelLines.push(longLabels[i]);
                    break;
                }
            }
        }

        if (labelField == db.conf.activeLabelField) {
            if (db.conf.topMarkers!==undefined) {
                labelLines.push("Top enriched/depleted markers: "+db.conf.topMarkers[clusterName].join(", "));
            }
            labelLines.push("");

            if (db.conf.markers!==undefined)
                labelLines.push("Click to show full marker gene list.");

            if (db.conf.clusterPngDir!==undefined) {
                var fullPath = cbUtil.joinPaths([db.name, db.conf.clusterPngDir, clusterName+".png"]);
                labelLines.push("<img src='"+fullPath+"'>");
            }
        }

        labelLines.push("Alt/Option-Click to select cluster; Shift-Click to add cluster to selection");
        showTooltip(ev.clientX+15, ev.clientY, labelLines.join("<br>"));
    }

    function onNoClusterNameHover(ev) {
        hideTooltip();
    }

    function sanitizeName(name) {
        /* ported from cellbrowser.py: remove non-alpha, allow underscores */
        var newName = name.replace(/\+/g, "Plus").replace(/-/g, "Minus").replace(/%/g, "Perc");
        var newName = newName.replace(/[^a-zA-Z_0-9+]/g, "");
        return newName;
    }

    function onlyAlphaNum(name) {
        /* only allow alphanumeric characters */
        var newName = name.replace(/[^a-zA-Z0-9+]/g, "");
        return newName;
    }

    function onActRendChange(otherRend) {
        /* called after the user has activated a view with a click */
        renderer.legend = gLegend;
        renderer = otherRend;
        gLegend = otherRend.legend;
        $("#tpLayoutCombo").val( otherRend.coordIdx ).trigger('chosen:updated');
        buildLegendBar();
    }

    function removeSplit(renderer) {
        /* stop split screen mode */
        if (!renderer)
            return;
        if (!renderer.childPlot && !renderer.parentPlot)
            return;
        if (!renderer.isMain) {
            // make sure the left renderer is the active one
            renderer.childPlot.activatePlot();
        }
        renderer.unsplit();
        $("#tpSplitMenuEntry").text("Split Screen");
    }

    function onSplitClick() {
        /* user clicked on View > Split Screen */
        if (!renderer.childPlot && !renderer.parentPlot) {
            // nothing is split yet -> start the split
            renderer.onActiveChange = onActRendChange;

            var currCoordIdx = $("#tpLayoutCombo").val();
            renderer.legend = gLegend;
            renderer.coordIdx = currCoordIdx; // keep for onActRendChange
            renderer.isMain = true;

            renderer.split();

            renderer.childPlot.legend = gLegend;
            renderer.childPlot.coordIdx = currCoordIdx; // keep for onActRendChange

            $("#tpSplitMenuEntry").text("Unsplit Screen");
            $("#mpCloseButton").click(renderer.unsplit);
            //$("#mpCloseButton").click(onSplitClick);

        } else {
            removeSplit(renderer);
        }
        renderer.drawDots();
    }

    function groupAverages(geneArrs, arrGroups, groupCount) {
        /* given an array of gene expression vectors (ints), and a 2nd array that assigns these to groups,
        return an array of the arrays with the averages for the groups (as integers)
        */
        let geneAvgs = [];

        for (let geneIdx=0; geneIdx < geneArrs.length; geneIdx++) {
            let geneArr = geneArrs[geneIdx];

            let groupSums = new Uint32Array(groupCount);
            let groupCounts = new Uint32Array(groupCount);
            //for (var groupIdx=0; groupIdx < groupCount; groupIdx++) {
            //    groupSums.push(0);
            //    groupCounts.push(0);
            //}

            for (let i=0; i<geneArr.length; i++) {
                let group = arrGroups[i];
                groupSums[group] += geneArr[i];
                groupCounts[group]++;
            }

            let groupAvgs = [];
            for (var groupIdx=0; groupIdx < groupCount; groupIdx++)
                groupAvgs.push(Math.round(groupSums[groupIdx]/groupCounts[groupIdx]));
            geneAvgs.push(groupAvgs);

        }
        return geneAvgs;
    }

    function onHeatCellClick(geneName, clusterName) {
        /* color by gene and select all cells in cluster */
        colorByLocus(geneName);
        // clusterName?
        selectByColor
    }

    function onHeatCellHover(rowIdx, colIdx, rowName, colName, value, ev) {
        /* user hovers over a cell on the heatmap */
        let htmls = [];
        if (rowName)
            htmls.push(rowName);
        if (colName)
            htmls.push(colName)
        if (value!==null)
            htmls.push(" "+(value*10)+"-"+((value+1)*10)+"%");
        showTooltip(ev.clientX+15, ev.clientY, htmls.join(" "));
    }

    function plotHeatmap(clusterMetaInfo, exprVecs, geneSyms) {
        /* Create the heatmap from exprVecs.
        */
        if (!geneSyms || geneSyms.length===0) {
            alert("No quick genes are defined. Heatmaps currently only work on pre-defined gene sets.");
            return;
        }

        var clusterCount = clusterMetaInfo.valCounts.length;

        var clusterNames = [];
        for (let valInfo of clusterMetaInfo.valCounts) {
            clusterNames.push(valInfo[0]); // 0=name, 1=count
        }

        var clusterArr = clusterMetaInfo.arr;
        var geneAvgs = groupAverages(exprVecs, clusterArr, clusterCount);

        var div = document.createElement("div");
        //let heatHeight = Math.min(150, 16*exprVecs.length);
        let heatHeight = parseInt(renderer.height*0.5);
        div.id = "tpHeat";
        div.style.height = heatHeight+"px";

        renderer.setSize(renderer.getWidth(), renderer.height-heatHeight, true);

        var canvLeft = metaBarWidth+metaBarMargin;
        var heatWidth = window.innerWidth - canvLeft - legendBarWidth;
        // create the div for the heat map view
        div.style.width = heatWidth+"px";
        div.style.left = metaBarWidth+"px";
        div.style.top = (menuBarHeight+toolBarHeight+renderer.height)+"px";
        div.style.position = "absolute";
        document.body.appendChild(div);

        var heatmap = new MaxHeat(div, {mainRenderer:renderer});
        //var colors = getFieldColors(clusterMetaInfo)
        var colors = makeColorPalette(cDefGradPaletteHeat, 10);

        heatmap.loadData(geneSyms, clusterNames, geneAvgs, colors);
        heatmap.draw();
        heatmap.onCellHover = onHeatCellHover;
        heatmap.onClick = onHeatCellClick;
        db.heatmap = heatmap;
    }


    function removeHeatmap() {
        /* remove the heatmap */
        let heatHeight = db.heatmap.height;
        document.getElementById("tpHeat").remove();
        delete db.heatmap;
        renderer.setSize(renderer.getWidth(), renderer.height+heatHeight, true);
        changeUrl({'heat':null});
    }

    function onHeatClick() {
        // TODO: rewrite this one day with promises...
        let resultCount = 0;
        let exprVecs = [];
        let geneSyms = [];
        let metaInfo = null;

        function partDone() {
            resultCount++;
            if (resultCount===2)
                plotHeatmap(metaInfo, exprVecs, geneSyms);
        }

        function onClusterMetaDone(metaArr, metaInfo) {
            metaInfo.arr = metaArr;
            partDone();
        }

        function onGenesDone(geneVecs) {
            /* */
            for (var geneInfo of geneVecs) {
                geneSyms.push(geneInfo[0]); // gene symbol
                exprVecs.push(geneInfo[1]); // binned expression vector
            }
            partDone();
        }

        /* user clicked on View > Heatmap */
        if (db && db.heatmap) {
            removeHeatmap();
        }
        else {
            if (!db.conf.quickGenes) {
                alert("No quick genes defined for this dataset. Heatmaps currently only work if "+
                    "a list of dataset-specific genes is defined. " +
                    "Add a statement quickGenesFile to cellbrowser.conf and put a few gene symbols "+
                    "into the file, one per line.");
                return;
            }
            db.loadGeneSetExpr(onGenesDone);
            metaInfo = getClusterFieldInfo();
            db.loadMetaVec(metaInfo, onClusterMetaDone, onProgress);
            changeUrl({"heat":"1"});
        }
    }

    function onClusterNameClick(clusterName, _, event) {
        /* build and open the dialog with the marker genes table for a given cluster */
        var metaInfo = getClusterFieldInfo();
        var isNumber = false;
        var nameIdx = null;
        if (metaInfo.type == "int" || metaInfo.type == "float") {
            isNumber = true;
        } else {
            nameIdx = metaInfo.ui.shortLabels.indexOf(clusterName);
        }
        if (event.altKey || event.shiftKey) {
            db.loadMetaVec(metaInfo, function(values) {
                var clusterCells = [];
                for (var i = 0, I = values.length; i < I; i++) {
                    if (isNumber && metaInfo.origVals[i].toFixed(2) == clusterName) {
                        clusterCells.push(i);
                    } else if (!isNumber && values[i] == nameIdx) {
                        clusterCells.push(i);
                    }
                }
                if (event.altKey) {
                    renderer.selectSet(clusterCells);
                } else if (event.shiftKey) {
                    var selection = renderer.getSelection();
                    selection = selection.concat(clusterCells);
                    renderer.selectSet(selection);
                }
                renderer.drawDots();
            });
            return;
        }
        if (db.conf.labelField != db.conf.activeLabelField) {
            return;
        }

        var tabInfo = db.conf.markers; // list with (label, subdirectory)

        console.log("building marker genes window for "+clusterName);
        var htmls = [];
        htmls.push("<div id='tpPaneHeader' style='padding:0.4em 1em'>");

        var buttons = [];

        if (tabInfo===undefined || tabInfo.length===0) {
            tabInfo = [];
            buttons.push({
                text:"Close",
                click: function() { $( this ).dialog( "close" ) }
            });
            htmls.push("No marker genes are available in this dataset. " +
                "To add marker genes, contact the original authors of the dataset and ask them to add " +
                " them to the cell browser.");
        } else {
            htmls.push("Click gene symbols below to color plot by gene<br>");
            buttons.push({
                text:"Download as file",
                click: function() {
                    document.location.href = markerTsvUrl;
                }
            });
        }
        htmls.push("</div>");

        var doTabs = (tabInfo.length>1);

        if (doTabs) {
            htmls.push("<div id='tabs'>");
            htmls.push("<ul>");
            for (var tabIdx = 0; tabIdx < tabInfo.length; tabIdx++) {
                var tabLabel = tabInfo[tabIdx].shortLabel;
                htmls.push("<li><a href='#tabs-"+tabIdx+"'>"+tabLabel+"</a>");
            }
            htmls.push("</ul>");
        }

        for (let tabIdx = 0; tabIdx < tabInfo.length; tabIdx++) {
            var divName = "tabs-"+tabIdx;
            var tabDir = tabInfo[tabIdx].name;
            var sanName = sanitizeName(clusterName);
            var markerTsvUrl = cbUtil.joinPaths([db.name, "markers", tabDir, sanName+".tsv.gz"]);
            htmls.push("<div id='"+divName+"'>");
            htmls.push("Loading...");
            htmls.push("</div>");

            loadClusterTsv(markerTsvUrl, loadMarkersFromTsv, divName, clusterName);
        }

        htmls.push("</div>"); // tabs

        var winWidth = window.innerWidth - 0.10*window.innerWidth;
        var winHeight = window.innerHeight - 0.10*window.innerHeight;
        var title = "Cluster markers for &quot;"+clusterName+"&quot;";

        var metaInfo = getClusterFieldInfo();
        if (metaInfo.ui.longLabels) {
            //var nameIdx = cbUtil.findIdxWhereEq(metaInfo.ui.shortLabels, 0, clusterName);
            //var acronyms = db.conf.acronyms;
            //title += " - "+acronyms[clusterName];
            var longLabel = metaInfo.ui.longLabels[nameIdx];
            if (clusterName!==longLabel)
                title += " - "+metaInfo.ui.longLabels[nameIdx];
        }

        //if (acronyms!==undefined && clusterName in acronyms)
            //title += " - "+acronyms[clusterName];
        showDialogBox(htmls, title, {width: winWidth, height:winHeight, "buttons":buttons});
        $(".ui-widget-content").css("padding", "0");
        $("#tabs").tabs();
    }

    function geneListFormat(htmls, s, symbol) {
    /* transform a string in the format dbName|linkId|mouseOver;... to html and push these to the htmls array */
        var dbParts = s.split(";");
        for (var i = 0; i < dbParts.length; i++) {
            var dbPart = dbParts[i];
            var idParts = dbPart.split("|");

            var dbName = idParts[0];
            var linkId = null;
            var mouseOver = "";

            // linkId and mouseOver are optional
            if (idParts.length>1) {
                linkId = idParts[1];
            }
            if (idParts.length>2) {
                mouseOver = idParts[2];
            }

            var dbUrl = dbLinks[dbName];
            if (dbUrl===undefined)
                htmls.push(dbName);
            else {
                if (linkId==="" || linkId===null)
                    linkId = symbol;
                htmls.push("<a target=_blank title='"+mouseOver+"' data-placement='auto left' class='link' href='"+dbUrl+linkId+"'>"+dbName+"</a>");
            }

            if (i!==dbParts.length-1)
                htmls.push(", ");
        }
    }

    function loadMarkersFromTsv(papaResults, url, divId, clusterName) {
        /* construct a table from a marker tsv file and write as html to the DIV with divID */
        console.log("got coordinate TSV rows, parsing...");
        var rows = papaResults.data;

        var headerRow = rows[0];

        var htmls = [];

        var markerListIdx = parseInt(divId.split("-")[1]);
        var selectOnClick = db.conf.markers[markerListIdx].selectOnClick;

        htmls.push("<table class='table' id='tpMarkerTable'>");
        htmls.push("<thead>");
        var hprdCol = null;
        var geneListCol = null;
        var exprCol = null;
        var pValCol = null
        //var doDescSort = false;
        for (var i = 1; i < headerRow.length; i++) {
            var colLabel = headerRow[i];
            var isNumber = false;

            if (colLabel.indexOf('|') > -1) {
                var parts = colLabel.split("|");
                colLabel = parts[0];
                var colType = parts[1];
                if (colType==="int" || colType==="float")
                    isNumber = true;
            }

            var width = null;
            if (colLabel==="_geneLists") {
                colLabel = "Gene Lists";
                geneListCol = i;
            }
            else if (colLabel==="P_value" || colLabel==="p_val" || colLabel==="pVal") {
                colLabel = "P-value";
                pValCol = i;
                //if (i===2)
                    //doDescSort = true;
            }

            else if (colLabel==="_expr") {
                colLabel = "Expression";
                exprCol = i;
            }
            else if (colLabel==="_hprdClass") {
                hprdCol = i;
                colLabel = "Protein Class (HPRD)";
                width = "200px";
            }

            var addStr = "";
            if (isNumber)
                addStr = " data-sort-method='number'";

            if (width===null)
                htmls.push("<th"+addStr+">");
            else
                htmls.push("<th style='width:"+width+"'"+addStr+">");
            colLabel = colLabel.replace(/_/g, " ");
            htmls.push(colLabel);
            htmls.push("</th>");
        }
        htmls.push("</thead>");

        var hubUrl = makeHubUrl();

        htmls.push("<tbody>");
        for (let i = 1; i < rows.length; i++) {
            var row = rows[i];
            if ((row.length===1) && row[0]==="") // papaparse sometimes adds empty lines to files
                continue;

            htmls.push("<tr>");
            var geneId = row[0];
            var geneSym = row[1];
            htmls.push("<td><a data-gene='"+geneSym+"' class='link tpLoadGeneLink'>"+geneSym+"</a>");
            if (hubUrl!==null) {
                var fullHubUrl = hubUrl+"&position="+geneSym+"&singleSearch=knownCanonical";
                htmls.push("<a target=_blank class='link' style='margin-left: 10px; font-size:80%; color:#AAA' title='link to UCSC Genome Browser' href='"+fullHubUrl+"'>Genome</a>");
            }
            htmls.push("</td>");

            for (var j = 2; j < row.length; j++) {
                var val = row[j];
                htmls.push("<td>");
                // added for the autism dataset, allows to add mouse overs with images
                // field has to start with ./
                if (val.startsWith("./")) {
                    var imgUrl = val.replace("./", db.url+"/");
                    var imgHtml = '<img width="100px" src="'+imgUrl+'">';
                    val = "<a data-toggle='tooltip' data-placement='auto' class='tpPlots link' target=_blank title='"+imgHtml+"' href='"+ imgUrl + "'>plot</a>";
                }
                if (j===geneListCol || j===exprCol)
                    geneListFormat(htmls, val, geneSym);
                else if (j===pValCol)
                    htmls.push(parseFloat(val).toPrecision(5)); // five digits ought to be enough for everyone
                else
                    htmls.push(val);
                htmls.push("</td>");
            }
            htmls.push("</tr>");
        }

        htmls.push("</tbody>");
        htmls.push("</table>");

        // sub function ----
        function onMarkerGeneClick(ev) {
            /* user clicks onto a gene in the table of the marker gene dialog window */
            var geneSym = ev.target.getAttribute("data-gene");
            $(".ui-dialog").remove(); // close marker dialog box
            if (selectOnClick) {
                clusterField = db.conf.labelField;
                var queryList = [{'m':clusterField, 'eq':clusterName}];
                findCellsMatchingQueryList(queryList, function(cellIds) {
                        renderer.selectSet(cellIds);
                        //changeUrl({'select':JSON.stringify(queryList)});
                });
            }
            colorByLocus(geneSym);
        }
        // ----

        $("#"+divId).html(htmls.join(""));
        var sortOpt = {};
        //if (doDescSort)
            //sortOpt.descending=true;
        new Tablesort(document.getElementById('tpMarkerTable'));
        $(".tpLoadGeneLink").on("click", onMarkerGeneClick);
        activateTooltip(".link");

        var ttOpt = {"html": true, "animation": false, "delay":{"show":100, "hide":100} };
        $(".tpPlots").bsTooltip(ttOpt);
    }

    var digitTest = /^\d+$/,
        keyBreaker = /([^\[\]]+)|(\[\])/g,
        plus = /\+/g,
        paramTest = /([^?#]*)(#.*)?$/;

    function deparam(params){
    /* https://github.com/jupiterjs/jquerymx/blob/master/lang/string/deparam/deparam.js */
        if(! params || ! paramTest.test(params) ) {
            return {};
        }

        var data = {},
            pairs = params.split('&'),
            current;

        for(var i=0; i < pairs.length; i++){
            current = data;
            var pair = pairs[i].split('=');

            // if we find foo=1+1=2
            if(pair.length !== 2) {
                pair = [pair[0], pair.slice(1).join("=")];
            }

            var key = decodeURIComponent(pair[0].replace(plus, " ")),
                value = decodeURIComponent(pair[1].replace(plus, " ")),
                parts = key.match(keyBreaker);

            for ( var j = 0; j < parts.length - 1; j++ ) {
                var part = parts[j];
                if (!current[part] ) {
                    // if what we are pointing to looks like an array
                    current[part] = digitTest.test(parts[j+1]) || parts[j+1] === "[]" ? [] : {};
                }
                current = current[part];

            }
            var lastPart = parts[parts.length - 1];
            if(lastPart === "[]"){
                current.push(value);
            }else{
                current[lastPart] = value;
            }
        }
        return data;
    }

    function changeUrl(vars, oldVars) {
    /* push the variables (object) into the history as the current URL. key=null deletes a variable. */
       // first get the current variables from the URL of the window
       var myUrl = window.location.href;
       myUrl = myUrl.replace("#", "");
       var urlParts = myUrl.split("?");
       var baseUrl = urlParts[0];

       let urlVars;
       if (oldVars===undefined) {
           var queryStr = urlParts[1];
           urlVars = deparam(queryStr); // parse key=val&... string to object
       } else {
           urlVars = oldVars;
       }

       // overwrite everthing that we got
       for (var key in vars) {
           var val = vars[key];
           if (val===null || val==="") {
               if (key in urlVars)
                   delete urlVars[key];
            } else
               urlVars[key] = val;
       }

       var argStr = jQuery.param(urlVars); // convert to query-like string
       argStr = argStr.replace(/%20/g, "+");

       var dsName = "noname";
       if (db!==null)
            dsName = db.getName();

       if (argStr.length > 1000)
           warn("Cannot save current changes to the URL, the URL would be too long. "+
               "You can try to shorten some cluster labels to work around the problem temporarily. "+
               "But please contact us at cells@ucsc.edu and tell us about the error. Thanks!");
        else
           history.pushState({}, dsName, baseUrl+"?"+argStr);
    }

    function delVars(varNames) {
        /* remove a CGI variable from the URL */
        var o = {};
        for (var varName of varNames)
            o[varName] = null;
        changeUrl(o);
    }

    function delState(varName) {
        /* remove a CGI variable from the URL */
        var o = {};
        o[varName] = null;
        changeUrl(o);
    }

    function addStateVar(varName, varVal) {
        /* add a CGI variable to the URL */
        var o = {};
        o[varName] = varVal;
        changeUrl(o);
    }

    function getVar(name, defVal) {
        /* get query variable from current URL or default value if undefined */
       var myUrl = window.location.href;
       myUrl = myUrl.split("#")[0]; // remove anchor from URL
       var urlParts = myUrl.split("?");
       var queryStr = urlParts[1];
       var varDict = deparam(queryStr); // parse key=val&... string to object
       if (varDict[name]===undefined)
           return defVal;
       else
           return varDict[name];
    }

    function getVarSafe(name, defVal) {
        let val = getVar(name, defVal);
        if (val)
            val = val.replace(/\W/g, '');
        return val;
    }

    function pushZoomState(zoomRange) {
        /* write the current zoom range to the URL. Null to remove it from the URL. */
       if (zoomRange===null)
            changeUrl({zoom:null});
       else
           changeUrl({zoom:zoomRange.minX.toFixed(5)+"_"+zoomRange.maxX.toFixed(5)+"_"+zoomRange.minY.toFixed(5)+"_"+zoomRange.maxY.toFixed(5)});
    }

    function getZoomRangeFromUrl() {
        /* return a zoomRange object based on current URL */
        var zoomStr = getVar("zoom", null);
        if (zoomStr===null)
            return null;
        var zs = zoomStr.split("_");
        if (zs.length!==4)
            return null;
        var zoomRange = {};
        zoomRange.minX = parseFloat(zs[0]);
        zoomRange.maxX = parseFloat(zs[1]);
        zoomRange.minY = parseFloat(zs[2]);
        zoomRange.maxY = parseFloat(zs[3]);
        return zoomRange;
    }

    function redirectIfSubdomain() {
        /* rewrite the URL if at ucsc and subdomain is specified
         * e.g. autism.cells.ucsc.edu -> cells.ucsc.edu?ds=autism */
        /* we cannot run in the subdomain, as otherwise localStorage and
         * cookies are not shared */

        // at UCSC, the dataset can be part of the hostname
        // we got a "* CNAME" in the campus DNS server for this.
        // it's easier to type, and pretty in manuscripts e.g.
        // autism.cells.ucsc.edu instead of cells.ucsc.edu?ds=autism
        var myUrl = new URL(window.location.href);
        var hostName = myUrl.hostname;
        if (hostName.endsWith("cells.ucsc.edu")) {
            var hostParts = hostName.split(".");
            if (hostParts.length===4) {
                var datasetName = hostParts[0];
                hostParts.shift();
                myUrl.hostname = hostParts.join(".");
                var newUrl = myUrl+"?ds="+datasetName;
                window.location.replace(newUrl);
                return true;
            }
        return false;
        }
    }

    function getDatasetNameFromUrl() {
        /* search for the "ds" parameter or a DNS hostname that indicates the dataset */
        // if ds=xxx was found in the URL, load the respective dataset
        var datasetName = getVar("ds");
        if (datasetName)
            datasetName = datasetName.replace(/ /g, "/"); // + is easier to type than %23

        if (datasetName===undefined)
            datasetName = "";
        // hacks for July 2018 and for backwards compatibility with previous version
        if (datasetName==="autism10X" || datasetName==="autism10x")
            datasetName = "autism";
        if (datasetName==="aparna")
            datasetName = "cortex-dev";

        // adult pancreas is the only dataset with an uppercase letter
        // make sure that at least at UCSC, dataset names are always lowercased.
        // The reason is that at UCSC, we the datasetname can be part of the URL,
        // e.g. cortex-dev.cells.ucsc.edu, which the user could enter as CoRTex-dev.cells.ucsc.edu
        // On all other servers, this is on an issue
        if (datasetName && datasetName!=="adultPancreas" && pageAtUcsc())
            datasetName = datasetName.toLowerCase();
        return datasetName;
    }

    /* ==== MAIN ==== ENTRY FUNCTION */
    function main(rootMd5) {
        /* start the data loaders, show first dataset. If in  */
        if (redirectIfSubdomain())
            return;

        setupKeyboard();
        buildMenuBar();

        var datasetName = getDatasetNameFromUrl()
        // pre-load dataset.json here?
        menuBarHeight = $('#tpMenuBar').outerHeight(true);

        var canvLeft = metaBarWidth+metaBarMargin;
        var canvTop  = menuBarHeight+toolBarHeight;
        var canvWidth = window.innerWidth - canvLeft - legendBarWidth;
        var canvHeight = window.innerHeight - menuBarHeight - toolBarHeight;

        if (renderer===null) {
           var div = document.createElement('div');
           div.id = "tpMaxPlot";
           renderer = new MaxPlot(div, canvTop, canvLeft, canvWidth, canvHeight);
           document.body.appendChild(div);
           activateTooltip(".mpButton"); // tpMaxPlot has no special tooltip support itself

           self.tooltipDiv = makeTooltipCont();
           document.body.appendChild(self.tooltipDiv);
       }

        buildEmptyLegendBar(metaBarWidth+metaBarMargin+renderer.width, toolBarHeight);

        renderer.setupMouse();
        $(window).resize(onWindowResize);

        renderer.onLabelClick = onClusterNameClick;
        renderer.onLabelHover = onClusterNameHover;
        renderer.onNoLabelHover = onNoClusterNameHover;
        renderer.onCellClick = onCellClickOrHover;
        renderer.onCellHover = onCellClickOrHover;
        renderer.onNoCellHover = clearMetaAndGene;
        renderer.onZoom100Click = onZoom100Click;
        renderer.onSelChange = onSelChange;
        renderer.canvas.addEventListener("mouseleave", hideTooltip);

        loadDataset(datasetName, false, rootMd5);
    }

    // only export these functions
    return {
        "main":main
    }

}();



function _tpReset() {
/* for debugging: reset the intro setting */
    localStorage.removeItem("introShown");
}
