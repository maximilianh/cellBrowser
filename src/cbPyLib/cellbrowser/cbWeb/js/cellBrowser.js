// A viewer for (x,y) scatter plots of small circles
// shows associated meta data (usually key->val attributes, one mapping per circle) 
// and expression data (string -> float, one mapping per circle)

/* jshint -W097 */
"use strict";

var tsnePlot = function() {
    var db = null; // the cbData object from cbData.js. Loads coords,
                   // annotations and gene expression vectors
    
    var gDatasetList = null; // array of dataset descriptions (objects)

    var gVersion = "0.3";
    var gCurrentCoordName = null; // currently shown coordinates

    // object with all information needed to map to the legend colors
    var gLegend = null;
    // all info about the current legend. gLegend.rows is:
    // [ colorHexStr, defaultColorHexStr, label, count, internalKey, uniqueKey ]
    // internalKey can be int or str, depending on current coloring mode. 
    // E.g. in meta coloring, it's the metadata string value.
    // When doing expression coloring, it's the expression bin index.
    // uniqueKey is used to save manually defined colors to localStorage

    var renderer = null;
    var gWinInfo = null; // .width and .height of the PIXI canvas

    // -- CONSTANTS
    var gTitle = "UCSC Cell Browser";
    const COL_PREFIX = "col_";

    // depending on the type of data, single cell or bulk RNA-seq, we call a circle a 
    // "sample" or a "cell". This will adapt help menus, menus, etc.
    var gSampleDesc = "cell";

    // width of left meta bar in pixels
    var metaBarWidth = 200;
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
    var datasetComboWidth = 200;
    var layoutComboWidth = 150;

    // height of bottom gene bar
    var geneBarHeight = 100;
    var geneBarMargin = 5;
    // color for missing value when coloring by expression value
    //var cNullColor = "CCCCCC";
    //const cNullColor = "DDDDDD";
    const cNullColor = "87CEFA";

    const cDefGradPalette = "reds";  // default legend gradient palette for numeric ranges
    const cDefQualPalette  = "rainbow"; // default legend palette for categorical values

    const exprBinCount = 10; //number of expression bins for genes

    var HIDELABELSNAME = "Hide labels";
    var SHOWLABELSNAME = "Show labels";
    var METABOXTITLE   = "Cell Annotations";

    // maximum number of distinct values that one can color on
    const MAXCOLORCOUNT = 200;

    // histograms show only the top X values and summarize the rest into "other"
    var HISTOCOUNT = 12;
    // the sparkline is a bit shorter
    var SPARKHISTOCOUNT = 12;

    // links to various external databases
    var dbLinks = {
        "HPO" : "http://compbio.charite.de/hpoweb/showterm?gene=", // entrez ID
        "OMIM" : "https://omim.org/entry/", // OMIM ID
        "COSMIC" : "http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=", // gene symbol
        "SFARI" : "https://gene.sfari.org/database/human-gene/", // gene symbol
        "BrainSpLMD" : "http://www.brainspan.org/lcm/search?exact_match=true&search_type=gene&search_term=", // entrez
        "BrainSpMouseDev" : "http://developingmouse.brain-map.org/gene/show/", // internal Brainspan ID
        "Eurexp" : "http://www.eurexpress.org/ee/databases/assay.jsp?assayID=", // internal ID
        "LMD" : "http://www.brainspan.org/lcm/search?exact_match=true&search_type=gene&search_term=" // entrez
    }

    var DEBUG = true;

    function _dump(o) {
    /* for debugging */
        console.log(JSON.stringify(o));
    }

    function debug(msg, args) {
        if (DEBUG) {
            console.log(formatString(msg, args));
        }
    }

    function warn(msg) {
        alert(msg);
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

    function keys(o) {
    /* return all keys of object */
        var allKeys = [];
        for(var k in o) allKeys.push(k);
        return allKeys;
    }

    function capitalize(s) {
        return s[0].toUpperCase() + s.slice(1);
    }

    function cartSave(key, value, defaultValue) {
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

    function cartGet(key, defaultValue) {
    /* get a value from localStorage or the current URL or return the default if not defined in either place.
     * The URL overrides localStorage. */
        var val = localStorage.getItem(key);
        if (val===null && defaultValue!==undefined)
            val = defaultValue;
        val = getVar(key, val);
        return val
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
        var ttOpt = {"html": true, "animation": false, "delay":{"show":350, "hide":100}, container:"body"}; 
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

    function updateToolbar() {
    /* update the toolbar with the current dataset info */
        var geneList = [];

        for(var key in gCurrentDataset.matrixOffsets) {
            geneList.push({id:key, text:key});
        }

        //$('#tpGeneCombo').select2({
            //placeholder : "Gene Symbol",
            //ajax : {},
            //dataAdapter : RefAdapter,
        //});
    }

    function updateMenu() {
    /* deactivate menu options based on current variables */
     // the "hide selected" etc menu options are only shown if some cells are selected
     if (renderer.selCells===null) {
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

    function prettyNumber(/*int*/ count) /*str*/ {
        /* convert a number to a shorter string, e.g. 1200 -> 1.2k, 1200000 -> 1.2M, etc */
        if (count>1000000) {
            var f = (count / 1000000);
            return f.toFixed(1)+"M";
        }
        if (count>10000) {
            var f = (count / 1000);
            return f.toFixed(0)+"k";
        }
        if (count>1000) {
            var f = (count / 1000);
            return f.toFixed(1)+"k";
        }

        return count;
    }

    function openDatasetLoadPane(selDatasetIdx) {
        /* open dataset dialog: load html into the three panes  */
        var datasetName = gDatasetList[selDatasetIdx].name;
        var descUrl = joinPaths([datasetName, "summary.html"]);
        $("#pane1").load(descUrl, function( response, status, xhr ) {
            if ( status === "error" ) {
                $( "#pane1" ).html("File "+descUrl+" was not found");
            }
            $("#tabLink1").tab("show");
        });

        var methodsUrl = joinPaths([datasetName, "methods.html"]);
        $("#pane2").load(methodsUrl, function( response, status, xhr ) {
            if ( status === "error" ) {
                $( "#pane2" ).html("File "+methodsUrl+" was not found");
            }
        });

        var downloadUrl = joinPaths([datasetName, "downloads.html"]);
        $("#pane3").load(downloadUrl, function( response, status, xhr ) {
            if ( status === "error" ) {
                $( "#pane3" ).html("File "+downloadUrl+" was not found");
            }
        });
    }

    function openDatasetDialog() {
    /* build dataset open dialog */

        var winWidth = window.innerWidth - 0.05*window.innerWidth;
        var winHeight = window.innerHeight - 0.05*window.innerHeight;
        var buttonWidth = 300;
        var tabsWidth = winWidth - buttonWidth - 50;


        var htmls = [];
        var activeIdx = 0;
        htmls.push("<div class='list-group' style='width:"+buttonWidth+"px'>");
        for (var i = 0; i < gDatasetList.length; i++) {
            var dataset = gDatasetList[i];
            var line = "<button id='tpDatasetButton_"+i+"' type='button' class='list-group-item' data-datasetid='"+i+"'>"; // bootstrap seems to remove the id
            htmls.push(line);

            if (dataset.sampleCount!==undefined) {
                var countDesc = prettyNumber(dataset.sampleCount);
                htmls.push("<span class='badge'>"+countDesc+"</span>");
            }

            if (dataset.tags!==undefined) {
                for (var tagI = 0; tagI < dataset.tags.length; tagI++) {
                var tag = dataset.tags[tagI];
                htmls.push("<span class='badge'>"+tag+"</span>");
                }
            }
            htmls.push(dataset.shortLabel+"</button>");
            if (db!==null && db.name===dataset.name)
                activeIdx = i;
        }
        htmls.push("</div>"); // list-group

        htmls.push("<div id='tpOpenDialogLabel' style='width:"+tabsWidth+"px; position:absolute; left: 340px; top: 10px;'>");
        htmls.push("<div id='tpOpenDialogTabs'>");
        htmls.push("<ul class='nav nav-tabs'>");
        htmls.push("<li class='active'><a id='tabLink1' data-toggle='tab' href='#pane1'>Abstract</a></li>");
        htmls.push("<li><a id='tabLink2' data-toggle='tab' href='#pane2'>Methods</a></li>");
        htmls.push("<li><a id='tabLink3' data-toggle='tab' href='#pane3'>Data Download</a></li>");
        htmls.push("</ul>");
        htmls.push("</div>");

        htmls.push("<div class='tab-content'>");

        htmls.push("<div id='pane1' class='tab-pane'>");
        htmls.push("<p>Loading abstract...</p>");
        htmls.push("</div>");

        htmls.push("<div id='pane2' class='tab-pane'>");
        htmls.push("<p>Loading methods...</p>");
        htmls.push("</div>");

        htmls.push("<div id='pane3' class='tab-pane'>");
        htmls.push("<p>Loading data download...</p>");
        htmls.push("</div>");

        htmls.push("</div>"); // tab-content

        htmls.push("</div>"); // tpOpenDialogLabel

        //htmls.push("<div id='tpSelectedId' data-selectedid='0'>"); // store the currently selected datasetId in the DOM
        var selDatasetIdx = 0;

        var buttons = {
        "Open Dataset" :
            function(event) {
                $( this ).dialog( "close" );
                var datasetName = gDatasetList[selDatasetIdx].name;
                loadDataset(datasetName, true);
            }
        };

        if (db!==null)
            buttons["Cancel"] = function() { $( this ).dialog( "close" ); };

        showDialogBox(htmls, "Open Cell Browser Dataset", {width: winWidth, height:winHeight, "buttons":buttons});

        $("button.list-group-item").eq(selDatasetIdx).css("z-index", "1000"); // fix up first overlap
        $("button.list-group-item").keypress(function(e) {
            // load the current dataset when the user presses Return
            if (e.which == '13') {
                var datasetName = gDatasetList[selDatasetIdx].name;
                loadDataset(datasetName, true);
                $(".ui-dialog-content").dialog("close");
            }
        });

        $(".list-group-item").click( function (ev) {
            selDatasetIdx = parseInt($(ev.target).data('datasetid')); // index of clicked dataset
            $(".list-group-item").removeClass("active");
            $('#tpDatasetButton_'+selDatasetIdx).bsButton("toggle"); // had to rename .button() in .html
            openDatasetLoadPane(selDatasetIdx);
        });

        $("#tabLink1").tab("show");

        $(".list-group-item").focus( function (event) {
            selDatasetIdx = parseInt($(event.target).data('datasetid')); // index of clicked dataset
            // bootstrap has a bug where the blue selection frame is hidden by neighboring buttons
            // we're working around this here by bumping up the current z-index.
            $("button.list-group-item").css("z-index", "0");
            $("button.list-group-item").eq(selDatasetIdx).css("z-index", "1000");
        });

        if (activeIdx!==null)
            $('#tpDatasetButton_'+activeIdx).bsButton("toggle"); // had to rename .button() in .html to bsButton
        
        // this is weird, but I have not found a better way to make the tab show up
        $("#tpOpenDialogTabs a:last").tab("show");
        $("#tpOpenDialogTabs a:first").tab("show");

        // finally, activate the default pane and load its html
        $("button.list-group-item").eq(activeIdx).trigger("focus");
        openDatasetLoadPane(activeIdx);
    }

    function drawLayoutMenu() {
       /* Add the View / Layout menu to the DOM */
      var htmls = [];  
      var coordFiles = gCurrentDataset.coordFiles;
      for (var i=0; i<coordFiles.length; i++) {
           var label = coordFiles[i].shortLabel;
           if (label===undefined) {
               label = "Dataset "+i;
               console.log("Dataset "+i+" has no 'shortLabel' attribute");
           }
           htmls.push('<li><a href="#" class="tpDataset" id="'+i+'">'+label+'</a></li>');
      }
      $("#tpLayoutMenu").empty();
      $("#tpLayoutMenu").append(htmls.join(""));
    }

    function onSelChange(cellIds) {
    /* called each time when the selection has been changed */
        updateGeneTableColors(cellIds);
        if ("geneSym" in gLegend)
            buildViolinPlot(gLegend.geneSym);

        //if (cellIds.length===0)
            //clearMetaAndGene();
        //else if (cellIds.length===1)
            //updateMetaBarOneCell(cellIds[0]);
        //else
            //updateMetaBarManyCells(cellIds);
    }

    function onSaveAsClick() {
    /* File - Save Image as ... */
        var canvas = $("canvas")[0];
        canvas.toBlob(function(blob) { saveAs( blob , "cellBrowser.png"); } , "image/png");
    }

    function onDownloadClick() {
    /* user clicks one of the "download data" submenu links: matrix, meta or coords  */
        var name = event.target.id.split("_")[1]; // the name of the dataset
        var baseUrl = gCurrentDataset.baseUrl;
        var url = null;
        if (name==="matrix") {
            url = joinPaths([baseUrl,"geneMatrix.tsv"]);
            // hack for aparna, to send the original matrix
            if (baseUrl==="aparna/" || baseUrl==="aparna") 
                url = joinPaths([baseUrl,"matrix.tsv.gz"]);
            document.location.href = url;
        }
        if (name==="meta") {
            url = joinPaths([baseUrl,"meta.tsv"]);
            document.location.href = url;
        }
    }

    function onSelectAllClick() {
    /* Edit - select all visible*/
        renderer.selectClear();
        renderer.selectVisible();
        //updateSelection();
        renderer.drawDots();
    }

    function onSelectNoneClick() {
    /* Edit - Select None */
        renderer.selectClear();
        renderer.drawDots();
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

    function onSelectByIdClick() {
        /* Edit - select cells by ID */
        var dlgHeight = 500;
        var dlgWidth = 500;
        var htmls = [];
        var buttons = {
        "OK" :
            function() {
                var idListStr = $("#tpIdList").val();
                idListStr = idListStr.trim().replace(/\r\n/g,"\n");
                gSelCellIds = {};
                if (idListStr==="") 
                    return;

                // first check the IDs
                var allCellIds = getAllCellIdsAsDict();
                var notFoundIds = [];
                var idList = idListStr.split("\n");
                for (var i = 0; i < idList.length; i++) {
                    var cellId = idList[i];
                    if (cellId==="")
                        continue;
                    if (!(cellId in allCellIds))
                        notFoundIds.push(cellId);
                }

                if (notFoundIds.length!==0) {
                    //alert("Could not find these "+gSampleDesc+" IDs:"+ notFoundIds.join(", "));
                    $('#tpNotFoundIds').text("Could not find these IDs: "+notFoundIds.join(", "));
                    $('#tpNotFoundHint').text("Please fix them and click the OK button to try again.");
                    }
                else {
                    for (i = 0; i < idList.length; i++) {
                        var cellId = idList[i];
                        if (cellId==="")
                            continue;
                        gSelCellIds[cellId] = true;
                        }
                    $( this ).dialog( "close" );
                    //updateSelection();
                    plotDots();
                    renderer.render(stage);
                    //alert(idList.length+ " " + gSampleDesc+"s are selected");
                }
            }
        };

        htmls.push("<textarea id='tpIdList' style='height:320px;width:350px;display:block'>");
        htmls.push("</textarea><div id='tpNotFoundIds'></div><div id='tpNotFoundHint'></div>");
        var title = "Paste a list of IDs (one per line) to select "+gSampleDesc+"s";
        showDialogBox(htmls, title, {showClose:true, height:dlgHeight, width:dlgWidth, buttons:buttons});
    }

    function onExportIdsClick() {
        /* Edit - Export cell IDs */
        var idList = [];
        for (var cellId in gSelCellIds) {
            idList.push(cellId);
        }

        var dlgHeight = 500;

        var htmls = [];
        if (idList.length===0)
            {
            htmls.push("No cells are selected. Shown below are the identifiers of all cells visible on the screen.<p>");
            for (var i = 0; i < pixelCoords.length; i++)
                idList.push(shownCoords[i][0]);
            }

        var idListEnc = encodeURIComponent(idList.join("\n"));
        htmls.push("<textarea style='height:320px;width:350px;display:block'>");
        htmls.push(idList.join("\n"));
        htmls.push("</textarea>");
        
        //htmls.push('<a style="text-decoration:underline" href="data:text/plain;charset=utf-8,'+idListEnc+'" download="identifiers.txt">Download as text file</a><br>');

        //htmls.push('<a style="text-decoration:underline" href="data:text/plain;charset=utf-8,'+idListEnc+'" download="identifiers.txt">Download as text file</a><br>');

        var buttons = {
        "Download as file" :
            function() {
                var blob = new Blob([idList.join("\n")], {type: "text/plain;charset=utf-8"});
                saveAs(blob, "identifiers.txt");
            },
        "Copy to clipboard" :
            function() {
                $("textarea").select();
                document.execCommand('copy');
                $( this ).dialog( "close" );
            }
        };


        showDialogBox(htmls, "List of "+idList.length+" IDs", {showClose:true, height:dlgHeight, width:400, buttons:buttons});
    }

    function buildMenuBar() {
        /* draw the menubar at the top */
       var htmls = [];
       htmls.push("<div id='tpMenuBar'>");
       htmls.push('<nav class="navbar navbar-default navbar-xs">');

       htmls.push('<div class="container-fluid">');

         htmls.push('<div class="navbar-header">');
           htmls.push('<a class="navbar-brand" href="#">'+gTitle+'</a>');
         htmls.push('</div>');

         htmls.push('<ul class="nav navbar-nav">');

         htmls.push('<li class="dropdown">');
           htmls.push('<a href="#" class="dropdown-toggle" data-toggle="dropdown" data-submenu role="button" aria-haspopup="true" aria-expanded="false">File</a>');
           htmls.push('<ul class="dropdown-menu">');
             htmls.push('<li><a href="#" id="tpOpenDatasetLink"><span class="dropmenu-item-label">Open Dataset...</span><span class="dropmenu-item-content">o</span></a></li>');
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
         htmls.push('<li><a id="tpSelectAll" href="#"><span class="dropmenu-item-label">Select all visible</span><span class="dropmenu-item-content">a</span></a></li>');
         htmls.push('<li><a id="tpSelectNone" href="#"><span class="dropmenu-item-label">Select none</span><span class="dropmenu-item-content">n</span></a></li>');
         //htmls.push('<li><a id="tpMark" href="#"><span class="dropmenu-item-label">Mark selected</span><span class="dropmenu-item-content">h m</span></a></li>');
         //htmls.push('<li><a id="tpMarkClear" href="#"><span class="dropmenu-item-label">Clear marks</span><span class="dropmenu-item-content">c m</span></a></li>');
         //htmls.push('<li><a id="tpSelectById" href="#">Search for ID...</a></li>');
         //htmls.push('<li><a id="tpExportIds" href="#">Export selected IDs...</a></li>');
         htmls.push('</ul>'); // View dropdown
         htmls.push('</li>'); // View dropdown

         htmls.push('<li class="dropdown">');
         htmls.push('<a href="#" class="dropdown-toggle" data-toggle="dropdown" data-submenu role="button" aria-haspopup="true" aria-expanded="false">View</a>');
         htmls.push('<ul class="dropdown-menu">');

         htmls.push('<li><a href="#" id="tpZoomPlus"><span class="dropmenu-item-label">Zoom in</span><span class="dropmenu-item-content">+</span></a></li>');
         htmls.push('<li><a href="#" id="tpZoomMinus"><span class="dropmenu-item-label">Zoom out</span><span class="dropmenu-item-content">-</span></a></li>');
         htmls.push('<li><a href="#" id="tpZoom100Menu"><span class="dropmenu-item-label">Zoom to 100%</span><span class="dropmenu-item-content">spc</span></a></li>');

         htmls.push('<li><hr class="half-rule"></li>');

         //htmls.push('<li><a href="#" id="tpOnlySelectedButton">Show only selected</a></li>');
         //htmls.push('<li><a href="#" id="tpFilterButton">Hide selected '+gSampleDesc+'s</a></li>');
         //htmls.push('<li><a href="#" id="tpShowAllButton">Show all '+gSampleDesc+'</a></li>');
         htmls.push('<li><a href="#" id="tpHideShowLabels">Hide labels<span class="dropmenu-item-content">c l</span></a></li>');
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
         htmls.push('<a href="#" class="dropdown-toggle" data-toggle="dropdown" data-submenu role="button" aria-haspopup="true" aria-expanded="false">Help</a>');
         htmls.push('<ul class="dropdown-menu">');
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
       $('#tpZoomPlus').click( onZoomInClick );
       $('#tpZoomMinus').click( onZoomOutClick );
       //$('#tpShowAllButton').click( onShowAllClick );
       $('#tpHideShowLabels').click( onHideShowLabelsClick );
       $('#tpExportIds').click( onExportIdsClick );
       $('#tpSelectById').click( onSelectByIdClick );
       $('#tpMark').click( onMarkClick );
       $('#tpMarkClear').click( onMarkClearClick );
       $('#tpTutorialButton').click( function()  { showIntro(false); } );
       $('#tpOpenDatasetLink').click( openDatasetDialog );
       $('#tpSaveImage').click( onSaveAsClick );
       $('#tpSelectAll').click( onSelectAllClick );
       $('#tpSelectNone').click( onSelectNoneClick );
       $('#tpDownloadMenu li a').click( onDownloadClick );

       // This version is more like OSX/Windows:
       // - menus only open when you click on them
       // - once you have clicked, they start to open on hover
       // - a click anywhere else will stop the hovering
       var doHover = false;
       $(".nav > .dropdown").click( function(){ doHover = true;} );
       $(".nav > .dropdown").hover(
           function(event) { 
               if (doHover) {  
                   $(".dropdown-submenu").removeClass("open"); $(".dropdown").removeClass("open"); $(this).addClass('open'); 
               } 
           });
       $(document).click ( function() { doHover= false; });

       $('[data-submenu]').submenupicker();

    }

    function resizeDivs() {
       /* resize all divs and the renderer to current window size */
       var rendererLeft = metaBarWidth+metaBarMargin;
       var rendererHeight  = window.innerHeight - menuBarHeight - toolBarHeight;

       var rendererWidth = window.innerWidth - legendBarWidth - rendererLeft;
       var legendBarLeft = rendererWidth+metaBarMargin+metaBarWidth;

       $("#tpToolBar").css("width", rendererWidth);
       $("#tpToolBar").css("height", toolBarHeight);
       $("#tpLeftSidebar").css("height", window.innerHeight - menuBarHeight);
       $("#tpLegendBar").css("height", window.innerHeight - menuBarHeight);
       $('#tpLegendBar').css('left', legendBarLeft+"px");

       renderer.setSize(rendererWidth, rendererHeight);
    }

    var progressUrls = [];

    function onProgressConsole(ev) {
        console.log(ev);
    }

    function onProgress(ev) {
        /* show progress bars */
        var url = ev.currentTarget.responseURL;
        if (url.search("exprMatrix.bin")!==-1)
            return;

        var index = progressUrls.indexOf(url);
        if (index===-1) {
            progressUrls.push(url);
            index = progressUrls.length-1;
        }

        var label = url;
        if (url.endsWith("coords.bin"))
            label = "Loading Coordinates";
        else if (url.endsWith(".bin"))
            label = "Loading cell annotations";
        var labelId = "#tpProgressLabel"+index;
        $(labelId).html(label);

        var percent = Math.round(100 * (ev.loaded / ev.total));

        if (percent===100) {
            $("#tpProgress"+index).css("width", percent+"%");
            $("#tpProgress"+index).show(0);
            progressUrls.splice(index, 1);
            $("#tpProgressDiv"+index).css("display", "none");
        }
        else {
            $("#tpProgress"+index).css("width", percent+"%");
            $("#tpProgressDiv"+index).css("display", "inherit");
        }
    }

    function colorByMetaField(fieldName, doneLoad) {
       /* load the meta data for a field, setup the colors, send it all to the renderer and call doneLoad */
       if (doneLoad===undefined)
           doneLoad = function() { renderer.drawDots() };

       if (fieldName===null) {
           // obscure hacky option: you can set the default color field to "None"
           // so there is no coloring at all on startup
           renderer.setColors(["black"]);
           var cellCount = db.conf.sampleCount;
           renderer.setColorArr(new Uint8Array(cellCount)); 
           gLegend.rows = [];
           gLegend.rows.push( [ "000000", null, "No Value", cellCount, 0, null] );
           doneLoad();
           return;
       }

       var fieldIdx  = db.fieldNameToIndex(fieldName);
       console.log("Color by meta field "+fieldName);

       // internal field names cannot contain non-alpha chars, so tolerate user errors here
       // otherwise throw an error
       if (fieldIdx === null) {
           fieldIdx = db.fieldNameToIndex(fieldName.replace(/[^0-9a-z]/gi, ''));
           if (fieldIdx === null) {
               alert("The field "+fieldName+" does not exist in the sample/cell annotations. Cannot color on it.")
               return;
           }
       }

       var fieldInfo = db.getMetaFields()[fieldIdx];

       if (fieldInfo.diffValCount > MAXCOLORCOUNT && fieldInfo.binMethod===undefined) {
           warn("This field has "+fieldInfo.diffValCount+" different values. Coloring on a field that has more than "+MAXCOLORCOUNT+" different values is not supported.");
           return null;
       }

       var defaultMetaField = db.getDefaultColorField()[1];
       if (fieldName!==defaultMetaField)
           changeUrl({"meta":fieldName, "gene":null});


       buildLegendForMetaIdx(fieldIdx);
       var renderColors = legendGetColors(gLegend.rows);
       renderer.setColors(renderColors);

       db.loadMetaVec(fieldIdx, function(carr) {renderer.setColorArr(carr); doneLoad(); } , onProgress);

       //changeUrl({"gene":null, "meta":fieldName});
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

    function splitExpr(exprVec, selCells) {
        /* split the expression vector into two vectors, one for selected and one for unselected cells */
        var selMap = {};
        for (var i = 0; i < selCells.length; i++) {
            selMap[i] = null;
        }

        var sel = [];
        var unsel = [];
        for (i = 0; i < exprVec.length; i++) {
            if (i in selMap)
                sel.push(exprVec[i]);
            else
                unsel.push(exprVec[i]);
        }

        return [sel, unsel]
    }

    function buildViolinPlot(geneSym) {
        /* create the violin plot */
        console.time("violin");

        var quickExpr = db.quickExpr[geneSym];
        // quickExpr format is: [exprVec, geneDesc, binInfo]
        var exprVec = quickExpr[0];

        $('#tpViolinCanvas').remove();
        var htmls = [];
        htmls.push("<canvas style='height:200px; padding-top: 10px; padding-bottom:30px' id='tpViolinCanvas'></canvas>");
        $('#tpViolin').append(htmls.join(""));

        const ctx = document.getElementById("tpViolinCanvas").getContext("2d");

        var dataList = [];
        var labelList = [];
        var selCells = renderer.getSelection();
        if (selCells===null) {
            dataList = [Array.prototype.slice.call(exprVec)];
            labelList = ['All cells'];
        } else {
            dataList = splitExpr(exprVec, selCells);
            labelList = ['Selected', 'Others'];
        }

 	var boxplotData = {
          labels : labelList,
	  datasets: [{
            data : dataList,
	    label: 'Dataset 1',
	    backgroundColor: 'rgba(255,0,0,0.5)',
	    borderColor: 'red',
	    borderWidth: 1,
	    outlierColor: '#999999',
	    padding: 7,
	    itemRadius: 0,
	    outlierColor: '#999999'
	}]
	};
	
        if ("violinChart" in window)
            window.violinChart.destroy();

	window.violinChart = new Chart(ctx, {
	    type: 'violin',
	    data: boxplotData,
	    options: {
	      //scales: {
		//xAxes: [{
		    //ticks: {
			//autoSkip: false,
			//maxRotation: 40,
			//minRotation: 40
		    //}
		//}],
	      //},
	      //responsive: true,
              maintainAspectRatio: false,
	      //legend: {
		//position: 'top',
	      //},
              legend: { display: false },
	      title: { display: false }
	    }
	  });
        console.timeEnd("violin");
    }

    function loadGeneAndColor(geneSym, onDone) {
        /* color by a gene, load the array into the renderer and call onDone or just redraw */
        if (onDone===undefined)
            onDone = function() { renderer.drawDots() };

        function gotGeneVec(exprArr, decArr, geneSym, geneDesc, binInfo) {
            /* called when the expression vector has been loaded */
            if (decArr===null)
                return;
            console.log("Received expression vector, gene "+geneSym+", geneId "+geneDesc);
            _dump(binInfo);
            makeLegendExpr(geneSym, geneDesc, binInfo);
            db.quickExpr[geneSym] = [exprArr, geneDesc, binInfo];
            buildLegendBar();
            renderer.setColors(legendGetColors(gLegend.rows));
            renderer.setColorArr(decArr);
            buildViolinPlot(geneSym);
            onDone();
        }

        changeUrl({"gene":geneSym, "meta":null, "pal":null});
        console.log("Loading gene expression vector for "+geneSym);
        db.loadExprVec(geneSym, gotGeneVec, onProgress, exprBinCount);

        // clear the meta combo
        $('#tpMetaCombo').val(0).trigger('chosen:updated');
    }

    function preloadAllMeta() {
        /* start loading the full meta value vectors and add them to db.quickMeta */
        var metaFieldInfo = db.getMetaFields();
        db.quickMeta = {};
        for (var fieldIdx = 0; fieldIdx < metaFieldInfo.length; fieldIdx++) {
           var fieldInfo = metaFieldInfo[fieldIdx];
           if (fieldInfo.type==="uniqueString")
               continue
           db.loadMetaVec(fieldIdx, function(carr) {db.quickMeta[fieldIdx]=carr} , null);
        }
    }

    function preloadQuickGenes() {
       /* start loading the quick gene expression vectors in the background now 
        * add them to db.quickExpr */
       var quickGenes = db.conf.quickGenes;
       db.quickExpr = {};
       var validGenes = db.getGenes();
       var loadCounter = 0;
       if (quickGenes) {
           for (var i=0; i<quickGenes.length; i++) {
               var sym = quickGenes[i][0];
               if (! (sym in validGenes)) {
                  alert("Error: "+sym+" is in quick genes list but is not a valid gene");
                  continue;
               }

               db.loadExprVec(
                   sym, 
                   function(exprVec, geneSym, geneDesc, binInfo) {
                      db.quickExpr[geneSym] = [exprVec, geneDesc, binInfo];
                      loadCounter++;
                      if (loadCounter===quickGenes.length)
                        updateGeneTableColors(null);
                   },
                   onProgressConsole, exprBinCount);
           }
       }
    }

   function gotCoords(coords, info, clusterMids) {
       /* called when the coordinates have been loaded */
       if (coords.length===0)
           alert("cellBrowser.js/gotCoords: coords.bin seems to be empty");
       renderer.setCoords(coords, clusterMids, info.minX, info.maxX, info.minY, info.maxY);
   }

    function renderData() {
    /* init the renderer, start loading and draw data when ready
     */
       var loadsDone = 0;

       function doneOnePart() {
       /* make sure renderer only draws when both coords and other data have loaded */
           loadsDone +=1;
           if (loadsDone===2) {
               buildLegendBar();
               renderer.setColors(legendGetColors(gLegend.rows));
               renderer.setTitle(db.conf.shortLabel);
               renderer.drawDots();
           }
       }

       function gotFirstCoords(coords, info, clusterMids) {
           /* very ugly way to implement promises. Need a better approach one day. */
           gotCoords(coords, info, clusterMids);
           doneOnePart();
       }

       renderer.initPlot(db.conf);

       buildLeftSidebar(db.getMetaFields());
       buildToolBar(db.conf.coords, db.conf.name, metaBarWidth+metaBarMargin, toolBarHeight);
       activateMode("move");

       db.loadCoords(0, gotFirstCoords, onProgress);

       var colorByInfo = db.getDefaultColorField();
       var colorType = colorByInfo[0];
       var colorBy = colorByInfo[1];

       // allow to override coloring by URL args
       if (getVar("gene", null)!==null) {
           colorType = "gene";
           colorBy = getVar("gene");
           activateTab("gene");
       }
       else if (getVar("meta", null)!==null) {
           colorType = "meta";
           colorBy = getVar("meta");
           activateTab("meta");
       }

       var forcePalName = getVar("pal", null);

       gLegend = {};
       if (colorType==="meta") {
           colorByMetaField(colorBy, doneOnePart);
           // update the meta field combo box
           var fieldIdx  = db.fieldNameToIndex(colorBy);
           $('#tpMetaCombo').val(fieldIdx).trigger('chosen:updated');
       }
       else {
           loadGeneAndColor(colorBy, doneOnePart);
           // update the gene combo box
           var sel = $('#tpGeneCombo')[0].selectize;
           sel.addOption({text: colorBy, value: colorBy});
           sel.refreshOptions();
           sel.setTextboxValue(colorBy);
       }

       if (forcePalName!==null) {
           legendChangePaletteAndRebuild(forcePalName);
           renderer.drawDots();
        }

       preloadQuickGenes();
       preloadAllMeta();
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
            rows.sort(function(a, b) { return naturalSort(a[2], b[2]); });
        }
        else {
            // sort this list by count = index 3
            rows.sort(function(a, b) { return b[3] - a[3]; }); // reverse-sort by count
        }
        buildLegendBar();
    }
        
    //function makeLegendObject(sortBy) {
    /* create the gLegend object */
        //if (gLegend.type=="meta")
            //{
            //// color by meta attribute
            //gLegend = makeLegendMeta(gLegend.metaFieldIdx, sortBy);
            //}
        //else {
            //// color by gene
            //var geneIdx = gLegend.geneIdx;
            //var geneInfo = gCurrentDataset.preloadExpr.genes[geneIdx];
            //var geneId = geneInfo[0];
            //var geneSym = geneInfo[1];
            //var deciles = gCurrentDataset.preloadExpr.deciles[geneId];
            //var cellExpr = gCurrentDataset.preloadExpr.cellExpr;
            //makeLegendExpr(geneIdx, geneSym, deciles, cellExpr);
        //}
    //}

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
        if ($("#tpHideShowLabels").text()===SHOWLABELSNAME) {
            renderer.setShowLabels(true);
            $("#tpHideShowLabels").text(HIDELABELSNAME);
        }
        else {
            renderer.setShowLabels(false);
            $("#tpHideShowLabels").text(SHOWLABELSNAME);
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
        //#$("#tpZoom100Button").tooltip('hide');
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
        resizeDivs();
    }

    function onColorPaletteClick(ev) {
        /* called when users clicks a color palette */
        var palName = ev.target.getAttribute("data-palette");
        legendChangePaletteAndRebuild(palName);
        renderer.drawDots();
    }

    function buildEmptyLegendBar(fromLeft, fromTop) {
        // create an empty right side legend bar
        var htmls = [];
        htmls.push("<div id='tpLegendBar' style='position:absolute;top:"+fromTop+"px;left:"+fromLeft+"px; width:"+legendBarWidth+"px'>");
        htmls.push("<div class='tpSidebarHeader'>Legend");

        //htmls.push("<div id='tpToolbarButtons' style='padding-bottom: 2px'>");
        htmls.push("<div style='float:right' class='btn-group btn-group-xs'>");
            htmls.push("<button type='button' class='btn btn-default dropdown-toggle' data-toggle='dropdown' aria-haspopup='true' aria-expanded='false' id='tpChangeColorScheme'>Colors&nbsp;<span class='caret'> </span></button>");
            htmls.push('<ul class="dropdown-menu pull-right">');
            htmls.push('<li><a class="tpColorLink" data-palette="default" href="#">Default</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="rainbow" href="#">Rainbow Qualitative</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="tol-dv" href="#">Paul Tol&#39;s Qualitative</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="blues" href="#">Shades of Blues</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="reds" href="#">Shades of Reds</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="tol-sq" href="#">Beige to red</a></li>');
            htmls.push('<li><a class="tpColorLink" data-palette="tol-rainbow" href="#">Blue to red</a></li>');
            htmls.push('</ul>');
        htmls.push("</div>"); // btn-group
        //htmls.push("</div>"); // tpToolbarButtons

        htmls.push("</div>"); // tpSidebarHeader

        htmls.push("<div id='tpLegendTitleBox' style='position:relative; width:100%; height:1.5em; font-weight: bold'>");
                    htmls.push("<div id='tpLegendContent'>");
                    htmls.push("</div>"); // content 
                    htmls.push("<div id='tpViolin'>");
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
        var colArr = []
        for (var i = 0; i < rows.length; i++) {
            var row = rows[i];
            var col = row[0];
            if (col===null)
                col = row[1]; // only use default color if nothing else set

            var idx = row[4]; // 4 = meta val index or expression bin XX
            colArr[idx] = col; // 0 = color
        }

        return colArr;
    }

    function legendChangePaletteAndRebuild(palName) {
        /* change the legend color palette and put it into the URL and redraw */
        var rows = gLegend.rows;
        var success = legendSetPalette(gLegend, palName, gLegend.metaFieldIdx);
        if (success) {
            if (palName==="default")
                changeUrl({"pal":null});
            else
                changeUrl({"pal":palName});
            buildLegendBar();
            var colors = legendGetColors(gLegend.rows);
            renderer.setColors(colors);
        }
    }

    function legendSetPalette(legend, origPalName, metaFieldIndex) {
    /* update the defColor [1] attribute of all legend rows. pal is an array of hex colors. 
     * metaFieldIndex is optional. If it is set, will use the predefined colors that are
     * in the field configuration.
     * */
        var palName = origPalName;

        //palName = getVar("pal", palName);

        if (origPalName==="default") {
            legendResetColors();
            //changeUrl({"pal":null});
            if (legend.rowType==="category")
                palName = cDefQualPalette;
            else
                palName = cDefGradPalette;
        }

        var rows = legend.rows;
        var n = rows.length;
        var pal = null;
        var usePredefined = false;
        // if this is a field for which colors were defined manually, use them
        if (metaFieldIndex!==undefined && db.conf.metaFields[metaFieldIndex].colors!==undefined && origPalName==="default") {
            pal = db.conf.metaFields[metaFieldIndex].colors;
            usePredefined = true;
        } else
            pal = makeColorPalette(palName, n);

        if (pal===null) {
            alert("Sorry, this palette does not have "+rows.length+" different colors");
            return false;
        }

        for (var i = 0; i < rows.length; i++) {
            rows[i][1] = pal[i];
        }
        legend.palName = palName;

        // update the dropdown menu
        $('.tpColorLink').parent().removeClass("active");
        // force the menu to the "defaults" entry if we're using predefined colors
        if (usePredefined)
            palName = "default";
        $('.tpColorLink[data-palette="'+palName+'"]').parent().addClass("active");
        return true;
    }

    function makeLegendRowsNumeric(binInfo) {
        /* return an array of legend lines given bin info from gene expression or a numeric meta field  */
        var legendRows = [];
        //var zeroIsGrey = false;

        // special case for the first "0" element = no value, make this always grey
        //var bin0Min = binInfo[0][0];
        //var bin0Max = binInfo[0][1];
        //if (bin0Min===0 && bin0Max===0)
            //zeroIsGrey = true;

        //var defColors = makeColorPalette(10, true, zeroIsGrey);
        var colIdx = 0;
        for (var binIdx = 0; binIdx < binInfo.length; binIdx++) {
            var binMin = binInfo[binIdx][0];
            var binMax = binInfo[binIdx][1];

            var count  = binInfo[binIdx][2];
            var legendId = binIdx;

            // pretty print the numbers
            var minDig = 2;
            //if (binMin % 1 === 0) // % 1 = fractional part
                //minDig = 0

            var maxDig = 2;
            //if (binMin % 1 === 0)
             //   maxDig = 0

            var legLabel = null;
            if (binMin!==binMax)
                legLabel = binMin.toFixed(minDig)+' - '+binMax.toFixed(maxDig);
            else
                legLabel = binMin.toFixed(minDig);

            var uniqueKey = legLabel;
            var legColor = null;

            // override any color with the color specified in the current URL
            var savKey = COL_PREFIX+legLabel;
            var legColor = getVar(savKey, null);

            if (binMin===0 && binMax===0) {
                legLabel = "No Value";
                uniqueKey = "noExpr";
                legColor = cNullColor;
            }
            else
                colIdx++;

            legendRows.push( [ legColor, null, legLabel, count, binIdx, uniqueKey] );
        }
        return legendRows;
    }

    function makeLegendExpr(geneSym, geneId, binInfo) {
        /* build gLegend object for coloring by expression
         * return the colors as an array of hex codes */

        activateTooltip("#tpGeneSym");
        
        var legendRows = makeLegendRowsNumeric(binInfo);
        var colors = legendGetColors(legendRows);

        gLegend.rows = legendRows;
        gLegend.title = "Gene: "+geneSym;
        gLegend.titleHover = geneId;
        gLegend.geneSym = geneSym;
        gLegend.rowType = "range";
        legendSetPalette(gLegend, "default");
        return colors;
    }

    function onGeneClick (event) {
    /* user clicked on a gene in the gene table */
        var geneIdx = parseInt(event.target.id.split("_")[1]); // the index of the gene
        $('.tpMetaBox').removeClass('tpMetaSelect');
        $('.tpGeneBarCell').removeClass("tpGeneBarCellSelected");
        $('#tpGeneBarCell_'+geneIdx).addClass("tpGeneBarCellSelected");
        var geneSym = db.conf.quickGenes[geneIdx][0];
        loadGeneAndColor(geneSym);
    }

    function showDialogBox(htmlLines, title, options) {
        /* show a dialog box with html in it */
        $('#tpDialog').remove();

        var addStr = "";
        if (options.width!==undefined)
            addStr = "max-width:"+options.width+"px;";
        var maxHeight = $(window).height()-200;
        // unshift = insert at pos 0
        htmlLines.unshift("<div style='display:none;"+addStr+"max-height:"+maxHeight+"px' id='tpDialog' title='"+title+"'>");
        htmlLines.push("</div>");
        $(document.body).append(htmlLines.join(""));

        var dialogOpts = {modal:true, closeOnEscape:true};
        if (options.width!==undefined)
            dialogOpts["width"] = options.width;
        if (options.height!==undefined)
            dialogOpts["height"] = options.height;
        dialogOpts["maxHeight"] = maxHeight;
        if (options.buttons!==undefined)
            dialogOpts["buttons"] =  options.buttons;
        else
            dialogOpts["buttons"] =  {};

        if (options.showOk!==undefined)
            dialogOpts["buttons"].OK = function() { $( this ).dialog( "close" ); };
        if (options.showClose!=undefined)
            dialogOpts["buttons"].Cancel = function() { $( this ).dialog( "close" ); };
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
        var buttons = {"OK" : onGeneDialogOkClick };
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
            cellId = cellIds[i];
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
                cellId = cellIds[cellI];
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

    /**
     from https://stackoverflow.com/questions/29855098/is-there-a-built-in-javascript-function-similar-to-os-path-join:
     * Joins 2 paths together and makes sure there aren't any duplicate seperators
     * @param parts the parts of the url to join. eg: ['http://google.com/', '/my-custom/path/']
     * @param separator The separator for the path, defaults to '/'
     * @returns {string} The combined path
     */
    function joinPaths(parts, separator) {
      return parts.map(function(part) { return part.trim().replace(/(^[\/]*|[\/]*$)/g, ''); }).join(separator || '/');
    }

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

    //    makeLegendExpr(0, sExpr.symbol, sExpr.deciles, cellToExpr);
    //    buildLegendBar();
    //    plotDots();
    //    renderer.render(stage);

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
        var url = joinPaths([baseUrl, "geneMatrix.tsv"]);
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

        if (notFoundGenes.length!=0) {
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

        var showOk = (notFoundGenes.length!=0);
        showDialogBox(htmls, "Downloading expression data", {width:350, showOk:true});

        var progressLabel = $( "#tpProgressLabel" );
        $("#tpGeneProgress").progressbar( {
              value: false,
              max  : gLoad_geneList.length
              });

    }

    function buildGeneTable(htmls, tableWidth, cellWidth, geneInfos) {
    /* create gene expression info table */
        $('#tpGenes').remove();

        if (db.conf.quickGenes===undefined)
            return;


        htmls.push("<div style='margin-top:8px' id='tpGenes'>");
        htmls.push("<div style='padding-left:3px; font-weight:bold'>Quick Genes<div>");
        htmls.push('<div style="margin-top:6px" class="tpHint">');
        htmls.push('Hover or select cells to update colors.');
        htmls.push('</div>');
        htmls.push('<table style="margin-top:10px" id="tpGeneTable"><tr>');
        //htmls.push('<td><button id="tpChangeGenes" title="Change the list of genes that are displayed in this table" class = "ui-button ui-widget ui-corner-all" style="width:95%">Change</button></td>');

        var colsPerRow = Math.round(tableWidth / cellWidth);
        var cellWidth = Math.round(tableWidth/colsPerRow); 

        var currWidth = 1;
        for (var i = 0; i < geneInfos.length; i++) {
            var geneInfo = geneInfos[i];
            var geneId   = geneInfo[0];
            var geneDesc = geneInfo[1];
            var pubDesc = geneInfo[2];
            if (((i % colsPerRow) == 0) && (i!=0)) {
                htmls.push("</tr><tr>");
            }
            htmls.push('<td title="'+geneDesc+'" id="tpGeneBarCell_'+i+'" class="tpGeneBarCell">'+geneId+'</td>');
        }
        htmls.push("</tr></table>");

    }

    function likeEmptyString(label) {
    /* some special values like "undefined" and empty string get colored in grey  */
        return (label===null || label.trim()==="" || label==="none" || label==="None" || label==="unknown" 
                || label==="nd" || label==="n.d."
                || label==="Unknown" || label==="NaN" || label==="NA" || label==="undefined" || label==="Na")
    }

    function numMetaToBinInfo(fieldInfo) {
        /* convert a numeric meta field info to look like gene expression info for the legend:
         * an array of [start, end, count] */
        var binInfo = [];
        if (fieldInfo.binMethod==="uniform") {
            var binMin = fieldInfo.minVal;
            var stepSize = fieldInfo.stepSize;
            var binCounts = fieldInfo.binCounts;
            var binCount = fieldInfo.binCounts.length;
            for (var i=0; i<binCount; i++) {
                binInfo.push( [binMin, binMin+stepSize, binCounts[i]] );
                binMin+=stepSize;
            }
        } else if (fieldInfo.binMethod==="quantiles") {
            var binMin = fieldInfo.minVal;
            var breaks = fieldInfo.breaks;
            var binCounts = fieldInfo.binCounts;
            var binCount = fieldInfo.binCounts.length;
            for (var i=0; i<binCount; i++) {
                binInfo.push( [breaks[i], breaks[i+1], binCounts[i]] );
            }
        }
            
        return binInfo;
    }

    function makeLegendMeta(metaIndex, sortBy) {
    /* Build a new gLegend object and return it */
        var legend = {};
        legend.type = "meta";
        legend.metaFieldIdx = metaIndex;
        legend.titleHover = null;

        var fieldInfo = db.getMetaFields()[metaIndex];
        legend.fieldName = fieldInfo.label;
        legend.title = legend.fieldName.replace(/_/g, " ");

        // numeric meta fields are a special case
        if (fieldInfo.type==="int" || fieldInfo.type==="float") {
            var binInfo = numMetaToBinInfo(fieldInfo);
            legend.rows = makeLegendRowsNumeric(binInfo);
            legend.rowType = "range";
            legendSetPalette(legend, "default", metaIndex);
            return legend;
        }

        if (fieldInfo.diffValCount > MAXCOLORCOUNT) {
            warn("This field has "+fieldInfo.diffValCount+" different values. Coloring on a field that has more than "+MAXCOLORCOUNT+" different values is not supported.");
            return null;
        }

        var metaCounts = fieldInfo.valCounts;

        // we are going to sort this list later, so we need to keep track of what the original
        // index in the list was, as every meta data value is stored by its index, not
        // its label. We store the index as [2] of the metaCounts array.
        if (metaCounts[0].length!==3)
            for (var valIdx=0; valIdx < metaCounts.length; valIdx++)
                metaCounts[valIdx].push(valIdx);

        var oldSortBy = cartGet("SORT", legend.fieldName);
        if (sortBy===undefined && oldSortBy!==undefined)
            sortBy = oldSortBy;

        var fieldName = fieldInfo.label;

        // force field names that look like "cluster" to a rainbow palette
        // even if they look like numbers
        if (sortBy===undefined) {
            // should cluster fields be sorted by their name
            if (fieldName.indexOf("luster"))
                sortBy = "count";
            else if (fieldInfo.type==="float" || fieldInfo.type==="int")
                sortBy = "name";
            else
                sortBy = "count";
        }

        // sort like numbers if the strings are mostly numbers, otherwise sort like strings
        var sortResult = sortPairsBy(metaCounts, sortBy);
        var countListSorted = sortResult.list;

        var useGradient = (fieldInfo.type==="float" || fieldInfo.type==="int");
        //var defaultColors = makeColorPalette(countListSorted.length, useGradient);

        var rows = [];
        //var defColIdx = 0;
        for (var legRowIdx = 0; legRowIdx < countListSorted.length; legRowIdx++) {
            var legRowInfo = countListSorted[legRowIdx];
            var label = legRowInfo[0];
            var count = legRowInfo[1];
            var valIdx = legRowInfo[2];
            var uniqueKey = label;

            //var defColor = defaultColors[legRowIdx];
            if (likeEmptyString(label))
                color = cNullColor;
            //else {
                //defColor = defaultColors[defColIdx];
                //defColIdx++;
            //    }
            // override any color with the color specified in the current URL
            var savKey = COL_PREFIX+label;
            var color = cartGet(savKey, null);

            rows.push( [ color, null, label, count, valIdx, uniqueKey] );
        }

        legend.rows = rows;
        legend.isSortedByName = sortResult.isSortedByName;
        legend.rowType = "category";
        legendSetPalette(legend, "default", metaIndex);
        return legend;
    }

    function findCellIdsForLegendIds (cellIdToLegendId, legendIds, cellIds) {
    /* given a legendId, return an object with cellIds that are associated to a given class */
        if (cellIds==undefined)
            cellIds = {};
        for (var i = 0; i < legendIds.length; i++) {
            var legendId = legendIds[i];
            for (var cellId in cellIdToLegendId) {
                if (gClasses[cellId] == legendId)
                    cellIds[cellId] = true;
                }
        }
        return cellIds;
    }

    function buildLegendForMetaIdx(fieldId) {
    /* rebuild the legend */
        var legend = makeLegendMeta(fieldId);
        if (legend==null)
            return;

        $('.tpMetaBox').removeClass('tpMetaSelect');
        $('.tpMetaValue').removeClass('tpMetaValueSelect');
        $('#tpMetaBox_'+fieldId).addClass('tpMetaSelect');
        $('#tpMeta_'+fieldId).addClass('tpMetaValueSelect');
        $('.tpGeneBarCell').removeClass('tpGeneBarCellSelected');

        gLegend = legend;
        buildLegendBar();
        $('#tpLegendTitle').text(legend.fieldName.replace(/_/g, " "));
    }

    function onMetaClick (event) {
    /* called when user clicks a meta data field or label */
        var fieldId = parseInt(event.target.id.split("_")[1]);
        if (isNaN(fieldId)) {
            // try up one level in the DOM tree
            fieldId = parseInt(event.target.parentElement.id.split("_")[1]);
        }
        var fieldName = db.getMetaFields()[fieldId].name;
        colorByMetaField(fieldName, );
    }

    function addMetaTipBar(htmls, valFrac, valStr) {
        /* add another bar to a simple histogram built from divs */
        htmls.push("<div>&nbsp;");
        htmls.push("<div class='tpMetaTipPerc'>"+(100*valFrac).toFixed(1)+"%</div>");
        htmls.push("<div class='tpMetaTipName'>"+valStr+"</div>");
        //htmls.push("<span class='tpMetaTipCount'>"+valCount+"</span>");
        var pxSize = (valFrac * metaTipWidth).toFixed(0);
        htmls.push("<div style='width:"+pxSize+"px' class='tpMetaTipBar'>&nbsp</div>");
        htmls.push("</div>");
    }

    function onMetaMouseOver (event) {
    /* called when user hovers over meta element: shows the histogram of selected values */
        var metaHist = gCurrentDataset.metaHist;
        if (metaHist===undefined)
            return;

        // mouseover over spans or divs will not find the id, so look at their parent, which is the main DIV
        var target = event.target;
        if (target.id==="")
            target = event.target.parentNode;
        if (target.id==="")
            target = event.target.parentNode;
        var targetId = target.id;

        var fieldId = parseInt(targetId.split("_")[1]);

        var valHist = metaHist[fieldId];
        var htmls = [];
        var otherCount = 0;
        var totalSum = 0;
        for (var i = 0; i < valHist.length; i++) {
            var valInfo  = valHist[i];
            var valCount = valInfo[0];
            var valFrac  = valInfo[1];
            var valStr   = valInfo[2];

            totalSum += valCount;
            // only show the top values, summarize everything else into "other"
            if (i > HISTOCOUNT) {
                otherCount += valCount;
                continue;
            }

            if (valStr==="")
                valStr = "<span style='color:indigo'>(empty)</span>";
            addMetaTipBar(htmls, valFrac, valStr);
        }

        if (otherCount!==0) {
            var otherFrac = (otherCount / totalSum);
            addMetaTipBar(htmls, otherFrac, "<span style='color:indigo'>(other)</span>");
        }
        
        $('#tpMetaTip').html(htmls.join(""));
        $('#tpMetaTip').css({top: event.target.offsetTop+"px", left: metaBarWidth+"px", width:metaTipWidth+"px"});
        $('#tpMetaTip').show();
    }

    function buildComboBox(htmls, id, entries, selIdx, placeholder, width) {
    /* make html for a combo box and add lines to htmls list */
        htmls.push('<select style="'+width+'px" id="'+id+'" data-placeholder="'+placeholder+'" class="tpCombo">');
        for (var i = 0; i < entries.length; i++) {
            var isSelStr = "";
            var entry = entries[i];
            if ((selIdx!==undefined && i===selIdx))
                isSelStr = " selected"
            htmls.push('<option value="'+entry[0]+'"'+isSelStr+'>'+entry[1]+'</option>');
        }
        htmls.push('</select>');
    }

    function loadCoordSet(coordIdx) {
        //var coordUrl = gCurrentDataset.coordFiles[coordIdx].url;
        //console.log("Loading coordinates from "+coordUrl);
        //startLoadTsv("coords", coordUrl, loadCoordsFromTsv);
        db.loadCoords(coordIdx, 
                function(a,b,c) { gotCoords(a,b,c); renderer.drawDots();},
                onProgress);
    }

    function onLayoutChange(ev, params) {
        /* user changed the layout in the combobox */
        var coordIdx = parseInt(params.selected);
        loadCoordSet(coordIdx);
        changeUrl({"layout":coordIdx, "zoom":null});
        // remove the focus from the combo box
        removeFocus();
    }

    function onGeneChange(ev) {
        /* user changed the gene in the combobox */
        var geneSym = ev.target.value;
        if (geneSym==="") 
            return; // do nothing if user just deleted the current gene
        loadGeneAndColor(geneSym);
        $(this).blur(); // remove focus
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

    function loadDataset(datasetName, resetVars) {
        /* load a dataset and optionally reset all the URL variables.
         * When a dataset is opened through the UI, the variables have to
         * be reset, as their values (gene or meta data) may not exist
         * there. If it's opened via a URL, the variables must stay. */
        db = new CbDbFile(datasetName); 

        var vars = undefined;
        if (resetVars)
            vars = {};

        changeUrl({"ds":datasetName}, vars);
        db.loadConfig(function() { renderData() });
        
        // start the tutorial after a while
        var introShownBefore = localStorage.getItem("introShown");
        if (introShownBefore==undefined)
           setTimeout(function(){ showIntro(true); }, 5000); // show after 5 secs

    }

    function onDatasetChange(ev, params) {
        /* user changed the dataset in the dropbox */
        var datasetIdx = parseInt(params.selected);
        var datasetName = gDatasetList[datasetIdx].name;
        $(this).blur();
        removeFocus();
        loadDataset(datasetName, true);
    }

    function buildLayoutCombo(htmls, files, id, width, left, top) {
        /* files is a list of elements with a shortLabel attribute. Build combobox for them. */
        htmls.push('<div class="tpToolBarItem" style="position:absolute;left:'+left+'px;top:'+top+'px"><label for="'+id+'">Layout</label>');
        var entries = [];
        for (var i = 0; i < files.length; i++) {
            var coordFiles = files[i];
            var label = coordFiles.shortLabel;
            if (label===undefined)
                warn("Layout coordinate file "+i+" has no .shortLabel attribute")
            entries.push([i, label]);
        }
        buildComboBox(htmls, id, entries, 0, "Select a layout algorithm...", width);
        htmls.push('</div>');
    }

    function buildDatasetCombo(htmls, datasets, id, width, left, top) {
        /* datasets with a list of elements with a shortLabel attribute. Build combobox for them. */
        htmls.push('<div class="tpToolBarItem" style="position:absolute;width:150px;left:'+left+'px;top:'+top+'px"><label for="'+id+'">Dataset</label>');
        var entries = [];
        for (var i = 0; i < datasets.length; i++) {
            var dataset = datasets[i];
            var label = dataset.shortLabel;
            if (label===undefined)
                label = "Dataset "+i;
            entries.push([i, label]);
        }
        buildComboBox(htmls, id, entries, 0, "Select a dataset...", width);
        htmls.push('</div>');
    }

    function buildMetaFieldCombo(htmls, metaFieldInfo, id, left) {
        htmls.push('<div id="tpMetaFieldComboBox" style="padding-left:2px">');
        var entries = [["_none", ""]];
        for (var i = 0; i < metaFieldInfo.length; i++) {
            var field = metaFieldInfo[i];
            var fieldName = field.label;
            var hasTooManyVals = (field.diffValCount>100);
            if (hasTooManyVals)
                continue;
            entries.push( ["tpMetaVal_"+i, fieldName] );
        }

        buildComboBox(htmls, id, entries, 0, "select a field...", metaBarWidth);
        htmls.push('</div>');
    }

    function buildGeneCombo(htmls, id, left, width) {
        /* datasets with a list of elements with a shortLabel attribute. Build combobox for them. */
        //htmls.push('<div class="tpToolBarItem" style="position:absolute;left:'+left+'px;top:'+toolBarComboTop+'px">');
        htmls.push('<div class="tpToolBarItem" style="padding-left: 3px">');
        htmls.push('<label style="display:block; margin-bottom:8px; padding-top: 8px;" for="'+id+'">Color by Gene</label>');
        htmls.push('<select style="width:'+width+'px" id="'+id+'" placeholder="search for gene..." class="tpCombo">');
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

    function geneComboSearch(query, callback) {
        /* called when the user types something into the gene box, returns matching gene symbols */
        if (!query.length) 
            return callback();

        //var geneList = [];

        //for (var geneSym in db.getGenes()) {
            //if (geneSym.toLowerCase().startsWith(query.toLowerCase()))
                //geneList.push({"id":geneSym, "text":geneSym});
        //}
        db.searchGenes(query.toLowerCase(), callback);
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


    function makeHubUrl(geneSym) {
        /* return URL of the hub.txt file, possibly jumping to a given gene  */
            var hubUrl = db.conf.hubUrl;
            if (hubUrl==undefined)
                return null;
            var ucscDb = db.conf.ucscDb;
            if (ucscDb===undefined) {
                alert("Internal error: ucscDb is not defined in cellbrowser.conf. Example values: hg19, hg38, mm10, etc. You have to set this variable to make track hubs work.");
                return "";
            }
            //var fullUrl = "https://genome.ucsc.edu/cgi-bin/hgTracks?hubUrl="+hubUrl+"&genome="+ucscDb;
            var fullUrl = "https://genome-test.gi.ucsc.edu/cgi-bin/hgTracks?hubUrl="+hubUrl+"&genome="+ucscDb;

            if (geneSym!==undefined)
                fullUrl += "&position="+geneSym+"&singleSearch=knownCanonical";

            return fullUrl;
    }

    function buildToolBar (coordInfo, datasetName, fromLeft, fromTop) {
    /* add the tool bar with icons of tools and add under body to the DOM */
        $("#tpToolBar").remove();

        var htmls = [];
        htmls.push("<div id='tpToolBar' style='position:absolute;left:"+fromLeft+"px;top:"+fromTop+"px'>");
        //htmls.push('<div id="tpIcons" style="display:inline-block">');
        //htmls.push('<div class="btn-group" role="group" style="vertical-align:top">');
        //htmls.push('<button data-placement="bottom" data-toggle="tooltip" title="Zoom-to-rectangle mode.<br>Keyboard: Windows/Command or z" id="tpIconModeZoom" class="ui-button tpIconButton" style="margin-right:0"><img src="img/zoom.png"></button>');
        //htmls.push('<button data-placement="bottom" title="Move mode. Keyboard: Alt or m" id="tpIconModeMove" data-toggle="tooltip" class="ui-button tpIconButton" style="margin-right:0"><img src="img/move.png"></button>');
        //htmls.push('<button data-placement="bottom" title="Select mode.<br>Keyboard: shift or s" id="tpIconModeSelect" class="ui-button tpIconButton" style="margin-right:0"><img src="img/select.png"></button>');
        //htmls.push('</div>');

        //htmls.push('&emsp;');
        //htmls.push('<button title="Zoom in" id="tpIconZoomIn" type="button" class="btn-small btn-outline-primary noPad"><i class="material-icons">zoom_in</i></button>');
        //htmls.push('<button title="Zoom out" id="tpIconZoomOut" type="button" class="btn-small btn-outline-primary noPad"><i class="material-icons">zoom_out</i></button>');
        //htmls.push('<button title="Zoom to 100%, showing all data, keyboard: space" data-placement="bottom" data-toggle="tooltip" id="tpZoom100Button" class="ui-button tpIconButton" style="margin-right:0"><img src="img/center.png"></button>');

        //htmls.push("&emsp;");

        //htmls.push('<div class="btn-group" role="group" style="vertical-align:top">');
        //htmls.push('<button title="More info about this dataset" id="tpIconDatasetInfo" type="button" class="ui-button tpIconButton"><img title="More info about this dataset" src="img/info.png"></button>');
        //htmls.push('</div>');
        //htmls.push('<img class="tpIconButton" id="tpIconDatasetInfo" data-placement="bottom" data-toggle="tooltip" title="More info about this dataset" src="img/info.png" style="height:18px;position:absolute;top:4px; left:'+(toolBarComboLeft+datasetComboWidth+60)+'px">');

        //htmls.push("&emsp;");
        buildLayoutCombo(htmls, coordInfo, "tpLayoutCombo", 250, 290, 2);
        //buildDatasetCombo(htmls, gDatasetList, "tpDatasetCombo", 100, 220, 0);
        
        htmls.push('<button id="tpOpenDatasetButton" class="gradientBackground ui-button ui-widget ui-corner-all" style="margin-top:3px; height: 24px; border-radius:3px; padding-top:3px">Open Dataset...</button>');

        var hubUrl = db.conf.hubUrl;
        if (hubUrl!==undefined) {
            var fullUrl = makeHubUrl();
            htmls.push('<a target=_blank href="'+fullUrl+'" id="tpOpenUcsc" class="gradientBackground ui-button ui-widget ui-corner-all" style="margin-left: 10px; margin-top:3px; height: 24px; border-radius:3px; padding-top:3px">Genome Browser</a>');
        }

        htmls.push("</div>");

        $(document.body).append(htmls.join(""));

        //var el = document.getElementById('tpOpenUcsc');
        //el.addEventListener("click"

        activateTooltip('.tpIconButton');

        //$('#tpIconModeMove').click( function() { activateMode("move")} );
        //$('#tpIconModeZoom').click( function() { activateMode("zoom")} );  
        //$('#tpIconModeSelect').click( function() { activateMode("select")} );
        //$('#tpZoom100Button').click( onZoom100Click );
        $('#tpIconDatasetInfo').click( function() { openDatasetDialog()});

        //$('#tpIconZoomIn').click( onZoomInClick );
        //$('#tpIconZoomOut').click( onZoomOutClick );
        //$('#tpIconZoom100').click( onZoom100Click );

        activateCombobox("tpDatasetCombo", datasetComboWidth);
        activateCombobox("tpLayoutCombo", layoutComboWidth);

        //$('#tpGeneCombo').select2({ placeholder : "Gene Symbol"}); // preliminary setup, real setup will be in updateToolbar()
        $('#tpGeneCombo').selectize({
                labelField : 'text',
                valueField : 'id',
                searchField : 'text',
                load : geneComboSearch
        });

        $('#tpDatasetCombo').change(onDatasetChange);
        // update the combobox, select the right dataset
        var datasetIdx = cbUtil.findIdxWhereEq(gDatasetList, "name", datasetName);
        $("#tpDatasetCombo").val(datasetIdx).trigger("chosen:updated");
        $('#tpLayoutCombo').change(onLayoutChange);
        $('#tpGeneCombo').change(onGeneChange);
        $('#tpOpenDatasetButton').click(openDatasetDialog);
    }

    function metaFieldToLabel(fieldName) {
    /* convert the internal meta field string to the label shown in the UI. Fix underscores, _id, etc */
        if (fieldName==="_id")
            fieldName = capitalize(gSampleDesc)+" identifier";
        else
            fieldName = fieldName.replace(/_/g, " ");
        return fieldName;
    }

    function buildMetaPanel(htmls, metaFieldInfo) {
        htmls.push("<div id='tpMetaPanel'>");
        for (var i = 0; i < metaFieldInfo.length; i++) {
            var field = metaFieldInfo[i];
            var fieldName = field.label;

            // fields without binning and with too many unique values are greyed out
            var isGrey = (field.diffValCount>100 && field.binMethod!==undefined);

            var addClass = "";
            var addTitle="";
            htmls.push("<div class='tpMetaBox' id='tpMetaBox_"+i+"'>");
            if (isGrey) {
                addClass=" tpMetaLabelGrey";
                addTitle=" title='This field contains too many different values. You cannot click it to color on it.'";
            }
            htmls.push("<div id='tpMetaLabel_"+i+"' class='tpMetaLabel"+addClass+"'"+addTitle+">"+fieldName+"</div>");

            var styleAdd="";
            if (field.opt!==undefined) {
                var opt = field.opt;
                if (opt.fontSize!==undefined)
                    styleAdd = ";font-size:"+field.opt.fontSize;
            }

            htmls.push("<div class='tpMetaValue' style='width:"+(metaBarWidth-2*metaBarMargin)+"px"+styleAdd+"' id='tpMeta_"+i+"'>&nbsp;</div>");
            htmls.push("</div>"); // tpMetaBox
        }
        htmls.push("<div style='background-color:white; float:right' id='tpMetaNote' style='display:none; height:1em'></div>");
        htmls.push("</div>"); // tpMetaPanel
    }

    function buildLeftSidebar (metaFieldInfo) {
    /* add the left sidebar with the meta data fields. db.loadConf 
     * must have completed before this can be run, we need the meta field info. */
        $("#tpLeftSidebar").remove();
        // setup the tabs
        var tabsWidth = metaBarWidth;

        var htmls = [];
        htmls.push("<div id='tpLeftSidebar' style='position:absolute;left:0px;top:"+menuBarHeight+"px;width:"+metaBarWidth+"px'>");

        htmls.push("<div class='tpSidebarHeader'>Color Control</div>");

        // a bar with the tabs
        htmls.push("<div id='tpLeftTabs'>");
        htmls.push("<ul>");
        htmls.push("<li><a href='#tpAnnotTab'>Cell Annotations</a></li>");
        htmls.push("<li><a href='#tpGeneTab'>Genes</a></li>");
        htmls.push("</ul>");

        htmls.push("<div id='tpAnnotTab'>");
        //htmls.push("<div id='tpSideMeta'>");
        htmls.push('<label style="padding-left: 2px; margin-bottom:8px; padding-top:8px" for="'+"tpMetaCombo"+'">Color by Annotation</label>');
        buildMetaFieldCombo(htmls, metaFieldInfo, "tpMetaCombo", 0);
        htmls.push('<div style="padding-top:4px; padding-bottom: 4px; padding-left:2px" class="tpHint">Mouse-over a '+gSampleDesc+' to show it below</div>');
        buildMetaPanel(htmls, metaFieldInfo);
        //htmls.push("</div>"); // tpSideMeta
        htmls.push("</div>"); // tpAnnotTab

        htmls.push("<div id='tpGeneTab'>");

        buildGeneCombo(htmls, "tpGeneCombo", 0, metaBarWidth-10);
        buildGeneTable(htmls, metaBarWidth, 66, db.conf.quickGenes);

        htmls.push("</div>"); // tpGeneTab

        htmls.push("</div>"); // tpLeftSidebar

        $(document.body).append(htmls.join(""));

        $("#tpLeftTabs").tabs();
        $('#tpLeftTabs').tabs("option", "active", 0); // open the first tab

        $('.tpGeneBarCell').click( onGeneClick );
        $('#tpChangeGenes').click( onChangeGenesClick );
        //$("#tpGenePaneLink").tab("show");
        //$("#tpMetaPaneLink").tab("show");

        /* htmls.push("<div id='tpMetaPanes' style='margin-top:2px; margin-left:1px'>");

                htmls.push("<div id='tpAnnotPane' class='tab-pane'>");
                buildMetaPanel(htmls, metaFieldInfo);
                htmls.push("</div>");

                htmls.push("<div id='tpGenePane' class='tab-pane' style='height:200px'>");
                buildGeneCombo(htmls, "tpGeneCombo", 0);
                htmls.push("</div>");

            htmls.push("</div>"); // tab-content

        htmls.push("</div>"); // opendialogtabs

        htmls.push("</div>"); // tpLeftSidebar

        $(document.body).append(htmls.join(""));
        $(document.body).append("<div id='tpMetaTip'></div>");

        // this is weird, but I have not found a better way to make the tab show up

        //clearMetaAndGene();

        */
        activateCombobox("tpMetaCombo", metaBarWidth-10);
        $("#tpMetaCombo").change( onMetaComboChange );
        $(".tpMetaLabel").click( onMetaClick );
        $(".tpMetaValue").click( onMetaClick );
        //$(".tpMetaValue").mouseover( onMetaMouseOver );
        $(".tpMetaValue").mouseleave ( function() { $('#tpMetaTip').hide()} );
        
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
        // setup the tooltips
        //$('[title!=""]').tooltip();
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
                intro: "Click into the legene to select "+gSampleDesc+"s.<br>Click a color to change it or select a palette from the 'Colors' menu.<br>To setup your own cell browser, see 'Help - Setup your own'",
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
    function hslToRgb(h, s, l) {
      var r, g, b;

      if (s == 0) {
	r = g = b = l; // achromatic
      } else {
	function hue2rgb(p, q, t) {
	  if (t < 0) t += 1;
	  if (t > 1) t -= 1;
	  if (t < 1/6) return p + (q - p) * 6 * t;
	  if (t < 1/2) return q;
	  if (t < 2/3) return p + (q - p) * (2/3 - t) * 6;
	  return p;
	}

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
        
    function makeColorPalette(palName, n) {
    /* return an array with n color hex strings */
    // Use Google's palette functions for now, first Paul Tol's colors, if that fails, use the usual HSV rainbow
        var pal = [];
        if (palName==="blues")
            pal = makeHslPalette(0.6, n);
        else if (palName==="reds")
            pal = makeHslPalette(0.0, n);
        else
            pal = palette(palName, n);

        return pal;
    }

    function colorByCluster() {
    /* called when meta and coordinates have been loaded: scale data and color by meta field  */
        //setZoomRange();
    }

    function startLoadTsv(fullUrl, func, addInfo) {
    /* load a tsv file relative to baseUrl and call a function when done */
        function conversionDone(data) {
            Papa.parse(data, {
                    complete: function(results, localFile) {
                                func(results, localFile, addInfo);
                            },
                    error: function(err, file) {
                                if (addInfo!==undefined)
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
        Mousetrap.bind('o', openDatasetDialog);
        Mousetrap.bind('c m', onMarkClearClick);
        Mousetrap.bind('h m', onMarkClick);

        Mousetrap.bind('space', onZoom100Click);

        Mousetrap.bind('z', function() { activateMode("zoom"); });
        Mousetrap.bind('m', function() { activateMode("move"); });
        Mousetrap.bind('s', function() { activateMode("select"); });

        Mousetrap.bind('-', onZoomOutClick);
        Mousetrap.bind('+', onZoomInClick);
        Mousetrap.bind('n', onSelectNoneClick);
        Mousetrap.bind('a', onSelectAllClick);
        Mousetrap.bind('d', function() {$('#tpDatasetCombo').trigger("chosen:open"); return false;});
        Mousetrap.bind('l', function() {$('#tpLayoutCombo').trigger("chosen:open"); return false;});
        Mousetrap.bind('g', function() {$("#tpGeneCombo").selectize()[0].selectize.focus(); return false;});
        Mousetrap.bind('c l', onHideShowLabelsClick );

    }

    // https://stackoverflow.com/a/33861088/233871
    function isInt(x) {
        return (typeof x==='number') && x % 1 == 0;
    }

    function naturalSort (a, b) {
    /* copied from https://github.com/Bill4Time/javascript-natural-sort/blob/master/naturalSort.js */
    "use strict";
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
};

    function sortPairsBy(countList, sortBy) {
    /* sort an array in the format [name, count] by either name (using naturalSort) or count */
    // XX no need anymore to return a dict, just return the list
        // convert the dict to a list of (count, key)
        //var countList = [];
        //var numCount = 0;
        //var count = null;
        //var useGradient = false; // use a rainbow or gradient palette?
        //for (var key in dict) {
            //count = dict[key];
            //if (!isNaN(key)) // key looks like a number
                //numCount++;
            //countList.push( [count, key] );
        //}

        // auto-detect: sort list by name if most names are numbers
        //if ((countList.length >= 4 && numCount >= countList.length-1)) {
            //if (sortBy===undefined)
                //sortBy = "name";
            //useGradient = true;
        //}

        var isSortedByName = null;

        if (sortBy==="name") {
            countList.sort(function(a, b) { return naturalSort(a[0], b[0]); });  // I have a feeling that this is very slow...
            isSortedByName = true;
        }
        else {
            // sort this list by count
            countList = countList.sort(function(a, b){ return b[1] - a[1]; }); // reverse-sort by count
            isSortedByName = false;
        }

        var ret = {};
        ret.list = countList;
        ret.isSortedByName = isSortedByName;
        // pallette should be a gradient for data types where this makes sense
        return ret;
    }

    function assignColors(fieldIdx, countList, doReset) {
    /* given an array of (count, value), assign a color to every value.
     Returns an dict of value : (color as integer, color as hex string without #) 
     Will check if user has a manually defined color for this field before and use it, if present.
     If doReset is true, will delete of manual definitions.
     */
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

    function findDotsWithMeta(metaIdx, hlValue) {
    /* return array of cellIds with a given meta value */
        var metaData = gCurrentDataset.metaData;
        ret = [];
        for (var i = 0; i < pixelCoords.length; i++) {
            var cellId = pixelCoords[i][0];
            var x = pixelCoords[i][1];
            var y = pixelCoords[i][2];
            var metaRow = metaData[cellId];
            var metaVal = metaRow[metaIdx];
            if (metaVal === hlValue) {
                ret.push(cellId);
            }
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

    function onLegendLabelClick(ev) {
    /* called when user clicks on legend entry. */
        var legendId = parseInt(ev.target.id.split("_")[1]);
        if (("lastClicked" in gLegend) && gLegend.lastClicked==legendId) {
            // user clicked the same entry as before: clear selection
            gLegend.lastClicked = null;
            renderer.selectClear();
            $('#tpLegend_'+legendId).removeClass('tpLegendSelect');
            menuBarHide("#tpFilterButton");
            menuBarHide("#tpOnlySelectedButton");
            updateGeneTableColors(null);
        }
        else {
            if (!ev.shiftKey && !ev.ctrlKey && !ev.metaKey) {
                renderer.selectClear();
                $('.tpLegend').removeClass('tpLegendSelect');
            }
            console.log(ev);
            var colIdx = gLegend.rows[legendId][4];
            //gSelCellIds = findCellIdsForLegendIds(gClasses, [legendId], gSelCellIds);
            renderer.selectByColor(colIdx);
            menuBarShow("#tpFilterButton");
            menuBarShow("#tpOnlySelectedButton");
            $('#tpLegend_'+legendId).addClass('tpLegendSelect');
            gLegend.lastClicked=legendId;
            //updateSelection();
        }
        renderer.drawDots();
    }

    function legendResetColors () {
    /* reset all manual colors in legend to defaults */
        var rows = gLegend.rows;
        if (rows===undefined)
            return;
        for (var i = 0; i < rows.length; i++) {
            var colorHex = rows[i][0];
            var defColor = rows[i][1];
            gLegend.rows[i][0] = null;
            var saveKey = rows[i][5];
            // also clear localStorage
            localStorage.removeItem(saveKey);
        }
    }

    function onSortByClick (ev) {
    /* flip the current legend sorting */
        var sortBy = null;
        if (ev.target.id.endsWith("Col1")) // column 1 is the Name
            sortBy = "name"
        else
            sortBy = "freq";

        cartSave("s_"+gLegend.fieldName,sortBy, gLegend.defaultSortBy);
        legendSort(sortBy);
        buildLegendBar();
    }


    function onMetaRightClick (key, options) {
    /* user right-clicks on a meta field */
        var metaIdx = parseInt(options.$trigger[0].id.split("_")[1]);

        //if (key==0) {
            //gCurrentDataset.labelField = gCurrentDataset.metaFields[metaIdx];
            //gClusterMids = null; // force recalc
            //plotDots();
            //renderer.render(stage);
            //updateMenu();
        //}
        if (key==0) {
            copyToClipboard("#tpMeta_"+metaIdx);
            //$("textarea").select();
            //document.execCommand('copy');
            //console.log(val);
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
            $('#tpLegendCol1').html('Name<span class="caret"></span>');
            $('#tpLegendCol2').html('Frequency<span class="caret"></span>');
        }
        else {
            $('#tpLegendCol1').html('Range');
            $('#tpLegendCol2').html('Frequency');
        }
    }

    function buildLegendBar(sortBy) {
    /* draws current legend as specified by gLegend.rows 
     * sortBy can be "name" or "count" or "undefined" (=auto-detect)
     * */
        if (gLegend.rows==undefined)
            return;

        $('#tpLegendContent').empty();

        var htmls = [];

        var colors = [];
        var rows = gLegend.rows;

        // add the "sort by" div
        //var sortLabel = null;
        //if (gLegend.isSortedByName===true)
            //sortLabel = "Sort by freq." // add the link to sort by the other possibility
        //else
            //sortLabel = "Sort by name"
            //$('#tpLegendTitle').html('<span id="tpGeneSym" title="' +geneId+'">'+geneSym+" expression</span>");
        //htmls.push("<span id='tpSortBy' style='color: #888; cursor:pointer; font-size:13px'>"+sortLabel+"</span>");
        //htmls.push("<div class='btn-group btn-group-xs' role='group'>");
        //htmls.push("<button type='button' class='btn btn-default' id='tpSortBy'>"+sortLabel+"</button>");
        //htmls.push("</div>");
        htmls.push('<span id="tpLegendTitle" title="' +gLegend.titleHover+'">'+gLegend.title+"</span>");
        htmls.push('<div class="tpHint">Click below to select '+gSampleDesc+'s</small></div>');
        htmls.push("</div>"); // title
        htmls.push('<div id="tpLegendHeader"><span id="tpLegendCol1"></span><span id="tpLegendCol2"></span></div>');

        // get the sum of all, to calculate frequency
        var sum = 0;
        for (var i = 0; i < rows.length; i++) {
            var count = rows[i][3];
            sum += count;
        }

        var acronyms = db.conf.acronyms;
        if (acronyms===undefined)
            acronyms = {};

        for (var i = 0; i < rows.length; i++) {
            var row = rows[i];
            var colorHex = row[0]; // manual color
            if (colorHex===null)
                colorHex = row[1]; // default color

            var label = row[2];
            var count = row[3];
            var freq  = 100*count/sum;

            colors.push(colorHex); // save for later

            var labelClass = "tpLegendLabel";
            label = label.replace(/_/g, " ").replace(/'/g, "&#39;").trim();
            if (label==="") {
                label = "(empty)";
                labelClass += " tpGrey";
            }
            else if (likeEmptyString(label))
                labelClass += " tpGrey";

            var labelDesc = label;
            var labelDesc = acronyms[label] || null;
            if (labelDesc===null) {
                // only show the full value on mouse over if the label is long, "" suppresses mouse over
                if (label.length > 20 || labelDesc===null)
                    labelDesc = label;
                else
                    labelDesc = "";
            }


            var classStr = "tpLegend";
            var line = "<div id='tpLegend_" +i+ "' class='" +classStr+ "'>";
            htmls.push(line);
            htmls.push("<input class='tpColorPicker' id='tpLegendColorPicker_"+i+"' />");

            htmls.push("<span class='"+labelClass+"' id='tpLegendLabel_"+i+"' data-placement='auto top' title='"+labelDesc+"'>");
            htmls.push(label);
            htmls.push("</span>");
            //htmls.push("<span class='tpLegendCount'>"+count+"</div>");
            var prec = 1;            
            if (freq<1)
                prec = 2;
            htmls.push("<span class='tpLegendCount' title='"+count+" of "+sum+"'>"+freq.toFixed(prec)+"%</div>");
            htmls.push("</span>");

            htmls.push("</div>");
            //htmls.push("<input class='tpLegendCheckbox' id='tpLegendCheckbox_"+i+"' type='checkbox' checked style='float:right; margin-right: 5px'>");
        }

        var htmlStr = htmls.join("");
        $('#tpLegendContent').append(htmlStr);
        setLegendHeaders(gLegend.rowType);

        activateTooltip("#tpResetColors");
        activateTooltip("#tpSortBy");
        $("#tpLegendCol1").click( onSortByClick );
        $("#tpLegendCol2").click( onSortByClick );

        $('.tpLegend').click( onLegendLabelClick );
        //$('.tpLegendLabel').attr( "title", "Click to select samples with this value. Shift click to select multiple values.");
        //$('#tpResetColors').click( onResetColorsClick );
        //$('#tpSortBy').click( onSortByClick );
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
        for (var i = 0; i < colors.length; i++) {
            var opt = {
                hideAfterPaletteSelect : true,
                color : colors[i],
                showPalette: true,
                //allowEmpty : true,
                showInput: true,
                preferredFormat: "hex",
                change: onColorPickerChange
                }
            $("#tpLegendColorPicker_"+i).spectrum(opt);
        }

    }

    function onColorPickerChange(color) {
        /* called when user manually selects a color in the legend with the color picker */
        var valueIdx = parseInt(this.id.split("_")[1]);
        var rows = gLegend.rows;
        var clickedRow = rows[valueIdx];
        var oldColorHex = clickedRow[0];
        var defColorHex = clickedRow[1];

        var newCol = $(this).spectrum('get');

        var newColHex = "";
        if (newCol===null)
            newColHex = oldColorHex; // if user clicked abort, revert to default color
        else
            newColHex = newCol.toHex();
        clickedRow[0] = newColHex;

        var saveKey = COL_PREFIX+clickedRow[5];
        // save color to cart if necessary
        cartSave(saveKey, newColHex, defColorHex);

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

        var pal = makeColorPalette(cDefGradPalette, exprBinCount);

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
            if (cellIds!==null && cellIds.length!=0) {
                for (var ci=0; ci<cellIds.length; ci++) {
                    sum += vec[cellIds[ci]];
                }
                avg = Math.round(sum / cellIds.length);
                //console.log("sum "+sum+" avg "+avg);
            }
            var color = pal[avg];
            //console.log("color "+color);
            $("#tpGeneBarCell_"+i).css("background-color", "#"+color);
	    var fontColor = "#333333";
	    if (isDark(color))
		fontColor = "white";
	    $("#tpGeneBarCell_"+i).css("color", fontColor);
        }
        console.timeEnd("avgCalc");
    }

    function countValues(arr) {
        /* count values in array, return an array of [value, count], sorted by count */
        // first make a list of all possible meta values and their counts
        var metaData = gCurrentDataset.metaData;
        var metaCounts = {};
        for (var i = 0; i < pixelCoords.length; i++) {
            var cellId = pixelCoords[i][0];
            var metaRow = metaData[cellId];
            var metaVal = metaRow[metaIndex];
            metaCounts[metaVal] = 1 + (metaCounts[metaVal] || 0);
        }

        var counts = {};
        for (var i = 0; i < arr.length; i++) {
            var val = arr[i];
            counts[val] = (counts[num] || 0) + 1;
            var metaVal = metaRow[metaIndex];
        }
    }

    function updateMetaBarManyCells(cellIds) {
    /* update the meta fields on the left to reflect/summarize a list of cellIds */
        var metaFields = gCurrentDataset.metaFields;
        var metaData = gCurrentDataset.metaData;
        var cellCount = cellIds.length;
        $('#tpMetaTitle').text("Meta data of "+cellCount+" "+gSampleDesc+"s");

        // for every field...
        var metaHist = {};
        for (var metaIdx = 0; metaIdx < metaFields.length; metaIdx++) {
            var metaCounts = {};
            // make an object of value -> count in the cells
            for (var i = 0; i < cellCount; i++) {
                var cellId = cellIds[i];
                var metaVal = metaData[cellId][metaIdx];
                metaCounts[metaVal] = 1 + (metaCounts[metaVal] || 0);
            }
            // convert the object to an array (count, percent, value) and sort it by count
            var histoList = [];
            for (var key in metaCounts) {
                var count = metaCounts[key];
                var frac  = (count / cellCount);
                histoList.push( [count, frac, key] );
            }
            histoList = histoList.sort(function(a, b){ return b[0] - a[0]; }); // reverse-sort by count
            metaHist[metaIdx] = histoList;

            // make a quick list of the top values for the sparklines, ignore the rest
            var countList = [];
            var otherCount = 0;
            for (var i=0; i < histoList.length; i++) {
                var count = histoList[i][0];
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
            var topVal   = histoList[0][2];
            var percStr = (100*topPerc).toFixed(1)+"%";

            if (topVal.length > 14)
                topVal = topVal.substring(0, 14)+"...";

            var label = "";
            if (histoList.length!==1){
                if (histoList[0][0]===1)
                    label = "<span class='tpMetaMultiVal'>" + histoList.length + " unique values</span>";
                else
                //label = "<span class='tpMetaMultiVal'>" + histoList.length + " values</span><span class='tpMetaHistLabel'> Histogram </span>&nbsp;&nbsp;<span id='tpMetaSpark_"+metaIdx+"'></span>";
                //label = "<span class='tpMetaMultiVal'>" + histoList.length + " values</span>&nbsp;<span id='tpMetaSpark_"+metaIdx+"'></span>";
                label = "<span class='tpMetaMultiVal'>" + percStr + " " +topVal+"</span>&nbsp;<span class='tpMetaSpark' id='tpMetaSpark_"+metaIdx+"'></span>";
            }
            else
                label = topVal;

            $('#tpMeta_'+metaIdx).html(label);
            $('#tpMetaSpark_'+metaIdx).sparkline(countList, {type:"bar", barSpacing:0, disableTooltips:true});

            //var topCount = countList[0][0];
            //var topPerc  = countList[0][1];
            //var topVal   = countList[0][2];
            //$('#tpMeta_'+metaIdx).text(topVal);
            //var label = metaFieldToLabel(metaFields[metaIdx]);
            //$('#tpMetaLabel_'+metaIdx).html(label+"<span class='tpMetaPerc'>"+(100*topPerc).toFixed(1)+"%</span>");
        }
        gCurrentDataset.metaHist = metaHist;
    }

    function clearMetaAndGene() {
        /* clear the meta and gene field field info */
        var fieldCount = db.getMetaFields().length;
        for (var i = 0; i < fieldCount; i++) {
            $('#tpMeta_'+i).attr('title', "");
            $('#tpMeta_'+i).html("");
        }
        updateGeneTableColors(null);
    }

    function updateMetaBarOneCell(cellInfo) {
        /* update the meta bar with meta data from a single cellId */
        $('#tpMetaTitle').text(METABOXTITLE);
        for (var i = 0; i < cellInfo.length; i++) {
            var fieldValue = cellInfo[i];
            if (fieldValue.startsWith("http") && fieldValue.endsWith(".png")) {
                $('#tpMeta_'+i).css('height', "40px");
                $('#tpMeta_'+i).html("<img src='"+fieldValue+"'></img>");
            } else
                $('#tpMeta_'+i).html(fieldValue);
            $('#tpMeta_'+i).attr('title', cellInfo[i]);
        }
    }

    function onDotMouseOver (mouseData) {
        /* user hovers over a circle with the mouse */
        if (mouseDownX!=null) // user is currently zooming
            return;
        if (! (gSelCellIds===null || jQuery.isEmptyObject(gSelCellIds))) // some cells are selected
            return;
        var cellId = mouseData.target.cellId;
        this.alpha = 1.0;
        //renderer.render(stage);

        //db.loadMetaForCell(cellId, function(ci) { updateMetaBarOneCell(ci);} onProgress);
        updateGeneTableColors([cellId]);
    }

    function onCellClickOrHover (cellIds, ev) {
        /* user clicks onto a circle with the mouse or hovers over one. 
         * ev is undefined if not a click. */

        // do nothing if only hover but we already have a selection
        var selCells = renderer.getSelection();
        if (ev===undefined && selCells!==null && selCells.length!==0)
            return;

        if (cellIds===null || cellIds.length===0)
            clearMetaAndGene();
        else {
            var cellId = cellIds[0];
            db.loadMetaForCell(cellId, function(ci) { updateMetaBarOneCell(ci);}, onProgress);
            if (cellIds.length===1)
                $("#tpMetaNote").hide();
            else {
                $("#tpMetaNote").html("...and "+(cellIds.length-1)+" other "+gSampleDesc+"s underneath");
                $("#tpMetaNote").show();
            }
        }

        updateGeneTableColors(cellIds);

        if (ev!==undefined) {
            // it was a click
            if (!ev.shiftKey && !ev.ctrlKey && !ev.metaKey)
                renderer.selectClear(); 
            renderer.selectAdd(cellId);
            renderer.drawDots();
            event.stopPropagation();
        }
    }

    function onClusterNameHover(clusterName, ev) {
       /* user hovers over cluster label */
       //var htmls = [];
       //htmls.push("<div class='tpHover'>"+clusterName+"</div>");
       //$(document.body).append(htmls.join(""));
       console.log("Hover over "+clusterName);
    }

    function sanitizeName(name) {
        /* ported from cellbrowser.py: remove non-alpha, allow underscores */
        var newName = name.replace(/[^a-zA-Z_0-9+]/g, "");
        return newName;
    }

    function onClusterNameClick(clusterName) {
        /* build and open the dialog with the marker genes table for a given cluster */
        console.log("building marker genes window for "+clusterName);
        var tabInfo = db.conf.markers; // list with (label, subdirectory)
        if (tabInfo===undefined)
            return;
        var doTabs = (tabInfo != undefined && tabInfo.length>1);

        var htmls = [];

        htmls.push("<div id='tpPaneHeader' style='padding:8px'>");
        //var hubUrl = db.conf.hubUrl;
        //if (hubUrl!==undefined) {
            //htmls.push("<p>");
            //htmls.push("<a target=_blank class='link' href='"+hubUrl+"'>Show Sequencing Reads on UCSC Genome Browser</a><p>");
        //}

        htmls.push("Click gene symbols below to color plot by gene<br>");

        htmls.push("</div>");

        if (doTabs) {
            htmls.push("<div id='tabs'>");
            htmls.push("<ul>");
            for (var tabIdx = 0; tabIdx < tabInfo.length; tabIdx++) {
                var tabLabel = tabInfo[tabIdx].shortLabel;
                htmls.push("<li><a href='#tabs-"+tabIdx+"'>"+tabLabel+"</a>");
            }
            htmls.push("</ul>");
        }

        for (var tabIdx = 0; tabIdx < tabInfo.length; tabIdx++) {
            var divName = "tabs-"+tabIdx;
            var tabDir = tabInfo[tabIdx].name;
            var sanName = sanitizeName(clusterName);
            var markerTsvUrl = joinPaths([db.name, "markers", tabDir, sanName+".tsv.gz"]);
            htmls.push("<div id='"+divName+"'>");
            htmls.push("Loading...");
            htmls.push("</div>");

            startLoadTsv(markerTsvUrl, loadMarkersFromTsv, divName);
        }

        htmls.push("</div>"); // tabs

        var buttons = {
        "Download as file" :
            function() {
                //url = joinPaths([baseUrl,"geneMatrix.tsv"]);
                document.location.href = markerTsvUrl;
            },
        };

        var winWidth = window.innerWidth - 0.10*window.innerWidth;
        var winHeight = window.innerHeight - 0.10*window.innerHeight;
        var title = "Cluster markers for &quot;"+clusterName+"&quot;";
        var acronyms = db.conf.acronyms;
        if (acronyms!==undefined && clusterName in acronyms)
            title += " - "+acronyms[clusterName];
        showDialogBox(htmls, title, {width: winWidth, height:winHeight, "buttons":buttons});
        //$(".tpLoadGeneLink").on("click", onMarkerGeneClick);
        //activateTooltip(".link");
        $(".ui-widget-content").css("padding", "0");
        $("#tabs").tabs();
        //$("table").focus();
        //removeFocus();
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

    function onMarkerGeneClick(ev) {
        /* user clicks onto a gene in the table of the marker gene dialog window */
        var geneSym = ev.target.getAttribute("data-gene");
        $(".ui-dialog").remove(); // close marker dialog box
        //loadSingleGeneFromMatrix(geneSym);
        loadGeneAndColor(geneSym);
    }

    function loadMarkersFromTsv(papaResults, url, divId) {
        /* construct a table from a marker tsv file and write as html to the DIV with divID */
        console.log("got coordinate TSV rows, parsing...");
        var rows = papaResults.data;
        var headerRow = rows[0];

        var htmls = [];


        htmls.push("<table class='table' id='tpMarkerTable'>");
        htmls.push("<thead>");
        var hprdCol = null;
        var geneListCol = null;
        var exprCol = null;
        var pValCol = null
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
            else if (colLabel==="P_value" || colLabel==="p_val" || colLabel=="pVal") {
                colLabel = "P-value";
                pValCol = i;
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
            htmls.push(colLabel);
            htmls.push("</th>");
        }
        htmls.push("</thead>");

        var hubUrl = makeHubUrl();

        htmls.push("<tbody>");
        for (var i = 1; i < rows.length; i++) {
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
                    htmls.push(parseFloat(val).toFixed(4)); // four digits ought to be enough for everyone
                else
                    htmls.push(val);
                htmls.push("</td>");
            }
            htmls.push("</tr>");
        }

        htmls.push("</tbody>");
        htmls.push("</table>");

        $("#"+divId).html(htmls.join(""));
        new Tablesort(document.getElementById('tpMarkerTable'));
        $(".tpLoadGeneLink").on("click", onMarkerGeneClick);
        activateTooltip(".link");

        var ttOpt = {"html": true, "animation": false, "delay":{"show":100, "hide":100} }; 
        $(".tpPlots").bsTooltip(ttOpt);

        removeFocus();
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
            if(pair.length != 2) { 
                pair = [pair[0], pair.slice(1).join("=")]
            }
              
            var key = decodeURIComponent(pair[0].replace(plus, " ")), 
                value = decodeURIComponent(pair[1].replace(plus, " ")),
                parts = key.match(keyBreaker);
    
            for ( var j = 0; j < parts.length - 1; j++ ) {
                var part = parts[j];
                if (!current[part] ) {
                    // if what we are pointing to looks like an array
                    current[part] = digitTest.test(parts[j+1]) || parts[j+1] == "[]" ? [] : {}
                }
                current = current[part];
            }
            var lastPart = parts[parts.length - 1];
            if(lastPart == "[]"){
                current.push(value)
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

       if (oldVars===undefined) {
           var queryStr = urlParts[1];
           var urlVars = deparam(queryStr); // parse key=val&... string to object
       } else {
           urlVars = oldVars
       }

       // overwrite everthing that we got
       for (var key in vars) {
           var val = vars[key];
           if (val===null || val in urlVars)
               delete urlVars[key];
           else
               urlVars[key] = val;
       }

       var argStr = jQuery.param(urlVars); // convert to query-like string
       history.pushState({}, db.getName(), baseUrl+"?"+argStr);
    }

    function delState(varName) {
        /* remove a CGI variable from the URL */
        var o = {};
        o[varName] = null;
        changeUrl(o);
    }

    function addStateVar(varName, varVal) {
        /* add a CGI variable from the URL */
        var o = {};
        o[varName] = varVal;
        changeUrl(o);
    }

    function getVar(name, defVal) {
        /* get query variable from current URL or default value if undefined */
       var myUrl = window.location.href;
       myUrl = myUrl.replace("#", "");
       var urlParts = myUrl.split("?");
       var queryStr = urlParts[1];
       var varDict = deparam(queryStr); // parse key=val&... string to object
       if (varDict[name]===undefined)
           return defVal;
       else
           return varDict[name];
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
        if (zs.length!=4)
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

    function extractDatasetFromUrl() {
        /* search for the "ds" parameter or a DNS hostname that indicates the dataset */
        // if ds=xxx was found in the URL, load the respective dataset
        var datasetName = getVar("ds");

        //if (datasetName===undefined)
            //datasetName = datasetList[0].name;
        // hacks for July 2018 and for backwards compatibility with previous version
        if (datasetName==="autism10X" || datasetName==="autism10x")
            datasetName = "autism";
        if (datasetName==="aparna")
            datasetName = "cortex-dev";
        return datasetName;
    }

    /* ==== MAIN ==== ENTRY FUNCTION */
    function loadData(datasetList, globalOpts) {
        /* start the data loaders, show first dataset */
        if (redirectIfSubdomain())
            return;
        gDatasetList = datasetList;

        if (globalOpts!=undefined) {
            if ("sampleType" in globalOpts)
                gSampleDesc = globalOpts["sampleType"];
            if ("title" in globalOpts)
                gTitle = globalOpts["title"];
        }
            
        setupKeyboard();
        buildMenuBar();

        var datasetName = extractDatasetFromUrl(datasetList)

        //menuBarHide("#tpShowAllButton");

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
       }

        buildEmptyLegendBar(metaBarWidth+metaBarMargin+renderer.width, toolBarHeight);

        renderer.setupMouse();
        $(window).resize(onWindowResize);

        renderer.onLabelClick = onClusterNameClick;
        renderer.onLabelHover = onClusterNameHover;
        renderer.onCellClick = onCellClickOrHover;
        renderer.onCellHover = onCellClickOrHover;
        renderer.onNoCellHover = clearMetaAndGene;
        renderer.onZoom100Click = onZoom100Click;
        renderer.onSelChange = onSelChange;

        if (datasetName===undefined)
            openDatasetDialog();
        else 
            loadDataset(datasetName, false);
       
    }

    // only export these functions 
    return {
        "loadData":loadData
    }

}();

function _tpReset() {
/* for debugging: reset the intro setting */
    localStorage.removeItem("introShown");
}
