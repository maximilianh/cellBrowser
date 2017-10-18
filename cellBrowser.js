// A viewer for (x,y) scatter plots 
// shows associated meta data (usually key->val attributes, one mapping per circle) 
// and expression data (string -> float, one mapping per circle)

/* jshint -W097 */
"use strict";

var tsnePlot = function() {
    var gOptions = null; // the whole object given to the loadData function:
    // .datasets: list of objects, with:
    // - .baseUrl (reqd) = URL where files are loaded from
    // - .label (reqd) = label for the dataset shown
    // - .labelField (reqd) = label for the dataset shown
    // - .coordFiles: optional list of (URLOrFilename, label)
    // - .metaColors: options, dict of metaValue -> hexcolor (six chars)

    var gCurrentDataset = null; // pointer to one of gOptions.datasets

    var jsonBaseUrl = null; // current base url, has to end with '/'

    var seuratCoords = null; // array of [cellId, coord1, coord2]. All coords are floats
    var pcaCoords = null;   // array of [cellId, coord1, coord2, coord3, ... All coords are floats.

    var allCoords = null;     // pointer to either seuratCoords or pcaCoords
    var shownCoords = null;   // currently shown float coords of allCoords

    var coord1Idx = 1; // index of X-axis value in the tuples of allCoords
    var coord2Idx = 2; // index of Y-axis value in the tuples of allCoords

    var metaFields = null; // array of meta field labels
    var metaData = null;   // dict cellId -> array of meta fields (strings)
    var gLegendColors = null; // assignment of meta value -> color, specified in config.json

    var matrixOffsets = null; // dict gene symbol -> [offset in matrix, line length in matrix]
    var geneFields = null; // array of (geneSym, geneDesc)
    var geneExpr = null; // dict cellId -> array of gene expr fields (->geneFields, floats)
    var gDeciles = null; // geneId -> Array of 11 values, the ranges for all deciles of a gene
    var gDecileColors = null; // global to avoid recalculting a palette of 10 colors

    var gLabelMetaIdx = null; // meta field to use for labels, default is no label
    var gClusterMids = null; // array of (clusterLabel, x, y), recalculated if null

    // only temporarily used to store data while it is loaded in parallel ajax requests
    var gLoad_geneList = []; // list of genes that were pasted into the gene box
    var gLoad_geneExpr = {}; // map of geneSymbol -> [geneId, array of floats with expression vectors]

    // object with all information needed to map to the legend colors
    var gLegend = null;
    // all info about the current legend. gLegend.rows is:
    // [ colorHexStr, defaultColorHexStr, label, count, internalKey, uniqueKey ]
    // internalKey can be int or str, depending on current coloring mode. 
    // E.g. in meta coloring, it's the metadata string value.
    // When doing expression coloring, it's the expression bin index.
    // uniqueKey is used to save manually defined colors to localStorage

    // mapping from dots to legend entries (and therefore colors)
    // dict with cellId -> index of legend in gLegend.rows or -1 if cell is hidden
    var gClasses = null; 

    // current zoom range: currently displayed min/max range of allCoords
    var zoomRange = null;   // dict of minX, maxX, minY, maxY -> float
    var zoomFullRange = null; // range to show 100%

    // An object of cellIds that are highlighted/selected
    var gSelCellIds = {};

    // the list of cellIds currently used for gene bar coloring
    var gGeneBarCells = null;

    var pixelCoords = null;  // array of [cellId, x, y] with x and y being integers on the screen

    var renderer = null; // the PIXI canvas renderer
    var stage = null; // the PIXI stage
    var winInfo = null; // .width and .height of the PIXI canvas

    // -- coloring-related

    // zooming
    var mouseDownX = null;
    var mouseDownY = null;

    // -- CONSTANTS
    // product name is always a variable, management always has the last word
    var gTitle = "UCSC Cluster Browser";

    // depending on the type of data, single cell or bulk RNA-seq, we call a circle a 
    // "sample" or a "cell". This will adapt help menus, menus, etc.
    var gSampleDesc = "cell";

    // width of left meta bar in pixels
    var metaBarWidth = 200;
    var metaBarMargin = 5;
    // width of legend
    var legendBarWidth = 200;
    var legendBarMargin = 10;
    // height of button bar at top
    var menuBarHeight = null; // defined at runtime by div.height
    // height of toolbar
    var toolBarHeight = 33;

    // height of bottom gene bar
    var geneBarHeight = 100;
    var geneBarMargin = 5;
    // current transparency and circle size
    var transparency = 0.6;
    var circleSize = 4;
    // color for missing value when coloring by expression value
    var cNullColor = "c7c7c7";
    var cNullForeground = "#AAAAAA";

    function _dump(o) {
    /* for debugging */
        console.log(JSON.stringify(o));
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
     if (jQuery.isEmptyObject(gSelCellIds)) {
         menuBarHide("#tpOnlySelectedButton");
         menuBarHide("#tpFilterButton");
     }
     else {
         menuBarShow("#tpOnlySelectedButton");
         menuBarShow("#tpFilterButton");
     }
     
     // the "show all" menu entry is only shown if some dots are actually hidden
     if ((pixelCoords!==null && shownCoords!==null) && pixelCoords.length===shownCoords.length)
         menuBarHide("#tpShowAllButton");
     else
         menuBarShow("#tpShowAllButton");

     $("#tpTrans"+(transparency*100)).addClass("active");
     $("#tpSize"+circleSize).addClass("active");

     // the "hide labels" menu entry is only shown if labels are visible
     if (gLabelMetaIdx === null)
         menuBarHide("#tpHideLabels");
    }

    function onOpenDatasetLink() {
    /* user clicks on File - Open Dataset */
        var htmls = [];

        htmls.push("<div class='list-group' style='width:300px'>");
        for (var i = 0; i < gOptions.datasets.length; i++) {
            var dataset = gOptions.datasets[i];
            var line = "<button type='button' class='list-group-item' data-datasetid='"+i+"'>"; // bootstrap seems to remove the id
            htmls.push(line);
            if (dataset.sampleCount!=undefined) {
                line = "<span class='badge'>"+dataset.sampleCount+" cells</span>";
                htmls.push(line);
            }
            htmls.push(dataset.label+"</button>");
        }
        htmls.push("</div>"); // list-group
        htmls.push("<div id='tpOpenDialogLabel' style='width:450px; position:absolute; left: 340px; top: 10px;'>");
        htmls.push("</div>");
        //htmls.push("<div id='tpSelectedId' data-selectedid='0'>"); // store the currently selected datasetId in the DOM
        var selDatasetIdx = 0;

        var buttons = {
        "Cancel" :
            function() {
                $( this ).dialog( "close" );
            },
        "Open Dataset" :
            function(event) {
                loadDataset(selDatasetIdx);
                $( this ).dialog( "close" );
            }
        };

        showDialogBox(htmls, "Open Dataset", {width: 800, height:500, "buttons":buttons});
        $("#tpOpenDialogLabel").html(gOptions.datasets[0].longLabel);
        $("button.list-group-item").eq(selDatasetIdx).css("z-index", "1000"); // fix up first overlap

        $(".list-group-item").click( function (ev) {
            /* user selects a dataset from the list */
            selDatasetIdx = parseInt($(event.target).data('datasetid')); // index of clicked dataset
            $("#tpOpenDialogLabel").html(gOptions.datasets[selDatasetIdx].longLabel);
            // bootstrap has a bug where the blue selection frame is hidden by neighboring buttons
            // we're working around this here by bumping up the current z-index.
            $("button.list-group-item").css("z-index", "0");
            $("button.list-group-item").eq(selDatasetIdx).css("z-index", "1000");
        });
    }

    function drawLayoutMenu() {
       /* Add the View / Layout menu to the DOM */
      var htmls = [];  
      var coordFiles = gCurrentDataset.coordFiles;
      for (var i=0; i<coordFiles.length; i++) {
            htmls.push('<li><a href="#" class="tpDataset" id="'+i+'">'+coordFiles[i][0]+'</a></li>');
      }
      $("#tpLayoutMenu").empty();
      $("#tpLayoutMenu").append(htmls.join(""));
    }

    function updateSelection() {
    /* has to be called each time the selection has been changed */
        plotDots();
        renderer.render(stage);
        updateMenu();
        var cellIds = keys(gSelCellIds);
        updateGeneBarColors(cellIds);
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
        var url = none;
        if (name==="matrix") {
            url = joinPaths([baseUrl,"geneMatrix.tsv"]);
            document.location.href = url;
        }
        if (name==="meta") {
            url = joinPaths([baseUrl,"meta.tsv"]);
            document.location.href = url;
        }
    }

    function onSelectAllClick() {
    /* Edit - select all */
        gSelCellIds = {};
        for (var i=0; i<shownCoords.length; i++) {
            var cellId = shownCoords[i][0];
            gSelCellIds[cellId] = true;
        }
        updateSelection();
    }

    function onSelectNoneClick() {
    /* Edit - Select None */
        gSelCellIds = {};
        updateSelection();
    }

    function onExportIdsClick() {
        /* Edit - Export cell IDs */
        var idList = [];
        for (var cellId in gSelCellIds) {
            idList.push(cellId);
        }

        var dlgHeight = 500;

        var htmls = [];
        if (idList.length==0)
            {
            //showDialogBox(["Please select cells first, e.g. by clicking into the legend"], "No selection", {showOk:true});
            htmls.push("No cells are selected. Shown below are the identifiers of all cells visible on the screen.<p>");
            for (var i = 0; i < shownCoords.length; i++)
                idList.push(shownCoords[i][0]);
            }

        var idListEnc = encodeURIComponent(idList.join("\n"));
        htmls.push("<textarea style='height:320px;width:250px;display:block'>");
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

    function drawMenuBar() {
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
             htmls.push('<li><a href="#" id="tpOpenDatasetLink"><span class="dropmenu-item-label">Open Dataset...</span><span class="dropmenu-item-content">O</span></a></li>');
             htmls.push('<li class="dropdown-submenu"><a tabindex="0" href="#">Download Data</a>');
               htmls.push('<ul class="dropdown-menu" id="tpDownloadMenu">');
                 htmls.push('<li><a href="#" id="tpDownload_matrix">Gene Expression Matrix</a></li>');
                 htmls.push('<li><a href="#" id="tpDownload_meta">Cell Metadata</a></li>');
                 //htmls.push('<li><a href="#" id="tpDownload_coords">Visible coordinates</a></li>');
               htmls.push('</ul>'); // Download sub-menu
             htmls.push('<li><a href="#" id="tpSaveImage">Download current image</a></li>');
             htmls.push('</li>');   // sub-menu container

           htmls.push('</ul>'); // File menu
         htmls.push('</li>'); // File dropdown

         htmls.push('<li class="dropdown">');
         htmls.push('<a href="#" class="dropdown-toggle" data-toggle="dropdown" data-submenu role="button" aria-haspopup="true" aria-expanded="false">Edit</a>');
         htmls.push('<ul class="dropdown-menu">');
         htmls.push('<li><a id="tpSelectAll" href="#"><span class="dropmenu-item-label">Select all</span><span class="dropmenu-item-content">A</span></a></li>');
         htmls.push('<li><a id="tpSelectNone" href="#"><span class="dropmenu-item-label">Select none</span><span class="dropmenu-item-content">N</span></a></li>');
         htmls.push('<li><a id="tpExportIds" href="#">Export '+gSampleDesc+' IDs</a></li>');
         htmls.push('</ul>'); // View dropdown
         htmls.push('</li>'); // View dropdown

         htmls.push('<li class="dropdown">');
         htmls.push('<a href="#" class="dropdown-toggle" data-toggle="dropdown" data-submenu role="button" aria-haspopup="true" aria-expanded="false">View</a>');
         htmls.push('<ul class="dropdown-menu">');

         htmls.push('<li><a href="#" id="tpZoom100Button"><span class="dropmenu-item-label">Zoom to 100%</span><span class="dropmenu-item-content">Z</span></a></li>');


         htmls.push('<li class="dropdown-submenu"><a tabindex="0" href="#">Layout Method</a>');
           htmls.push('<ul class="dropdown-menu" id="tpLayoutMenu">');
           htmls.push('</ul>'); // Layout sub-menu

         htmls.push('</li>');   // sub-menu container

         htmls.push('<li><hr class="half-rule"></li>');

         htmls.push('<li><a href="#" id="tpOnlySelectedButton">Show only selected</a></li>');
         htmls.push('<li><a href="#" id="tpFilterButton">Hide selected '+gSampleDesc+'s</a></li>');
         htmls.push('<li><a href="#" id="tpShowAllButton">Show all '+gSampleDesc+'</a></li>');
         htmls.push('<li><a href="#" id="tpHideLabels">Hide cluster labels</a></li>');
         htmls.push('<li><hr class="half-rule"></li>');

         htmls.push('<li class="dropdown-submenu"><a tabindex="0" href="#">Transparency</a>');
           htmls.push('<ul class="dropdown-menu" id="tpTransMenu">');
             htmls.push('<li id="tpTrans0"><a href="#">0%</a></li>');
             htmls.push('<li id="tpTrans40"><a href="#">40%</a></li>');
             htmls.push('<li id="tpTrans60"><a href="#">60%</a></li>');
             htmls.push('<li id="tpTrans80"><a href="#">80%</a></li>');
           htmls.push('</ul>'); // Transparency sub-menu
         htmls.push('</li>');   // sub-menu container

         htmls.push('<li class="dropdown-submenu"><a tabindex="0" href="#">Circle size</a>');
           htmls.push('<ul class="dropdown-menu" id="tpSizeMenu">');
             htmls.push('<li id="tpSize1"><a href="#">1 px</a></li>');
             htmls.push('<li id="tpSize2"><a href="#">2 px</a></li>');
             htmls.push('<li id="tpSize3"><a href="#">3 px</a></li>');
             htmls.push('<li id="tpSize4"><a href="#">4 px</a></li>');
             htmls.push('<li id="tpSize5"><a href="#">5 px</a></li>');
             htmls.push('<li id="tpSize6"><a href="#">6 px</a></li>');
             htmls.push('<li id="tpSize7"><a href="#">7 px</a></li>');
             htmls.push('<li id="tpSize8"><a href="#">8 px</a></li>');
           htmls.push('</ul>'); // Circle size sub-menu
         htmls.push('</li>');   // sub-menu container

         htmls.push('</ul>'); // View dropdown-menu
         htmls.push('</li>'); // View dropdown container

         htmls.push('<li class="dropdown">');
         htmls.push('<a href="#" class="dropdown-toggle" data-toggle="dropdown" data-submenu role="button" aria-haspopup="true" aria-expanded="false">Help</a>');
         htmls.push('<ul class="dropdown-menu">');
         htmls.push('<li><a href="#" id="tpTutorialButton">Tutorial</a></li>');
         htmls.push('</ul>'); // Help dropdown-menu
         htmls.push('</li>'); // Help dropdown container

       htmls.push('</ul>'); // navbar-nav

       //htmls.push('<ul class="nav navbar-nav navbar-right">');
           //htmls.push('<li><a href="#">Zoom out 100%</a></li>');
       //htmls.push('</ul>'); // navbar-nav

       htmls.push('</div>'); // container
       htmls.push('</nav>'); // navbar
       htmls.push('</div>'); // tpMenuBar

       $(document.body).append(htmls.join(""));
      
       $('#tpTransMenu li a').click( onTransClick );
       $('#tpSizeMenu li a').click( onSizeClick );
       $('#tpFilterButton').click( onHideSelectedClick );
       $('#tpOnlySelectedButton').click( onShowOnlySelectedClick );
       $('#tpZoom100Button').click( onZoom100Click );
       $('#tpShowAllButton').click( onShowAllClick );
       $('#tpHideLabels').click( onHideLabelsClick );
       $('#tpExportIds').click( onExportIdsClick );
       $('#tpTutorialButton').click( function()  { showIntro(false); } );
       $('#tpOpenDatasetLink').click( onOpenDatasetLink );
       $('#tpSaveImage').click( onSaveAsClick );
       $('#tpSelectAll').click( onSelectAllClick );
       $('#tpSelectNone').click( onSelectNoneClick );
       $('#tpDownloadMenu li a').click( onDownloadClick );

       updateMenu();

       // make sure that anchors under disabled li-elements are not clickable
       //$('tpMenuBar').on('click', '.disabled > a', function(e) {
           //e.preventDefault();
           //return false;
       //});

       menuBarHeight = $('#tpMenuBar').outerHeight(true);
       
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

       // This version is like the existing genome browser menu hover effect
       // activate mouse hovers for the navbar http://stackoverflow.com/a/21486327/233871
       //$(".dropdown").hover(
       //    function(){ $(this).addClass('open') },
       //    function(){ $(this).removeClass('open') }
       //);
       // do the same thing for submenus
       //$(".dropdown-submenu").hover(
       //    function(){ $(this).addClass('open') },
       //    function(){ $(this).removeClass('open') }
       //);
       // activate the submenus
       $('[data-submenu]').submenupicker();
    }

    function drawMenu() {
       var htmls = [];
       htmls.push("<div id='tpMenuBar'>");

       htmls.push("<div id='tpZoom100Button' class='tpButton'>Zoom 100%</div>");

       htmls.push("<div id='tpDatasetButton' class='tpButton'>Dataset <span style='font-size:0.9em'>&#9662;</span>");
       htmls.push("<ul id='tpDatasetMenu' class='tpMenu'>");
       for (var i = 0; i < gOptions.datasets.length; i++) {
           var dataset = gOptions.datasets[i];
           var line = "<li><div id='tpDataset_"+i+"'>"+dataset.label+"</div></li>";
           htmls.push(line);
       }
       htmls.push("</ul>");
       htmls.push("</div>");

       htmls.push("<div id='tpMethodButton' class='tpButton'>Method <span style='font-size:0.9em'>&#9662;</span>");
       htmls.push("<ul id='tpMethodMenu' class='tpMenu'>");
       htmls.push("<li><div id='seurat'>Seurat T-SNE</div></li>");
       htmls.push("<li><div id='pca12'>PCA: PC1 / PC2</div></li>");
       htmls.push("<li><div id='pca23'>PCA: PC2 / PC3</div></li>");
       htmls.push("<li><div id='pca34'>PCA: PC3 / PC4</div></li>");
       htmls.push("</ul>");
       htmls.push("</div>");

       htmls.push("<div id='tpSizeButton' class='tpButton'>Size <span style='font-size:0.9em'>&#9662;</span>");
       htmls.push("<ul id='tpSizeMenu' class='tpMenu'>");
       htmls.push("<li><div>1 px</div></li>");
       htmls.push("<li><div>2 px</div></li>");
       htmls.push("<li><div>3 px</div></li>");
       htmls.push("<li><div>4 px</div></li>");
       htmls.push("<li><div>5 px</div></li>");
       htmls.push("<li><div>6 px</div></li>");
       htmls.push("</ul>");
       htmls.push("</div>");

       htmls.push("<div id='tpTransButton' class='tpButton'>Transparency <span style='font-size:0.9em'>&#9662;</span>");
       htmls.push("<ul id='tpTransMenu' class='tpMenu'>");
       htmls.push("<li><div>0%</div></li>");
       htmls.push("<li><div>10%</div></li>");
       htmls.push("<li><div>20%</div></li>");
       htmls.push("<li><div>30%</div></li>");
       htmls.push("<li><div>40%</div></li>");
       htmls.push("<li><div>50%</div></li>");
       htmls.push("<li><div>60%</div></li>");
       htmls.push("<li><div>70%</div></li>");
       htmls.push("<li><div>80%</div></li>");
       htmls.push("<li><div>90%</div></li>");
       htmls.push("</ul>");
       htmls.push("</div>");

       htmls.push("<div id='tpFilterButton' class='tpButton tpButtonInactive'>Hide selected</div>");
       htmls.push("<div id='tpOnlySelectedButton' class='tpButton tpButtonInactive'>Show only selected</div>");
       htmls.push("<div id='tpShowAllButton' class='tpButton tpButtonInactive'>Show all</div>");
       htmls.push("<div id='tpHelpButton' class='tpButton'>Tutorial</div>");

       htmls.push("</div>"); // tpMenuBar

       $(document.body).append(htmls.join(""));
       $( ".tpMenu" ).menu();
       //jQuery( function() { $('#tpSizeMenu').width($('#tpSizeButton').width()); } );
       jQuery( function() { $('#tpTransMenu').width($('#tpTransButton').width()); } );

       $('#tpZoom100Button').click( onZoom100Click );

       $('#tpSizeButton').click( function()   { $("#tpSizeMenu").toggle(); } );
       $('#tpMethodButton').click( function() { $("#tpMethodMenu").toggle(); } );
       $('#tpTransButton').click( function()  { $("#tpTransMenu").toggle(); } );
       $('#tpDatasetButton').click( function()  { $("#tpDatasetMenu").toggle(); } );
       $('#tpHelpButton').click( function()  { showIntro(false); } );

       $('#tpMethodMenu li div').click( onMethodClick );
       $('#tpTransMenu li div').click( onTransClick );
       $('#tpSizeMenu li div').click( onSizeClick );
       $('#tpDatasetMenu li div').click( onDatasetClick );

       $('#tpFilterButton').click( onHideSelectedClick );
       $('#tpOnlySelectedButton').click( onShowOnlySelectedClick );
       $('#tpShowAllButton').click( onShowAllClick );
    }

    function resizeRenderer() {
       /* resize the canvas based on current window size */
       var fromLeft = legendBarWidth+legendBarMargin;
       var width = window.innerWidth - metaBarWidth - fromLeft;
       var height  = window.innerHeight - geneBarHeight - menuBarHeight;
       renderer.resize(width, height);
       return {"renderer":renderer, "width" : width, "height" : height};
    }

    function onBackgroundMouseMove(ev) {
        /* called when the mouse is moved over the Canvas */
           if (mouseDownX === null)
               return;
           //$("#tpSelectBox").css({left:mouseDownX+legendBarWidth+legendMargin, top:mouseDownY+menuBarHeight}).show();
           var x = ev.data.global.x;
           var y = ev.data.global.y;
           var selectWidth = Math.abs(x - mouseDownX);
           var selectHeight = Math.abs(y - mouseDownY);
           var minX = Math.min(ev.data.global.x, mouseDownX);
           var minY = Math.min(ev.data.global.y, mouseDownY);
           var posCss = {
                "left":minX + metaBarWidth+metaBarMargin, 
                "top": minY + menuBarHeight,
                "width":selectWidth,
                "height":selectHeight
           };
           $("#tpSelectBox").css(posCss).show();
    }

    function onBackgroundMouseClick(ev) {
       //if (mouseDownX !== null) // currently zooming? do not select
        //   return;
       var x = ev.data.global.x;
       var y = ev.data.global.y;
       // 
       if ((x-mouseDownX)!==0 || (y-mouseDownY)!==0) {
            //mouseDownX = null;
            //mouseDownY = null;
           return;
       }

       console.log("click onto free space, clearing selection");
       mouseDownX = null;
       mouseDownY = null;
       if (!jQuery.isEmptyObject(gSelCellIds)) {
           gSelCellIds = {};
           updateSelection();
       }
    }

    function onBackgroundMouseDown (ev) {
    /* click on background */
        // just keep the coordinates, do nothing else
        mouseDownX = ev.data.global.x;
        mouseDownY = ev.data.global.y;
    }

    function onBackgroundMouseUp(ev) {
       if (mouseDownX === null) // user did not start the click on the canvas
           return;
       $("#tpSelectBox").hide();
       $("#tpSelectBox").css({"width":0, "height":0});

       var x = ev.data.global.x;
       var y = ev.data.global.y;

       // do nothing if it was just a single click
       if ((x-mouseDownX)===0 && (y-mouseDownY)===0) {
            mouseDownX = null;
            mouseDownY = null;
           return;
       }

       // address case when the mouse movement was upwards
       var pxMinX = Math.min(x, mouseDownX);
       var pxMaxX = Math.max(x, mouseDownX);
       var pxMinY = Math.min(y, mouseDownY);
       var pxMaxY = Math.max(y, mouseDownY);

       // convert screen pixel click coords to data float values
       var spanX = zoomRange.maxX - zoomRange.minX;
       var spanY = zoomRange.maxY - zoomRange.minY;
       var xMult = spanX / winInfo.width;
       var yMult = spanY / winInfo.height;

       var oldMinX = zoomRange.minX;
       var oldMinY = zoomRange.minY;
       zoomRange.maxX = oldMinX + (pxMaxX * xMult);
       zoomRange.minX = oldMinX + (pxMinX * xMult);
       zoomRange.maxY = oldMinY + (pxMaxY * yMult);
       zoomRange.minY = oldMinY + (pxMinY * yMult);
       console.log(zoomRange);

       mouseDownX = null;
       mouseDownY = null;

       pixelCoords = scaleData(shownCoords);
       plotDots();
       renderer.render(stage);
    }

    function createRenderer() {
    /* Create the renderer, sideeffect: sets the "stage" global. 
     * menuBarHeight needs to be set before calling this.
     */
       var opt = { antialias: true, transparent: true };

       //renderer = new PIXI.WebGLRenderer(0, 0, opt);
       //renderer = PIXI.autoDetectRenderer(0, 0, opt);
       renderer = new PIXI.CanvasRenderer(0, 0, opt);

       // setup the mouse zooming callbacks
       renderer.plugins.interaction.on('mousedown', onBackgroundMouseDown);
       renderer.plugins.interaction.on('mousemove', onBackgroundMouseMove );
       renderer.plugins.interaction.on('mouseup', onBackgroundMouseUp);

       // add the divs for the selection rectangle
       var htmls = [];
       htmls.push("<div style='position:absolute; display:none' id='tpSelectBox'></div>");
       //htmls.push("<div style='position:absolute; left:"+fromLeft+"px;' id='tpMiddle'>");
       $(document.body).append(htmls.join(""));

       $(window).resize(onResize);
       
       renderer.view.style.border = "1px solid #AAAAAA";
       renderer.backgroundColor = 0xFFFFFF;

       // fill the whole window
       renderer.view.style.position = "absolute";
       renderer.view.style.display = "block";
       renderer.view.style.left = metaBarWidth+metaBarMargin+"px";
       renderer.view.style.top = menuBarHeight+"px";
       renderer.autoResize = true;

       if (renderer.type instanceof PIXI.CanvasRenderer) {
           console.log('Using Canvas');
       } else {
           console.log('Using WebGL');
       } 

       $(document.body).append(renderer.view);
       return resizeRenderer();

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

    function buildLegend(sortBy) {
    /* build the gLegend and gClasses globals */
        if (gLegend.type=="meta")
            {
            gLegend = buildLegendForMeta(gLegend.metaFieldIdx, sortBy);
            gClasses = assignCellClasses();
            }
        else
            buildLegendForExprAndClasses(sortBy);
    }

    function filterCoordsAndUpdate(cellIds, mode) {
    /* hide/show currently selected cell IDs or "show all". Rebuild the legend and the coloring. */
        if (mode=="hide")
            shownCoords = removeCellIds(shownCoords, cellIds);
        else if (mode=="showOnly")
            shownCoords = showOnlyCellIds(shownCoords, cellIds);
        else
            shownCoords = allCoords.slice();

        pixelCoords = scaleData(shownCoords);

        buildLegend();
        drawLegend();
        gSelCellIds = {};
        plotDots();
        renderer.render(stage);
        menuBarShow("#tpShowAllButton");
    }

    function onHideSelectedClick(ev) {
    /* user clicked the hide selected button */
        filterCoordsAndUpdate(gSelCellIds, "hide");
        menuBarHide("#tpFilterButton");
        menuBarHide("#tpOnlySelectedButton");
        menuBarShow("#tpShowAllButton");
        ev.preventDefault();
    }

    function onShowOnlySelectedClick(ev) {
    /* user clicked the only selected button */
        filterCoordsAndUpdate(gSelCellIds, "showOnly");
        menuBarHide("#tpFilterButton");
        menuBarHide("#tpOnlySelectedButton");
        menuBarShow("#tpShowAllButton");
        ev.preventDefault();
    }

    function onShowAllClick(ev) {
    /* user clicked the show all menu entry */
        //gSelCellIds = {};
        filterCoordsAndUpdate(gSelCellIds, "showAll");
        shownCoords = allCoords.slice(); // complete copy of list, fastest in Blink
        pixelCoords = scaleData(shownCoords);
        buildLegend();
        drawLegend();
        gClasses = assignCellClasses();
        plotDots();
        menuBarHide("#tpFilterButton");
        menuBarHide("#tpOnlySelectedButton");
        menuBarHide("#tpShowAllButton");
        gLegend.lastClicked = null;
        renderer.render(stage);
    }

    function onHideLabelsClick(ev) {
    /* user clicked the hide labels button */
        //gSelCellIds = {};
        gLabelMetaIdx = null;
        plotDots();
        menuBarHide("#tpHideLabels");
        renderer.render(stage);
    }

    function onMethodClick(ev) {
    /* user clicked a Method menu entry */
        var methodType = ev.target.id;
        switch (methodType) {
            case "pca12" :
                coord1Idx = 1;
                coord2Idx = 2;
                allCoords = pcaCoords;
                break;
            case "pca23" :
                coord1Idx = 2;
                coord2Idx = 3;
                allCoords = pcaCoords;
                break;
            case "pca34" :
                coord1Idx = 3;
                coord2Idx = 4;
                allCoords = pcaCoords;
                break;
            case "seurat" :
                coord1Idx = 1;
                coord2Idx = 2;
                allCoords = seuratCoords;
                break;
        }

        scaleDataAndColorByCluster();
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
        zoomRange = cloneObj(zoomFullRange);
        pixelCoords = scaleData(shownCoords);
        plotDots();
        renderer.render(stage);
        ev.preventDefault();
    }

    function onResize(ev) {
        /* called when window is resized by user */
        winInfo = resizeRenderer();
        pixelCoords = scaleData(shownCoords);
        var legendBarLeft = winInfo.width+metaBarMargin+metaBarWidth;
        $('#tpLegendBar').css('left', legendBarLeft+"px");
        var geneBarTop = winInfo.height+menuBarHeight+geneBarMargin;
        drawGeneBar();
        plotDots();
        renderer.render(stage);
    }

    function drawLegendBar() {
        // create an empty right side legend bar
        var legendLeft = metaBarWidth+metaBarMargin+winInfo.width;
        var htmls = [];
        htmls.push("<div id='tpLegendBar' style='position:absolute;top:"+menuBarHeight+"px;left:"+legendLeft+"px; width:"+legendBarWidth+"px'>");
        htmls.push("<div id='tpLegendTitle' style='position:relative; width:100%; height:1.5em; font-weight: bold'>");
        htmls.push("Legend:");
        htmls.push("</div>"); // title
            //htmls.push("<div id='tpLegendTabs'>");
                //htmls.push("<ul>");
                //htmls.push("<li><a href='#tpLegendTab1'>Fill</a></li>");
                //htmls.push("<li><a href='#tpLegendTab2'>Outline</a></li>");
                //htmls.push("</ul>");
                //htmls.push("<div id='tpLegendTab1'>");
                    htmls.push("<div id='tpLegendContent'>");
                    htmls.push("</div>"); // content 
                //htmls.push("</div>"); // tab1 
                //htmls.push("<div id='tpLegendTab2'>");
                //htmls.push("</div>"); // tab2 
            //htmls.push("</div>"); // tabs 
        htmls.push("</div>"); // bar 
        $(document.body).append(htmls.join(""));
        //$('#tpLegendTabs').tabs();
        //$('#tpLegendTabs > ul').hide();
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

    function makeColorMap(legendRows) {
    /* go over the legend lines, create a map colorMapKey : color */
        var colorMap = {};
        for (var i = 0; i < legendRows.length; i++) {
            var colorMapKey = legendRows[i][4];
            var color = legendRows[i][0];
            colorMap[colorMapKey] = color;
        }
        return colorMap;
    }

    function loadColors(legendRows) {
    /* go over the legend lines, use the unique key to find the manually
     * defined colors in localStorage or in gOptions.metaColors */
        _dump(gOptions);
        var legendColors = gCurrentDataset.metaColors;
        for (var i = 0; i < legendRows.length; i++) {
            var legKey = legendRows[i][5];
            var defColor = legendRows[i][1];
            var legLabel = legendRows[i][2];
            // first check localstorage
            var color = localStorage.getItem(legKey);
            if (legendColors !== undefined && ((color===undefined) || (color===null))) {
                    color = legendColors[legLabel];
            }

            if (color===undefined || color === null) {
                color=defColor;
            }
            legendRows[i][0] = color;
        }
        return legendRows;
    }

    function assignCellClasses() {
    /* return a map cellId:colorInteger, based on the currently active meta field */
        var metaIdx = gLegend.metaFieldIdx;

        if (gLegend.rows===undefined)
            return;

        // build map meta val -> legendId
        var metaToLegend = {};
        var rows = gLegend.rows;
        var metaVal = null;
        for (var i = 0; i < rows.length; i++) {
            metaVal = rows[i][4];
            metaToLegend[metaVal] = i;
        }

        var cellIdToLegendId = {};
        for (var j = 0; j < allCoords.length; j++) {
            var cellId = allCoords[j][0];
            metaVal = metaData[cellId][metaIdx];
            cellIdToLegendId[cellId] = metaToLegend[metaVal];
        }
        return cellIdToLegendId;
    }

    function newArray(n, val) {
        /* return an array initialized with val */
        var arr = new Array(n);
        for (var i = 0; i < n; i++)
            arr[i] = val;
        return arr;
    }

    function findBin(ranges, val) {
    /* given an array of values, find the index i where ranges[i] < val < ranges[i+1] */
    /* XX performance - use binary search */
    for (var i = 0; i < ranges.length-1; i++) {
        if ((val >= ranges[i]) && (val <= ranges[i+1]))
            return i;
        }
    }

    function buildLegendForExprAndClasses() {
        /* build legend for expression and update the gClasses global var */
        var geneIdx = gLegend.geneIdx;
        var geneId = geneFields[geneIdx][0];
        var deciles = gDeciles[geneId];

        $('#tpLegendTitle').html(geneId+" expression");
        
        // get count for each decile and the special bin for null values
        var decileCounts = newArray(10, 0);
        var nullCount = 0;
        var cellToClass = {};
        for (var i = 0; i < shownCoords.length; i++) {
            var cellId = shownCoords[i][0];
            var exprVal = geneExpr[cellId][geneIdx]; 
            if (exprVal===null)
                {
                nullCount++;
                cellToClass[cellId] = 0;
                }
            else {
                var binIdx = findBin(deciles, exprVal);
                decileCounts[binIdx]++;
                cellToClass[cellId] = binIdx+1; // binIdx 0 is reserved for the null count
            }
        }

        var legendRows = [];
        // "no value" 
        legendRows.push( [ cNullColor, cNullColor, "No Value", nullCount, 0, geneId+"|null"] );

        if ( deciles!==null ) {
            var defColors = makeColorPalette(10, true);
            for (var j = 0; j < (deciles.length-1); j++) {
                var count = decileCounts[j];
                var binMin = deciles[j];
                var binMax = deciles[j+1];
                var legendId = j+1; // bin 0 is reserved for null, so legendId = binIdx+1
                var legLabel = binMin.toFixed(2)+' to '+binMax.toFixed(2);

                var defColor = defColors[j];
                var uniqueKey = geneId+"|"+legLabel;
                legendRows.push( [ defColor, defColor, legLabel, count, legendId, uniqueKey] );
            }
        }

    gLegend.rows = loadColors(legendRows);
    gClasses = cellToClass;
    }

    function buildLegendForExpr_old() {
        // get expression values, keep null values separate
        var geneIdx = gLegend.geneIdx;
        var geneId = geneFields[geneIdx][0];
        var exprVals = [];
        var nullIds = [];
        for (var i = 0; i < shownCoords.length; i++) {
            var cellId = shownCoords[i][0];
            var exprVal = geneExpr[cellId][geneIdx]; 
            if (exprVal==null)
                nullIds.push(cellId);
            else
                exprVals.push( [exprVal, cellId] );
        }

        // sort them
        exprVals = exprVals.sort(function(a, b){ return a[0] - b[0]; }); // sort by exprVal

        // go over the values and assign cellIds to bins
        var cell2bin = {};
        var binSize = exprVals.length / 10;
        var binMaxArr = newArray(11, -Number.MAX_VALUE); // 10 bins + bin0 for null values
        var binCounts = newArray(11, 0);
        for (i = 0; i < exprVals.length; i++) {
            var exprRow = exprVals[i];
            var exprVal = exprRow[0];
            var cellId = exprRow[1];

            var binIdx = 1+Math.floor(i / binSize); // binIdx=0 is reserved
            binMaxArr[binIdx] = Math.max(binMaxArr[binIdx], exprVal);
            binCounts[binIdx]++;
            cell2bin[cellId] = binIdx;
        }

        // add the cellIds with null values
        for (i = 0; i < nullIds.length; i++) {
            cell2bin[nullIds[i]] = 0;
        }

        gClasses = cell2bin;

        // contruct the legend
        var geneSym = geneFields[geneIdx][0];
        var legendRows = [];

        // bin0 is "no value" 
        legendRows.push( [ cNullColor, cNullColor, "No Value", nullIds.length, 0, geneSym+"|null"] );

        if (exprVals.length>0) {
            var defColors = makeColorPalette(binMaxArr.length, true);
            var lastMax = exprVals[0][0];
            for (i = 1; i < binCounts.length; i++) {
                var count = binCounts[i];
                var binMin = lastMax;
                var binMax = binMaxArr[i];
                var binId = i;
                var legLabel = binMin.toFixed(2)+' to '+binMax.toFixed(2);
                var defColor = defColors[i];
                var uniqueKey = geneSym+"|"+legLabel;
                legendRows.push( [ defColor, defColor, legLabel, count, binId, uniqueKey] );
                lastMax = binMax;
            }
        }

        gLegend.rows = loadColors(legendRows);

    }

    function onGeneClick (event) {
    /* user clicked on a gene in the gene bar: sort expression vals into deciles and color by them  */
        var geneIdx = parseInt(event.target.id.split("_")[1]); // the index of the gene
        $('.tpMetaBox').removeClass('tpMetaSelect');
        $('.tpGeneBarCell').removeClass("tpGeneBarCellSelected");
        $('#tpGeneBarCell_'+geneIdx).addClass("tpGeneBarCellSelected");
        gSelCellIds = {};
        gLegend = {};
        gLegend.type = "gene";
        gLegend.geneIdx = geneIdx;
        buildLegend();
        drawLegend();
        plotDots();
        renderer.render(stage);
    }

    function showDialogBox(htmlLines, title, options) {
        /* show a dialog box with html in it */
        $('#tpDialog').remove();
        var addStr = "";
        if (options.width!=undefined)
            addStr = "max-width:"+options.width+"px;";
        var maxHeight = $(window).height()-200;
        // unshift = insert at pos 0
        htmlLines.unshift("<div style='display:none;"+addStr+"max-height:"+maxHeight+"px' id='tpDialog' title='"+title+"'>");
        htmlLines.push("</div>");
        $(document.body).append(htmlLines.join(""));
        var dialogOpts = {modal:true};
        if (options.width!=undefined)
            dialogOpts["width"] = options["width"];
        if (options.height!=undefined)
            dialogOpts["height"] = options["height"];
        dialogOpts["maxHeight"] = maxHeight;
        if (options.buttons!=undefined)
            dialogOpts["buttons"] =  options.buttons;
        else
            dialogOpts["buttons"] =  {};

        if (options.showOk!=undefined)
            dialogOpts["buttons"]["OK"] = function() { $( this ).dialog( "close" ); };
        if (options.showClose!=undefined)
            dialogOpts["buttons"]["Close"] = function() { $( this ).dialog( "close" ); };
        //dialogOpts["position"] = "center";
        //dialogOpts["height"] = "auto";
        //dialogOpts["width"] = "auto";

        $( "#tpDialog" ).dialog(dialogOpts);
    }

    function onChangeGenesClick(ev) {
    /* called when user clicks the "change" button in the gene list */
        var htmls = [];
        //htmls.push("<div id='tpGeneDialog'>");
        htmls.push("<p style='padding-bottom: 5px'>Enter a list of gene symbols, one per line:</p>");
        //htmls.push("<div>");
        htmls.push("<textarea id='tpGeneListBox' class='form-control' style='height:320px'>");


        for (var i = 0; i < geneFields.length; i++) {
            htmls.push(geneFields[i][0]+"\r\n");
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

        // get the sorted list of cell IDs - XX is there a risk that the matrix does not have a sorted order? 
        // XX Risk of python sorting a different way than javascript?
        var cellIds = [];
        for (var cellId in metaData) {
            cellIds.push(cellId);
        }
        cellIds.sort();

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
        for (var geneI = 0; geneI < gLoad_geneList.length; geneI++) {
            var geneSym = gLoad_geneList[geneI];
            var geneInfo = gLoad_geneExpr[geneSym];

            var geneId = geneInfo[0];
            var geneVec = geneInfo[1];
            newGeneFields.push( [geneSym, geneId] );

            for (var cellI = 0; cellI < cellIdCount; cellI++) {
                cellId = cellIds[cellI];
                newExprData[cellId][geneI] = geneVec[cellI];
            }
        }

        gLoad_geneList = null;
        gLoad_geneExpr = null;
        geneExpr = newExprData;
        geneFields = newGeneFields;
        drawGeneBar(); 

    }

    function onReceiveExprLine(line) {
        /* called when a line of the expression matrix has been loaded */
        var symbol = this.geneSymbol;
        console.log("Got gene "+symbol);
        var fields = line.split("\t");
        var geneId = fields[0];
        // convert all but the first field to a new vector of floats
        var exprVec = [];
        for (var i = 1; i < fields.length; i++) {
            exprVec.push( parseFloat(fields[i]));
        }
        var progressbar = $( "#tpGeneProgress" );
        var val = progressbar.progressbar( "value" ) || 0 ;
        val++;
        progressbar.progressbar( "value", val );
        $( "#tpProgressLabel" ).text(symbol);
        gLoad_geneExpr[symbol] = [geneId, exprVec];

        var progrMax = progressbar.progressbar("option", "max");
        if (val >= progrMax)
            onGeneLoadComplete();
    }

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
                //progressInc : progressPerGene,
                success: onReceiveExprLine
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

    function drawGeneBar() {
    /* create bottom gene expression info bar */
        var canvasWidth = winInfo.width;
        var canvasHeight = winInfo.height;
        $('#tpGeneBar').remove();

        $(document.body).append("<div id='tpGeneBar' style='position:absolute;top:"+(menuBarHeight+canvasHeight+geneBarMargin)+"px;left:"+(metaBarWidth+metaBarMargin)+"px;width:"+canvasWidth+"px; height:"+geneBarHeight+"px; background-color:light-grey'></div>");

        var html = [];
        html.push('<table id="tpGeneTable"><tr>');
        html.push('<td><button id="tpChangeGenes" title="Change the list of genes that are displayed in this table" class = "ui-button ui-widget ui-corner-all" style="width:95%">Change</button></td>');

        var cellWidth = 83; // -> must correspond to CSS class tpGeneBarCell + padding + border
        var colsPerRow = Math.floor(canvasWidth / cellWidth);

        var currWidth = 1;
        for (var i = 0; i < geneFields.length; i++) {
            var geneName = geneFields[i][0];
            var geneDesc = geneFields[i][1];
            if (((i % colsPerRow) == 0) && (i!=0)) {
                html.push("</tr><tr>");
                currWidth=0;
            }
            html.push('<td title="'+geneDesc+'" id="tpGeneBarCell_'+i+'" class="tpGeneBarCell">'+geneName+'</td>');
            currWidth += cellWidth;
        }
        html.push("</tr></table>");

        var geneBarEl = $('#tpGeneBar');
        geneBarEl.append(html.join("")); // much faster than creating indiv. DOM elements.
        $('.tpGeneBarCell').click( onGeneClick );
        $('#tpChangeGenes').click( onChangeGenesClick );
    }

    function addMenus() {
        // construct the link for the menu
        //$('#tpMetaBar').append('<a href="#" id="tpColorByLink">Color By Meta &#9662;</a>&nbsp;');

        // construct the link for the menu
        //$('#tpMetaBar').append('<a href="#" id="tpColorByGeneLink">Color By Gene &#9662;</a><br>');


        // build the color-by menu
        $('#tpMetaBar').append('<ul class="tpMenu" id="tpColorByMenu" style="display:none">');
        for (var i = 0; i < fields.length; i++) {
            var fieldName = fields[i];
            $('#tpColorByMenu').append('<li><div>'+fieldName+'</div></li>');
        }
        $('#tpColorByMenu').append('</ul>');
        $( "#tpColorByMenu" ).menu();

        var tpColorByLink = $("#tpColorByLink");
        tpColorByLink.click(function() {
          var menu = $('#tpColorByMenu');
          if (menu.is(":visible"))
              {
              console.log("hide it");
              menu.hide();
              }
          else
              {
              console.log("show it");
              menu.show();
              }
        });

        // close menu when document is clicked anywhere else -- ? performance ?
        var tpMenuEl = $('.tpMenu');
        $(document).click(function(e) {
            clicked = $(e.target);
            if (clicked[0]===tpColorByLink[0])
                return;
            if(!$.contains(tpMenuEl, clicked)) {
                console.log("clicked outside, now hiding");
                $('.tpMenu').hide();
            }
        });
    }

    function likeEmptyString(label) {
    /* some special values like "undefined" and empty string get colored in grey  */
        return (label===null || label.trim()==="" || label==="none" || label==="None" 
                || label==="Unknown" || label==="NaN" || label==="NA" || label==="undefined" || label=="Na")
    }

    function buildLegendForMeta(metaIndex, sortBy) {
    /* Build a new gLegend object and return it */
        var legend = {};
        legend.type = "meta";
        legend.metaFieldIdx = metaIndex;

        // first make a list of all possible meta values and their counts
        var metaCounts = {};
        for (var i = 0; i < pixelCoords.length; i++) {
            var cellId = pixelCoords[i][0];
            var metaRow = metaData[cellId];
            var metaVal = metaRow[metaIndex];
            metaCounts[metaVal] = 1 + (metaCounts[metaVal] || 0);
        }

        if (keys(metaCounts).length > 100) {
            alert("Cannot color on a field that has more than 100 different values");
            return null;
        }

        // sort like numbers if the strings are mostly numbers, otherwise sort like strings
        var sortResult = legendSort(metaCounts, sortBy);
        var countListSorted = sortResult.list;
        
        if (fieldName!==undefined)
            $("#tpLegendTitle").html(fieldName.replace("_", " "));

        var fieldName = metaFields[metaIndex];
        // force field names that look like "cluster" to a rainbow palette
        // even if they look like numbers
        if (fieldName.indexOf("luster") != -1 && sortBy===undefined)
            sortBy = "count";

        var defaultColors = makeColorPalette(countListSorted.length, sortResult.useGradient);

        var rows = [];
        for (i = 0; i < countListSorted.length; i++) {
            var count = countListSorted[i][0];
            var label = countListSorted[i][1];
            var color = defaultColors[i];
            var uniqueKey = fieldName+"|"+label;
            var colorMapKey = label;
            if (likeEmptyString(label))
                color = cNullColor;

            rows.push( [ color, color, label, count, colorMapKey, uniqueKey] );
        }

        legend.rows = loadColors(rows);
        legend.isSortedByName = sortResult.isSortedByName;
        return legend;
    }

    function findCellIdsForLegendId (cellIdToLegendId, legendId, cellIds) {
    /* given a legendId, return an object with cellIds that are associated to a given class */
        if (cellIds==undefined)
            cellIds = {};
        for (var cellId in cellIdToLegendId) {
            if (gClasses[cellId] == legendId)
                cellIds[cellId] = true;
            }
        return cellIds;
    }

    function colorByMetaField(fieldId) {
    /* color by a meta field and rebuild/redraw everything */
        var legend = buildLegendForMeta(fieldId);
        if (legend==null)
            return;

        $('.tpMetaBox').removeClass('tpMetaSelect');
        $('#tpMetaBox_'+fieldId).addClass('tpMetaSelect');
        $('.tpGeneBarCell').removeClass('tpGeneBarCellSelected');

        gLegend = legend;
        drawLegend();
        gClasses = assignCellClasses();
        plotDots();
        renderer.render(stage);
    }

    function onMetaClick (event) {
    /* called when user clicks a meta data field or label */
        var fieldId = parseInt(event.target.id.split("_")[1]);
        gSelCellIds = {};
        colorByMetaField(fieldId);
    }

    function drawToolBar () {
    /* add the tool bar with icons of tools */
        $("#tpToolBar").remove();
        $(document.body).append("<div id='tpToolBar' style='position:absolute;left:0px;top:"+menuBarHeight+"px'></div>");

        var htmls = [];
        htmls.push('<button title="Activate zoom" id="tpIconZoom" type="button" class="btn-small btn-outline-primary"><i class="material-icons">zoom_in</i></button>');
        htmls.push('<button title="Zoom to 100%" id="tpIconZoom100" type="button" class="btn-small btn-outline-primary"><i style="height:20px" class="material-icons">crop_free</i></button>');
        $("#tpToolBar").append(htmls.join(""));

       $('#tpIconZoom100').click( onZoom100Click );
    }

    function drawMetaBar () {
    /* add the left sidebar with the meta data fields */
        $("#tpMetaBar").remove();
        $(document.body).append("<div id='tpMetaBar' style='position:absolute;left:0px;top:"+(menuBarHeight+toolBarHeight)+"px'></div>");

        var htmls = [];
        for (var i = 0; i < metaFields.length; i++) {
            var fieldName = metaFields[i].replace(/_/g, " ");
            htmls.push("<div class='tpMetaBox' id='tpMetaBox_"+i+"'><div id='tpMetaLabel_"+i+"' class='tpMetaLabel'>"+fieldName+":</div><div class='tpMetaValue' style='width:"+(metaBarWidth-2*metaBarMargin)+"px' id='tpMeta_"+i+"'>&nbsp;</div></div>");
        }
        $("#tpMetaBar").append(htmls.join(""));
        $(".tpMetaLabel").click( onMetaClick );
        $(".tpMetaValue").click( onMetaClick );
        
        // setup the right-click menu
        var menuItems = [{name: "Use as cluster label"}];
        var menuOpt = {
            selector: ".tpMetaBox",
            items: menuItems, 
            className: 'contextmenu-customwidth',
            callback: onMetaRightClick
        };
        $.contextMenu( menuOpt );
    }


    function addInfoBars() {
    /* add the info-panels with the meta fields and legend on the sides of the screen */
        drawGeneBar();
        drawMetaBar();
    }

    function findMinMax(coords) {
    /* given list of (x,y) tuples, find min/max range of values */
        var minX = 9999999;
        var maxX = -9999999;
        var minY = 9999999;
        var maxY = -9999999;
 
        for (var i = 0; i < coords.length; i++) {
            var x = coords[i][1];
            var y = coords[i][2];
                
            minX = Math.min(minX, x);
            maxX = Math.max(maxX, x);
            minY = Math.min(minY, y);
            maxY = Math.max(maxY, y);
        }
        
        var ret = {};
        ret.minX = minX;
        ret.maxX = maxX;
        ret.minY = minY;
        ret.maxY = maxY;
        return ret;
    }
        
    function scaleData(coords) {
    /* scale coords (x,y float coordinates) to integer pixels on screen and return new list of (cellId, x, y) */
        var winWidth = winInfo.width;
        var winHeight = winInfo.height;
        console.log("Starting coordinate scaling");

        // var access is faster than object access
        var minX = zoomRange.minX;
        var maxX = zoomRange.maxX;
        var minY = zoomRange.minY;
        var maxY = zoomRange.maxY;
        
        var spanX = maxX - minX;
        var spanY = maxY - minY;

        // transform from data floats to screen pixel coordinates
        var newCoords = [];
        for (var i = 0; i < shownCoords.length; i++) {
            var cellId = shownCoords[i][0];
            var x = shownCoords[i][coord1Idx];
            var y = shownCoords[i][coord2Idx];
            // ignore anything outside of current zoom range. Performance: Quad-tree?
            if ((x < minX) || (x > maxX) || (y < minY) || (y > maxY))
                continue;
            var newX = Math.round(((x-minX)/spanX)*winWidth);
            var newY = Math.round(((y-minY)/spanY)*winHeight);
            var newRow = [cellId, newX, newY];
            newCoords.push(newRow);
        }
        console.log("Coordinate scaling done");
        gClusterMids = null; // make sure the cluster midpoints get re-calculated
        return newCoords;
    }

    function drawCircle(x, y, fill, line) {
    /* create a circle object */
       var c = new PIXI.Graphics();
       c.beginFill(fill);
       if (line != undefined)
           c.lineStyle(1, line, 1);
       c.drawCircle(x, y, circleSize);
       c.endFill();
       return c;
    }

    function drawRect(x, y, fill) {
    /* create a rectangle object */
       var rect = new PIXI.Graphics();
       rect.beginFill(fill);
       //rect.lineStyle(4, 0xFF3300, 1);
       rect.drawRect(x, y, dotWidth, dotHeight);
       rect.endFill();
       return rect;
    }

    function showIntro(addFirst) {
        /* add the intro.js data */
        //var htmls = [];
        //htmls.push("<a href='http://example.com/' data-intro='Hello step one!'></a>");
        //$(document.body).append(htmls.join(""));

        localStorage.setItem("introShown", "true");
        var intro = introJs();
        intro.setOption("hintAnimation", false);
        if (addFirst)
            intro.addStep({
                element: document.querySelector('#tpHelpButton'),
                intro: "Are you here for the first time and wondering what this is?<br>The tutorial takes only 1 minute. To skip the tutorial now, click 'Skip'.<br>You can always show it again by clicking 'Tutorial'.",
              });

        intro.addSteps(
            [
              {
                intro: "You can see a scatter plot of circles in the middle of the screen. Each circle is a "+gSampleDesc+".",
              },
              {
                element: document.querySelector('#tpMetaBar'),
                intro: "Meta data: when you move the mouse over a circle, its meta data fields will be shown here.<br>Click on a field to color the circles by the values of a field.<br>Right-click to select a field for the text labels.",
                position: 'right'
              },
              {
                element: document.querySelector('#tpGeneBar'),
                intro: "Expression data: when you move the mouse, expression values will be shown here.<br>Click on a gene to color the circles by gene expression level (log'ed).",
                position: 'top'
              },
              {
                element: document.querySelector('#tpLegendBar'),
                intro: "The legend shows the mapping from colors to meta or gene expression values.<br>Click on a color to change it.<br>Right-click to hide or show-only "+gSampleDesc+"s with certain values.<br>Click on the legend text to highlight "+gSampleDesc+"s with this value.",
                position: 'left'
              },
            ]);
        intro.start();
    }

    function loadConfig(jsonDict) {
    /* accepts dict with 'metaFields' and 'geneFields', loads into global vars and sets up window */
        metaFields = jsonDict.metaFields;
        geneFields = jsonDict.geneFields;
        gLabelMetaIdx = metaFields.indexOf(gOptions.labelField);
        if (gLabelMetaIdx==undefined) {
            gLabelMetaIdx = null;
            menuBarHide("#tpHideLabels");
        } else
            menuBarShow("#tpHideLabels");
        addInfoBars(metaFields);
        resizeRenderer();
    }

    function loadMeta(jsonDict) {
    /* accepts dict with 'meta', loads into global var 'meta' */
        metaData = jsonDict.meta;
        //if (currDataset.clusterLabelField!==undefined)
            //gLabelMetaIdx = metaFields.indexOf(gOptions.labelField);
    }

    function makeColorPalette(n, isGradient) {
    /* return an array with n color hex strings */
    // Use Google's palette functions for now, first Paul Tol's colors, if that fails, use the usual HSV rainbow
    if (!isGradient)
        return palette(["tol", "rainbow"], n);
    else
        return palette(["tol-sq"], n);
    }

    function cloneObj(d) {
    /* returns a copy of an object */
        // see http://stackoverflow.com/questions/122102/what-is-the-most-efficient-way-to-deep-clone-an-object-in-javascript
        return JSON.parse(JSON.stringify(d));
    }

    function keys(o) {
    /* return all keys of object */
        var allKeys = [];
        for(var k in o) allKeys.push(k);
        return allKeys;
    }

    function scaleDataAndColorByCluster() {
    /* called when meta and seurat coords have been loaded: scale data and color by meta field  */
        // find range of data
        zoomFullRange = findMinMax(allCoords);
        zoomRange = cloneObj(zoomFullRange);
        shownCoords = allCoords.slice();
        pixelCoords = scaleData(shownCoords);

        var clusterFieldName = gCurrentDataset.clusterField;
        if (clusterFieldName!=undefined) {
            var clusterFieldIdx = metaFields.indexOf(clusterFieldName);
            //if (clusterFieldIdx===undefined)
                //console.log("Option clusterField is not defined for dataset "+gCurrentDatasetIdx);
            if (clusterFieldIdx===undefined)
                alert("Internal error: labelling field (option 'clusterField') "+clusterFieldName+" is not a valid field");

            gLabelMetaIdx = clusterFieldIdx;
            gClusterMids = null; // force recalc
            colorByMetaField(clusterFieldIdx);
        }
        else {
            gLabelMetaIdx = null;
            gClusterMids = null; // force recalc
            colorByMetaField(0);
        }

        var introShownBefore = localStorage.getItem("introShown");
        if (introShownBefore==undefined)
            setTimeout(function(){ showIntro(true); }, 3000); // first show user the UI and let them think 3 secs
    }

    function loadTsneCoords(jsonDict) {
    /* accepts dict with 'shownCoords' loads into global var 'coords' */
        seuratCoords = jsonDict["seuratCoords"];
        allCoords = seuratCoords;
    }

    function loadPcaCoords(jsonDict) {
    /* accepts dict with 'shownCoords' loads into global var 'pcaCoords' */
        pcaCoords = jsonDict["pcaCoords"];
    }

    function loadGeneExpr(jsonDict) {
    /* accepts dict with 'geneExpr' : cellId : array of floats */
        geneExpr = jsonDict["geneExpr"];
        gDeciles = jsonDict["deciles"];
    }

    function indexLoadDone(jsonDict) {
        /* called when the matrix index is loaded */
        matrixOffsets = jsonDict;
    }

    function startLoadJson(path, func, showError) {
    /* load a json file relative to baseUrl and call a function when done */
    var fullUrl = jsonBaseUrl + path;
    return jQuery.getJSON(fullUrl, function(data){ func(data); } )
        .fail(function() { 
            if (showError) 
                alert("Could not load "+fullUrl);
        });
    }

    function loadDataset(datasetIdx) {
        /* reset the view and load a whole new dataset */
        gClusterMids = null;
        gSelCellIds = {};
        pcaCoords = null;
        seuratCoords = null;
        gLegend = null;

        gCurrentDataset = gOptions.datasets[datasetIdx];
        var coordFiles = gCurrentDataset.coordFiles;
        if (coordFiles==undefined)
            gCurrentDataset.coordFiles = [
                ["Seurat T-SNE", "seuratCoords.json"],
                ["PCA: PC1 / PC2", "pc12.json"],
                ["PCA: PC2 / PC3", "pc23.json"],
                ["PCA: PC3 / PC4", "pc34.json"],
                ["PCA: PC4 / PC5", "pc45.json"],
            ];
       drawLayoutMenu();

        jsonBaseUrl = gCurrentDataset.baseUrl;
        startLoadJson("geneExpr.json", loadGeneExpr, false );
        // make sure that we plot only when both the config, the meta data and the tsne coords
        // have been loaded
        jQuery.when(
            startLoadJson("config.json", loadConfig, true),
            startLoadJson("sampleMeta.json", loadMeta, true ),
            startLoadJson("seuratCoords.json", loadTsneCoords, true ) )
        .then( scaleDataAndColorByCluster );
        startLoadJson("geneMatrixOffsets.json", indexLoadDone, false );
        //startLoadJson("seuratCoords.json", loadTsneCoords);
        //startLoadJson("sampleMeta.json", loadMeta );
            //startLoadJson("seuratCoords.json", loadTsneCoords ) )

        startLoadJson("pcaCoords.json", loadPcaCoords, true );
    }

    function setupKeyboard() {
    /* bind the keyboard shortcut keys */
        Mousetrap.bind('o', onOpenDatasetLink);
        Mousetrap.bind('z', onZoom100Click);
        Mousetrap.bind('n', onSelectNoneClick);
        Mousetrap.bind('a', onSelectAllClick);
    }

    function loadData(opt, globalOpts) {
        /* start the JSON data loading, uses first dataset */
        if (globalOpts!=undefined) {
            if ("sampleType" in globalOpts)
                gSampleDesc = globalOpts["sampleType"];
            if ("title" in globalOpts)
                gTitle = globalOpts["title"];
        }
            
        gOptions = opt;
        drawMenuBar();
        drawToolBar();
        setupKeyboard();
        //drawMenu();
        winInfo = createRenderer();
        drawLegendBar();
        loadDataset(0);
        menuBarHide("#tpShowAllButton");
    }

    function legendSort(dict, sortBy) {
    /* return a dictionary as reverse-sorted list (count, fieldValue).
    If there are more than 4 entries and all but one key are numbers, 
    sorts by key instead.
    */
        // convert the dict to a list of (count, key)
        var countList = [];
        var numCount = 0;
        var count = null;
        var useGradient = false; // use a rainbow or gradient palette?
        for (var key in dict) {
            count = dict[key];
            if (!isNaN(key)) // key looks like a number
                numCount++;
            countList.push( [count, key] );
        }

        // auto-detect: treaet labels as numbers if most values are numbers
        if ((sortBy===undefined && countList.length >= 4 && numCount >= countList.length-1)) {
            sortBy = "name";
            useGradient = true;
        }

        var newList = [];
        var isSortedByName = false;

        if (sortBy==="name") {
            isSortedByName = true;
            
            // change the keys in the list from strings to floats, if possible
            // this makes sure that the subsequent sorting is really by number
            for (var i = 0; i < countList.length; i++) {
                var el = countList[i];
                count = el[0];
                var label = el[1];
                if (label!=="" && !isNaN(label)) // label string looks like a number
                    label = parseFloat(label);  // -> convert label to float
                newList.push( [count, label] );
            }

            // sort by name/label
            newList.sort(function(a, b) { 
                //return a[1] > b[1];  // this does not work. Why?
                return a[1].localeCompare(b[1]); 
            }); 

            // change values back to strings
            var newList2 = [];
            for (i = 0; i < newList.length; i++) {
                var el = newList[i];
                count = el[0];
                var label = el[1].toString();
                newList2.push( [count, label] );
            }
            newList = newList2;
        }
        else {
            // sort this list by count
            newList = countList.sort(function(a, b){ return b[0] - a[0]; }); // reverse-sort by count
        }

        var ret = {};
        ret.list = newList;
        ret.isSortedByName = isSortedByName;
        ret.useGradient = useGradient;
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
        }
    }

    function findDotsWithMeta(metaIdx, hlValue) {
    /* return array of cellIds with a given meta value */
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

    function onLegendLabelClick(event) {
    /* called when user clicks on legend entry. updates gLegend and gSelCellIds */
        var legendId = parseInt(event.target.id.split("_")[1]);
        if (("lastClicked" in gLegend) && gLegend.lastClicked==legendId) {
            // user clicked the same entry as before: remove all highlights
            gLegend.lastClicked = null;
            gSelCellIds = {};
            $('#tpLegend_'+legendId).removeClass('tpLegendSelect');
            menuBarHide("#tpFilterButton");
            menuBarHide("#tpOnlySelectedButton");
            updateGeneBarColors([]);
        }
        else {
            var domEvent = event.originalEvent;
            if (!domEvent.shiftKey && !domEvent.ctrlKey && !domEvent.metaKey) {
                gSelCellIds = {}; // remove highlight
                $('.tpLegend').removeClass('tpLegendSelect');
            }
            var metaVal = gLegend.rows[legendId][4];
            gSelCellIds = findCellIdsForLegendId(gClasses, legendId, gSelCellIds);
            menuBarShow("#tpFilterButton");
            menuBarShow("#tpOnlySelectedButton");
            $('#tpLegend_'+legendId).addClass('tpLegendSelect');
            gLegend.lastClicked=legendId;
            updateGeneBarColors(keys(gSelCellIds));
        }
        plotDots();
        renderer.render(stage);
    }

    function onResetColorsClick (ev) {
    /* reset all colors in legend to defaults */
        // clear localStorage
        var rows = gLegend.rows;
        for (var i = 0; i < rows.length; i++) {
            var colorHex = rows[i][0];
            var defColor = rows[i][1];
            gLegend.rows[i][0] = defColor;
            var saveKey = rows[i][5];
            localStorage.removeItem(saveKey);
        }
        
        drawLegend();
        updateGeneBarColors(null);
        plotDots();
        renderer.render(stage);
    }

    function onSortByClick (ev) {
    /* flip the current legend sorting */
        var sortBy = null;
        var nextSortBy = null;
        if (gLegend.isSortedByName) {
            sortBy = "count";
            nextSortBy = "name";
        }
        else {
            sortBy = "name";
            nextSortBy = "count";
        }
        buildLegend(sortBy);
        drawLegend();
        jQuery("#tpSortBy").text("Sort by "+nextSortBy);
    }


    function onMetaRightClick (key, options) {
    /* user right-clicks on a meta field */
        var metaIdx = parseInt(options.$trigger[0].id.split("_")[1]);

        if (key==0) {
            gLabelMetaIdx = metaIdx;
            gClusterMids = null; // force recalc
            menuBarShow("#tpHideLabels");
            plotDots();
            renderer.render(stage);
        }
            
    }

    function onLegendRightClick (key, options) {
    /* user right-clicks on a legend label */
        var legendId = parseInt(options.$trigger[0].id.split("_")[1]);
        var cellIds = findCellIdsForLegendId(gClasses, legendId);
        var mode;
        if (key==0) // hide
            mode = "hide";
        else
            mode = "showOnly";
        filterCoordsAndUpdate(cellIds, mode);
    }

    function drawLegend(sortBy) {
    /* draws current legend as specified by gLegend.rows 
     * sortBy can be "name" or "count" or "undefined" (=auto-detect)
     * */
        if (gLegend.rows==undefined)
            return;

        $('#tpLegendContent').empty();

        var htmls = [];
        var colors = [];
        var rows = gLegend.rows;

        for (var i = 0; i < rows.length; i++) {
            var colorHex = rows[i][0];
            var label = rows[i][2];
            var count = rows[i][3];

            colors.push(colorHex); // save for later
            var classStr = "tpLegend";
            var line = "<div id='tpLegend_" +i+ "' class='" +classStr+ "'>";
            htmls.push(line);
            htmls.push("<input class='tpColorPicker' id='tpLegendColorPicker_"+i+"' />");

            classStr = "tpLegendLabel";
            label = label.replace(/_/g, " ");
            if (label.trim()==="") {
                label = "(empty)";
                classStr += " tpGrey";
            }

            htmls.push("<div class='"+classStr+"' id='tpLegendLabel_"+i+"'>");
            htmls.push(label);
            htmls.push("<div class='tpLegendCount'>"+count+"</div>");
            htmls.push("</div>");
            //htmls.push("<input class='tpLegendCheckbox' id='tpLegendCheckbox_"+i+"' type='checkbox' checked style='float:right; margin-right: 5px'>");
            htmls.push("</div>");
        }
        htmls.push("<span id='tpResetColors' style='color: #888; cursor:pointer; font-size:13px'>Reset colors</span>&emsp;");

        // add the "sort by" div
        var sortLabel = null;
        if (gLegend.isSortedByName===true)
            sortLabel = "Sort by count" // add the link to sort by the other possibility
        else
            sortLabel = "Sort by name"
        htmls.push("<span id='tpSortBy' style='color: #888; cursor:pointer; font-size:13px'>"+sortLabel+"</span>");

        var htmlStr = htmls.join("");
        $('#tpLegendContent').append(htmlStr);
        $('.tpLegend').click( onLegendLabelClick );
        $('.tpLegendLabel').attr( "title", "Click to select samples with this value. Shift click to select multiple values.");
        $('#tpResetColors').click( onResetColorsClick );
        $('#tpSortBy').click( onSortByClick );

        // setup the right-click menu
        var menuItems = [{name: "Hide "+gSampleDesc+"s with this value"}, {name:"Show only "+gSampleDesc+"s with this value"}];
        var menuOpt = {
            selector: ".tpLegend",
            items: menuItems, 
            className: 'contextmenu-customwidth',
            callback: onLegendRightClick 
        };
        $.contextMenu( menuOpt );

        // activate the color pickers
        for (var i = 0; i < colors.length; i++) {
            var opt = {
                hideAfterPaletteSelect : true,
                color : colors[i],
                showPalette: true,
                allowEmpty : true,
                showInput: true,
                preferredFormat: "hex",
                change: onColorPickerChange
                }
            $("#tpLegendColorPicker_"+i).spectrum(opt);
        }

        //$('#tpLegendTabs').tabs();
        $('#tpLegendTabs > ul').hide();
    }

    function onColorPickerChange(color) {
        /* called when user manually selects a color in the legend with the color picker */
        var rows = gLegend.rows;
        for (var i = 0; i < rows.length; i++) {
            var oldColorHex = rows[i][0];
            var defColorHex = rows[i][1];
            var newCol = $("#tpLegendColorPicker_"+i).spectrum("get");
            var newColHex = null;
            if (newCol==null)
                newColHex = oldColorHex.toHex();
            else
                newColHex = newCol.toHex();
            rows[i][0] = newColHex;

            // delete from storage if not needed anymore
            var saveKey = rows[i][5];
            if (newColHex == defColorHex)
                localStorage.removeItem(saveKey);
            else {
                // save the new color to localStorage
                localStorage.setItem(saveKey, newColHex);
                }
        }

        plotDots();
        renderer.render(stage);
        updateGeneBarColors(null);
    }

    function updateGeneBarColors(cellIds) {
    /* change the colors of the bottom gene bar to reflect the expression in cellIds */
        if (cellIds==null)
            cellIds = gGeneBarCells;
        else
            gGeneBarCells = cellIds;
        // translate the gene expression values to intensities
        var geneCount = geneFields.length;
        var intensities = newArray(geneCount, 0);

        for (var i = 0; i < cellIds.length; i++) {
            var cellId = cellIds[i];
            var exprRow = geneExpr[cellId];
            for (var j = 0; j < geneCount; j++) {
                var val = exprRow[j];
                if (val != null) {
                    var geneId = geneFields[j][0];
                    var decileRanges = gDeciles[geneId];
                    if (decileRanges==undefined) {
                        //console.log("Warning: No decile in json file for "+geneId);
                        continue;
                    }
                    var binIdx = findBin(decileRanges, val);
                    intensities[j] += binIdx;
                }
            }
        }

        // use the average binIdx for each gene
        if (cellIds.length != 0)
            for (var i = 0; i < intensities.length; i++) {
                intensities[i] = Math.round(intensities[i] / cellIds.length);
            }

        if (gDecileColors === null)
            gDecileColors = makeColorPalette(10, true);

        // get the colors from the current legend
        var nullColor = cNullColor;
        var colors = gDecileColors.slice();
        // if right now coloring by expression, the user may have defined colors manually
        if (gLegend.type=="gene") {
            nullColor = gLegend.rows[0][0];
            var rows = gLegend.rows;
            for (var i = 1; i < rows.length; i++)
                colors[i-1] = rows[i][0];
        }

        // use the average intensities to pick the colors
        for (var i = 0; i < geneFields.length; i++) {
            var fieldName = "#tpGeneBarCell_"+i;
            var foreColor = "black";
            var binIdx = intensities[i];
            var backColor;
            if (binIdx==0) {
                backColor = nullColor;
                foreColor = cNullForeground;
            }
            else
                backColor = colors[binIdx-0];

            if (binIdx >= 6)
                foreColor = "white";
            $(fieldName).css("background-color", "#"+backColor);
            $(fieldName).css("color", foreColor);
        }
    }

    function updateMetaBar(cellId) {
        // update the meta bar with meta data from cellId
        var cellInfo = metaData[cellId];
        for (var i = 0; i < metaFields.length; i++) {
            var fieldName = metaFields[i];
            $('#tpMeta_'+i).text(cellInfo[i]);
        }

    }

    function onDotMouseOver (mouseData) {
        /* user hovers over a circle with the mouse */
        if (mouseDownX!=null) // user is currently zooming
            return;
        if (! (gSelCellIds===null || jQuery.isEmptyObject(gSelCellIds))) // some cells are selected
            return;
        var cellId = mouseData.target.cellId;
        //this.alpha = 0.5;
        //renderer.render(stage);

        updateMetaBar(cellId)
        updateGeneBarColors([cellId]);
    }

    function onDotMouseClick (event) {
        /* user clicks onto a circle with the mouse */
        if (mouseDownX!=null) // user is currently zooming
            return;
        var cellId = event.target.cellId;
        var domEvent = event.data.originalEvent;
        if (!domEvent.shiftKey && !domEvent.ctrlKey && !domEvent.metaKey)
            gSelCellIds = {};
        gSelCellIds[cellId] = true;
        updateMetaBar(cellId)
        updateGeneBarColors([cellId]);
        plotDots();
        renderer.render(stage);
        event.stopPropagation();
    }

    function plotClusterLabels (metaFieldIdx) {
    /* plot semi-transparent labels for clusters */
        if (metaFieldIdx === null)
            return;

        if (gClusterMids === null) {
            // arrange the current coordinates to format cluster -> array of coords
            var clusterCoords = {};
            for (var i = 0; i < pixelCoords.length; i++) {
                var cellId = pixelCoords[i][0];
                var x = pixelCoords[i][1];
                var y = pixelCoords[i][2];
                var meta = metaData[cellId];
                var clusterLabel = meta[metaFieldIdx];
                if (clusterCoords[clusterLabel] == undefined)
                    clusterCoords[clusterLabel] = [];
                clusterCoords[clusterLabel].push( [x, y] );
            }

            // for each cluster, calculate the midpoints of all coords and append to array as (label, x, y)
            gClusterMids = [];
            for (var label in clusterCoords) {
                var coords = clusterCoords[label];

                var xSum = 0;
                var ySum = 0;
                for (var i = 0; i < coords.length; i++) {
                    xSum += coords[i][0];
                    ySum += coords[i][1];
                }

                var dotCount = clusterCoords[label].length;
                var midX = Math.floor(xSum / dotCount);
                var midY = Math.floor(ySum / dotCount);
                gClusterMids.push( [label, midX, midY] );            
            }
        }

        var txtStyle = {
            font : 'bold 13px Arial', 
            fill : 'black',
            stroke : '#ffffff',
            strokeThickness: 2,
            //dropShadow : true,
            //dropShadowColor : '#ffffff',
            //dropShadowDistance : 2,
        };
        var anchor = new PIXI.Point(0.5, 0.5); // center the text

        for (var i = 0; i < gClusterMids.length; i++) {
            var mid = gClusterMids[i];
            var label = mid[0].replace(/_/g, " ");
            var midX = mid[1];
            var midY = mid[2];

            var t = new PIXI.Text(label, txtStyle);
            t.x = midX;
            t.y = midY;
            t.anchor = anchor;
            t.alpha=0.8;
            stage.addChild(t);
        }
    }

    function plotDots() {
        /* plot the data onto the canvas. colorByIndx is the index of the meta field to color on. */
        if (stage!=null)
            {
            console.time("clear");
            stage.destroy(true); // free all children in OpenGL's GPU memory
            console.timeEnd("clear");
            }

        // a stage is a collection of drawables
        stage = new PIXI.Container();

       //var backgroundLayer = new PIXI.DisplayObjectContainer();
       var background = new PIXI.Graphics();
       background.beginFill(0x000000, 0.0);
       background.drawRect(0, 0, winInfo.width, winInfo.height);
       background.interactive = true;
       background.on("click", onBackgroundMouseClick);

       stage.addChild(background);

       //renderer.plugins.interaction.on('mousedown', onBackgroundMouseDown);
       //renderer.plugins.interaction.on('mousemove', onMouseMove );
       //renderer.plugins.interaction.on('mouseup', onBackgroundMouseUp);
       
        // write grey cell/sample count in upper left corner
        var selCount = Object.keys(gSelCellIds).length;
        var hideCount = allCoords.length - shownCoords.length;

        var textLines = [];
        var line1 = pixelCoords.length + ' '+gSampleDesc+'s';
        if (hideCount!=0)
            line1 += " shown";
        textLines.push(line1);

        if (selCount!=0)
            textLines.push(selCount+" selected");
        if (hideCount!=0)
            textLines.push(hideCount+" hidden");
        var topText = textLines.join("\n");

        var t = new PIXI.Text(topText, {'fontSize': '12px'});
        t.x = 1;
        t.y = 0;
        t.alpha=0.5;
        stage.addChild(t);

        var greyCol = parseInt("888888", 16);
        var classColors = null;

        if (gLegend!=null) {
            // need integer colors for pixi, one per legend row = bin
            classColors = {};
            var rows = gLegend.rows;
            for (var i = 0; i < rows.length; i++) {
                var colorInt = parseInt(rows[i][0], 16);
                classColors[i] = colorInt;
            }
        }

        var hlCoords = [];
        console.time("drawing");

        for (var i = 0; i < pixelCoords.length; i++) {
            var cellId = pixelCoords[i][0];
            var x = pixelCoords[i][1];
            var y = pixelCoords[i][2];
            var fill = null;
            if (classColors==null) {
                fill = greyCol;
                }
            else
                {
                var gClass = gClasses[cellId];
                fill = classColors[gClass];
                }
            //if (fill==-1) // -1 == hidden, do not draw
                //continue
            //var dot = drawRect(x, y, fill);

            if (cellId in gSelCellIds)
                {
                // don't draw selected cells now, keep them for later
                hlCoords.push([x, y, fill]);
                continue;
                }

            var dot = drawCircle(x, y, fill);
            dot.alpha = transparency;

            // setup mouse over 
            dot.interactive=true;
            dot.hitArea = new PIXI.Rectangle(x-(circleSize/2), y-(circleSize/2), circleSize, circleSize);
            dot.cellId = cellId;
            dot.mouseover = onDotMouseOver;
            //dot.mouseout = function(mouseData) {
              //this.alpha = 0.5;
            //}
            dot.on('mousedown', onDotMouseClick);
            stage.addChild(dot);
        }
        
        // finally, draw the currently highlighted dots on top
        plotSelection(hlCoords);

        plotClusterLabels(gLabelMetaIdx);

        console.timeEnd("drawing");
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
