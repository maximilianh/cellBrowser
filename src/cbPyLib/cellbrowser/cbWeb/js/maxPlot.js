'use strict';
// maxPlot: a fast scatter plot class
/*jshint globalstrict: true*/
/*jshint -W069 */
/*jshint -W104 */
/*jshint -W117 */

// TODO:
// fix mouseout into body -> marquee stays

function getAttr(obj, attrName, def) {
    var val = obj[attrName];
    if (val===undefined)
        return def;
    else
        return val;
}

function cloneObj(d) {
/* returns a deep copy of an object, wasteful and destroys old references */
    // see http://stackoverflow.com/questions/122102/what-is-the-most-efficient-way-to-deep-clone-an-object-in-javascript
    return JSON.parse(JSON.stringify(d));
}

function copyObj(src, trg) {
/* object copying: copies all values from src to trg */
    var key;
    for (key in src) {
        trg[key] = src[key]; // copies each property to the objCopy object
  }
}


function MaxPlot(div, top, left, width, height, args) {
    // a class that draws circles onto a canvas, like a scatter plot
    // div is a div DOM element under which the canvas will be created
    // top, left: position in pixels, integers
    // width and height: integers, in pixels, includes the status line

    const HIDCOORD = 12345; // magic value for missing coordinates
      // In rare instances, coordinates are saved but should not be shown. This way of implementing hiding
      // may look hacky, but it simplifies the logic and improves performance.

    var self = this; // 'this' has two conflicting meanings in javascript.
    // I use 'self' to refer to object variables, so I can use 'this' to refer to the caller context

    const gTextSize = 16; // size of cluster labels
    const gTitleSize = 18; // size of title text
    const gStatusHeight = 14; // height of status bar
    const gZoomButtonSize = 30; // size of zoom buttons
    const gZoomFromLeft = 10;  // position of zoom buttons from left
    const gZoomFromBottom = 140;  // position of zoom buttons from bottom
    const gButtonBackground = "rgb(230, 230, 230, 0.85)" // grey level of buttons
    const gButtonBackgroundClicked = "rgb(180, 180, 180, 0.6)"; // grey of buttons when clicked
    const gCloseButtonFromRight = 60; // distance of "close" button from right edge

    // the rest of the initialization is done at the end of this file,
    // because the init involves many functions that are not defined yet here

    this.initCanvas = function (div, top, left, width, height) {
        /* initialize a new Canvas */

        div.style.top = top+"px";
        div.style.left = left+"px";
        div.style.position = "absolute";
        div.style.display = "block";
        self.div = div;

        self.gSampleDescription = "cell";
        self.ctx = null; // the canvas context
        self.canvas = addCanvasToDiv(div, top, left, width, height-gStatusHeight );

        self.interact = false;

        if (args && args.showClose===true) {
            self.closeButton = addChildControls(10, width-gCloseButtonFromRight);
        }

        if (args===undefined || (args["interact"]!==false)) {
            self.interact = true;

            addZoomButtons(height-gZoomFromBottom, gZoomFromLeft, self);
            addModeButtons(10, 10, self);
            addStatusLine(height-gStatusHeight, left, width, gStatusHeight);
            addTitleDiv(height-gTitleSize-gStatusHeight, left+5);

            /* add the div used for the mouse selection/zoom rectangle to the DOM */
            var selectDiv = document.createElement('div');
            selectDiv.id = "mpSelectBox";
            selectDiv.style.border = "1px dotted black";
            selectDiv.style.position = "absolute";
            selectDiv.style.display  = "none";
            selectDiv.style.pointerEvents = "none";
            self.div.appendChild(selectDiv);

            // callbacks when user clicks or hovers over label or cell
            self.onLabelClick = null; // called on label click, args: text of label and event
            self.onCellClick = null; // called on cell click, args: array of cellIds and event
            self.onCellHover = null; // called on cell hover, arg: array of cellIds
            self.onNoCellHover = null; // called on hover over empty background
            self.onSelChange = null; // called when the selection has been changed, arg: array of cell Ids
            self.onLabelHover = null; // called when mouse hovers over a label
            self.onNoLabelHover = null; // called when mouse does not hover over a label
            // self.onZoom100Click: called when user clicks the zoom100 button. Implemented below.
            self.selectBox = selectDiv; // we need this later
            self.setupMouse();

            // connected plots
            self.childPlot = null;    // plot that is syncing from us, see split()
            self.parentPlot = null;   // plot that syncs to us, see split()
        }

        addProgressBars(top+Math.round(height*0.3), left+30);



        // timer that is reset on every mouse move
        self.timer = null;

        // all other object variables are added by the "initPlot(args)" function below

        // when the user starts to select but lifts the mouse button outside
        // the canvas, the current selection must be reset
        //self.canvas.addEventListener("mouseleave", function() {self.;});

        //document.body.addEventListener("mouseup", function(evt) {
            //var targetId = evt.target.id;
            //console.log("mouseup on body");
            //if(targetId !== "mpCanvas" && targetId!=="mpSelectBox") {  // do nothing on the canvas
                //console.log("mouseup outside of canvas, resetting marquee");
                //self.resetMarquee();
                //evt.stopPropagation();
            //}
        //});
        // for this to work, the body has to really cover the whole page
        //document.body.style.height = "100vh";
    }

    function isHidden(x, y) {
        /* special coords are used for circles that are off-screen or otherwise not visible */
       return ((x===HIDCOORD && y===HIDCOORD)) // not shown (e.g. no coordinate or off-screen)
    }

    this.initPort = function(args) {
        /* init all viewport related state (zoom, radius, alpha) */
        self.port = {};
        self.port.zoomRange = {}; // object with keys minX, , maxX, minY, maxY
        self.port.radius     = getAttr(args, "radius", null);    // current radius of the circles, 0=one pixel dots

        // we keep a copy of the 'initial' arguments at 100% zoom
        self.port.initZoom   = {};
        self.port.initRadius = self.port.radius;                      // circle radius at full zoom
        self.port.initAlpha   = getAttr(args, "alpha", 0.3);
    };

    this.initPlot = function(args) {
        /* create a new scatter plot on the canvas */
        if (args===undefined)
            args = {};

        self.globalOpts = args;

        self.mode = 1;   // drawing mode


        // everything related to circle coordinates
        self.coords = {};
        self.coords.orig = null;   // coordinates of cells in original coordinates
        self.coords.labels    = null;   // cluster label positions in pixels, array of [x,y,text]

        self.coords.px   = null;   // coordinates of cells and labels as screen pixels or (HIDCOORD,HIDCOORD) if not shown
        self.coords.labelBbox = null;   // cluster label bounding boxes, array of [x1,x2,x2,y2]


        self.col = {};
        self.col.pal = null;        // list of six-digit hex codes
        self.col.arr = null;        // length is coords.px/2, one byte per cell = index into self.col.pal

        self.selCells = new Set();  // IDs of cells that are selected (drawn in black)

        self.doDrawLabels = true;  // should cluster labels be drawn?
        self.initPort(args);

        // mouse drag is modal: can be "select", "move" or "zoom"
        self.dragMode = "select";

        // for zooming and panning
        self.mouseDownX = null;
        self.mouseDownY = null;
        self.panCopy    = null;

        // to detect if user just clicked on a dot
        self.dotClickX = null;
        self.dotClickY = null;

        self.activateMode(getAttr(args, "mode", "move"));
    };

    // call the constructor
    //self.newObject(div, top, left, width, height, args);

    this.clear = function() {
        clearCanvas(self.ctx, self.canvas.width, self.canvas.height);
    };

    //this.setPos = function(left, top) {
       /* position the canvas on the page */
       //self.div.style.left = left+"px";
       //self.div.style.top = top+"px";
    //};

    this.setTitle = function (text) {
        self.title = text;
        self.titleDiv.innerHTML = text;
    };


    // -- (private) helper functions
    // -- these are normal functions, not methods, they do not access "self"

    function gebi(idStr) {
        return document.getElementById(idStr);
    }

    function removeElById(idStr) {
        var el = gebi(idStr);
        if (el!==null) {
            el.parentNode.removeChild(el);
        }
    }

    function activateTooltip(selector) {
        /* uses bootstrap tooltip. Use noconflict in html, I had to rename BS's tooltip to avoid overwrite by jquery
         */
        if (window.jQuery && $.fn.bsTooltip!==undefined) {
            var ttOpt = {"html": true, "animation": false, "delay":{"show":400, "hide":100}, container:"body"};
            $(selector).bsTooltip(ttOpt);
        }
    }

    function guessRadius(coordCount) {
        /* a few rules to find a good initial radius, depending on the number of dots */
        if (coordCount > 50000)
            return 0;
        else if (coordCount > 10000)
            return 2;
        else if (coordCount > 4000)
            return 4;
        else
            return 5;
    }

    function createButton(width, height, id, title, text, imgFname, paddingTop, paddingBottom, addSep, addThickSep) {
        /* make a light-grey div that behaves like a button, with text and/or an image on it
         * Images are hard to vertically center, so padding top can be specified.
         * */
        var div = document.createElement('div');
        div.id = id;
        div.className = "mpButton";
        div.style.backgroundColor = gButtonBackground;
        div.style.width = width+"px";
        div.style.height = height+"px";
        div.style["text-align"]="center";
        div.style["vertical-align"]="middle";
        div.style["line-height"]=height+"px";
        if (text!==null)
            if (text.length>3)
                div.style["font-size"]="11px";
            else
                div.style["font-size"]="14px";
        div.style["font-weight"]="bold";
        div.style["font-family"]="sans-serif";

        if (title!==null)
            div.title = title;
        if (text!==null)
            div.textContent = text;
        if (imgFname!==null && imgFname!==undefined) {
            var img = document.createElement('img');
            img.src = imgFname;
            if (paddingTop!==null && paddingTop!==undefined) {
                img.style.paddingTop = paddingTop+"px";
            if (paddingBottom)
                img.style.paddingBottom = paddingBottom+"px";
            }
            div.appendChild(img);
        }
        if (addSep===true)
            div.style["border-bottom"] = "1px solid #D7D7D7";
        if (addThickSep===true)
            div.style["border-bottom"] = "2px solid #C7C7C7";


        // make color dark grey when mouse is pressed
        div.addEventListener("mousedown", function() {
                this.style.backgroundColor = gButtonBackgroundClicked;
        });

        div.addEventListener("mouseup", function() {
                this.style.backgroundColor = gButtonBackground;
        });
        return div;
    }

    function makeCtrlContainer(top, left) {
        /* make a container for half-transprent ctrl buttons over the canvas */
        var ctrlDiv = document.createElement('div');
        ctrlDiv.id = "mpCtrls";
        ctrlDiv.style.position = "absolute";
        ctrlDiv.style.left = left+"px";
        ctrlDiv.style.top = top+"px";
        ctrlDiv.style["border-radius"]="2px";
        ctrlDiv.style["cursor"]="pointer";
        ctrlDiv.style["box-shadow"]="0px 2px 4px rgba(0,0,0,0.3)";
        ctrlDiv.style["border-top-left-radius"]="2px";
        ctrlDiv.style["border-top-right-radius"]="2px";
        ctrlDiv.style["user-select"]="none";
        return ctrlDiv;
    }

    function addZoomButtons(top, left, self) {
        /* add the plus/minus buttons to the DOM and place at position x,y on the screen */
        var width = gZoomButtonSize;
        var height = gZoomButtonSize;

        var plusDiv = createButton(width, height, "mpCtrlZoomPlus", "Zoom in. Keyboard: +", "+", null, null, null, true);
        //plusDiv.style["border-bottom"] = "1px solid #D7D7D7";

        var fullDiv = createButton(width, height, "mpCtrlZoom100", "Zoom in. Keyboard: space", "100%", null, null, null, true);
        //full.style["border-bottom"] = "1px solid #D7D7D7";

        var minusDiv = createButton(width, height, "mpCtrlZoomMinus", "Zoom out. Keyboard: -", "-");

        var ctrlDiv = makeCtrlContainer(top, left);
        ctrlDiv.appendChild(plusDiv);
        ctrlDiv.appendChild(fullDiv);
        ctrlDiv.appendChild(minusDiv);
        self.zoomDiv = ctrlDiv;

        self.div.appendChild(ctrlDiv);

        minusDiv.addEventListener('click', function() { self.zoomBy(0.75); self.drawDots(); });
        fullDiv.addEventListener('click', function() { self.zoom100(); self.drawDots()});
        plusDiv.addEventListener('click', function() { self.zoomBy(1.333); self.drawDots(); });
    }

    function addTitleDiv(top, left) {
        var div = document.createElement('div');
        div.style.cursor = "default";
        div.style.left = "4px";
        div.style.width = "max-content"; // not in MS Edge
        div.style.top = top+"px";
        div.style.display = "block";
        div.style.position = "absolute";
        div.style.fontSize = gTitleSize;
        div.id = 'mpTitle';
        div.style['color'] = "#B0B0B0";
        self.div.appendChild(div);
        self.titleDiv = div;
    }

    function addCloseButton(top, left) {
        /* add close button and sync checkbox */
        var div = document.createElement('div');
        div.style.cursor = "default";
        div.style.left = left+"px";
        div.style.top = top+"px";
        div.style.display = "block";
        div.style.position = "absolute";
        div.style.fontSize = gTitleSize;
        div.style.padding = "3px";
        div.style.borderRadius = "3px";
        div.style.border = "1px solid #c5c5c5";
        div.style.backgroundColor = "#f6f6f6";
        div.style.color = "#454545";
        div.id = 'mpCloseButton';
        div.textContent = "Close";
        self.div.appendChild(div);
        return div;
    }

    function addChildControls(top, left) {
        addCloseButton(top, left);
    }

    function appendButton(parentDiv, id, title, imgName) {
        /* add a div styled like a button under div */
        var div = document.createElement('div');
        div.title = title;
        div.id = id;

    }

    function addModeButtons(top, left, self) {
        /* add the zoom/move/select control buttons to the DOM */
        var ctrlDiv = makeCtrlContainer(top, left);

        var bSize = gZoomButtonSize;

        var selectButton = createButton(bSize, bSize, "mpIconModeSelect", "Select mode. Keyboard: shift or s", null, "img/select.png", 0, 4, true, true);
        selectButton.addEventListener ('click',  function() { self.activateMode("select")}, false);

        var zoomButton = createButton(bSize, bSize, "mpIconModeZoom", "Zoom-to-rectangle mode. Keyboard: Windows/Command or z", null, "img/zoom.png", 4, 4, true);
        zoomButton.addEventListener ('click', function() { self.activateMode("zoom")}, false);

        var moveButton = createButton(bSize, bSize, "mpIconModeMove", "Move mode. Keyboard: Alt or m", null, "img/move.png", 4, 4);
        moveButton.addEventListener('click', function() { self.activateMode("move");}, false);

        self.icons = {};
        self.icons["move"] = moveButton;
        self.icons["select"] = selectButton;
        self.icons["zoom"] = zoomButton;

        //ctrlDiv.innerHTML = htmls.join("");
        ctrlDiv.appendChild(moveButton);
        ctrlDiv.appendChild(selectButton);
        ctrlDiv.appendChild(zoomButton);

        self.div.appendChild(ctrlDiv);
        self.toolDiv = ctrlDiv;

        activateTooltip('.mpIconButton');
    }

    function setStatus(text) {
        self.statusLine.innerHTML = text;
    }

    function addStatusLine(top, left, width, height) {
        /* add a status line div */
        var div = document.createElement('div');
        div.id = "mpStatus";
        div.style.backgroundColor = "rgb(240, 240, 240)";
        div.style.position = "absolute";
        div.style.top = top+"px";
        //div.style.left = left+"px";
        div.style.width = width+"px";
        div.style.height = height+"px";
        div.style["border-left"]="1px solid #DDD";
        div.style["border-right"]="1px solid #DDD";
        div.style["border-top"]="1px solid #DDD";
        div.style["font-size"]=(gStatusHeight-1)+"px";
        div.style["cursor"]="pointer";
        div.style["font-family"] = "sans-serif";
        self.div.appendChild(div);
        self.statusLine = div;
    }

    function addProgressBars(top, left) {
       /* add the progress bar DIVs to the DOM */
       var div = document.createElement('div');
       div.id = "mpProgressBars";
       div.style.top = top+"px";
       div.style.left = left+"px";
       div.style.position = "absolute";

       var htmls = [];
       for (var i=0; i<3; i++) {
           htmls.push('<div id="mpProgressDiv'+i+'" style="display:none; height:17px; width:300px; background-color: rgba(180, 180, 180, 0.3)" style="">');
           htmls.push('<div id="mpProgress'+i+'" style="background-color:#666; height:17px; width:10%"></div>');
           htmls.push('<div id="mpProgressLabel'+i+'" style="color:white; line-height:17px; position:absolute; top:'+(i*17)+'px;left:100px">Loading...</div>');
           htmls.push('</div>');
       }

       div.innerHTML = htmls.join("");
       self.div.appendChild(div);
    }

    function addCanvasToDiv(div, top, left, width, height) {
        /* add a canvas element to the body element of the current page and keep left/top/width/eight in self */
        var canv = document.createElement('canvas');
        canv.id = 'mpCanvas';
        //canv.style.border = "1px solid #AAAAAA";
        canv.style.backgroundColor = "white";
        canv.style.position = "absolute";
        canv.style.display = "block";
        canv.style.width = width+"px";
        canv.style.height = height+"px";
        //canv.style.top = top+"px";
        //canv.style.left = left+"px";
        // No scaling = one unit on screen is one pixel. Essential for speed.
        canv.width = width;
        canv.height = height;

        // need to keep these as ints, need them all the time
        self.width = width;
        self.height = height;
        self.top = top; // location of the canvas in pixels
        self.left = left;

        div.appendChild(canv); // adds the canvas to the div element
        self.canvas = canv;
        // alpha:false recommended by https://developer.mozilla.org/en-US/docs/Web/API/Canvas_API/Tutorial/Optimizing_canvas
        self.ctx = self.canvas.getContext("2d", { alpha: false });
        // by default, the canvas background is transparent+black
        // we use alpha=false, so we need to initialize the canvas with white pixels
        clearCanvas(self.ctx, width, height);

        return canv;
    }

    function scaleLabels(labels, zoomRange, borderSize, winWidth, winHeight) {
        /* scale cluster label position to pixel coordinates */
        winWidth = winWidth-(2*borderSize);
        winHeight = winHeight-(2*borderSize);

        var minX = zoomRange.minX;
        var maxX = zoomRange.maxX;
        var minY = zoomRange.minY;
        var maxY = zoomRange.maxY;

        var spanX = maxX - minX;
        var spanY = maxY - minY;
        var xMult = winWidth / spanX;
        var yMult = winHeight / spanY;

        // scale the label coords
        var pxLabels = [];
        for (var i = 0; i < labels.length; i++) {
            var annot = labels[i];
            var x = annot[0];
            var y = annot[1];
            var text = annot[2];
            // XX ignore anything outside of current zoom range. Performance?
            if (isHidden(x,y) || (x < minX) || (x > maxX) || (y < minY) || (y > maxY)) {
                pxLabels.push(null);
            }
            else {
                var xPx = Math.round((x-minX)*xMult)+borderSize;
                var yPx = winHeight - Math.round((y-minY)*yMult)+borderSize;
                pxLabels.push([xPx, yPx, text]);
            }
        }
        return pxLabels;
    }

    function constrainVal(x, min, max) {
        /* if x is not in range min, max, limit to min or max */
        if (x < min)
            return min;
        if (x > max)
            return max;
        return x;
    }

    function scaleLines(lines, zoomRange, winWidth, winHeight) {
        /* scale an array of (x1, y1, x2, y2), cutting lines at the screen edges */
        var minX = zoomRange.minX;
        var maxX = zoomRange.maxX;
        var minY = zoomRange.minY;
        var maxY = zoomRange.maxY;

        var spanX = maxX - minX;
        var spanY = maxY - minY;
        var xMult = winWidth / spanX;
        var yMult = winHeight / spanY;

        // transform from data floats to screen pixel coordinates
        var pxLines = [];
        for (var lineIdx = 0; lineIdx < lines.length; lineIdx++) {
            var line = lines[lineIdx];
            var x1 = line[0];
            var y1 = line[1];
            var x2 = line[2];
            var y2 = line[3];

            var startInvis = ((x1 < minX) || (x1 > maxX) || (y1 < minY) || (y1 > maxY));
            var endInvis = ((x2 < minX) || (x2 > maxX) || (y2 < minY) || (y2 > maxY));

            // line is entirely hidden
            if (startInvis && endInvis)
                continue
            if (startInvis) {
                x1 = constrainVal(x1, minX, maxX);
                y1 = constrainVal(y1, minY, maxY);
            }
            if (endInvis) {
                x2 = constrainVal(x2, minX, maxX);
                y2 = constrainVal(y2, minY, maxY);
            }

            var x1Px = Math.round((x1-minX)*xMult);
            var y1Px = winHeight - Math.round((y1-minY)*yMult);
            var x2Px = Math.round((x2-minX)*xMult);
            var y2Px = winHeight - Math.round((y2-minY)*yMult);
            pxLines.push( [x1Px, y1Px, x2Px, y2Px] );
        }
        return pxLines;
    }

    function scaleCoords(coords, borderSize, zoomRange, winWidth, winHeight, annots) {
    /* scale list of [x (float),y (float)] to integer pixels on screen and
     * annots is an array with on-screen annotations in the format (x, y,
     * otherInfo) that is also scaled.  return [array of (x (int), y (int)),
     * scaled annots array]. Take into account the current zoom range.      *
     * Canvas origin is top-left, but usually plotting origin is bottom-left,
     * so also flip the Y axis. sets invisible coords to HIDCOORD
     * */
        console.time("scale");
        var minX = zoomRange.minX;
        var maxX = zoomRange.maxX;
        var minY = zoomRange.minY;
        var maxY = zoomRange.maxY;

        winWidth = winWidth-(2*borderSize);
        winHeight = winHeight-(2*borderSize);

        var spanX = maxX - minX;
        var spanY = maxY - minY;
        var xMult = winWidth / spanX;
        var yMult = winHeight / spanY;

        // transform from data floats to screen pixel coordinates
        var pixelCoords = new Uint16Array(coords.length);
        for (var i = 0; i < coords.length/2; i++) {
            var x = coords[i*2];
            var y = coords[i*2+1];
            // set everything outside of current zoom range as hidden
            if ((x < minX) || (x > maxX) || (y < minY) || (y > maxY)) {
                pixelCoords[2*i] = HIDCOORD; // see isHidden()
                pixelCoords[2*i+1] = HIDCOORD;
            }
            else {
                var xPx = Math.round((x-minX)*xMult)+borderSize;
                // flipY: y-axis is flipped, so we do winHeight - pixel value
                var yPx = winHeight - Math.round((y-minY)*yMult)+borderSize;
                pixelCoords[2*i] = xPx;
                pixelCoords[2*i+1] = yPx;
            }
        }

        console.timeEnd("scale");
        return pixelCoords;
    }

    function drawRect(ctx, pxCoords, coordColors, colors, radius, alpha, selCells) {
        /* draw not circles but tiny rectangles. Maybe good enough for 2pixels sizes */
       console.log("Drawing "+coordColors.length+" rectangles, with fillRect");
       ctx.save();
       ctx.globalAlpha = alpha;
       var dblSize = 2*radius;
       var count = 0;
       for (var i = 0; i < pxCoords.length/2; i++) {
           var pxX = pxCoords[2*i];
           var pxY = pxCoords[2*i+1];
           if (isHidden(pxX, pxY))
               continue;
           var col = colors[coordColors[i]];
           ctx.fillStyle="#"+col;
           ctx.fillRect(pxX-radius, pxY-radius, dblSize, dblSize);
           count++;
       }

       // draw the selection as black rectangles
       ctx.globalAlpha = 0.7;
       ctx.fillStyle="black";
        selCells.forEach(function(cellId) {
           let pxX = pxCoords[2*cellId];
           let pxY = pxCoords[2*cellId+1];
           ctx.fillRect(pxX-radius, pxY-radius, dblSize, dblSize);
           count += 1;
        })
       console.log(count+" rectangles drawn (including selection)");
       ctx.restore();
       return count;
    }

    function drawCirclesStupid(ctx, pxCoords, coordColors, colors, radius, alpha, selCells) {
    /* draw little circles onto canvas. pxCoords are the centers.  */
       console.log("Drawing "+coordColors.length+" circles with stupid renderer");
       ctx.globalAlpha = alpha;
       var dblSize = 2*radius;
       var count = 0;
       for (var i = 0; i < pxCoords.length/2; i++) {
           var pxX = pxCoords[2*i];
           var pxY = pxCoords[2*i+1];
           if (isHidden(pxX, pxY))
               continue;
           var col = colors[coordColors[i]];
           ctx.fillStyle="#"+col;
           //ctx.fillRect(pxX-size, pxY-size, dblSize, dblSize);
           ctx.beginPath();
           ctx.arc(pxX, pxY, radius, 0, 2 * Math.PI);
           ctx.closePath();
           ctx.fill();
           count++;
       }
       return count;
    }

    function intersectRect(r1left, r1right, r1top, r1bottom, r2left, r2right, r2top, r2bottom) {
      /* return true if two rectangles overlap,
       https://stackoverflow.com/questions/2752349/fast-rectangle-to-rectangle-intersection
	*/
      return !(r2left > r1right || r2right < r1left || r2top > r1bottom || r2bottom < r1top);
    }

    function drawLines(ctx, pxLines, width, height, attrs) {
        /* draw lines defined by array with (x1, y1, x2, y2) arrays.
         * color is a CSS name, so usually prefixed by # if a hexcode
         * width is the width in pixels.
         * */
        ctx.save();
        //ctx.globalAlpha = 1.0;

        ctx.strokeStyle = attrs.lineColor || "#888888";
        ctx.lineWidth = attrs.lineWidth || 3;
        ctx.globalAlpha = attrs.lineAlpha || 0.5;
        //ctx.miterLimit =2;
        //ctx.strokeStyle = "rgba(200, 200, 200, 0.3)";

        for (var i=0; i < pxLines.length; i++) {
            var line = pxLines[i];
            var x1 = line[0];
            var y1 = line[1];
            var x2 = line[2];
            var y2 = line[3];

            ctx.beginPath();
            ctx.moveTo(x1, y1);
            ctx.lineTo(x2, y2);
            ctx.stroke();
        }
        ctx.restore();
    }

    function drawLabels(ctx, labelCoords, winWidth, winHeight, zoomFact) {
        /* given an array of [x, y, text], draw the text. returns bounding
         * boxes as array of [x1, y1, x2, y2]  */

        console.time("labels");
        ctx.save();
        ctx.font = "bold "+gTextSize+"px Sans-serif"
        ctx.globalAlpha = 1.0;

        ctx.strokeStyle = '#EEEEEE';
        ctx.lineWidth = 5;
        ctx.miterLimit =2;
        ctx.strokeStyle = "rgba(200, 200, 200, 0.3)";
        ctx.textBaseline = "top";

        ctx.shadowBlur=6;
        ctx.shadowColor="white";
        ctx.fillStyle = "rgba(0,0,0,0.8)";
        ctx.textAlign = "left";

        var addMargin = 1; // how many pixels to extend the bbox around the text, make clicking easier
        var bboxArr = []; // array of click hit boxes

        for (var i=0; i < labelCoords.length; i++) {
            var coord = labelCoords[i];
            if (coord===null) { // outside of view range, push a null to avoid messing up the order of bboxArr
                bboxArr.push( null );
                continue;
            }

            var x = coord[0];
            var y = coord[1];
            var text = coord[2];

            var textWidth = Math.round(ctx.measureText(text).width);
            // move x to the left, so text is centered on x
            x = x - Math.round(textWidth*0.5);

            var textX1 = x;
            var textY1 = y;
            var textX2 = Math.round(x+textWidth);
            var textY2 = y+gTextSize;

            //if (zoomFact===1.0) {
                // at 100% zoom, make some minimal effort to keep labels on screen
                //if (x < 0)
                    //x = 0;
                //if ((x + textWidth) > winWidth)
                    //x = winWidth - textWidth;
                //if (y+gTextSize > winHeight)
                    //y = winHeight-gTextSize;

                // also only at 100% zoom, make a minimal effort to avoid label overlaps
                // a perfect solution would take much more time
                //for (var j=0; j < bboxArr.length; j++) {
                    //var bbox = bboxArr[j];
                    //if (bbox===null) // = outside of screen
                        //continue;
                    //var bx1 = bbox[0];
                    //var by1 = bbox[1];
                    //var bx2 = bbox[2];
                    //var by2 = bbox[3];
                    //if (intersectRect(textX1, textX2, textY1, textY2, bx1, bx2, by1, by2)) {
                            // push the overlapping label away a little
                            //var diff = Math.round(0.75*gTextSize);
                            //if (textY1 < by1)
                                //y -= diff;
                            //else
                                //y += diff;
                        //}
                //}
            //}
            // don't draw labels where the midpoint is off-screen
            if (x<0 || y<0 || x>winWidth || y>winHeight) {
                bboxArr.push( null );
                continue;
            }


            ctx.strokeText(text,x,y);
            ctx.fillText(text,x,y);

            bboxArr.push( [textX1-addMargin, textY1-addMargin, textX2+addMargin, textY2+addMargin] );
        }
        ctx.restore();
        console.timeEnd("labels");
        return bboxArr;
    }

    // function drawLabels_dom(ctx, labelCoords, isFull) {
    //     /* given an array of [x, y, text], draw the text. returns bounding boxes as array of [x1, y1, x2, y2]  */
    //     for (var i=0; i < labelCoords.length; i++) {
    //         var coord = labelCoords[i];
    //         var x = coord[0];
    //         var y = coord[1];
    //         var text = coord[2];

    //         var div = document.createElement('div');
    //         div.id = id;
    //         //div.style.border = "1px solid #DDDDDD";
    //         div.style.backgroundColor = "rgb(230, 230, 230, 0.6)";
    //         div.style.width = width+"px";
    //         div.style.height = height+"px";
    //         div.style["text-align"]="center";
    //         div.style["vertical-align"]="middle";
    //         div.style["line-height"]=height+"px";
    //     }
    //     ctx.restore();
    //     return bboxArr;
    // }

    // https://stackoverflow.com/questions/5560248/programmatically-lighten-or-darken-a-hex-color-or-rgb-and-blend-colors
    function shadeColor(color, percent) {
        var f=parseInt(color,16),t=percent<0?0:255,p=percent<0?percent*-1:percent,R=f>>16,G=f>>8&0x00FF,B=f&0x0000FF;
        return (0x1000000+(Math.round((t-R)*p)+R)*0x10000+(Math.round((t-G)*p)+G)*0x100+(Math.round((t-B)*p)+B)).toString(16).slice(1);
    }

    function drawCirclesDrawImage(ctx, pxCoords, coordColors, colors, radius, alpha, selCells) {
    /* predraw and copy circles into canvas. pxCoords are the centers.  */
       // almost copied from by https://stackoverflow.com/questions/13916066/speed-up-the-drawing-of-many-points-on-a-html5-canvas-element
       // around 2x faster than drawing full circles
       // create an off-screen canvas

       ctx.save();
       console.log("Drawing "+coordColors.length+" coords with drawImg renderer, radius="+radius);
       var off = document.createElement('canvas'); // not added to DOM, will be gc'ed
       var diam = Math.round(2 * radius);
       var tileWidth = diam + 2; // must add one pixel on each side, space for antialising
       var tileHeight = tileWidth; // otherwise circles look cut off
       off.width = (colors.length+1) * tileWidth;
       off.height = tileHeight;
       var ctxOff = off.getContext('2d');

       //pre-render circles into the off-screen canvas.
       for (var i = 0; i < colors.length; ++i) {
           //ctxOff.lineWidth=1;
           ctxOff.fillStyle = "#"+colors[i];
           ctxOff.beginPath();
           // arc(x, y, r, 0, 2*pi)
           ctxOff.arc(i * tileWidth + radius +1, radius+1, radius, 0, 2 * Math.PI);
           ctxOff.closePath();
           ctxOff.fill();

           // only draw outline for big circles
           ctxOff.lineWidth=1.0;
           if (radius>5) {
               var strokeCol = "#"+shadeColor(colors[i], 0.9);
               ctxOff.strokeStyle=strokeCol;

               ctxOff.beginPath();
               ctxOff.arc(i * tileWidth + radius + 1, radius +1, radius, 0, 2 * Math.PI);
               ctxOff.closePath();
               ctxOff.stroke();
           }

       }

       // pre-render a black circle outline for the selection, quality of anti-aliasing?
       var selImgId = colors.length;
       ctxOff.lineWidth=2;
       ctxOff.strokeStyle="black";
       ctxOff.beginPath();
       // args: arc(x, y, r, 0, 2*pi)
       ctxOff.arc((selImgId * tileWidth) + radius +1, radius+1, radius-1, 0, 2 * Math.PI);
       ctxOff.stroke();

       if (alpha!==undefined)
           ctx.globalAlpha = alpha;

       // blit the circles onto the main canvas
       var count = 0;
       for (let i = 0; i < pxCoords.length/2; i++) {
           var pxX = pxCoords[2*i];
           var pxY = pxCoords[2*i+1];
           if (isHidden(pxX, pxY))
               continue;
           var col = coordColors[i];
           count++;
           // drawImage(img,sx,sy,swidth,sheight,x,y,width,height);
           ctx.drawImage(off, col * tileWidth, 0, tileWidth, tileHeight, pxX - radius - 1, pxY - radius - 1, tileWidth, tileHeight);
       }

       // overdraw the selection as solid black circle outlines
       ctx.globalAlpha = 0.7;
        selCells.forEach(function(cellId) {
           let pxX = pxCoords[2*cellId];
           let pxY = pxCoords[2*cellId+1];
           if (isHidden(pxX, pxY))
                return;
           // make sure that old leftover overlapping black circles don't shine through
           let col = coordColors[cellId];
           ctx.drawImage(off, col * tileWidth, 0, tileWidth, tileHeight, pxX - radius -1, pxY - radius-1, tileWidth, tileHeight);

           ctx.drawImage(off, selImgId * tileWidth, 0, tileWidth, tileHeight, pxX - radius -1, pxY - radius-1, tileWidth, tileHeight);
        })

       console.log(count +" circles drawn");
       ctx.restore();
       return count;
    }

    function hexToInt(colors) {
    /* convert a list of hex values to ints */
        var intList = [];
        for (var i = 0; i < colors.length; i++) {
            var colHex = colors[i];
            var colInt = parseInt(colHex, 16);
            if (colInt===undefined) {
                alert("Illegal color value, not a six-digit hex code: "+colHex);
                intList.push(0);
            }
            else
                intList.push(colInt);
        }
        return intList;
    }

    function drawRectBuffer(ctx, width, height, pxCoords, colorArr, colors, alpha, selCells) {
        /* Draw little rectangles with size 3 using a memory buffer*/
       var canvasData = ctx.getImageData(0, 0, width, height);
       var cData = canvasData.data;

       var rgbColors = hexToInt(colors);
       var invAlpha = 1.0 - alpha;

       // alpha-blend pixels into array
       for (var i = 0; i < pxCoords.length/2; i++) {
           var pxX = pxCoords[2*i];
           var pxY = pxCoords[2*i+1];
           if (isHidden(pxX, pxY))
               continue;
           var p = 4 * (pxY*width+pxX); // pointer to red value of pixel at x,y

           var oldR = cData[p];
           var oldG = cData[p+1];
           var oldB = cData[p+2];

           var newRgb = rgbColors[colorArr[i]];
           var newR = (newRgb >>> 16) & 0xff;
           var newG = (newRgb >>> 8)  & 0xff;
           var newB = (newRgb)        & 0xff;

           var mixR = ~~(oldR * invAlpha + newR * alpha);
           var mixG = ~~(oldG * invAlpha + newG * alpha);
           var mixB = ~~(oldB * invAlpha + newB * alpha);

           cData[p] = mixR;
           cData[p+1] = mixG;
           cData[p+2] = mixB;
           cData[p+3] = 255; // no transparency... ever?
       }
    }

    function drawPixels(ctx, width, height, pxCoords, colorArr, colors, alpha, selCells) {
        /* draw single pixels into a pixel buffer and copy the buffer into a canvas */

       // by default the canvas has black pixels
       // so not doing: var canvasData = ctx.createImageData(width, height);
       // XX is this really faster than manually zero'ing the array?
       var canvasData = ctx.getImageData(0, 0, width, height);
       var cData = canvasData.data;

       var rgbColors = hexToInt(colors);
       var invAlpha = 1.0 - alpha;

       var count = 0;

       // alpha-blend pixels into array
       for (var i = 0; i < pxCoords.length/2; i++) {
           var pxX = pxCoords[2*i];
           var pxY = pxCoords[2*i+1];
           if (isHidden(pxX, pxY))
               continue;
           var p = 4 * (pxY*width+pxX); // pointer to red value of pixel at x,y

           var oldR = cData[p];
           var oldG = cData[p+1];
           var oldB = cData[p+2];

           var newRgb = rgbColors[colorArr[i]];
           var newR = (newRgb >>> 16) & 0xff;
           var newG = (newRgb >>> 8)  & 0xff;
           var newB = (newRgb)        & 0xff;

           var mixR = ~~(oldR * invAlpha + newR * alpha);
           var mixG = ~~(oldG * invAlpha + newG * alpha);
           var mixB = ~~(oldB * invAlpha + newB * alpha);

           cData[p] = mixR;
           cData[p+1] = mixG;
           cData[p+2] = mixB;
           cData[p+3] = 255; // no transparency... ever?
           count++;
       }

       // overdraw the selection as black pixels
        selCells.forEach(function(cellId) {
           let pxX = pxCoords[2*cellId];
           let pxY = pxCoords[2*cellId+1];
           if (isHidden(pxX, pxY))
                return;
           let p = 4 * (pxY*width+pxX); // pointer to red value of pixel at x,y
           cData[p] = 0;
           cData[p+1] = 0;
           cData[p+2] = 0;
        })

       self.ctx.putImageData(canvasData, 0, 0);
       return count;
    }

    function findRange(coords) {
    /* find range of pairs-array and return obj with attributes minX/maxX/minY/maxY */
        var minX = 9999999;
        var maxX = -9999999;
        var minY = 9999999;
        var maxY = -9999999;

        for (var i = 0; i < coords.length/2; i++) {
            var x = coords[i*2];
            var y = coords[i*2+1];

            minX = Math.min(minX, x);
            maxX = Math.max(maxX, x);

            minY = Math.min(minY, y);
            maxY = Math.max(maxY, y);
        }

        var obj = {};
        obj.minX = minX;
        obj.maxX = maxX;
        obj.minY = minY;
        obj.maxY = maxY;
        return obj; // not needed, but more explicit
    }

    function clearCanvas(ctx, width, height) {
    /* clear with a white background */
        // jsperf says this is fastest on Chrome, and still OK-ish in FF
        //console.time("clear");
        ctx.save();
        ctx.globalAlpha = 1.0;
        ctx.fillStyle = "rgb(255,255,255)";
        ctx.fillRect(0, 0, width, height);
        ctx.restore();
        //console.timeEnd("clear");
    }

    // -- object methods (=access the self object)

    this.onZoom100Click = function(ev) {
        self.zoom100();
        self.drawDots();
    };

    this.scaleData = function() {
       /* scale coords and labels to current zoom range, write results to pxCoords and pxLabels */
       if (!self.coords) // window resize can call this before coordinates are loaded.
           return;

       var borderMargin = self.port.radius;
       self.calcRadius();

       self.coords.px = scaleCoords(self.coords.orig, borderMargin, self.port.zoomRange, self.canvas.width, self.canvas.height);
       if (self.coords.lines)
           self.coords.pxLines = scaleLines(self.coords.lines, self.port.zoomRange, self.canvas.width, self.canvas.height);
    }

    this.setTopLeft = function(top, left) {
        /* set top and left position in pixels of the canvas */
        self.top = top;
        self.left = left; // keep an integer version of these numbers
        self.div.style.top = top+"px";
        self.div.style.left = left+"px";

        self.setSize(self.width, self.height, false); // resize the various buttons
    }

    this.quickResize = function(width, height) {
       /* resize the canvas and move the status line, don't rescale or draw  */
       self.div.style.width = width+"px";
       self.div.style.height = height+"px";

       if (self.childPlot) {
           width = width/2;
           //self.childPlot.left = self.left+width;
           //self.childPlot.canvas.style.left = self.childPlot.left+"px";
           self.childPlot.setPos(null, self.left+width);
           self.childPlot.setSize(width, height, true);
       }

       if (self.closeButton) {
           self.closeButton.style.left = width - gCloseButtonFromRight;
       }

       // css and actual canvas sizes: these must be identical, otherwise canvas gets super slow
       self.canvas.style.width = width+"px";
       self.width = width;
       self.height = height;
       //let canvHeight = height - gStatusHeight;

       let canvHeight = height - gStatusHeight;
       self.canvas.height = canvHeight;
       self.canvas.width = width;
       self.canvas.style.height = canvHeight+"px";
       self.zoomDiv.style.top = (height-gZoomFromBottom)+"px";
       self.zoomDiv.style.left = (gZoomFromLeft)+"px";

       var statusDiv = self.statusLine;
       statusDiv.style.top = (height-gStatusHeight)+"px";
       statusDiv.style.width = width+"px";

       self.titleDiv.style.top = (height-gStatusHeight-gTitleSize)+"px";

    }

    this.setPos = function(top, left) {
       /* position canvas. Does not affect child  */
       if (top) {
          self.top = top;
          self.div.style.top = top+"px";
       }
       if (left) {
          self.left = left;
          self.div.style.left = left+"px";
       }
    }

    this.setSize = function(width, height, doRedraw) {
       /* resize canvas on the page re-scale the data and re-draw, unless doRedraw is false */
       if (width===null)
           width = self.div.getBoundingClientRect().width;

       self.quickResize(width, height);

       if (self.coords)
           self.scaleData();
       //clearCanvas(self.ctx, width, height);
       if (doRedraw===undefined || doRedraw===true)
           self.drawDots();
    };

    this.setCoords = function(coords, clusterLabels, minX, maxX, minY, maxY, opts) {
       /* specify new coordinates of circles to draw, an array of (x,y) coordinates */
       /* Scale data to current screen dimensions */
       /* clusterLabels is optional: array of [x, y, labelString]*/
       /* minX, maxX, etc are optional
        * opts are optional arguments like radius, alpha etc, see initPlot/args */
       if (coords.length === 0)
           alert("cbDraw-setCoords called with no coordinates");

       var coordOpts = cloneObj(self.globalOpts);
       copyObj(opts, coordOpts);
       // XX
       var oldRadius = self.port.initRadius;
       var oldAlpha = self.port.initAlpha;
       var oldLabels = self.coords.pxLabels;
       self.port = {};
       self.initPort(coordOpts);
       if (oldRadius)
           self.port.initRadius = oldRadius;
       if (oldAlpha)
           self.port.initAlpha = oldAlpha;
       self.coords = {};
       // XX


       var newZr = {};
       if (minX===undefined || maxX===undefined || minY===undefined || maxY===undefined)
           newZr = findRange(coords);
       else {
           newZr = {minX:minX, maxX:maxX, minY:minY, maxY:maxY};
       }
       copyObj(newZr, self.port.initZoom);
       copyObj(newZr, self.port.zoomRange);

       self.coords.orig = coords;
       self.coords.labels = clusterLabels;

       var count = 0;
       for (var i = 0; i < coords.length/2; i++) {
           var cellX = coords[i*2];
           var cellY = coords[i*2+1];
           if (!(isHidden(cellX, cellY)))
               count++;
       }

       //setStatus((coords.length/2)+" "+self.gSampleDescription+"s loaded");
       setStatus(count+ " visible " + self.gSampleDescription+"s loaded");

       if (opts.lines)
           self._setLines(opts["lines"], opts);
       self.scaleData();
    };

    this.setLabelCoords = function(labelCoords) {
        var prevLabels = self.coords.labels.length > 0;
        self.coords.labels = labelCoords;
        return prevLabels;
    };

    this.setColorArr = function(colorArr) {
    /* set the color array, one array with one index per coordinate */
       self.col.arr = colorArr;
    };

    this.setColors = function(colors) {
    /* set the colors, one for each value of a in setColorArr(a). colors is an
     * array of six-digit hex strings. Not #-prefixed! */
       self.col.pal = colors;
    };

    this.calcRadius = function() {
        /* calculate the radius from current zoom factor and set radius, alpha and zoomFact in self.port */
        // make the circles a bit smaller than expected
        var zr = self.port.zoomRange;
        var iz = self.port.initZoom;
        var initAlpha = self.port.initAlpha;
        var initSpan = iz.maxX-iz.minX;
        var currentSpan = zr.maxX-zr.minX;
        var zoomFact = initSpan/currentSpan;

        var baseRadius = self.port.initRadius;
        if (baseRadius===0)
            baseRadius = 0.7;
        var radius = Math.floor(baseRadius * Math.sqrt(zoomFact));

        // the higher the zoom factor, the higher the alpha value
        var zoomFrac = Math.min(1.0, zoomFact/100.0); // zoom as fraction, max is 1.0
        var alpha = initAlpha + 3.0*zoomFrac*(1.0 - initAlpha);
        alpha = Math.min(0.8, alpha);
        console.log("Zoom factor: ", zoomFact, ", Radius: "+radius+", alpha: "+alpha);

        self.port.zoomFact = zoomFact;
        self.port.alpha = alpha;
        self.port.radius = radius;
    }

    this.drawDots = function() {
        /* draw coordinates to canvas with current colors */
        console.time("draw");

        self.clear();

        var radius = self.port.radius;
        var alpha = self.port.alpha;
        var zoomFact = self.port.zoomFact;
        var coords = self.coords.px;
        var pal = self.col.pal;
        var colArr = self.col.arr;
        var count = 0;

        if (alpha===undefined)
             alert("internal error: alpha is not defined");
        if (coords===null)
             alert("internal error: cannot draw if coordinates are not set yet");
        if (colArr.length !== (coords.length>>1))
            alert("internal error: cbDraw.drawDots - colorArr is not 1/2 of coords array. Got "+pal.length+" color values but coordinates for "+(coords.length/2)+" cells.");

        if (radius===0) {
            count = drawPixels(self.ctx, self.canvas.width, self.canvas.height, coords,
                colArr, pal, alpha, self.selCells);
        }

        else if (radius===1 || radius===2) {
            count = drawRect(self.ctx, coords, colArr, pal, radius, alpha, self.selCells);
        }
        else {
            switch (self.mode) {
                case 0:
                    count = drawCirclesStupid(self.ctx, coords, colArr, pal, radius, alpha, self.selCells);
                    break;
                case 1:
                    count = drawCirclesDrawImage(self.ctx, coords, colArr, pal, radius, alpha, self.selCells);
                    break;
                case 2:
                    break;
            }
        }

        self.count = count;

        console.timeEnd("draw");

        if (self.doDrawLabels===true && self.coords.labels!==null) {
            self.redrawLabels();
        }

        if (self.coords.pxLines) {
            console.time("draw lines");
            drawLines(self.ctx, self.coords.pxLines, self.canvas.width, self.canvas.height, self.coords.lineAttrs);
            console.timeEnd("draw lines");
        }

        if (self.childPlot)
            self.childPlot.drawDots();
    };

    this.redrawLabels = function() {
        self.coords.pxLabels = scaleLabels(
            self.coords.labels,
            self.port.zoomRange,
            self.port.radius,
            self.canvas.width,
            self.canvas.height
        );
        self.coords.labelBbox = drawLabels(
            self.ctx,
            self.coords.pxLabels,
            self.canvas.width,
            self.canvas.height,
            self.port.zoomFact
        );
    };

    this.cellsAtPixel = function(x, y) {
        /* return the Ids of all cells at a particular pixel */
        var res = [];
        var pxCoords = self.coords.px;
        for (var i = 0; i < pxCoords.length/2; i++) {
            var cellX = pxCoords[i*2];
            var cellY = pxCoords[i*2+1];
            if (cellX===x || cellY===y)
                res.push(i);
        }
        return res;
    };

    this.cellsInRect = function(x1, y1, x2, y2) {
        /* return the Ids of all cells within certain pixel boundaries */
        var res = [];
        var pxCoords = self.coords.px;
        for (var i = 0; i < pxCoords.length/2; i++) {
            var cellX = pxCoords[i*2];
            var cellY = pxCoords[i*2+1];
            if ((cellX >= x1) && (cellX <= x2) && (cellY >= y1) && (cellY <= y2))
                res.push(i);
        }
        return res;
    };

    this.zoom100 = function() {
       /* zoom to 100% and redraw */
       copyObj(self.port.initZoom, self.port.zoomRange);
       self.scaleData();
       //self.radius = self.initRadius;
       self.drawDots();
    };

    this.zoomTo = function(x1, y1, x2, y2) {
       /* zoom to rectangle defined by two points */
       // make sure that x1<x2 and y1<y2 - can happen if mouse movement was upwards
       console.log("Zooming to pixels: ", x1, y1, x2, y2);
       var pxMinX = Math.min(x1, x2);
       var pxMaxX = Math.max(x1, x2);

       var pxMinY = Math.min(y1, y2);
       var pxMaxY = Math.max(y1, y2);

       // force the zoom rectangle to have the same aspect ratio as our canvas
       // by adapting the height. This is what Microsoft software does
       // It may be better to fix the aspect ratio of the zoom rectangle while zooming?
       // We probably do not want to distort the ratio.
       //var aspectRatio = self.width / self.height;
       //var rectWidth  = (pxMaxX-pxMinX);
       //var newHeight = rectWidth/aspectRatio;
       //pxMaxY = pxMinY + newHeight;

       var zoomRange = self.port.zoomRange;
       // window size in data coordinates
       var spanX = zoomRange.maxX - zoomRange.minX;
       var spanY = zoomRange.maxY - zoomRange.minY;

       // multiplier to convert from pixels to data coordinates
       var xMult = spanX / self.canvas.width; // multiplier dataRange/pixel
       var yMult = spanY / self.canvas.height;

       var oldMinX = zoomRange.minX;
       var oldMinY = zoomRange.minY;

       zoomRange.minX = oldMinX + (pxMinX * xMult);
       zoomRange.minY = oldMinY + (pxMinY * yMult);

       zoomRange.maxX = oldMinX + (pxMaxX * xMult);
       zoomRange.maxY = oldMinY + (pxMaxY * yMult);

       self.port.zoomRange = zoomRange;
       console.log("Marquee zoom window: "+JSON.stringify(self.port.zoomRange));

       self.scaleData();
    };

    this.zoomBy = function(zoomFact, xPx, yPx) {
    /* zoom centered around xPx,yPx by a given factor. Returns new zoom range.
     * zoomFact = 1.2 means zoom +20%
     * zoomFact = 0.8 means zoom -20%
     * */
        var zr = self.port.zoomRange;
        var iz = self.port.initZoom;

        console.log("old zoomfact "+self.port.zoomFact);

        var xRange = Math.abs(zr.maxX-zr.minX);
        var yRange = Math.abs(zr.maxY-zr.minY);

        var minWeightX = 0.5; // how zooming should be distributed between min/max
        var minWeightY = 0.5;
        if (xPx!==undefined) {
            minWeightX = (xPx/self.width);
            minWeightY = (yPx/self.canvas.height);
        }
        var scale = (1.0-zoomFact);

        var newRange = {};
        newRange.minX = zr.minX - (xRange*scale*minWeightX);
        newRange.maxX = zr.maxX + (xRange*scale*(1-minWeightX));

        // inversed, because we flip the Y axis (flipY)
        newRange.minY = zr.minY - (yRange*scale*(1-minWeightY));
        newRange.maxY = zr.maxY + (yRange*scale*(minWeightY));

        // extreme zoom factors don't make sense, at some point we reach
        // the limit of the floating point numbers
        var newZoom = ((iz.maxX-iz.minX)/(newRange.maxX-newRange.minX));
        if (newZoom < 0.01 || newZoom > 1500)
            return zr;

        console.log("x min max "+zr.minX+" "+zr.maxX);
        console.log("y min max "+zr.minY+" "+zr.maxY);

        self.port.zoomRange = newRange;

        self.scaleData();

        // a special case for connected plots that are not sharing our pixel coordinates
        if (self.childPlot && self.coords!==self.childPlot.coords) {
            self.childPlot.zoomBy(zoomFact, xPx, yPx);
        }

        return newRange;
    };

    this.movePerc = function(xDiffFrac, yDiffFrac) {
        /* move a certain percentage of current view. xDiff/yDiff are floats, e.g. 0.1 is 10% up */
        var zr = self.port.zoomRange;
        var xRange = Math.abs(zr.maxX-zr.minX);
        var yRange = Math.abs(zr.maxY-zr.minY);

        var xDiffAbs = xRange*xDiffFrac;
        var yDiffAbs = yRange*yDiffFrac;

        var newRange = {};
        newRange.minX = zr.minX + xDiffAbs;
        newRange.maxX = zr.maxX + xDiffAbs;
        newRange.minY = zr.minY + yDiffAbs;
        newRange.maxY = zr.maxY + yDiffAbs;

        copyObj(newRange, self.port.zoomRange);
        self.scaleData();
    };

    this.panStart = function() {
       /* called when starting a panning sequence, makes a snapshop of the current image */
       self.panCopy = document.createElement('canvas'); // not added to DOM, will be gc'ed
       self.panCopy.width = self.canvas.width;
       self.panCopy.height = self.canvas.height;
       var destCtx = self.panCopy.getContext("2d", { alpha: false });
       destCtx.drawImage(self.canvas, 0, 0);
    }

    this.panBy = function(xDiff, yDiff) {
        /* pan current image by x/y pixels */
        console.log('panning by '+xDiff+' '+yDiff);

       //var srcCtx = self.panCopy.getContext("2d", { alpha: false });
       clearCanvas(self.ctx, self.canvas.width, self.canvas.height);
       self.ctx.drawImage(self.panCopy, -xDiff, -yDiff);
       // keep these for panEnd
       self.panDiffX = xDiff;
       self.panDiffY = yDiff;
    }

    this.panEnd = function() {
        /* end a sequence of panBy calls, called when the mouse is released */
        self.moveBy(self.panDiffX, -self.panDiffY); // -1 because of flipY
        self.panCopy = null;
        self.panDiffX = null;
        self.panDiffY = null;
    }

    // BEGIN SELECTION METHODS (could be an object?)

    this.selectClear = function(skipNotify) {
        /* clear selection */
        self.selCells.clear();
        setStatus("");
        if (self.onSelChange!==null && skipNotify!==true)
            self.onSelChange(self.selCells);
    };

    this.selectSet = function(cellIds) {
        /* set selection to an array of integer cellIds */
        self.selCells.clear();
        //self.selCells.push(...cellIds); // "extend" = array spread syntax, https://davidwalsh.name/spread-operator
        let selCells = self.selCells;
        for (let i=0; i<cellIds.length; i++)
            selCells.add(cellIds[i]);
        self._selUpdate();
    };

    this.selectAdd = function(cellIdx) {
        /* add a single cell to the selection. If it already exists, remove it. */
        console.time("selectAdd");
        if (self.selCells.has(cellIdx))
            self.selCells.delete(cellIdx);
        else
            self.selCells.add(cellIdx);
        console.time("selectAdd");
        self._selUpdate();
    };

    this.selectAll = function(cellIdx) {
        /* add all cells to selection */
        var selCells = self.selCells;
        var pxCoords = self.coords.px;
        for (var i = 0, I = pxCoords.length / 2; i < I; i++) {
            selCells.add(i);
        }
        self.selCells = selCells;
        self._selUpdate();
    };

    this.selectVisible = function() {
        /* add all visible cells to selection */
        var selCells = self.selCells;
        var pxCoords = self.coords.px;
        for (var i = 0; i < pxCoords.length/2; i++) {
            var pxX = pxCoords[2*i];
            var pxY = pxCoords[2*i+1];
            if (isHidden(pxX, pxY))
               continue;
            selCells.add(i);
        }
        self.selCells = selCells;
        self._selUpdate();
    }

    this.selectByColor = function(colIdx) {
        /* add all cells with a given color to the selection */
        var colArr = self.col.arr;
        var selCells = self.selCells;
        var cnt = 0;
        for (var i = 0; i < colArr.length; i++) {
            if (colArr[i]===colIdx) {
                selCells.add(i);
                cnt++;
            }
        }
        self.selCells = selCells;
        console.log(cnt + " cells appended to selection, by color");
        self._selUpdate();
    };

    this.unselectByColor = function(colIdx) {
        /* remove all cells with a given color from the selection */
        var colArr = self.col.arr;
        var selCells = self.selCells;
        var cnt = 0;
        for (var i = 0; i < colArr.length; i++) {
            if (colArr[i]===colIdx) {
                selCells.delete(i);
                cnt++;
            }
        }

        self.selCells = selCells;
        console.log(cnt + " cells removed from selection, by color");
        self._selUpdate();
    };

    this.selectInRect = function(x1, y1, x2, y2) {
        /* find all cells within a rectangle and add them to the selection. */
        var minX = Math.min(x1, x2);
        var maxX = Math.max(x1, x2);

        var minY = Math.min(y1, y2);
        var maxY = Math.max(y1, y2);

        console.time("select");
        var pxCoords = self.coords.px;
        for (var i = 0; i < pxCoords.length/2; i++) {
            var pxX = pxCoords[2*i];
            var pxY = pxCoords[2*i+1];
            if (isHidden(pxX, pxY))
               continue;
            if ((minX <= pxX) && (pxX <= maxX) && (minY <= pxY) && (pxY <= maxY)) {
                self.selCells.add(i);
            }

        }
        console.timeEnd("select");
        self._selUpdate();
    };

    this.getSelection = function() {
        /* return selected cells as a list of ints */
        var cellIds = [];
        self.selCells.forEach(function(x) {cellIds.push(x)});
        return cellIds;
    };

    this.selectInvert = function() {
        /* invert selection */
        var selCells = self.selCells;
        var cellCount = self.getCount();
        for (let i = 0; i < cellCount; i++) {
            if (selCells.has(i)) {
                selCells.delete(i);
            } else {
                selCells.add(i);
            }
        }
        self.selCells = selCells;
        self._selUpdate();
    };

    this.getCount = function() {
        /* return maximum number of cells in dataset, may include hidden cells, see isHidden() */
        return self.coords.orig.length / 2;
    };

    // END SELECTION METHODS (could be an object?)

    this._selUpdate = function() {
        /* called after the selection has been updated, calls the onSelChange callback */
        setStatus(self.selCells.size + " " + self.gSampleDescription + "s selected");
        if (self.onSelChange!==null)
            self.onSelChange(self.selCells);
    }

    this.moveBy = function(xDiff, yDiff) {
        /* update the pxCoords by a certain x/y distance and redraw */

        // convert pixel range to data scale range
        var zr = self.port.zoomRange;
        var xDiffData = xDiff * ((zr.maxX - zr.minX) / self.canvas.width);
        var yDiffData = yDiff * ((zr.maxY - zr.minY) / self.canvas.height);

        // move zoom range
        zr.minX = zr.minX + xDiffData;
        zr.maxX = zr.maxX + xDiffData;
        zr.minY = zr.minY + yDiffData;
        zr.maxY = zr.maxY + yDiffData;

        self.scaleData();

        // a special case for connected plots that are not sharing our pixel coordinates
        if (self.childPlot && self.coords!==self.childPlot.coords) {
            self.childPlot.moveBy(xDiff, yDiff);
        }
    };

    this.labelAt = function(x, y) {
        /* return the index and the text of the label at position x,y or null if nothing there */
        //console.time("labelCheck");
        var clusterLabels = self.coords.labels;
        if (clusterLabels===null || clusterLabels===undefined)
            return null;
        var labelCoords = self.coords.labels;
        var boxes = self.coords.labelBbox;

        if (boxes==null) // no cluster labels
            return null;

        if (labelCoords.length!==clusterLabels.length)
            alert("internal error maxPLot.js: coordinates of labels are different from clusterLabels");

        for (var i=0; i < labelCoords.length; i++) {
            var box = boxes[i];
            if (box===null) // = outside of the screen
                continue;
            var x1 = box[0];
            var y1 = box[1];
            var x2 = box[2];
            var y2 = box[3];
            if ((x >= x1) && (x <= x2) && (y >= y1) && (y <= y2)) {
                //console.timeEnd("labelCheck");
                var labelText = clusterLabels[i][2];
                return [labelText, i];
            }
        }
        //console.timeEnd("labelCheck");
        return null;
    }

    this.cellsAt = function(x, y) {
        /* check which cell's bounding boxes contain (x, y), return a list of the cell IDs, sorted by distance */
        //console.time("cellSearch");
        var pxCoords = self.coords.px;
        if (pxCoords===null)
            return null;
        var possIds = [];
        var radius = self.port.radius;
        for (var i = 0; i < pxCoords.length/2; i++) {
           var pxX = pxCoords[2*i];
           var pxY = pxCoords[2*i+1];
            if (isHidden(pxX, pxY))
               continue;
            var x1 = pxX - radius;
            var y1 = pxY - radius;
            var x2 = pxX + radius;
            var y2 = pxY + radius;
            if ((x >= x1) && (x <= x2) && (y >= y1) && (y <= y2)) {
                var dist = Math.sqrt(Math.pow(x-pxX, 2) + Math.pow(y-pxY, 2));
                possIds.push([dist, i]);
            }
        }

        //console.timeEnd("cellSearch");
        if (possIds.length===0)
            return null;
        else {
            possIds.sort(function(a,b) { return a[0]-b[0]} ); // sort by distance

            // strip the distance information
            var ret = [];
            for (let i=0; i < possIds.length; i++) {
                ret.push(possIds[i][1]);
            }
            return ret;
        }
    };

    this.resetMarquee = function() {
       /* make the marquee disappear and reset its internal status */
       if (!self.interact)
           return;

       self.mouseDownX = null;
       self.mouseDownY = null;
       self.lastPanX = null;
       self.lastPanY = null;
       self.selectBox.style.display = "none";
       self.selectBox.style.width = 0;
       self.selectBox.style.height = 0;
    };

    this.drawMarquee = function(x1, y1, x2, y2, forceAspect) {
        /* draw the selection or zooming marquee using the DIVs created by setupMouse */
        var selectWidth = Math.abs(x1 - x2);
        var selectHeight = 0;
        if (forceAspect) {
            var aspectRatio = self.width / self.canvas.height;
            selectHeight = selectWidth/aspectRatio;
        } else
            selectHeight = Math.abs(y1 - y2);

        var minX = Math.min(x1, x2);
        var minY = Math.min(y1, y2);
        var div = self.selectBox;
        div.style.left = (minX-self.left)+"px";
        div.style.top = (minY-self.top)+"px";
        div.style.width = selectWidth+"px";
        div.style.height = selectHeight+"px";
        div.style.display = "block";
    };

    this.activatePlot = function() {
        /* draw black border around plot, remove black border from all connected plots, call onActive */
        if (self.parentPlot===null)
            return false;

        // only need to do something if we're not already the active plot
        self.canvas.style["border"] = "2px solid black";
        self.parentPlot.canvas.style["border"] = "2px solid white";

        // flip the parent/child relationship
        self.childPlot = self.parentPlot;
        self.childPlot.parentPlot = self;
        self.parentPlot = null;
        self.childPlot.childPlot = null;

        // hide/show the tool and zoom buttons
        self.childPlot.zoomDiv.style.display = "none";
        self.childPlot.toolDiv.style.display = "none";
        self.zoomDiv.style.display = "block";
        self.toolDiv.style.display = "block";

        // notify the UI
        self.onActiveChange(self);
        return true;
    }

    this.onMouseMove = function(ev) {
        /* called when the mouse is moved over the Canvas */

        // set a timer so we can get "hover" functionality without too much CPU
        if (self.timer!==null)
            clearTimeout(self.timer);
        self.timer = setTimeout(self.onNoMouseMove, 130);
        // save mouse pos for onNoMouseMove timer handler
        self.lastMouseX = ev.clientX;
        self.lastMouseY = ev.clientY;

        // label hit check requires canvas coordinates x/y
        var clientX = ev.clientX;
        var clientY = ev.clientY;
        var canvasTop = self.top;
        var canvasLeft = self.left;
        var xCanvas = clientX - canvasLeft;
        var yCanvas = clientY - canvasTop;

        // when the cursor is over a label, change it to a hand, but only when there is no marquee
        if (self.coords.labelBbox!==null && self.mouseDownX === null) {
            var labelInfo = self.labelAt(xCanvas, yCanvas);
            if (labelInfo===null) {
                self.canvas.style.cursor = self.canvasCursor;
                self.onNoLabelHover(ev);
            } else {
                self.canvas.style.cursor = 'pointer'; // not 'hand' anymore ! and not 'grab' yet!
                if (self.onLabelHover!==null)
                    self.onLabelHover(labelInfo[0], labelInfo[1], ev);
                }
        }

        if (self.mouseDownX!==null) {
            // we're panning
            if (((ev.altKey || self.dragMode==="move")) && self.panCopy!==null) {
                var xDiff = self.mouseDownX - clientX;
                var yDiff = self.mouseDownY - clientY;
                self.panBy(xDiff, yDiff);
            }
            else  {
               // zooming or selecting
               var forceAspect = false;
               var anyKey = (ev.metaKey || ev.altKey || ev.shiftKey);
               if ((self.dragMode==="zoom" && !anyKey) || ev.metaKey )
                   forceAspect = true;
               self.drawMarquee(self.mouseDownX, self.mouseDownY, clientX, clientY, forceAspect);
            }
        }
    };

    this.onNoMouseMove = function() {
        /* called after some time has elapsed and the mouse has not been moved */
        if (self.coords.px===null)
            return;
        var x = self.lastMouseX - self.left; // need canvas, not screen coordinates
        var y = self.lastMouseY - self.top;
        var cellIds = self.cellsAt(x, y);
        // only call onNoCellHover if callback exists and there is nothing selected
        if (cellIds===null && self.onNoCellHover!==null && self.selCells===null)
                self.onNoCellHover();
        else if (self.onCellHover!==null)
                self.onCellHover(cellIds);
    };


    this.onMouseDown = function(ev) {
    /* user clicks onto canvas */
       if (self.activatePlot())
           return; // ignore the first click into the plot, if it was the activating click
       console.log("background mouse down");
       var clientX = ev.clientX;
       var clientY = ev.clientY;
       if ((ev.altKey || self.dragMode==="move") && !ev.shiftKey && !ev.metaKey) {
           console.log("alt key or move mode: starting panning");
           self.panStart();
       }
       self.mouseDownX = clientX;
       self.mouseDownY = clientY;
    };

    this.onMouseUp = function(ev) {
       console.log("background mouse up");
       // these are screen coordinates
       var clientX = ev.clientX;
       var clientY = ev.clientY;
       var mouseDidNotMove = (self.mouseDownX === clientX && self.mouseDownY === clientY);

       if (self.panCopy!==null && !mouseDidNotMove) {
           console.log("ending panning operation");
           self.panEnd();
           self.mouseDownX = null;
           self.mouseDownY = null;
           self.drawDots();
           return;
       } else {
           // abort panning
           self.panCopy = null;
       }

       if (self.mouseDownX === null && self.lastPanX === null)  {
           // user started the click outside of the canvas: do nothing
           console.log("first click must have been outside of canvas");
           return;
       }

       // the subsequent operations require canvas coordinates x/y
       var canvasTop = self.top;
       var canvasLeft = self.left;
       var x1 = self.mouseDownX - canvasLeft;
       var y1 = self.mouseDownY - canvasTop;
       var x2 = clientX - canvasLeft;
       var y2 = clientY - canvasTop;

       // user did not move the mouse, so this is a click
       if (mouseDidNotMove) {
            // recognize a double click -> zoom
            if (self.lastClick!==undefined && x2===self.lastClick[0] && y2===self.lastClick[1]) {
                self.zoomBy(1.33);
                self.lastClick = [-1,-1];
            } else {
                self.lastClick = [x2, y2];
            }

            var labelInfo = self.labelAt(x2, y2);
            if (labelInfo!==null && self.doDrawLabels)
                self.onLabelClick(labelInfo[0], labelInfo[1], ev);
            else {
                var clickedCellIds = self.cellsAt(x2, y2);
                // click on a cell -> update selection and redraw
                if (clickedCellIds!==null && self.onCellClick!==null) {
                    self.selectClear(true);
                    for (var i = 0; i < clickedCellIds.length; i++) {
                        self.selCells.add(clickedCellIds[i]);
                    }
                    self.drawDots();
                    self.onCellClick(clickedCellIds, ev);

                }
                else {
                // user clicked onto background:
                // reset selection and redraw
                    console.log("not moved at all: reset "+clientX+" "+self.mouseDownX+" "+self.mouseDownY+" "+clientY);
                    self.selectClear();

                    self.drawDots();
                }

                self.lastPanX = null;
                self.lastPanY = null;
            }
            self.mouseDownX = null;
            self.mouseDownY = null;
            return;
       }
       //console.log("moved: reset "+x+" "+mouseDownX+" "+mouseDownY+" "+y);

       // it wasn't a click, so it was a drag
       var anyKey = (ev.metaKey || ev.altKey || ev.shiftKey);

       // zooming
       if ((self.dragMode==="zoom" && !anyKey) || ev.metaKey ) {
            // get current coords of the marquee in canvas pixels
            var div = self.selectBox;
            let zoomX1 = parseInt(div.style.left.replace("px",""));
            let zoomY1 = parseInt(div.style.top.replace("px",""));
            let zoomX2 = zoomX1+parseInt(div.style.width.replace("px",""));
            let zoomY2 = zoomY1+parseInt(div.style.height.replace("px",""));
            zoomY1 = self.canvas.height - zoomY1;
            zoomY2 = self.canvas.height - zoomY2;
            self.zoomTo(zoomX1, zoomY1, zoomX2, zoomY2);
            // switch back to the mode before zoom was clicked
            if (self.prevMode) {
                self.activateMode(self.prevMode);
                self.prevMode = null;
            }


       }
       // marquee select
       else if ((self.dragMode==="select" && !anyKey) || ev.shiftKey ) {
           if (! ev.shiftKey)
               self.selectClear(true);
           self.selectInRect(x1, y1, x2, y2);
       }
       else {
           console.log("Internal error: no mode?");
       }

       self.resetMarquee();

       self.drawDots();
    };

    this.onWheel = function(ev) {
        /* called when the user moves the mouse wheel */
        if (self.parentPlot!==null)
            return;
        console.log(ev);
        var normWheel = normalizeWheel(ev);
        console.log(normWheel);
        var pxX = ev.clientX - self.left;
        var pxY = ev.clientY - self.top;
        var spinFact = 0.1;
        if (ev.ctrlKey) // = OSX pinch and zoom gesture (and no other OS/mouse combination?)
            spinFact = 0.08;  // is too fast, so slow it down a little
        var zoomFact = 1-(spinFact*normWheel.spinY);
        console.log("Wheel Zoom by "+zoomFact);
        self.zoomBy(zoomFact, pxX, pxY);
        self.drawDots();
        ev.preventDefault();
        ev.stopPropagation();
    };

    this.setupMouse = function() {
       // setup the mouse callbacks
       self.canvas.addEventListener('mousedown', self.onMouseDown);
       self.canvas.addEventListener("mousemove", self.onMouseMove);
       self.canvas.addEventListener("mouseup", self.onMouseUp);
       // when the user moves the mouse, the mouse is often NOT on the canvas,
       // but on the marquee box, so connect this one, too.
       self.selectBox.addEventListener("mouseup", self.onMouseUp);

       self.canvas.addEventListener("wheel", self.onWheel);
    };

    this.setShowLabels = function(doShow) {
        self.doDrawLabels = doShow;
    };

    this.getLabels = function() {
        /* get current labels */
        var ret = [];
        var labels = self.coords.labels;
        for (var i = 0; i<labels.length; i++)
            ret.push(labels[i][2]);
        return ret;
    }

    this.setLabels = function(newLabels) {
        /* set new label text */
        if (newLabels.length!==self.coords.labels.length) {
            console.log("maxPlot:setLabels error: new labels have wrong length.");
            return;
        }

        for (var i = 0; i<newLabels.length; i++)
            self.coords.labels[i][2] = newLabels[i];

       self.coords.pxLabels = scaleLabels(self.coords.labels, self.port.zoomRange, self.port.radius,
                                           self.canvas.width, self.canvas.height);

        // a special case for connected plots that are not sharing our pixel coordinates
        if (self.childPlot && self.coords!==self.childPlot.coords) {
            self.childPlot.setLabels(newLabels);
        }
    };

    this._setLines = function(lines, attrs) {
        if (lines===undefined)
            return;
        self.coords.lines = lines;
        //self.coords.pxLines = scaleLines(self.coords.lines, self.port.zoomRange, self.canvas.width, self.canvas.height);
        if (!attrs)
            self.coords.lineAttrs = {};
        else
            self.coords.lineAttrs = attrs;
    }

    this.activateMode = function(modeName) {
    /* switch to one of the mouse drag modes: zoom, select or move */
        if (modeName==="zoom")
            self.prevMode = self.dragMode;
        else
            self.prevMode = null;

        self.dragMode=modeName;

        var cursor = null;

        if (modeName==="move")
            cursor = 'all-scroll';
        else if (modeName==="zoom")
            cursor = "zoom-in"
        else if (modeName=="select")
            cursor = 'crosshair';
        //else
            //cursor= 'default';

        self.canvas.style.cursor = cursor;
        self.canvasCursor = cursor;

        self.resetMarquee();

        if (self.interact) {
            self.icons["move"].style.backgroundColor = gButtonBackground;
            self.icons["zoom"].style.backgroundColor = gButtonBackground;
            self.icons["select"].style.backgroundColor = gButtonBackground;
            self.icons[modeName].style.backgroundColor = gButtonBackgroundClicked;
        }
        if (self.childPlot)
            self.childPlot.activateMode(modeName);
    }

    this.randomDots = function(n, radius, mode) {
        /* draw x random dots with x random colors*/
	function randomArray(ArrType, length, max) {
            /* make Array and fill it with random numbers up to max */
            var arr = new ArrType(length);
            for (var i = 0; i<length; i++) {
                arr[i] = Math.round(Math.random() * max);
            }
            return arr;
	}

        if (mode!==undefined)
            self.mode = mode;
        self.port.radius = radius;
	self.setCoords(randomArray(Uint16Array, 2*n, 65535));
	self.setColors(["FF0000", "00FF00", "0000FF", "CC00CC", "008800"]);
	self.setColorArr(randomArray(Uint8Array, n, 4));

        console.time("draw");
        self.drawDots();
        console.timeEnd("draw");
        return self;
    };

    this.split = function() {
        /* reduce width of renderer, create new renderer and place both side-by-side.
         * They initially share the .coords but setCoords() can break that relationship.
         * */
        var canvHeight  = self.canvas.height;
        var canvLeft    = self.left;
        var newWidth = self.width/2;
        var newTop   = self.top;
        var newLeft  = self.left+newWidth;
        var newHeight = canvHeight+gStatusHeight;

        var newDiv = document.createElement('div');
        newDiv.id = "mpPlot2";
        document.body.appendChild(newDiv);

        var opts = cloneObj(self.globalOpts);
        opts.showClose = true;

        var plot2 = new MaxPlot(newDiv, newTop, newLeft, newWidth, newHeight, {"showClose" : true});
        //plot2.canvas.style.borderLeft = "1px solid grey";

        plot2.statusLine.style.display = "none";

        plot2.port = self.port;
        plot2.selCells = self.selCells;

        //plot2.coords = Object.assign({}, self.coords); // = shallow copy
        plot2.coords = self.coords;

        plot2.col = {};
        plot2.col.pal = self.col.pal;
        plot2.col.arr = self.col.arr;

        self.setSize(newWidth, newHeight, false); // will call scaleData(), but not redraw.

        plot2.onLabelClick = self.onLabelClick;
        plot2.onCellClick = self.onCellClick;
        plot2.onCellHover = self.onCellHover;
        plot2.onNoCellHover = self.onNoCellHover;
        plot2.onSelChange = self.onSelChange;
        plot2.onLabelHover = self.onLabelHover;
        plot2.onNoLabelHover = self.onNoLabelHover;
        plot2.onActiveChange = self.onActiveChange;

        //var closeButton = gebi('mpCloseButton');
        //closeButton.addEventListener('click', self.unsplit);

        plot2.drawDots();

        self.childPlot = plot2;
        plot2.parentPlot = self;

        // add a thick border and hide the menus in the child
        self.canvas.style["border"] = "2px solid black";
        self.childPlot.zoomDiv.style.display = "none";
        self.childPlot.toolDiv.style.display = "none";

        return plot2;
    };

    this.unsplit = function() {
        /* remove the connected non-active renderer */
        //var canvWidth = window.innerWidth - canvLeft - legendBarWidth;
        var otherRend = self.childPlot;
        self.childPlot = undefined;
        if (!otherRend) {
            otherRend = self.parentPlot;
            self.parentPlot = undefined;
        }
        self.setSize(self.width*2, self.height, false);


        otherRend.div.remove();
        self.canvas.style["border"] = "none";
        return;
    }

    this.getWidth = function() {
        /* return total size of renderer, including any split child renderers */
        if (self.childPlot)
            return self.width + self.childPlot.width;
        else
            return self.width;
    }

    // object constructor code
    self.initCanvas(div, top, left, width, height);
    self.initPlot(args);

}
