'use strict';
// maxPlot: a fast scatter plot class
/*jshint globalstrict: true*/

// TODO:
// remove onNoHover
// fix mouseout into body -> marquee stays
// fix the aspect ratio of the zoom marquee

function getAttr(obj, attrName, def) {
    var val = obj[attrName];
    if (val===undefined)
        return def;
    else
        return val;
}

function cloneObj(d) {
/* returns a copy of an object, wasteful */
    // see http://stackoverflow.com/questions/122102/what-is-the-most-efficient-way-to-deep-clone-an-object-in-javascript
    return JSON.parse(JSON.stringify(d));
}


function MaxPlot(div, top, left, width, height, args) {
    // a class that draws circles onto a canvas, like a scatter plot
    // div is a div DOM element under which the canvas will be created
    // top, left: position in pixels, integers
    // width and height: integers, in pixels, includes the status line
    
    var self = this; // 'this' has two conflicting meanings in javascript. 
    // I use 'self' to refer to object variables, so I can use 'this' to refer to the caller context
    
    const gTextSize = 16; // size of cluster labels
    const gStatusHeight = 12; // height of status bar
    const gZoomButtonSize = 30; // size of zoom buttons
    const gZoomFromRight = 60;  // position of zoom buttons from right
    const gZoomFromBottom = 120;  // position of zoom buttons from bottom
    const gButtonBackground = "rgb(230, 230, 230, 0.6)" // grey level of buttons
    const gButtonBackgroundClicked = "rgb(180, 180, 180, 0.6)"; // grey of buttons when clicked

    // the rest of the initialization is done at the end of this file,
    // because the init involves many functions that are not defined yet here

    this.initCanvas = function (div, top, left, width, height) {
        /* initialize a new Canvas */

        self.div = div;
        self.gSampleDescription = "cell"; 
        self.ctx = null; // the canvas context
        self.canvas = addCanvasToDiv(div, top, left, width, height-gStatusHeight );

        addZoomButtons(top+height-gZoomFromBottom, left+width-gZoomFromRight, self);
        addModeButtons(top+10, left+10, self);
        addStatusLine(top+height-gStatusHeight, left, width, gStatusHeight);
        addProgressBars(top+Math.round(height*0.3), left+30);

        /* add the div used for the mouse selection/zoom rectangle to the DOM */
        var selectDiv = document.createElement('div');
        selectDiv.id = "mpSelectBox";
        selectDiv.style.border = "1px dotted black";
        selectDiv.style.position = "absolute";
        selectDiv.style.display  = "none";
        selectDiv.pointerEvents = "none";
        self.div.appendChild(selectDiv);
        self.selectBox = selectDiv; // we need this later

        self.setupMouse();
        
        // callbacks when user clicks or hovers over label or cell 
        self.onLabelClick = null; // called on label click, args: text of label and event
        self.onCellClick = null; // called on cell click, args: array of cellIds and event
        self.onCellHover = null; // called on cell hover, arg: array of cellIds
        self.onNoCellHover = null; // called on hover over empty background
        self.onSelChange = null; // called when the selection has been changed, arg: array of cell Ids
        // self.onZoom100Click: called when user clicks the zoom100 button. Implemented below.

        // timer that is reset on every mouse move
        self.timer = null;

        // all other object variables are added by the "initPlot(args)" function below

        // when the user starts to select but lifts the mouse button outside
        // the canvas, the current selection must be reset
        //self.canvas.addEventListener("mouseleave", function() {self.resetMarquee();});

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
    
    this.initPlot = function(args) {
        /* create a new scatter plot on the canvas */
        if (args===undefined)
            args = {};

        self.mode = 1;   // drawing mode

        self.zoomRange = null; // object with keys minX, , maxX, minY, maxY
        self.coords     = null;
        self.clusterLabels = null; // cluster labels, array of [x,y,text]

        self.pxCoords   = null;   // coordinates of cells as screen pixels or 0,0 if not shown
        self.pxLabels   = null;   // cluster labels in pixels, array of [x,y,text] 
        self.pxLabelBbox = null;   // cluster label bounding boxes, array of [x1,x2,x2,y2]

        self.doDrawLabels = true;  // should cluster labels be drawn?

        self.colors     = null;   // list of six-digit hex codes
        self.colorArr   = null;   // length is pxCoords/2, one byte per cell = index into self.colors
        self.radius     = getAttr(args, "radius", null);    // current radius of the circles, 0=one pixel dots
        self.alpha      = getAttr(args, "alpha", 0.3);
        self.selCells   = null;  // IDs of cells that are "selected" and as such highlighted in some way

        // we keep a copy of the 'initial' arguments at 100% zoom
        self.initZoom   = cloneObj(self.zoomRange);
        self.initRadius = self.radius;                      // circle radius at full zoom

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
        clearCanvas(self.ctx, self.width, self.height);
    };

    this.setPos = function(left, top) {
       /* position the canvas on the page */
       self.canvas.style.left = left+"px";
       self.canvas.style.top = top+"px";
    };

    this.setTitle = function (text) {
        self.title = text;
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

    function createButton(width, height, id, title, text, imgFname, paddingTop, addSep) {
        /* make a light-grey div that behaves like a button, with text and/or an image on it 
         * Images are hard to vertically center, so padding top can be specified.
         * */
        var div = document.createElement('div');
        div.id = id;
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
        if (imgFname!==null && imgFname!=undefined) {
            var img = document.createElement('img');
            img.src = imgFname;
            if (paddingTop!==null && paddingTop!=undefined)
                img.style.paddingTop = paddingTop+"px";
            div.appendChild(img);
        }
        if (addSep===true)
            div.style["border-bottom"] = "1px solid #D7D7D7";

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

        var plusDiv = createButton(width, height, "mpCtrlZoomPlus", "Zoom in", "+", null, null, true);
        //plusDiv.style["border-bottom"] = "1px solid #D7D7D7";

        var fullDiv = createButton(width, height, "mpCtrlZoom100", "Zoom in", "100%", null, null, true);
        //full.style["border-bottom"] = "1px solid #D7D7D7";

        var minusDiv = createButton(width, height, "mpCtrlZoomMinus", "Zoom out", "-");

        var ctrlDiv = makeCtrlContainer(top, left);
        ctrlDiv.appendChild(plusDiv);
        ctrlDiv.appendChild(fullDiv);
        ctrlDiv.appendChild(minusDiv);

        self.div.appendChild(ctrlDiv);

        minusDiv.addEventListener('click', function() { self.zoomBy(0.75); self.drawDots(); });
        fullDiv.addEventListener('click', function() { self.zoom100(); self.drawDots()});
        plusDiv.addEventListener('click', function() { self.zoomBy(1.333); self.drawDots(); });
    }
    
    function appendButton(parentDiv, id, title, imgName) {
        /* add a div styled like a button under div */
        var div = document.createElement('div');
        div.title = title;
        div.id = id;

    }

    function addModeButtons(top, left, self) {
        /* add the zoom/move/select control buttons to the DOM */
        //var htmls = [];
        //htmls.push('<div id="mpIcons" style="display:inline-block">');
        //htmls.push('<div style="vertical-align:top">');
        //htmls.push('<div title="Select mode.<br>Keyboard: shift or s" id="mpIconModeSelect"><img src="img/select.png"></button>');
        //htmls.push('<div title="Zoom-to-rectangle mode. Keyboard: Windows/Command or z" id="mpIconModeZoom" style="display: block; margin-right:0"><img src="img/zoom.png"></button>');
        //htmls.push('<div data-placement="bottom" title="Move mode. Keyboard: Alt or m" id="mpIconModeMove" data-toggle="tooltip" class="ui-button tpIconButton" class="mpModeButton" style="margin-right:0"><img src="img/move.png"></button>');
        //htmls.push('</div>');
        //htmls.push('<button title="Zoom to 100%, showing all data, keyboard: space" data-placement="bottom" data-toggle="tooltip" id="mpZoom100Button" class="ui-button tpIconButton" style="font-size:10px; font-weight: bold; margin-top: 4px; margin-right:0; display:block; padding:0">100%</button>');
        //<img style="width:22px; height:22px" src="img/center.png">

        var ctrlDiv = makeCtrlContainer(top, left);

        var bSize = gZoomButtonSize;

        var selectButton = createButton(bSize, bSize, "mpIconModeSelect", "Select mode. Keyboard: shift or s", null, "img/select.png", 4, true);
        selectButton.addEventListener ('click',  function() { self.activateMode("select")}, false);

        var zoomButton = createButton(bSize, bSize, "mpIconModeZoom", "Zoom-to-rectangle mode. Keyboard: Windows/Command or z", null, "img/zoom.png", 4, true);
        zoomButton.addEventListener ('click', function() { self.activateMode("zoom")}, false);  

        var moveButton = createButton(bSize, bSize, "mpIconModeMove", "Move mode. Keyboard: Alt or m", null, "img/move.png", 4);
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

        activateTooltip('.mpIconButton');

        //gebi('mpZoom100Button').addEventListener ('click', function(ev) { return self.onZoom100Click(ev) } );
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
        div.style.left = left+"px";
        div.style.width = width+"px";
        div.style.height = height+"px";
        div.style["border-left"]="1px solid #DDD";
        div.style["border-right"]="1px solid #DDD";
        div.style["border-top"]="1px solid #DDD";
        //div.style["border-bottom"]="1px solid #DDD";
        div.style["font-size"]=(gStatusHeight-1)+"px";
        //div.style["vertical-align"]="middle";
        //div.style["line-height"]=(height-2)+"px";
        //div.style["box-shadow"]="0px 2px 4px rgba(0,0,0,0.3)";
        //div.style["border-radius"]="2px";
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
        canv.style.top = top+"px";
        canv.style.left = left+"px";
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
            if ((x < minX) || (x > maxX) || (y < minY) || (y > maxY)) {
                pxLabels.push(null);
                continue;
            }
            var xPx = Math.round((x-minX)*xMult)+borderSize;
            var yPx = winHeight - Math.round((y-minY)*yMult)+borderSize;
            pxLabels.push([xPx, yPx, text]);
        }
        return pxLabels;
    }

    function scaleCoords(coords, borderSize, zoomRange, winWidth, winHeight, annots) {
    /* scale list of [x (float),y (float)] to integer pixels on screen and
     * annots is an array with on-screen annotations in the format (x, y,
     * otherInfo) that is also scaled.  return [array of (x (int), y (int)),
     * scaled annots array]. Take into account the current zoom range.      *
     * Canvas origin is top-left, but usually plotting origin is bottom-left,
     * so also flip the Y axis.
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
            // simply ignore anything outside of current zoom range.
            if ((x < minX) || (x > maxX) || (y < minY) || (y > maxY))
                continue;
            var xPx = Math.round((x-minX)*xMult)+borderSize;
            // flipY: y-axis is flipped, so we do winHeight - pixel value
            var yPx = winHeight - Math.round((y-minY)*yMult)+borderSize;
            pixelCoords[2*i] = xPx;
            pixelCoords[2*i+1] = yPx;
        }

        console.timeEnd("scale");
        return pixelCoords;
    }

    function drawRect(ctx, pxCoords, coordColors, colors, radius, alpha, selCells) {
       console.log("Drawing "+coordColors.length+" rectangles, with fillRect");
       ctx.save();
       ctx.globalAlpha = alpha;
       var dblSize = 2*radius;
       var count = 0;
       for (var i = 0; i < pxCoords.length/2; i++) {
           var pxX = pxCoords[2*i];
           var pxY = pxCoords[2*i+1];
           if (pxX===0 && pxY===0)
               continue;
           var col = colors[coordColors[i]];
           ctx.fillStyle="#"+col;
           ctx.fillRect(pxX-radius, pxY-radius, dblSize, dblSize);
           count++;
       }
       
       // draw the selection as black rectangles
       ctx.globalAlpha = 0.7;
       ctx.fillStyle="black";
       if (selCells!==null)
           for (i = 0; i < selCells.length; i++) {
               var cellId = selCells[i];
               var pxX = pxCoords[2*cellId];
               var pxY = pxCoords[2*cellId+1];
               ctx.fillRect(pxX-radius, pxY-radius, dblSize, dblSize);
               count += 1;
           }
       console.log(count+" rectangles drawn (including selection)");
       ctx.restore();
    }

    function drawCirclesStupid(ctx, pxCoords, coordColors, colors, radius, alpha, selCells) {
    /* draw little circles onto canvas. pxCoords are the centers.  */
       console.log("Drawing "+coordColors.length+" circles with stupid renderer");
       ctx.globalAlpha = alpha;
       var dblSize = 2*radius;
       for (var i = 0; i < pxCoords.length/2; i++) {
           var pxX = pxCoords[2*i];
           var pxY = pxCoords[2*i+1];
           var col = colors[coordColors[i]];
           ctx.fillStyle="#"+col;
           //ctx.fillRect(pxX-size, pxY-size, dblSize, dblSize);
           ctx.beginPath();
           ctx.arc(pxX, pxY, radius, 0, 2 * Math.PI);
           ctx.closePath();
           ctx.fill();
       }
    }

    function intersectRect(r1left, r1right, r1top, r1bottom, r2left, r2right, r2top, r2bottom) {
      /* return true if two rectangles overlap, 
       https://stackoverflow.com/questions/2752349/fast-rectangle-to-rectangle-intersection
	*/
      return !(r2left > r1right || r2right < r1left || r2top > r1bottom || r2bottom < r1top);
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
            if (coord===null) { // outside of view range
                bboxArr.push( null );
                continue
            }

            var x = coord[0];
            var y = coord[1];
            var text = coord[2];

            var textWidth = Math.round(ctx.measureText(text).width);
            // move x to the left, so text is centered on x
            x = x - Math.round(textWidth*0.5);

            // don't draw labels where the midpoint is off-screen
            if (x<0 || y<0 || x>winWidth || y>winWidth)
                continue;

            var textX1 = x;
            var textY1 = y;
            var textX2 = Math.round(x+textWidth);
            var textY2 = y+gTextSize;

            if (zoomFact===1.0) {
                // at 100% zoom, make some minimal effort to keep labels on screen
                if (x < 0)
                    x = 0;
                if ((x + textWidth) > winWidth)
                    x = winWidth - textWidth;
                if (y+gTextSize > winHeight)
                    y = winHeight-gTextSize;

                // also only at 100% zoom, make a minimal effort to avoid label overlaps
                // a perfect solution would take much more time
                for (var j=0; j < bboxArr.length; j++) {
                    var bbox = bboxArr[j];
                    if (bbox===null) // = outside of screen
                        continue;
                    var bx1 = bbox[0];
                    var by1 = bbox[1];
                    var bx2 = bbox[2];
                    var by2 = bbox[3];
                    if (intersectRect(textX1, textX2, textY1, textY2, bx1, bx2, by1, by2)) {
                            // push the overlapping label away a little
                            var diff = Math.round(0.75*gTextSize);
                            if (textY1 < by1)
                                y -= diff;
                            else
                                y += diff;
                        }
                }
            }

            ctx.strokeText(text,x,y); 
            ctx.fillText(text,x,y); 

            bboxArr.push( [textX1-addMargin, textY1-addMargin, textX2+addMargin, textY2+addMargin] );
        }
        ctx.restore();
        console.timeEnd("labels");
        return bboxArr;
    }

    function drawLabels_dom(ctx, labelCoords, isFull) {
        /* given an array of [x, y, text], draw the text. returns bounding boxes as array of [x1, y1, x2, y2]  */
        for (var i=0; i < labelCoords.length; i++) {
            var coord = labelCoords[i];
            var x = coord[0];
            var y = coord[1];
            var text = coord[2];

            var div = document.createElement('div');
            div.id = id;
            //div.style.border = "1px solid #DDDDDD";
            div.style.backgroundColor = "rgb(230, 230, 230, 0.6)";
            div.style.width = width+"px";
            div.style.height = height+"px";
            div.style["text-align"]="center";
            div.style["vertical-align"]="middle";
            div.style["line-height"]=height+"px";
        }
        ctx.restore();
        return bboxArr;
    }

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
           if (radius>6) {
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
       for (var i = 0; i < pxCoords.length/2; i++) {
           var pxX = pxCoords[2*i];
           var pxY = pxCoords[2*i+1];
           if (pxX===0 && pxY===0)
               continue;
           var col = coordColors[i];
           count++;
           // drawImage(img,sx,sy,swidth,sheight,x,y,width,height);
           ctx.drawImage(off, col * tileWidth, 0, tileWidth, tileHeight, pxX - radius - 1, pxY - radius - 1, tileWidth, tileHeight);
       }
       
       // overdraw the selection as solid black circle outlines
       ctx.globalAlpha = 0.7;
       if (selCells!==null) {
           for (i = 0; i < selCells.length; i++) {
               var cellId = selCells[i];
               var pxX = pxCoords[2*cellId];
               var pxY = pxCoords[2*cellId+1];
               if (pxX===0 && pxY===0)
                   continue;
               // make sure that old leftover overlapping black circles don't shine through
               var col = coordColors[cellId];
               ctx.drawImage(off, col * tileWidth, 0, tileWidth, tileHeight, pxX - radius -1, pxY - radius-1, tileWidth, tileHeight);

               ctx.drawImage(off, selImgId * tileWidth, 0, tileWidth, tileHeight, pxX - radius -1, pxY - radius-1, tileWidth, tileHeight);
           }
       }

       console.log(count +" circles drawn");
       ctx.restore();
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
           if (pxX===0 && pxY===0)
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

       // alpha-blend pixels into array
       for (var i = 0; i < pxCoords.length/2; i++) {
           var pxX = pxCoords[2*i];
           var pxY = pxCoords[2*i+1];
           if (pxX===0 && pxY===0)
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
       
       // overdraw the selection as black pixels
       if (selCells!==null)
           for (i = 0; i < selCells.length; i++) {
               var cellId = selCells[i];
               var pxX = pxCoords[2*cellId];
               var pxY = pxCoords[2*cellId+1];
               if (pxX===0 && pxY===0)
                   continue;
               var p = 4 * (pxY*width+pxX); // pointer to red value of pixel at x,y
               cData[p] = 0; 
               cData[p+1] = 0;
               cData[p+2] = 0;
           }

       self.ctx.putImageData(canvasData, 0, 0);
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
        console.time("clear");
        ctx.save();
        ctx.globalAlpha = 1.0;
        ctx.fillStyle = "rgb(255,255,255)";
        ctx.fillRect(0, 0, width, height);
        ctx.restore();
        console.timeEnd("clear");
    }

    // -- object methods (=access the self object)
 
    this.onZoom100Click = function(ev) {
        self.zoom100();
        self.drawDots();
    };

    this.scaleData = function() {
       /* scale coords and labels to current zoom range, write results to pxCoords and pxLabels */
       var borderMargin = self.radius;
       self.pxCoords = scaleCoords(self.coords, borderMargin, self.zoomRange, self.canvas.width, self.canvas.height);

       self.pxLabels = null;
       if (self.clusterLabels!==undefined && self.clusterLabels!==null)
           self.pxLabels = scaleLabels(self.clusterLabels, self.zoomRange, borderMargin, self.canvas.width, self.canvas.height);
    }

    this.setSize = function(width, height) {
       /* resize canvas on the page re-scale the data and re-draw */
       
       // css and canvas sizes: these must be identical, otherwise canvas gets super slow
       var canvWidth = width;
       var canvHeight = height-gStatusHeight;
       self.canvas.style.width = canvWidth+"px";
       self.canvas.style.height = canvHeight+"px";
       self.canvas.width = canvWidth;
       self.canvas.height = canvHeight;

       self.width = width;
       self.height = height;

       var zoomDiv = gebi('mpCtrls');
       zoomDiv.style.top = (self.top+height-gZoomFromBottom)+"px";
       zoomDiv.style.left = (self.left+width-gZoomFromRight)+"px";

       var statusDiv = self.statusLine;
       statusDiv.style.top = (self.top+height-gStatusHeight)+"px";
       statusDiv.style.width = width+"px";

       self.scaleData();
       //clearCanvas(self.ctx, width, height);
       self.drawDots();
    };

    this.setCoords = function(coords, clusterLabels, minX, maxX, minY, maxY) {
       /* specify new coordinates of circles to draw, an array of (x,y) coordinates */
       /* Scale data to current screen dimensions */
       /* clusterLabels is optional: array of [x, y, labelString]*/
       /* minX, maxX, etc are optional */
       // "==null" checks for both undefined and null
       if (coords.length === 0)
           alert("cbDraw-setCoords called with no coordinates");

       if (minX===undefined || maxX===undefined || minY===undefined || maxY===undefined)
           self.initZoom = findRange(coords);
       else
           self.initZoom = {minX:minX, maxX:maxX, minY:minY, maxY:maxY}
       self.zoomRange = cloneObj(self.initZoom);

       self.coords = coords;
       self.clusterLabels = clusterLabels;
       setStatus((coords.length/2)+" "+self.gSampleDescription+"s loaded");

       self.scaleData();
       if (self.radius===null || self.radius==undefined) {
           self.radius = guessRadius(coords.length);
           self.initRadius = self.radius;
           }
    };

    this.setColorArr = function(colorArr) {
    /* set the color array, one array with one index per coordinate */
       self.colorArr = colorArr;
    };

    this.setColors = function(colors) {
    /* set the colors, one for each value of a in setColorArr(a). colors is an
     * array of six-digit hex strings. Not #-prefixed! */
       self.colors = colors;
    };

    this.drawTitle = function() {
        var ctx = self.ctx;
        ctx.save();
        ctx.font = "bold "+gTextSize+"px Sans-serif"
        //ctx.globalAlpha = 1.0;

        //ctx.strokeStyle = '#EEEEEE'; 
        //ctx.lineWidth = 5; 
        ctx.fillStyle = "rgba(220, 220, 220)";
        ctx.textBaseline = "top";
        //ctx.textAlign = "left";
        ctx.fillText(self.title,5,self.height - gTextSize - 3); 
        ctx.restore();
    };

    this.drawDots = function() {
        /* draw coordinates to canvas with current colors */
        console.time("draw");

        self.clear();

        if (self.title!==undefined)
            self.drawTitle();

        if (self.alpha===undefined)
             alert("internal error: alpha is not defined");
        if (self.pxCoords===null)
             alert("internal error: cannot draw if coordinates are not set yet");
        if (self.colorArr.length !== (self.pxCoords.length>>1))
            alert("internal error: cbDraw.drawDots - colorArr is not 1/2 of coords array. Got "+self.colors.length+" color values but coordinates for "+(self.pxCoords.length/2)+" cells.");

        self.zoomFact = ((self.initZoom.maxX-self.initZoom.minX)/(self.zoomRange.maxX-self.zoomRange.minX));

        console.log("Zoom factor: "+self.zoomFact);
        // make the circles a bit smaller than expected
        var baseRadius = self.initRadius;
        if (baseRadius===0)
            baseRadius = 0.7;
        self.radius = Math.floor(baseRadius * Math.sqrt(self.zoomFact));

        // the higher the zoom factor, the higher the alpha value
        var zoomFrac = Math.min(1.0, self.zoomFact/100.0); // zoom as fraction, max is 1.0
        var alpha = self.alpha + 3.0*zoomFrac*(1.0 - self.alpha);
        alpha = Math.min(0.8, alpha);
        console.log("Radius: "+self.radius+", alpha: "+alpha);

        if (self.radius===0) {
            drawPixels(self.ctx, self.canvas.width, self.canvas.height, self.pxCoords, 
                self.colorArr, self.colors, alpha, self.selCells);
        }

        else if (self.radius===1 || self.radius===2) {
            drawRect(self.ctx, self.pxCoords, self.colorArr, self.colors, self.radius, alpha, self.selCells);
        }
        else {

            switch (self.mode) {
                case 0:
                    drawCirclesStupid(self.ctx, self.pxCoords, self.colorArr, self.colors, self.radius, alpha, self.selCells);
                    break;
                case 1:
                    drawCirclesDrawImage(self.ctx, self.pxCoords, self.colorArr, self.colors, self.radius, alpha, self.selCells);
                    break;
                case 2:
                    break;
            }
        }

        console.timeEnd("draw");

        if (self.doDrawLabels===true && self.pxLabels!==null) {
            self.pxLabelBbox = drawLabels(self.ctx, self.pxLabels, self.canvas.width, self.canvas.height, self.zoomFact);
        }
    };

    this.cellsAtPixel = function(x, y) {
        /* return the Ids of all cells at a particular pixel */
        var res = [];
        var pxCoords = self.pxCoords;
        for (var i = 0; i < self.pxCoords.length/2; i++) {
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
        var pxCoords = self.pxCoords;
        for (var i = 0; i < self.pxCoords.length/2; i++) {
            var cellX = pxCoords[i*2];
            var cellY = pxCoords[i*2+1];
            if ((cellX >= x1) && (cellX <= x2) && (cellY >= y1) && (cellY <= y2))
                res.push(i);
        }
        return res;
    };

    //this.setSelection = function (idList) {
        /* select some dots. idList is a list of integers. Specify null to clear. */
        //self.selCells = new FastBitSet(idList);
    //};

    this.zoom100 = function() {
       /* zoom to 100% and redraw */
       self.zoomRange = cloneObj(self.initZoom);
       self.scaleData();
       self.radius = self.initRadius;
       self.drawDots();
    };

    this.zoomTo = function(x1, y1, x2, y2) {
       /* zoom to rectangle defined by two points */
       // make sure that x1<x2 and y1<y2 - can happen if mouse movement was upwards
       var pxMinX = Math.min(x1, x2);
       var pxMaxX = Math.max(x1, x2);

       var pxMinY = Math.min(y1, y2);
       var pxMaxY = Math.max(y1, y2);

       // force the zoom rectangle to have the same aspect ratio as our canvas
       // by adapting the height. This is what Microsoft software does
       // It may be better to fix the aspect ratio of the zoom rectangle while zooming?
       // We probably do not want to distort the ratio.
       var aspectRatio = self.width / self.height;
       var rectWidth  = (pxMaxX-pxMinX);
       var newHeight = rectWidth/aspectRatio;
       pxMaxY = pxMinY + newHeight;

       var zoomRange = self.zoomRange;
       // window size in data coordinates
       var spanX = zoomRange.maxX - zoomRange.minX;
       var spanY = zoomRange.maxY - zoomRange.minY;

       // multiplier to convert from pixels to data coordinates
       var xMult = spanX / self.width; // multiplier dataRange/pixel
       var yMult = spanY / self.height;

       var oldMinX = zoomRange.minX;
       var oldMinY = zoomRange.minY;

       zoomRange.minX = oldMinX + (pxMinX * xMult);
       zoomRange.minY = oldMinY + (pxMinY * yMult);

       zoomRange.maxX = oldMinX + (pxMaxX * xMult);
       zoomRange.maxY = oldMinY + (pxMaxY * yMult);

       self.zoomRange = zoomRange;

       self.scaleData();
    };

    this.zoomBy = function(zoomFact, xPx, yPx) {
    /* zoom centered around xPx,yPx by a given factor. Returns new zoom range. 
     * zoomFact = 1.2 means zoom +20%
     * zoomFact = 0.8 means zoom -20%
     * */
        var zr = self.zoomRange;
        var iz = self.initZoom;

        console.log("zoomfact "+self.zoomFact);

        var xRange = Math.abs(zr.maxX-zr.minX);
        var yRange = Math.abs(zr.maxY-zr.minY);

        var minWeightX = 0.5; // how zooming should be distributed between min/max
        var minWeightY = 0.5;
        if (xPx!==undefined) {
            minWeightX = (xPx/self.width);
            minWeightY = (yPx/self.height);
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
       
        self.zoomRange = newRange;
        self.scaleData();
        return newRange;
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

    this.selectClear = function() {
        /* clear selection */
        self.selCells = null;
        setStatus("");
        if (self.onSelChange!==null)
            self.onSelChange([]);
    };

    this.selectAdd = function(cellIdx) {
        /* add a single cell to the selection. If it already exists, remove it. */
        if (self.selCells===null) {
            self.selCells = [cellIdx];
            self._selUpdate();
            return;
        }

        console.time("selectAdd");
        var foundIdx = self.selCells.indexOf(cellIdx);
        if (foundIdx===-1)
            self.selCells.push(cellIdx);
        else
            self.selCells.splice(foundIdx, 1);
        console.time("selectAdd");
        self._selUpdate();
    };

    this.selectVisible = function() {
        /* add all visible cells to selection */
        if (self.selCells === null)
            self.selCells = [];

        var selCells = self.selCells;
        var pxCoords = self.pxCoords;
        for (var i = 0; i < pxCoords.length/2; i++) {
            var pxX = pxCoords[2*i];
            var pxY = pxCoords[2*i+1];
            if (pxX===0 && pxY===0)
                continue;
            selCells.push(i);
        }
        self.selCells = selCells;
        self._selUpdate();
    }

    this.selectByColor = function(colIdx) {
        /* add all cells with a given color to the selection */
        var colArr = self.colorArr;
        var selCells = self.selCells;
        if (selCells===null)
            selCells = [];
        for (var i = 0; i < colArr.length; i++) {
            if (colArr[i]===colIdx)
                selCells.push(i);
        }
        self.selCells = selCells;
        console.log(self.selCells.length+" cells selected, by color");
        self._selUpdate();
    };

    this.selectInRect = function(x1, y1, x2, y2) {
        /* find all cells within a rectangle and add them to the selection. */
        var minX = Math.min(x1, x2);
        var maxX = Math.max(x1, x2);

        var minY = Math.min(y1, y2);
        var maxY = Math.max(y1, y2);

        console.time("select");
        if (self.selCells===null)
            self.selCells = [];

        var pxCoords = self.pxCoords;
        for (var i = 0; i < pxCoords.length/2; i++) {
            var pxX = pxCoords[2*i];
            var pxY = pxCoords[2*i+1];
            if (pxX===0 && pxY===0)
                continue;
            if ((minX <= pxX) && (pxX <= maxX) && (minY <= pxY) && (pxY <= maxY)) {
                self.selCells.push(i);
            }

        }
        console.timeEnd("select");
        self._selUpdate();
    };

    this._selUpdate = function() {
        /* called after the selection has been updated, calls the onSelChange callback */
        setStatus(self.selCells.length+" "+self.gSampleDescription+"s selected");
        if (self.onSelChange!==null)
            self.onSelChange(self.selCells);
    }

    this.moveBy = function(xDiff, yDiff) {
        /* update the pxCoords by a certain x/y distance and redraw */

        // convert pixel range to data scale range
        var zr = self.zoomRange;
        var xDiffData = xDiff * ((zr.maxX - zr.minX) / self.canvas.width);
        var yDiffData = yDiff * ((zr.maxY - zr.minY) / self.canvas.height);
        
        // move zoom range 
        zr.minX = zr.minX + xDiffData;
        zr.maxX = zr.maxX + xDiffData;
        zr.minY = zr.minY + yDiffData;
        zr.maxY = zr.maxY + yDiffData;

        self.scaleData();
    };

    this.labelAt = function(x, y) {
        /* return the text of the label at position x,y or null if nothing there */
        //console.time("labelCheck");
        var clusterLabels = self.clusterLabels;
        if (clusterLabels===null || clusterLabels===undefined)
            return null;
        var labelCoords = self.pxLabels;
        var boxes = self.pxLabelBbox;

        for (var i=0; i < labelCoords.length; i++) {
            var box = boxes[i];
            if (box===null) // = outside of the screen
                continue;
            var labelText = clusterLabels[i][2];
            var x1 = box[0];
            var y1 = box[1];
            var x2 = box[2];
            var y2 = box[3];
            if ((x >= x1) && (x <= x2) && (y >= y1) && (y <= y2)) {
                //console.timeEnd("labelCheck");
                return labelText;
            }
        }
        //console.timeEnd("labelCheck");
        return null;
    }

    this.cellsAt = function(x, y) {
        /* check which cell's bounding boxes contain (x, y), return a list of the cell IDs, sorted by distance */
        console.time("cellSearch");
        var pxCoords = self.pxCoords;
        var possIds = [];
        for (var i = 0; i < pxCoords.length/2; i++) {
           var pxX = pxCoords[2*i];
           var pxY = pxCoords[2*i+1];
           if (pxX===0 && pxY===0)
               continue;
            var x1 = pxX - self.radius;
            var y1 = pxY - self.radius;
            var x2 = pxX + self.radius;
            var y2 = pxY + self.radius;
            if ((x >= x1) && (x <= x2) && (y >= y1) && (y <= y2)) {
                var dist = Math.sqrt(Math.pow(x-pxX, 2) + Math.pow(y-pxY, 2));
                possIds.push([dist, i]);
            }
        }

        console.timeEnd("cellSearch");
        if (possIds.length===0)
            return null;
        else {
            possIds.sort(function(a,b) { return a[0]-b[0]} ); // sort by distance

            // strip the distance information
            var ret = [];
            for (var i=0; i < possIds.length; i++) {
                ret.push(possIds[i][1]);
            }
            return ret;
        }
    };

    this.getSelection = function() {
        /* return selected cells as a list of ints */
        return self.selCells;
    };

    this.resetMarquee = function() {
       /* make the marquee disappear and reset its internal status */
       self.mouseDownX = null;
       self.mouseDownY = null;
       self.lastPanX = null;
       self.lastPanY = null;
       self.selectBox.style.display = "none";
       self.selectBox.style.width = 0;
       self.selectBox.style.height = 0;
    };

    this.drawMarquee = function(x1, y1, x2, y2) {
        /* draw the selection or zooming marquee using the DIVs created by setupMouse */
        var selectWidth = Math.abs(x1 - x2);
        var selectHeight = Math.abs(y1 - y2);
        var minX = Math.min(x1, x2);
        var minY = Math.min(y1, y2);
        var div = self.selectBox;
        div.style.left = minX+"px";
        div.style.top = minY+"px";
        div.style.width = selectWidth+"px";
        div.style.height = selectHeight+"px";
        div.style.display = "block";
    };

    //this.onHamster = function(ev, delta, deltaX, deltaY) {
      //console.log(delta, deltaX, deltaY);
      //console.log(ev);
      //var pxX = ev.originalEvent.clientX - self.left;
      //var pxY = ev.originalEvent.clientY - self.top;
      //self.zoomBy(1+(0.01*delta), pxX, pxY);
      //self.drawDots();
      //ev.preventDefault();
    //};

    this.onMouseMove = function(ev) {
        /* called when the mouse is moved over the Canvas */

        // set a timer so we can get "hover" functionality without too much CPU
        if (self.timer!=null)
            clearTimeout(self.timer);
        self.timer = setTimeout(self.onNoMouseMove, 170);
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

        // when the cursor is over a label, change it to a hand
        if (self.pxLabelBbox!==null && self.pxLabels!==null) {
            var labelInfo = self.labelAt(xCanvas, yCanvas);
            if (labelInfo===null)
                self.canvas.style.cursor = self.canvasCursor;
            else {
                self.canvas.style.cursor = 'pointer'; // not 'hand' anymore ! and not 'grab' yet!
                if (self.onLabelHover!==null)
                    self.onLabelHover(labelInfo, ev);
                }
        }

        if (self.mouseDownX!==null) {
            if ((ev.altKey || self.dragMode==="move") && self.panCopy!==null) {
                var xDiff = self.mouseDownX - clientX;
                var yDiff = self.mouseDownY - clientY;
                self.panBy(xDiff, yDiff);
            }
            else 
                self.drawMarquee(self.mouseDownX, self.mouseDownY, clientX, clientY);
        }
    };

    this.onNoMouseMove = function() {
        /* called after some time has elapsed and the mouse has not been moved */
        var x = self.lastMouseX - self.left; // need canvas, not screen coordinates
        var y = self.lastMouseY - self.top;
        var cellIds = self.cellsAt(x, y);
        // only call onNoCellHover if callback exists and there is nothing selected
        if (cellIds===null && self.onNoCellHover!==null && self.selCells===null )
                self.onNoCellHover();
        else if (self.onCellHover!==null)
                self.onCellHover(cellIds);
    };

    this.onMouseDown = function(ev) {
    /* user clicks onto canvas */
       console.log("background mouse down");
       var clientX = ev.clientX;
       var clientY = ev.clientY;
       if (ev.altKey || self.dragMode==="move") {
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
            var clickedLabel = self.labelAt(x2, y2);
            if (clickedLabel!==null)
                self.onLabelClick(clickedLabel, ev);
            else {
                var clickedCellIds = self.cellsAt(x2, y2);
                if (clickedCellIds!==null && self.onCellClick!==null) {
                    self.selCells = clickedCellIds;
                    self.onCellClick(clickedCellIds, ev);
                    self.drawDots();
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
       if ((self.dragMode==="zoom" && !anyKey) || ev.metaKey )
           self.zoomTo(x1, y1, x2, y2);
       // marquee select 
       else if ((self.dragMode==="select" && !anyKey) || ev.shiftKey ) {
           if (! ev.shiftKey)
               self.selectClear();
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
        console.log(ev);
        var normWheel = normalizeWheel(ev);
        console.log(normWheel);
        var pxX = ev.clientX - self.left;
        var pxY = ev.clientY - self.top;
        var spinFact = 0.1;
        if (ev.ctrlKey) // = OSX pinch and zoom gesture
            spinFact = 0.02;  // is too fast, so slow it down a little
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

    this.activateMode = function(modeName) {
    /* switch to one of the mouse drag modes: zoom, select or move */
        self.dragMode=modeName; 

        var cursor = null;

        if (modeName==="move")
            cursor = 'all-scroll';
        else if (modeName==="zoom")
            cursor = 'crosshair';
        else 
            cursor= 'default';

        self.canvas.style.cursor=cursor;
        self.canvasCursor = cursor;

        self.resetMarquee();

        self.icons["move"].style.backgroundColor = gButtonBackground; 
        self.icons["zoom"].style.backgroundColor = gButtonBackground; 
        self.icons["select"].style.backgroundColor = gButtonBackground; 
        self.icons[modeName].style.backgroundColor = gButtonBackgroundClicked; 

        //var upModeName = modeName[0].toUpperCase() + modeName.slice(1);
        //var buttonId = "mpIconMode"+upModeName;
    }

    this.randomDots = function(n, radius, mode) {
        /* draw x random dots with x random colors*/
	function randomArray(arrType, length, max) {
            /* make Array and fill it with random numbers up to max */
            var arr = new arrType(length);
            for (var i = 0; i<length; i++) {
                arr[i] = Math.round(Math.random() * max);
            }
            return arr;
	}

        if (mode!==undefined)
            self.mode = mode;
        self.radius = radius;
	self.setCoords(randomArray(Uint16Array, 2*n, 65535));
	self.setColors(["FF0000", "00FF00", "0000FF", "CC00CC", "008800"]);
	self.setColorArr(randomArray(Uint8Array, n, 4));

        console.time("draw");
        self.drawDots();
        console.timeEnd("draw");
        return self;

    };

    // object constructor code
    self.initCanvas(div, top, left, width, height);
    self.initPlot(args);

}
