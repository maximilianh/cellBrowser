'use strict';
// CbCanvas: mostly a class for drawing circles onto a canvas
/*jshint globalstrict: true*/

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


function CbCanvas(top, left, width, height) {
    // a class that draws circles onto a canvas
    
    var self = this; // 'this' has two conflicting meanings in javascript. 
    // I use 'self' to refer to object variables and 'this' to refer to the caller context

    // --- object variables 
    
    const gTextSize = 15;

    self.ctx = null; // the canvas context
    self.canvas = addCanvasToBody( top, left, width, height );
    
    // callbacks when user clicks or hovers over label or cell 
    self.onLabelClick = null; // gets text of label and event
    self.onCellClick = null; // gets array of cellIds and event
    self.onCellHover = null; // gets array of cellIds

    // timer that is reset on every mouse move
    self.timer = null;

    // all other object variables are added by the "initDataset(args)" function below
    
    // -- (private) helper functions
    // -- these are normal functions, not methods, they do not access "self"
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

    function setupMouseWheel(canvas) {
      canvas.on( 'DOMMouseScroll mousewheel', function ( event ) {
          var delta = event.originalEvent.wheelDelta;
          // event.originalEvent.detail > 0 || event.originalEvent.wheelDelta < 0 
          if (delta > 0) {
              //scroll down
              console.log('Down');
              zoom(-0.03);
          } else {
              //scroll up
              console.log('Up');
              zoom(0.03);
          }
      return false;
    });
    }

    function addCanvasToBody(top, left, width, height) {
        /* add a canvas element to the body element of the current page and keep left/top/width/eight in self */
        var canv = document.createElement('canvas');
        canv.id = 'tpCanvas';
        canv.style.border = "1px solid #AAAAAA";
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

        document.body.appendChild(canv); // adds the canvas to the body element
        self.canvas = canv;
        self.ctx = self.canvas.getContext("2d", { alpha: false });
        // by default, the canvas background is transparent+black
        // we use alpha=false, so we need to initialize the canvas with white pixels
        clearCanvas(self.ctx, width, height);

        return canv;
    }

    function scaleData(coords, borderSize, zoomRange, winWidth, winHeight, annots) {
    /* scale list of [x (float),y (float)] to integer pixels on screen and
     * annots is an array with on-screen annotations in the format (x, y, otherInfo) that is also scaled.
     * return [array of (x (int), y (int)), scaled annots array]. Take into account the current zoom range.  
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
        //var pixelCoords = [];
        //var pixelCount = 0;
        //var borderAdd = 2*borderSize;
        for (var i = 0; i < coords.length/2; i++) {
            var x = coords[i*2];
            var y = coords[i*2+1];
            // XX ignore anything outside of current zoom range. Performance?
            if ((x < minX) || (x > maxX) || (y < minY) || (y > maxY))
                continue;
            var xPx = Math.round((x-minX)*xMult)+borderSize;
            var yPx = Math.round((y-minY)*yMult)+borderSize;
            pixelCoords[2*i] = xPx;
            pixelCoords[2*i+1] = yPx;
            //pixelCount++;
            //pixelCoords.push(xPx);
            //pixelCoords.push(yPx);
        }
        //if (pixelCount!==coords.length/2)
            //pixelCoords = pixelCoords.slice(0, pixelCount*2);

        // also transform the labels
        var newAddCoords = [];
        if (annots!==undefined && annots!==null) {
            for (var i = 0; i < annots.length; i++) {
                var annot = annots[i];
                var x = annot[0];
                var y = annot[1];
                var other = annot[2];
                // XX ignore anything outside of current zoom range. Performance?
                if ((x < minX) || (x > maxX) || (y < minY) || (y > maxY))
                    continue;
                var xPx = Math.round((x-minX)*xMult)+borderSize;
                var yPx = Math.round((y-minY)*yMult)+borderSize;
                newAddCoords.push([xPx, yPx, other]);
            }
        }

        console.timeEnd("scale");
        var ret = [];
        ret[0] = pixelCoords;
        ret[1] = newAddCoords;
        return ret;
    }

    function drawRect(ctx, pxCoords, coordColors, colors, radius, alpha, selCells) {
       clearCanvas(ctx, width, height);
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
       clearCanvas(ctx, width, height);
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

    function drawLabels(ctx, labelCoords) {
        /* given an array of [x, y, text], draw the text. returns bounding boxes as array of [x1, y1, x2, y2]  */
        ctx.save();
        ctx.font = "bold "+gTextSize+"px Sans-serif"
        ctx.globalAlpha = 1.0;

        ctx.strokeStyle = '#EEEEEE'; 
        ctx.lineWidth = 4; 
        ctx.miterLimit=2;
        ctx.strokeStyle = "rgba(200, 200, 200, 0.3)";
        ctx.textBaseline = "top";

        ctx.shadowBlur=6;
        ctx.shadowColor="white";
        ctx.fillStyle = "rgba(0,0,0,0.8)";

        var addMargin = 1; // how many pixels to extend the bbox around the text, make clicking easier
        var bboxArr = [];
        for (var i=0; i < labelCoords.length; i++) {
            var coord = labelCoords[i];
            var x = coord[0];
            var y = coord[1];
            var text = coord[2];

            ctx.strokeText(text,x,y); 
            ctx.fillText(text,x,y); 

            var width = ctx.measureText(text).width;
            bboxArr.push( [x-addMargin, y-addMargin, x+width+addMargin, y+gTextSize+addMargin] );
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
       clearCanvas(ctx, width, height);

       ctx.save();
       console.log("Drawing "+coordColors.length+" coords with drawImg renderer, radius="+radius);
       var off = document.createElement('canvas'); // not added to DOM, will be gc'ed
       var diam = 2*radius;
       var tileWidth = diam+2;
       var tileHeight = diam+2;
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
               var strokeCol = "#"+shadeColor(colors[i], 0.3);
               ctxOff.strokeStyle=strokeCol;

               //ctxOff.beginPath();
               //ctxOff.arc(i * diam + radius, radius, radius, 0, 2 * Math.PI);
               //ctxOff.closePath();
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
       
       clearCanvas(ctx, width, height);
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
       clearCanvas(ctx, width, height);
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
 
    this.initDataset = function(args) {
        if (args===undefined)
            args = {};

        self.mode = 1;   // drawing mode

        self.zoomRange = {};
        self.zoomRange.minX = getAttr(args, "minX", null);
        self.zoomRange.maxX = getAttr(args, "maxX", null);
        self.zoomRange.minY = getAttr(args, "minY", null);
        self.zoomRange.maxY = getAttr(args, "maxX", null);

        self.coords     = null;
        self.clusterLabels = null; // cluster labels, array of [x,y,text]

        self.pxCoords   = null;   // coordinates of cells as screen pixels or 0,0 if not shown
        self.pxLabels   = null;   // cluster labels in pixels, array of [x,y,text] 
        self.pxLabelBbox = null;   // cluster label bounding boxes, array of [x1,x2,x2,y2]

        self.colors     = null;   // list of six-digit hex codes
        self.colorArr   = null;   // length is pxCoords/2, one byte per cell = index into self.colors
        self.radius     = getAttr(args, "radius", null);    // current radius of the circles, 0=one pixel dots
        self.alpha      = getAttr(args, "alpha", 0.3);
        self.selCells   = null;  // IDs of cells that are "selected" and as such highlighted in some way

        // we keep a copy of the 'initial' arguments at 100% zoom
        self.initZoom   = cloneObj(self.zoomRange);
        self.initRadius = self.radius;                      // circle radius at full zoom

        // mouse drag is modal: can be "select", "move" or "zoom"
        self.dragMode = "zoom";
         
        // for zooming and panning
        self.mouseDownX = null;
        self.mouseDownY = null;
        self.panCopy    = null;

        // to detect if user just clicked on a dot
        self.dotClickX = null;
        self.dotClickY = null;

    };
    
    this.clear = function() {
        clearCanvas(self.ctx, self.width, self.height);
    };

    this.setPos = function(left, top) {
       /* position the canvas on the page */
       self.canvas.style.left = left+"px";
       self.canvas.style.top = top+"px";
    };

    this.scaleData = function() {
       /* scale data to current zoom range */
       var s = scaleData(self.coords, self.radius, self.zoomRange, self.width, self.height, self.clusterLabels);
       self.pxCoords = s[0];
       self.pxLabels = s[1];
    }

    this.setSize = function(width, height) {
       /* resize canvas on the page re-scale the data and re-draw */
       self.canvas.style.width = width+"px";
       self.canvas.style.height = height+"px";
       self.canvas.width = width;
       self.canvas.height = height;
       self.width = width;
       self.height = height;

       self.scaleData();
       clearCanvas(self.ctx, width, height);
       self.drawDots();
    };

    this.setCoords = function(coords, clusterLabels, minX, maxX, minY, maxY) {
       /* specify new coordinates of circles to draw, an array of (x,y) coordinates */
       /* Scale data to current screen dimensions */
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

       self.scaleData();
       if (self.radius===null || self.radius==undefined)
           self.radius = guessRadius(coords.length);
    };

    this.setColorArr = function(colorArr) {
    /* set the color array, one array with one index per coordinate */
       self.colorArr = colorArr;
    };

    this.setColors = function(colors) {
    /* set the colors, one for each value of the color array. Values are hex strings. */
       self.colors = colors;
    };

    this.drawDots = function() {
        /* draw coordinates to canvas with current colors */
        if (self.alpha===undefined)
             alert("internal error: alpha is not defined");
        if (self.pxCoords===null)
             alert("internal error: cannot draw if coordinates are not set yet");
        if (self.colorArr.length !== self.pxCoords.length*0.5)
            alert("internal error: cbDraw.drawDots - colorArr is not 1/2 of coords array");

        self.zoomFact = ((self.initZoom.maxX-self.initZoom.minX)/(self.zoomRange.maxX-self.zoomRange.minX));

        // we make the circles only half as wide as expected
        self.radius = Math.floor((Math.max(0.8, self.initRadius) * self.zoomFact) * 0.4);

        // at high zoom levels, make the circles a bit smaller, to reduce overlap
        var drawRadius = self.radius;
        if (self.radius > 30)
            drawRadius = Math.round(drawRadius*0.3);
        else if (self.radius > 18)
            drawRadius = Math.round(drawRadius*0.4);
        else if (self.radius > 5)
            drawRadius= Math.round(drawRadius*0.7);
        //if (drawRadius > 6)
            //drawRadius = Math.round(drawRadius*0.6);
        console.log("Theoretical radius "+self.radius+", corrected radius "+drawRadius);
        self.radius = drawRadius;
                
        console.log("drawing, zoom factor "+self.zoomFact+", radius "+self.radius+", alpha "+self.alpha);
        console.time("draw");

        if (self.radius===0) {
            drawPixels(self.ctx, self.width, self.height, self.pxCoords, 
                self.colorArr, self.colors, self.alpha, self.selCells);
        }

        else if (self.radius===1 || self.radius===2) {
            drawRect(self.ctx, self.pxCoords, self.colorArr, self.colors, self.radius, self.alpha, self.selCells);
        }
        else {
            // the higher the zoom factor, the higher the alpha value
            var alpha = self.alpha;
            if (self.radius>4)
                alpha *= 2.0;
            else if (self.radius>10)
                alpha = 0.8;
            alpha = Math.min(1.0, alpha);

            switch (self.mode) {
                case 0:
                    drawCirclesStupid(self.ctx, self.pxCoords, self.colorArr, self.colors, drawRadius, alpha, self.selCells);
                    break;
                case 1:
                    drawCirclesDrawImage(self.ctx, self.pxCoords, self.colorArr, self.colors, drawRadius, alpha, self.selCells);
                    break;
                case 2:
                    break;
            }
        }

        console.timeEnd("draw");
        self.pxLabelBbox = drawLabels(self.ctx, self.pxLabels);
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
       /* zoom to 100% */
       self.zoomRange = cloneObj(self.initZoom);
       self.scaleData();
       self.radius = self.initRadius;
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
       var aspectRatio = self.width / self.height;
       var rectWidth  = (pxMaxX-pxMinX);
       var newHeight = rectWidth/aspectRatio;
       pxMaxY = pxMinY + newHeight;

       // adapt the circle size
       var currRadius = self.radius;
       if (currRadius===0) // 0 = 1 pixel wide
           currRadius = 0.5;

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

    this.zoomBy = function(scale) {
    /* zoom to the center by a given factor. Returns new zoom range. */
        var zoomRange = self.zoomRange;
        var xRange = Math.abs(zoomRange.maxX-zoomRange.minX);
        zoomRange.minX = zoomRange.minX - (xRange*scale);
        zoomRange.maxX = zoomRange.maxX + (xRange*scale);

        var yRange = Math.abs(zoomRange.maxY-zoomRange.minY);
        zoomRange.minY = zoomRange.minY - (yRange*scale);
        zoomRange.maxY = zoomRange.maxY + (yRange*scale);
       
        self.scaleData();
        return zoomRange;
    };

    this.panStart = function() {
       /* called when starting a panning sequence, makes a snapshop of the current image */
       self.panCopy = document.createElement('canvas'); // not added to DOM, will be gc'ed
       self.panCopy.width = self.width;
       self.panCopy.height = self.height;
       var destCtx = self.panCopy.getContext("2d", { alpha: false });
       destCtx.drawImage(self.canvas, 0, 0);
    }

    this.panBy = function(xDiff, yDiff) {
        /* pan current image by x/y pixels */
        console.log('panning by '+xDiff+' '+yDiff);

       //var srcCtx = self.panCopy.getContext("2d", { alpha: false });
       clearCanvas(self.ctx, self.width, self.height);
       self.ctx.drawImage(self.panCopy, -xDiff, -yDiff);
       // keep these for panEnd
       self.panDiffX = xDiff;
       self.panDiffY = yDiff;
    }

    this.panEnd = function() {
        /* end a sequence of panBy calls, called when the mouse is released */
        self.moveBy(self.panDiffX, self.panDiffY);
        self.panCopy = null;
        self.panDiffX = null;
        self.panDiffY = null;
    }

    this.selectClear = function() {
        /* clear selection */
        self.selCells = null;
    };

    this.selectAdd = function(cellIdx) {
        /* add a single cell to the selection. If it already exists, remove it. */
        if (self.selCells===null) {
            self.selCells = [cellIdx];
            return;
        }

        console.time("selectAdd");
        var foundIdx = self.selCells.indexOf(cellIdx);
        if (foundIdx===-1)
            self.selCells.push(cellIdx);
        else
            self.selCells.splice(foundIdx, 1);
        console.time("selectAdd");
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
        console.log(self.selCells.length+" total cells selected");
        console.timeEnd("select");
    };

    this.moveBy = function(xDiff, yDiff) {
        /* update the pxCoords by a certain x/y distance and redraw */

        // convert pixel range to data scale range
        var zr = self.zoomRange;
        var xDiffData = xDiff * ((zr.maxX - zr.minX) / self.width);
        var yDiffData = yDiff * ((zr.maxY - zr.minY) / self.height);
        
        // move zoom range 
        zr.minX = zr.minX + xDiffData;
        zr.maxX = zr.maxX + xDiffData;
        zr.minY = zr.minY + yDiffData;
        zr.maxY = zr.maxY + yDiffData;

        var s = scaleData(self.coords, self.radius, self.zoomRange, self.width, self.height, self.clusterLabels);
        self.pxCoords = s[0];
        self.pxLabels = s[1];
    };

    this.labelAt = function(x, y) {
        /* return the text of the label at position x,y or null if nothing there */
        console.time("labelCheck");
        var labelCoords = self.pxLabels;
        var boxes = self.pxLabelBbox;
        for (var i=0; i < labelCoords.length; i++) {
            var label = labelCoords[2];
            var box = boxes[i];
            var x1 = box[0];
            var y1 = box[1];
            var x2 = box[2];
            var y2 = box[3];
            if ((x >= x1) && (x <= x2) && (y >= y1) && (y <= y2))
                return label;
        }
        console.timeEnd("labelCheck");
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
    }

    this.drawMarquee = function(x1, y1, x2, y2) {
        /* draw the selection or zooming marquee using the DIVs created by setupMouse */
        var selectWidth = Math.abs(x1 - x2);
        var selectHeight = Math.abs(y1 - y2);
        var minX = Math.min(x1, x2);
        var minY = Math.min(y1, y2);
        var posCss = {
             "left":minX, 
             "top": minY,
             "width":selectWidth,
             "height":selectHeight
        };
        $("#tpSelectBox").css(posCss).show();
    }

    this.onMouseMove = function(ev) {
        /* called when the mouse is moved over the Canvas */
        //console.log("background move");
        if (self.timer!=null)
            clearTimeout(self.timer);
        self.timer = setTimeout(self.onNoMouseMove, 170);

        // save mouse pos for onNoMouseMove
        self.lastMouseX = ev.clientX;
        self.lastMouseY = ev.clientY;

        //if (self.mouseDownX === null) {
            // can stop quickly if no button pressed and no panning active
            //ev.preventDefault();
            //return;
        //}

        var clientX = ev.clientX;
        var clientY = ev.clientY;
        if (self.mouseDownX!==null) {
            if ((ev.altKey || self.dragMode==="move") && self.panCopy!==null) {
                var xDiff = self.mouseDownX - clientX;
                var yDiff = self.mouseDownY - clientY;
                self.panBy(xDiff, yDiff);
            }
            else 
                self.drawMarquee(self.mouseDownX, self.mouseDownY, clientX, clientY);
        }
        ev.preventDefault(); 
    }

    this.onNoMouseMove = function() {
        /* called after some time has elapsed and the mouse has not been moved */
        var x = self.lastMouseX - self.left; // need canvas, not screen coordinates
        var y = self.lastMouseY - self.top;
        var cellIds = self.cellsAt(x, y);
        if (cellIds!==null)
            self.onCellHover(cellIds);
    };

    this.onMouseDown = function(ev) {
    /* user clicks onto canvas */
       console.log("background mouse down");
       var clientX = ev.clientX;
       var clientY = ev.clientY;
       if (ev.altKey || self.dragMode==="move") {
           console.log("background mouse down, with meta");
           self.panStart();
       } 
       self.mouseDownX = clientX;
       self.mouseDownY = clientY;
       ev.preventDefault(); 
    };

    this.onMouseUp = function(ev) {
       console.log("background mouse up");
       if (self.panCopy!==null) {
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
       // these are screen coordinates
       var clientX = ev.clientX;
       var clientY = ev.clientY;
       
       // the subsequent operations require canvas coordinates x/y
       var canvasTop = self.top;
       var canvasLeft = self.left;
       var x1 = self.mouseDownX - canvasLeft;
       var y1 = self.mouseDownY - canvasTop;
       var x2 = clientX - canvasLeft;
       var y2 = clientY - canvasTop;

       // user did not move the mouse, so this is a click
       if (self.mouseDownX === clientX && self.mouseDownY === clientY) {
            var clickedLabel = self.labelAt(x2, y2);
            if (clickedLabel!==null)
                self.onLabelClick(clickedLabel, ev);
            else {
                var clickedCellIds = self.cellsAt(x2, y2);
                if (clickedCellIds!==null)
                    self.onCellClick(clickedCellIds, ev);
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
            //$("#tpSelectBox").hide();
            self.mouseDownX = null;
            self.mouseDownY = null;
            return;
       }
       //console.log("moved: reset "+x+" "+mouseDownX+" "+mouseDownY+" "+y);

       var anyKey = (ev.metaKey || ev.altKey || ev.shiftKey);

       // panning
       //if ((self.dragMode==="move" && !anyKey) || ev.altKey ) {
           //self.moveTo(x1, y1, x2, y2);
       //}
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

       self.mouseDownX = null;
       self.mouseDownY = null;
       self.lastPanX = null;
       self.lastPanY = null;
       $("#tpSelectBox").hide();
       $("#tpSelectBox").css({"width":0, "height":0});

       self.drawDots();
       ev.preventDefault(); 
    }

    this.setupMouse = function() {
       // add the div used for the mouse selection/zoom rectangle to the DOM
       var htmls = [];
       htmls.push("<div style='position:absolute; display:none; pointer-events:none;' id='tpSelectBox'></div>");
       $(document.body).append(htmls.join(""));
       // setup the mouse zooming callbacks
       self.canvas.onmousedown = self.onMouseDown;
       self.canvas.onmousemove = self.onMouseMove;
       self.canvas.onmouseup = self.onMouseUp;
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
	self.setColors(randomArray(Uint8Array, n, 5), ["FF0000", "00FF00", "0000FF", "CC00CC", "008800"]);

        console.time("draw");
        self.drawDots();
        console.timeEnd("draw");
        return self;
    };


}
