"use strict";
// CbCanvas: mostly a class for drawing circles onto a canvas
/*jshint globalstrict: true*/

function CbCanvas(top, left, width, height, args) {
    // a class that draws circles onto a canvas
    
    var self = this; // 'this' has two conflicting meanings in javascript. 
    // I use 'self' to refer to object variables and 'this' to refer to the caller context

    // --- object variables 
    
    if (args===undefined)
        args = {};

    self.mode = 1;
    //self.db = cbDb;
    self.ctx = null; // the canvas context
    self.width = width; // size of the canvas in pixels
    self.height = height;

    self.canvas = addCanvasToBody( left, top, width, height );

    self.minX = (args.minX | null);
    self.maxX = (args.maxX | null);
    self.minY = (args.minY | null);
    self.maxY = (args.maxY | null);

    //self.colorType = (args.colorType | "meta"); // "meta" or "gene"
    //self.colorOn = (args.colorOn | null) // gene symbol or meta data field name
    //self.coordIdx = (args.coordIdx | 0) // index of current coordinates

    self.coords     = null;
    self.clusterLabels = null; // cluster labels, array of [x,y,text]

    self.pxCoords   = null;
    self.pxLabels   = null;   // cluster labels in pixels, array of [x,y,text] 

    self.colorArr   = null;
    self.radius       = null;    // radius of the circles
    self.alpha      = 0.2;
    
    // -- (private) helper functions
    // -- these are normal functions, not methods, they do not access "self"
    function guessRadius(coordCount) {
        /* a few rules to find the optimal initial radius, depending on the number of dots */
        if (coordCount > 10000)
            return 0;
        if (coordCount > 4000)
            return 4;
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
        /* add a canvas element to the body element of the current page */
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

        document.body.appendChild(canv); // adds the canvas to the body element
        self.canvas = canv;
        self.ctx = self.canvas.getContext("2d", { alpha: false });
        // by default, the background is transparent
        // we use alpha=false, so we need to initialize the canvas with white pixels
        clearCanvas(canv, width, height);

        return canv;
    }

    function scaleData(coords, borderSize, minX, maxX, minY, maxY, winWidth, winHeight, annots) {
    /* scale list of [x (float),y (float)] to integer pixels on screen and
     * annots is an array with on-screen annotations in the format (x, y, otherInfo) that is also scaled.
     * return [array of (x (int), y (int)), scaled annots array]. Take into account the current zoom range.  
     * */
        winWidth = winWidth-(2*borderSize);
        winHeight = winHeight-(2*borderSize);

        var spanX = maxX - minX;
        var spanY = maxY - minY;
        var xMult = winWidth / spanX;
        var yMult = winHeight / spanY;

        // transform from data floats to screen pixel coordinates
        var startIdx = 0;
        var pixelCoords = new Uint16Array(coords.length-startIdx);
        //var borderAdd = 2*borderSize;
        for (var i = startIdx; i < coords.length/2; i++) {
            var x = coords[i*2];
            var y = coords[i*2+1];
            // XX ignore anything outside of current zoom range. Performance?
            if ((x < minX) || (x > maxX) || (y < minY) || (y > maxY))
                continue;
            var xPx = Math.round((x-minX)*xMult)+borderSize;
            var yPx = Math.round((y-minY)*yMult)+borderSize;
            pixelCoords[2*(i-startIdx)] = xPx;
            pixelCoords[2*(i-startIdx)+1] = yPx;
        }

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

        var ret = [];
        ret[0] = pixelCoords;
        ret[1] = newAddCoords;
        return ret;
    }

    function drawCirclesStupid(ctx, pxCoords, coordColors, colors, radius, alpha) {
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

    function drawLabels(ctx, labelCoords) {
        /* given an array of [x, y, text], draw the text */
        ctx.font = "bold 14px Sans-serif"
        ctx.globalAlpha = 1.0;
        for (var i=0; i < labelCoords.length; i++) {
            var coord = labelCoords[i];
            var x = coord[0];
            var y = coord[1];
            var text = coord[2];

            ctx.strokeStyle = '#EEEEEE'; 
            ctx.lineWidth = 4; 
            ctx.miterLimit=2;
            ctx.strokeStyle = "rgba(200, 200, 200, 0.3)";
            ctx.strokeText(text,x,y); 

            ctx.shadowBlur=6;
            ctx.shadowColor="white";
            ctx.fillStyle = "rgba(0,0,0,0.8)";
            ctx.fillText(text,x,y); 
        }
    }

    function drawCirclesDrawImage(ctx, pxCoords, coordColors, colors, radius, alpha) {
    /* blit circles onto canvas. pxCoords are the centers.  */
       // almost copied from by https://stackoverflow.com/questions/13916066/speed-up-the-drawing-of-many-points-on-a-html5-canvas-element
       // around 2x faster than drawing full circles
       // create an off-screen canvas
       console.log("Drawing "+coordColors.length+" circles with drawImg renderer");
       var diam = 2*radius;
       var off = document.createElement('canvas'); // not added to DOM, will be gc'ed
       off.width = colors.length * (2*diam);
       off.height = diam;
       var ctxOff = off.getContext('2d');  

       //pre-render circles into the off-screen canvas.
       for (var i = 0; i < colors.length; ++i) {
           ctxOff.fillStyle = "#"+colors[i];
           ctxOff.beginPath();
           ctxOff.arc(i * diam + radius, radius, radius, 0, 2 * Math.PI);
           ctxOff.closePath();
           ctxOff.fill();
       }

       if (alpha!==undefined)
           ctx.globalAlpha = alpha;

       // blit the circles onto the main canvas
       for (var i = 0; i < pxCoords.length/2; i++) {
           var pxX = pxCoords[2*i];
           var pxY = pxCoords[2*i+1];
           var col = coordColors[i];
           ctx.drawImage(off, col * diam, 0, diam, diam, pxX - radius, pxY - radius, diam, diam);
       }
    }

    function hexToInt(colors) {
    /* convert a list of hex values to ints */
        var intList = [];
        for (var i = 0; i < colors.length; i++) {
            var colHex = colors[i];
            var colInt = parseInt(colHex, 16);
            intList.push(colInt);
        }
        return intList;
    }

    function drawPixels(ctx, width, height, pxCoords, colorArr, colors, alpha) {
        /* draw single pixels into a pixel buffer and copy the buffer into a canvas */
       var canvasData = ctx.createImageData(width, height);
       var cData = canvasData.data;

       var rgbColors = hexToInt(colors);
       var invAlpha = 1.0 - alpha;

       // alpha-blend pixels into array
       for (var i = 0; i < pxCoords.length/2; i++) {
           var pxX = pxCoords[2*i];
           var pxY = pxCoords[2*i+1];
           var p = 4 * (pxY*width+pxX); // pointer to red value of pixel at x,y

           var oldR = cData[p];
           var oldG = cData[p+1];
           var oldB = cData[p+2];

           var newRgb = rgbColors[colorArr[i]];
           var newR = (newRgb >>> 16) & 0xff; // big endian?
           var newG = (newRgb >>> 8)  & 0xff;
           var newB = (newRgb)        & 0xff;

           var mixR = Math.round(oldR * invAlpha + newR * alpha);
           var mixG = Math.round(oldG * invAlpha + newG * alpha);
           var mixB = Math.round(oldB * invAlpha + newB * alpha);

           cData[p] = mixR;
           cData[p+1] = mixG;
           cData[p+2] = mixB;
           cData[p+3] = 255;
       }
       
       self.ctx.putImageData(canvasData, 0, 0);
    }

    function findRange(coords, obj) {
    /* find range of data and add as attributes minX/maxX/minY/maxY to obj, return obj */
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

        obj.minX = minX;
        obj.maxX = maxX;
        obj.minY = minY;
        obj.maxY = maxY;
        return obj; // not needed, but more explicit
    }

    function clearCanvas(canvas, width, height) {
    /* clear with a white background */
        // jsperf says this is fastest on Chrome, and still OK-ish in FF
        var ctx = self.canvas.getContext("2d", { alpha: false });
        ctx.fillStyle = "rgba(255,255,255,1)";
        ctx.fillRect(0, 0, width, height);
    }

    // -- methods
 
    this.clear = function() {
        clearCanvas(self.canvas);
    };

    this.setPos = function(left, top) {
       /* position the canvas on the page */
       self.canvas.style.left = left+"px";
       self.canvas.style.top = top+"px";
    };

    this.setSize = function(width, height) {
       /* resize canvas on the page, re-scale the data and re-draw */
       self.canvas.style.width = width+"px";
       self.canvas.style.height = height+"px";
       self.canvas.width = width;
       self.canvas.height = height;
       self.width = width;
       self.height = height;

       var s = scaleData(self.coords, self.radius, self.minX, self.maxX, 
           self.minY, self.maxY, self.width, self.height, self.clusterLabels);
       self.pxCoords = s[0];
       self.pxLabels = s[1];

       clearCanvas(self.canvas, width, height);
       self.drawDots();
    };

    this.setCoords = function(coords, clusterLabels, minX, maxX, minY, maxY) {
       /* specify new coordinates of circles to draw, an array of (x,y) coordinates */
       /* Scale data to current screen dimensions */
       /* minX, maxX, etc are optional */
       // "==null" checks for both undefined and null
       if (minX===undefined || maxX===undefined || minY===undefined || maxY===undefined)
           self = findRange(coords, self);

       self.coords = coords;
       self.clusterLabels = clusterLabels;

       var s = scaleData(self.coords, self.radius, self.minX, self.maxX, 
           self.minY, self.maxY, self.width, self.height, self.clusterLabels);
       self.pxCoords = s[0];
       self.pxLabels = s[1];
       self.radius = guessRadius(coords.length);
    };

    this.setColors = function(colorArr, colors) {
    /* set the colors, as one array with one index per coordinate, and another array with 
     * the colors for each index */
       self.colorArr = colorArr;
       self.colors = colors;
    };

    this.drawDots = function() {
        /* draw coordinates to canvas with current colors */
       if (self.pxCoords===null)
            alert("internal error: cannot draw if coordinates are not set yet");
       if (self.colorArr.length !== self.pxCoords.length*0.5)
           alert("internal error: cbDraw.drawDots - colorArr is not 1/2 of coords array");

        if (self.radius===0) {
            drawPixels(self.ctx, self.width, self.height, self.pxCoords, 
                self.colorArr, self.colors, self.alpha);
            return;
        }
            
        switch (self.mode) {
            case 0:
                drawCirclesStupid(self.ctx, self.pxCoords, self.colorArr, self.colors, self.radius, self.alpha);
                break;
            case 1:
                drawCirclesDrawImage(self.ctx, self.pxCoords, self.colorArr, self.colors, self.radius, self.alpha);
                break;
        }

        drawLabels(self.ctx, self.pxLabels);
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


    this.moveBy = function(xDiff, yDiff) {
        /* update the pxCoords by a certain x/y distance and redraw */
        // XX SLOW: we could move the canvas data, recalc pxCoords, and redraw only new area ? 

        // convert pixel range to data scale range
        var xDiffData = xDiff * ((self.maxX - self.minX) / self.width);
        var yDiffData = yDiff * ((self.maxY - self.minY) / self.height);
        
        // move zoom range 
        self.minX = self.minX + xDiffData;
        self.maxX = self.maxX + xDiffData;
        self.minY = self.minY + yDiffData;
        self.maxY = self.maxY + yDiffData;

        var s = scaleData(self.coords, self.radius, self.minX, self.maxX, 
           self.minY, self.maxY, self.width, self.height, self.clusterLabels);
        self.pxCoords = s[0];
        self.pxLabels = s[1];
    };

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
