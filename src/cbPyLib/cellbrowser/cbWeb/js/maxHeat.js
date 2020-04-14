"use strict";
// maxHeat: a simple scatter plot class using canvas. Data has to be binned beforehand.

/* jshint -W097 */ // allow file-wide use strict
/* jshint -W104 */  // allow some es6 parts 

function MaxHeat(div, args) {
    // a class to draw a heatmap with canvas
    // div is a div DOM element under which the canvas will be created
    // width and height: integers, in pixels
    // old-style class as otherwise not sure where to put constants labelFontSize;
    // TODO: convert to a new-style class?
    var self = this;
 
    if (!args)
        args = {};

    let labelFontSize = args.fontSize || 11; // height of row labels, will decrease with increasing row count
    let rowLabelSize = args.rowLabelSize || 60; // space for row labels on the left 
    let colLabelSize = args.colLabelSize || 34; // empty space above the heatmap
    let drawMode = 2; // 1 = stupid simple, 2 = 2x faster

    // the rest of the constructor is done at the end of this file,
    // the init involves many functions that are not defined yet

    function clearCanvas(ctx, width, height) {
    /* clear with a white background */
        // jsperf says this is fastest on Chrome, and still OK-ish in FF
        ctx.save();
        ctx.globalAlpha = 1.0;
        ctx.fillStyle = "rgb(255,255,255)";
        ctx.fillRect(0, 0, width, height);
        ctx.restore();
    }

    function addCanvasToDom(self, div, id, otherRenderer) {
        /* add a canvas element under div using the div's width and height */
        //function sepMouseMove(ev)
        //    /* called when the user moves the mouse after they have clicked onto the separator */
        //    /* this seems overly complicated. I tried split.js but it didn't make it a lot easier.
        //     * Isn't there a better plugin to handle resizing ? */
        //    let rect = sep.getBoundingClientRect();
        //    let shiftY = ev.clientY - rect.top;
        //    console.log(shiftY);
        //    sep.style.top = ev.pageY - shiftY + 'px';
        //    console.log(ev);
        //    ev.stopPropagation();
        //    return false;

        var divRect = div.getBoundingClientRect();
        var width = divRect.width;
        var height = divRect.height;

        // the selection rectangle
        var sel = document.createElement("div");
        sel.id = "hmSel";
        sel.style.position = "absolute";
        sel.style.backgroundColor = "transparent";
        sel.style.boxSizing = "border-box";
        sel.style.border = "1px solid black";
        sel.style.display = "none";
        div.appendChild(sel);
        self.selEl = sel;

        //sep.style.cursor = "row-resize";
        //sep.style.userSelect = "none";
        //sep.style.zIndex = 1000;
        
        //var sep = document.createElement("div");
        //sep.id = "hmSep";
        //sep.style.position = "absolute";
        //var parEl = div.parentElement;
        //sep.style.top = divRect.top+"px";
        //sep.style.left = divRect.left+"px";
        //sep.style.width = width+"px";
        //sep.style.height = "3px";
        //sep.style.backgroundColor = "black";
        //sep.style.cursor = "row-resize";
        //sep.style.userSelect = "none";
        //sep.style.zIndex = 1000;
        //sep.style.backgroundColor = "#BBB";
        //sep.ondragstart = function() { return false; } // disable built-in drag/drop handling
        //sep.onmousedown = function() 
            //document.body.append(sep);
            //document.onmousemove = sepMouseMove;
            //document.onmouseup = function() { document.onmousemove = document.onmouseup = null };
        //document.body.appendChild(sep);

        var canv = document.createElement("canvas");
        //canv.style.width = width+"px";
        //canv.style.height = height+"px";
        canv.id = id;
        // No scaling = one unit on screen is one pixel. Essential for speed.
        //canv.width = divRect.width;
        //canv.height = divRect.height;

        // keep these as ints, need them all the time
        //self.width = divRect.width;
        //self.height = divRect.height;
        self.canvas = canv;
        // alpha:false recommended by https://developer.mozilla.org/en-US/docs/Web/API/Canvas_API/Tutorial/Optimizing_canvas
        self.ctx = self.canvas.getContext("2d", { alpha: false });
        // by default, the canvas background is transparent+black
        // we use alpha=false, so we need to initialize the canvas with white pixels
        self.setSize(width, height);
        self.clear(); // make sure that we clear the canvas before it's added to the DOM, as it's black by default
        div.appendChild(canv); // adds the canvas to the div element

        // handle resizing

        $(div).resizable({
            handles : "n",
            start : function(ev) {
                self.origHeight = self.height;
                self.origOtherHeight = otherRenderer.height;
            },
            resize : function(ev) {
                if (otherRenderer) {
                    var newHeight = self.div.getBoundingClientRect().height;
                    otherRenderer.quickResize(null, self.origOtherHeight + (self.origHeight - newHeight));
                }
                self.setSize(self.width, newHeight);
                self.draw();
            },
            stop : function(ev) {
                if (otherRenderer) {
                    otherRenderer.setSize(null, self.origOtherHeight + (self.origHeight - self.height));
                    otherRenderer.drawDots();
                }
            },
        });

    }

    function evToRowCol(ev) {
        /* given a click or hover event, return [rowName, colName] */
        var rect = ev.target.getBoundingClientRect();
        var x = ev.clientX - rect.left; //x position within the canvas.
        var y = ev.clientY - rect.top;  //y position within the canvas.
        
        var rowIdx = parseInt((y-colLabelSize)/self.rowHeight);
        var colIdx = parseInt((x-rowLabelSize)/self.colWidth);

        let colName;
        if (x<rowLabelSize)
            colIdx = null;
        else
            colName = self.colLabels[self.colOrder[colIdx]];

        let rowName;
        if (y<colLabelSize)
            rowIdx = null;
        else
            rowName = self.rowLabels[self.rowOrder[rowIdx]];
        //console.log("mouse over coords:", x, y, rowIdx, colIdx);
        return [rowName, colName];
    }

    function onClick(ev) {
        /* call self.click(rowName, colName) */
        let rowAndCol = evToRowCol(ev);
        let rowName = rowAndCol[0];
        let colName = rowAndCol[1];
        if (self.onClick)
            self.onClick(rowName, colName);
    }

    function onMouseMove(ev) {
        /* mouse hover functionality */
        let rowAndCol = evToRowCol(ev);
        let rowName = rowAndCol[0];
        let colName = rowAndCol[1];
        
        var value = null;

        if (colIdx===null | rowIdx===null)
            self.selEl.style.display="none"
        else {
            var rowStart = self.rowStartsSizes[rowIdx*2];
            var rowHeight  = self.rowStartsSizes[rowIdx*2+1];
            var colStart = self.colStartsSizes[colIdx*2];
            var colWidth = self.colStartsSizes[colIdx*2+1];
            var selTop = rowStart;
            var selLeft = colStart;
            var selEl = self.selEl;
            selEl.style.top=(colLabelSize+selTop)+"px";
            selEl.style.left=(rowLabelSize+selLeft)+"px";
            selEl.style.height=rowHeight+"px";
            selEl.style.width=colWidth+"px";
            selEl.style.display="block"; 
            value = self.rows[rowIdx][colIdx];
        }

        if (colName==="")
            colName = "(empty)";

        if (self.onCellHover)
            self.onCellHover(self.rowOrder[rowIdx], self.colOrder[colIdx], rowName, colName, value, ev);
    }

    this.initDrawing = function (div, opts) {
        /* initialize a new plot */
        self.div = div;
        let otherRenderer = opts.mainRenderer; // otherRenderer will be resized/redrawn when needed
        addCanvasToDom(self, div, "mhCanvas", otherRenderer);
        self.onCellHover = null; // called on cell hover, arg: rowIdx, colIdx, ev. 
        // all other object variables are added by the "initPlot(args)" function below
        self.canvas.addEventListener("mousemove", onMouseMove);
        self.canvas.addEventListener("click", onClick);
    };

    this.initPlot = function(args) {
        console.log(args);
    };

    this.setPalette = function(palette) {
        self.palette = palette;
        if (self.palette.length!==this.maxVal)
            alert("internal error heat map: palette has incorrect length");
    };

    this.calcBoundaries= function() {
        /* calculate the row and column sizes. They're not all the same, due to half pixels */
        var rowCount = self.rowLabels.length;
        var colCount = self.colLabels.length;
        self.rowStartsSizes = makeIntRanges(self.height-colLabelSize, rowCount);
        self.colStartsSizes = makeIntRanges(self.width-rowLabelSize, colCount);
    }

    this.setSize = function(width, height) {
        /* change size of canvas and div and ourselves */

        self.width = width;
        self.height = height;
        self.canvas.width = width;
        self.canvas.height = height;
        self.canvas.style.width = width+"px";
        self.canvas.style.height = height+"px";

        if (self.rowLabels)
            self.calcBoundaries();

    };

    this.loadData = function(rowLabels, colLabels, rows, palette, rowOrder, colOrder) {
        /* load data into object, rowOrder and colOrder are optional */
        self.rowLabels = rowLabels;
        self.colLabels = colLabels;
        self.maxVal = palette.length; // maximum value that ever appears in 'rows'. The minimum is 0.
        self.rows = rows;
        self.checkData();
        self.setPalette(palette);
        self.calcBoundaries();

        //self.rowOrder = []; // used to convert screen to data row index, e.g. rowOrder[0] = 3
        //self.colOrder = [];

        //if (rowOrder===undefined) {
            // default row order is simply 0...n
            //for (let i=0; i<rowLabels.length; i++)
                //self.rowOrder.push(rowLabels.length-i-1);
        //} else
            //self.rowOrder = rowOrder;


        //if (colOrder===undefined) {
            // default col order is simply 0...n
            //for (let i=0; i<colLabels.length; i++)
                //self.colOrder.push(colLabels.length-i-1);
                ////self.colOrder.push(i);
        //} else
            //self.colOrder = colOrder;

        self.optimalLeaf();
    };

    function makeIntRanges(max, count) {
        /* given a max and a count, create count roughly equally sized bins from 0 to max, but as integers */
        /* The sizes all sum up to max, but every bin is not exactly the same, due to fractions e.g. for
        /* (max=10, count=3) returns [3,7,9] */
        var binSize = max/count;
        var startSizes = [];
        for (var i=0; i<count; i++) {
            var start = Math.round(i*binSize);
            var end = Math.round((i+1)*binSize);
            startSizes.push(start);
            startSizes.push(end-start);
        }
        return startSizes;
    }

    function drawRectsSimple(ctx, rowStartsSizes, colStartsSizes, rows, pal) {
        /* a completely naive implementation of the rectangle drawing */
        var rowCount = rowStartsSizes.length/2;
        var colCount = colStartsSizes.length/2;
        for (var rowIdx=0; rowIdx<rowCount; rowIdx++) {
            var rowStart = rowStartsSizes[rowIdx*2];
            var rowSize  = rowStartsSizes[rowIdx*2+1];
            var row = rows[rowIdx];

            for (var colIdx=0; colIdx < colCount; colIdx++) {
                ctx.fillStyle = "#"+pal[row[colIdx]];
                var colStart = colStartsSizes[colIdx*2];
                var colSize = colStartsSizes[colIdx*2+1];
                ctx.fillRect(rowLabelSize+colStart, colLabelSize+rowStart, colSize, rowSize);           
            }
        }    
    }

    function drawRectsOpt1(ctx, rowStartsSizes, colStartsSizes, rows, pal, maxVal, rowOrder, colOrder) {
        /* an implementation of the rectangle drawing that reduces context switches */
        // array index -> target array of x,y coords in heatmap
        // x,y positions are located at (i,i+1)
        var rowCount = rowStartsSizes.length/2;
        var colCount = colStartsSizes.length/2;
        
        var valToCoords = [];
        for (var i=0; i<maxVal+1; i++) // why +1 ?
            valToCoords.push([]);
        
        // convert from rows (array of arrays) to arrays of of coords, one array per color
        for (var rowI=0; rowI<rowCount; rowI++) {
            var dataRowI = rowOrder[rowI];
            var row = rows[dataRowI];
            for (var colI=0; colI < colCount; colI++) {
                var realColI = colOrder[colI];
                var val = row[realColI];
                var valArr = valToCoords[val];
                valArr.push(rowI);
                valArr.push(colI); 
            }
        }

        // plot the arrays, one color at a time
        for (var valI=0; valI < maxVal; valI++) {
            ctx.fillStyle = "#"+pal[valI];
            var coords = valToCoords[valI];
            for (i=0; i<coords.length/2; i++) {
                var startIdx = 2*i;
                var rowIdx = coords[startIdx];
                var colIdx = coords[startIdx+1];
                var rowStart = rowStartsSizes[rowIdx*2];
                var rowSize  = rowStartsSizes[rowIdx*2+1];
                var colStart = colStartsSizes[colIdx*2];
                var colSize = colStartsSizes[colIdx*2+1];
                ctx.fillRect(rowLabelSize+colStart, colLabelSize+rowStart, colSize, rowSize); 
            }
        }
    }

    this.draw = function() {

        this.clear();

        var pal = self.palette;
        var rows = self.rows;
        var rowCount = self.rowLabels.length;
        var colCount = self.colLabels.length;
        var rowHeight = (self.height - colLabelSize) / rowCount;
        var colWidth  = (self.width - rowLabelSize) / colCount;

        var colStartsSizes = self.colStartsSizes;
        var rowStartsSizes = self.rowStartsSizes;
        
        var ctx = self.ctx;
        ctx.save();
        
        var fontSize = labelFontSize;
        if (rowHeight < labelFontSize)
            fontSize = rowHeight-1; 
        var rowEndToTextBase = (rowHeight - fontSize)/2;
        ctx.font = fontSize+"px sans-serif";
        
        if (rowHeight>4) {
            // draw the row labels
            console.time("draw column labels");
            var rowLabels = self.rowLabels;
            for (var labelI=0; labelI<rowCount; labelI++) {
                var realLabelI = self.rowOrder[labelI]; // rows can be reordered: real... is after reordering
                var textY = parseInt(colLabelSize+(labelI*rowHeight) + rowEndToTextBase + fontSize - (0.2*rowHeight));
                ctx.fillText(rowLabels[realLabelI], 4, textY);
            }
            console.timeEnd("draw labels");
        }

        if (colLabelSize>5) {
            // draw the col labels
            console.time("draw column labels");
            var rot = -Math.PI / 7;
            var colLabels = self.colLabels;
            for (var labelI=0; labelI<colCount; labelI++) {
                var realLabelI = self.colOrder[labelI]; // rows can be reordered: real... is after reordering
                var startIdx = 2*labelI;
                var textX = rowLabelSize+colStartsSizes[startIdx];
                var rowWidth = colStartsSizes[startIdx+1];
                ctx.save();
                ctx.translate(textX+(0.3*rowWidth), colLabelSize);
                ctx.rotate(rot);
                ctx.fillText(colLabels[realLabelI], 0, 0);
                ctx.restore();
            }
            console.timeEnd("draw labels");
        }

        //ctx.strokeStyle = "rgba(1, 1, 1, 0)"; // no stroke for rectangles

        var rowStartsSizes = self.rowStartsSizes;
        var colStartsSizes = self.colStartsSizes;

        console.time("draw rects");
        switch (drawMode) {
        case 1:
            //drawRectsSimple(ctx, rowStartsSizes, colStartsSizes, rows, pal);
            alert("simple drawing mode is not supported anymore");
            break;
        case 2:
            drawRectsOpt1(ctx, rowStartsSizes, colStartsSizes, rows, pal, self.maxVal, 
                    self.rowOrder, self.colOrder);
            break;
        }
        console.timeEnd("draw rects");

        //var rowOrder = self.rowOrder;
        //var rowScreenToData = [];
        //for (var rowI=0; rowI<rowCount; rowI++)
            //rowScreenToData

        ctx.restore();

        self.rowHeight = rowHeight;
        self.colWidth = colWidth;
    };

    this.optimalLeaf = function () {
        /* order rows and columns with the optimalLeaf algorithm */
        console.time("heatmap optimal leaf ordering");
         if (self.rows.length == 1) { // if there's just one row, fudge the ordering
            self.rowOrder = [0];
            let colCount = self.rows[0].length;
            self.colOrder = [];
            for (var colIdx=0; colIdx < colCount; colIdx++) {
                self.colOrder[colIdx] = colIdx;
            }
        } else {
            let orderFunc = reorder.optimal_leaf_order();
            self.rowOrder = orderFunc(self.rows);
            let cols = reorder.transpose(self.rows);
            self.colOrder = orderFunc(cols);
        }

        console.timeEnd("heatmap optimal leaf ordering");
    };

    this.loadRandomData = function(rowCount, colCount, maxVal) {
        var rows = [];
        var rowLabels = [];
        for (var rowIdx=0; rowIdx < rowCount; rowIdx++) {
            var row = [];
            rowLabels.push( Math.random().toString(36).substring(7) );
            for (var colIdx=0; colIdx < colCount; colIdx++) {
                var val = Math.floor(Math.random() * Math.floor(maxVal));
                row.push(val);
            }
            rows.push(row);
        }

        var colLabels = [];
        for (var i=0; i < colCount; i++) {
            var rndString = Math.random().toString(36).substring(7);
            colLabels.push(rndString);
        }
        self.rows = rows;
        self.colLabels = colLabels;
        self.rowLabels = rowLabels;
        self.maxVal = maxVal;
        self.checkData();
    };

    this.checkData = function() {
        if (self.colLabels.length > self.width - rowLabelSize)
            alert("You are trying to show more columns on the heatmap than the window has pixels. "+
                "The heatmap will not be very useful until you reduce the number of columns that are shown.");
        if (self.rowLabels.length > self.height - colLabelSize)
            alert("You are trying to show more rows on the heatmap than the window has pixels. "+
                "The heatmap will not be very useful until you reduce the number of rows that are shown.");

    };

    this.clear = function() {
        clearCanvas(self.ctx, self.width, self.height);
    };

    // constructor
    self.initDrawing(div, args);
    self.initPlot(args);
}
