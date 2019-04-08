"use strict";
// maxHeat: a simple scatter plot class using canvas. Data has to be binned beforehand.

function MaxHeat(div, args) {
    // a class to draw a heatmap with canvas
    // div is a div DOM element under which the canvas will be created
    // width and height: integers, in pixels
    // old-style class as otherwise not sure where to put constants labelFontSize;

    var self = this;
    var labelFontSize = 16; // height of row labels, will decrease with increasing row count
    const gRowLabelSize = 50; // width of left text bar
    const gColLabelSize = 30; // empty space above the heatmap
    const drawMode = 2; // 1 = stupid simple, 2 = 2x faster

    // the rest of the initialization is done at the end of this file,
    // because the init involves many functions that are not defined yet here

    function clearCanvas(ctx, width, height) {
    /* clear with a white background */
        // jsperf says this is fastest on Chrome, and still OK-ish in FF
        ctx.save();
        ctx.globalAlpha = 1.0;
        ctx.fillStyle = "rgb(255,255,255)";
        ctx.fillRect(0, 0, width, height);
        ctx.restore();
    }

    function addCanvasToDom(self, div, id) {
        /* add a canvas element under div using the div's width and height */
        var canv = document.createElement("canvas");
        canv.id = id;
        var rect = div.getBoundingClientRect();
        canv.style.width = rect.width+"px";
        canv.style.height = rect.height+"px";
        // No scaling = one unit on screen is one pixel. Essential for speed.
        canv.width = rect.width;
        canv.height = rect.height;

        // need to keep these as ints, need them all the time
        self.width = rect.width;
        self.height = rect.height;
        div.appendChild(canv); // adds the canvas to the div element
        self.canvas = canv;
        // alpha:false recommended by https://developer.mozilla.org/en-US/docs/Web/API/Canvas_API/Tutorial/Optimizing_canvas
        self.ctx = self.canvas.getContext("2d", { alpha: false });
        // by default, the canvas background is transparent+black
        // we use alpha=false, so we need to initialize the canvas with white pixels
        clearCanvas(self.ctx, self.width, self.height);
    }

    function onMouseMove(ev) {
        //console.log(ev);
        var rect = ev.target.getBoundingClientRect();
        var x = ev.clientX - rect.left; //x position within the element.
        var y = ev.clientY - rect.top;  //y position within the element.
        
        var rowIdx = parseInt((y-gColLabelSize)/self.rowHeight);
        var colIdx = parseInt((x-gRowLabelSize)/self.colWidth);
        if (x<gColLabelSize)
            colIdx = null;
        if (y<gRowLabelSize)
            rowIdx = null;
        console.log("mouse over coords:", x, y, rowIdx, colIdx);

        if (self.onCellHover)
            self.onCellHover(rowIdx, colIdx, ev);
    }

    this.initDrawing = function (div) {
        /* initialize a new plot */
        self.div = div;
        addCanvasToDom(self, div, "mhCanvas");
        self.onCellHover = null; // called on cell hover, arg: rowIdx, colIdx, ev. 
        // all other object variables are added by the "initPlot(args)" function below
        self.canvas.addEventListener("mousemove", onMouseMove);
    };

    this.initPlot = function(args) {
        console.log(args);
    };

    this.setPalette = function(palette) {
        self.palette = palette;
    };

    this.loadData = function(rowLabels, colLabels, rows, maxVal) {
        /* load data into object */
        self.rowLabels = rowLabels;
        self.colLabels = colLabels;
        self.maxVal = maxVal; // maximum value that ever appears in 'rows', we assume that the minimum is 0.
        self.rows = rows;
        self.checkData();
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
                ctx.fillStyle = pal[row[colIdx]];
                var colStart = colStartsSizes[colIdx*2];
                var colSize = colStartsSizes[colIdx*2+1];
                ctx.fillRect(gRowLabelSize+colStart, gColLabelSize+rowStart, colSize, rowSize);           
            }
        }    
    }

    function drawRectsOpt1(ctx, rowStartsSizes, colStartsSizes, rows, pal, maxVal) {
        /* an implementation of the rectangle drawing that reduces context switches */
        // array index -> target array of x,y coords in heatmap
        // x,y positions are located at (i,i+1)
        var rowCount = rowStartsSizes.length/2;
        var colCount = colStartsSizes.length/2;
        
        var valToCoords = [];
        for (var i=0; i<maxVal; i++)
            valToCoords.push([]);
        
        // fill the arrays
        for (var rowI=0; rowI<rowCount; rowI++) {
            var row = rows[rowI];
            for (var colI=0; colI < colCount; colI++) {
                var val = row[colI];
                var valArr = valToCoords[val];
                valArr.push(rowI);
                valArr.push(colI); 
            }
        }

        // plot the arrays, one color at a time
        for (var valI=0; valI < maxVal; valI++) {
            ctx.fillStyle = pal[valI];
            var coords = valToCoords[valI];
            for (i=0; i<coords.length/2; i++) {
                var startIdx = 2*i;
                var rowIdx = coords[startIdx];
                var colIdx = coords[startIdx+1];
                var rowStart = rowStartsSizes[rowIdx*2];
                var rowSize  = rowStartsSizes[rowIdx*2+1];
                var colStart = colStartsSizes[colIdx*2];
                var colSize = colStartsSizes[colIdx*2+1];
                ctx.fillRect(gRowLabelSize+colStart, gColLabelSize+rowStart, colSize, rowSize); 
            }
        }
    }

    this.draw = function() {

        this.clear();

        var pal = self.palette;
        var rows = self.rows;
        var rowCount = self.rowLabels.length;
        var colCount = self.colLabels.length;
        var rowHeight = (self.height - gColLabelSize) / rowCount;
        var colWidth  = (self.width - gRowLabelSize) / colCount;

        var ctx = self.ctx;
        ctx.save();
        
        if (rowHeight < labelFontSize)
            labelFontSize = rowHeight-2; 
        var rowEndToTextBase = (rowHeight - labelFontSize)/2;
        
        if (rowHeight>8) {
            // draw the row labels
            console.time("draw labels");
            var rowLabels = self.rowLabels;
            for (var labelI=0; labelI<rowCount; labelI++) {
                var textY = parseInt(gColLabelSize+(labelI*rowHeight) + rowEndToTextBase + labelFontSize);
                ctx.fillText(rowLabels[labelI], 4, textY);
            }
            console.timeEnd("draw labels");
        }

        //ctx.strokeStyle = "rgba(1, 1, 1, 0)"; // no stroke for rectangles

        console.time("draw rects");
        var rowStartsSizes = makeIntRanges(self.height, rowCount);
        var colStartsSizes = makeIntRanges(self.width-gRowLabelSize, colCount);

        switch (drawMode) {
        case 1:
            drawRectsSimple(ctx, rowStartsSizes, colStartsSizes, rows, pal);
            break;
        case 2:
            drawRectsOpt1(ctx, rowStartsSizes, colStartsSizes, rows, pal, self.maxVal);
        }
        console.timeEnd("draw rects");

        ctx.restore();

        self.rowHeight = rowHeight;
        self.colWidth = colWidth;
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
        if (self.colLabels.length > self.width - gRowLabelSize)
            alert("You are trying to show more columns on the heatmap than the window has pixels. "+
                "The heatmap will not be very useful until you reduce the number of columns that are shown.");
        if (self.rowLabels.length > self.height - gColLabelSize)
            alert("You are trying to show more rows on the heatmap than the window has pixels. "+
                "The heatmap will not be very useful until you reduce the number of rows that are shown.");

    };

    this.clear = function() {
        clearCanvas(self.ctx, self.width, self.height);
    };

    // constructor
    self.initDrawing(div);
    self.initPlot(args);
}
