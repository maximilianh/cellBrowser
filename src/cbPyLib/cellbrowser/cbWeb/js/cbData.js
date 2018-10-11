// scDb: a class accessing single cell data from a URL
'use strict';
/*jshint globalstrict: true*/

/* a module with some helper functions */
var cbUtil = (function () {
    var my = {};

    my.joinPaths = function joinPaths(parts, separator) {
    // join paths together, taking care of duplicated /s
      return parts.map(function(part) { return part.trim().replace(/(^[\/]*|[\/]*$)/g, ''); }).join(separator || '/');
    };

    my.dumpObj = function (o) {
    /* for debugging */
        console.log(JSON.stringify(o));
    };

    my.loadJson = function(url, onSuccess, silent) {
    /* load json file from url and run onSuccess when done. Alert if it doesn't work. */
        var req = jQuery.ajax({
            "url" : url, 
            type : "GET",
            dataType : "json",
            mimeType : 'text/plain; charset=x-user-defined', // for local files, avoids errors
            success: function(data) {
                onSuccess(data);
            },
            error : function() {
                if (!silent)
                    alert("Could not load "+url); 
                onSuccess(null);
            }
        });
    };
    
    my.loadFile = function(url, arrType, onDone, onProgress, otherInfo, start, end) {
        /* load text or binary file with HTTP GET into fileData variable and call function
         * onDone(fileData, otherInfo) when done. .gz URLs are automatically decompressed.
         * optional: return as an object of arrType. 
         * If it is a text file, specify arrType=null.
         * optional: byte range request for start-end (0-based, inclusive).
         * */
        var oReq = new XMLHttpRequest();
        oReq.open("GET", url, true);

        if (arrType!==null)
            oReq.responseType = "arraybuffer";
        else
            oReq.responseType = "text";

        oReq.onload = cbUtil.onDoneBinaryData;
        oReq.onprogress = onProgress;
        oReq.onerror = function(e) { alert("Could not load file "+url); };
        oReq._onDone = onDone; // keep this for the callback
        oReq._otherInfo = otherInfo; // keep this for the callback
        oReq._arrType = arrType; // keep this for the callback casting
        oReq._url = url; // keep this for the callback error messages
        if (start!==undefined) {
            oReq.setRequestHeader('Range', "bytes="+start+"-"+end); // the bytes (incl.) you request
            oReq._start = start;
            oReq._end = end;
        }
        oReq.send(null);
    };

    my.loadTsvFile = function(url, onDone, addInfo) {
    /* load a tsv file relative to baseUrl and call a function when done */
        Papa.parse(url, {
                delimiter : "\t",
                skipEmptyLines : true,
                fastMode : true,
                download: true,
                complete: function(results, localFile) {
                            onDone(results, localFile, addInfo);
                        },
                error: function(err, file) {
                            if (addInfo!==undefined)
                                alert("could not load "+url);
                        }
                });
    };

    my.onDoneBinaryData = function(oEvent) {
        /* called when binary file has been loaded. gunzips .gz URLs. */
        var url = this._url;
        // 200 = OK, 206 = Partial Content OK
        if (this.status !== 200 && this.status !== 206) {
            alert("Could not load "+url+", error " + oEvent.statusText);
            return;
        }
        
        var binData = this.response;

        if (!binData) {
          alert("internal error when loading "+url+": no reponse from server?");
          return;
        }

        if (url.endsWith(".gz"))
            binData = pako.ungzip(binData);

        // if the user wants only a byte range...
        if (this._start!==undefined) {
            // check if the expected length is OK. Some webservers don't support byte range 
            // requests.
            var expLength = this._end-this._start+1; // byte range is inclusive
            if (binData.byteLength < expLength)
                alert("internal error cbData.js: chunk is too small. Does the HTTP server really support byte range requests?");

            if (binData.byteLength > expLength) {
                console.log("Webserver does not support byte range requests, working around it, but this may be slow");
                if (this._arrType) {
                    binData = new Uint8Array(binData);
                    binData = binData.slice(this._start, this._end).buffer; // slice does not include the end
                }
                else // binData is actually a normal string
                    binData = binData.slice(this._start, this._end); // slice does not include the end
                }
        }

        if (this._arrType)
            binData = new this._arrType(binData);
        this._onDone(binData, this._otherInfo);
    };

    my.makeType = function(typeStr) {
        /* given a string, return the correct type array for it */
        typeStr = typeStr.toLowerCase();
        if ((typeStr==="double" || typeStr=="float64"))
            return Float64Array;
        if ((typeStr==="float" || typeStr=="float32"))
            return Float32Array;
        else if (typeStr==="uint32" || typeStr=="dword")
            return Uint32Array;
        else if (typeStr==="uint16" || typeStr=="word")
            return Uint16Array;
        else if (typeStr==="uint8" || typeStr==="byte")
            return Uint8Array;
        else
            alert("Not a valid array type: "+typeStr);
            return null;
    };

    my.findObjWhereEq = function(objArr, keyName, searchName) {
        // given a list of objects, return the one where keyName==searchName
        var found = null;
        for (var i = 0; i < objArr.length; i++) {
            var el = objArr[i];
            if (el[keyName]===searchName) {
                found = el;
                break;
            }
        }
        return found;
    };

    my.findIdxWhereEq = function(objArr, keyName, searchName) {
        // given a list of objects, return the index where keyName==searchName
        var found = null;
        for (var i = 0; i < objArr.length; i++) {
            var el = objArr[i];
            if (el[keyName]===searchName) {
                found = i;
                break;
            }
        }
        return found;
    };

    my.baReadOffset = function(ba, o) {
    /* given a byte array, return the long int (little endian) at offset o */
        var offset = ba[o] |  ba[o+1] << 8 | ba[o+2] << 16 | ba[o+3] << 24;
        return offset;
    };

    my.baReadUint16 = function(ba, o) {
    /* read 16 bits, little endian, from byte array */
        var num = ba[o] | ba[o+1] << 16;
        return num;
    };


    return my;
}());


function CbDbFile(url) {
    // a class that loads all data from binary files loading for the cell browser: 
    // loading coordinates, loading meta data, resolve sample names to sample
    // indices, resolve sample indices to sample names load meta for a given
    // cell index
    var self = this; // this has two conflicting meanings in javascript. 
    // we use 'self' to refer to object variables and 'this' to refer to the calling object

    self.name = url;
    self.url = url;
    self.geneOffsets = null;

    self.exprCache = {}; // cached compressed expression vectors

    this.conf = null;

    this.loadConfig = function(onDone) {
    /* load config and gene offsets from URL and call func when done */

        var doneCount = 0;
        function gotOneFile() {
            doneCount++;
            if (doneCount===2)
                onDone();
        }

        // load config and call onDone
        var dsUrl = cbUtil.joinPaths([this.url, "dataset.json"]);
        // should I deactivate the cache here?
        // this is a file that users typicaly change often
        dsUrl = dsUrl+"?"+Math.floor(Math.random()*100000000);
        cbUtil.loadJson(dsUrl, function(data) { self.conf = data; gotOneFile();});

        // start loading gene offsets, this takes a while
        var osUrl = cbUtil.joinPaths([this.url, "exprMatrix.json"]);
        cbUtil.loadJson(osUrl, function(data) { self.geneOffsets = data; gotOneFile();});
    };

    this.loadCoords = function(coordIdx, onDone, onProgress) {
    /* load coordinates from URL and call onDone(array of (x,y), coordInfoObj, labelMids) when done */
        //var coordInfo = cbUtil.findObjWhereEq(self.conf.coords, "name", this.coordName);
        var i = 0;
        var binData = null;
        var meta = null;
        var labelMids = undefined; // null means: no json file found

        function binDone(data, other) {
            binData = data;
            meta = other;
            if (labelMids!==undefined)
                onDone(binData, meta, labelMids);
        }
        function jsonDone(data) {
            labelMids = data;
            if (binData!==null)
                onDone(binData, meta, labelMids);
        }

        var coordInfo = self.conf.coords[coordIdx];
        if (coordInfo===null) {
           alert("Could not find coordinates with name "+this.coordName);
           return;
        }
        var binUrl = cbUtil.joinPaths([self.url, "coords", coordInfo.name, "coords.bin"]);
        var arrType = cbUtil.makeType(coordInfo.type || "float64");
        cbUtil.loadFile(binUrl, arrType, binDone, onProgress, coordInfo);

        var jsonUrl = cbUtil.joinPaths([self.url, "coords", coordInfo.name, "clusterLabels.json"]);
        cbUtil.loadJson(jsonUrl, jsonDone, true);
    };

    this.getDefaultColorField = function() {
        /* return a pair: [0] is "meta" or "gene" and [1] is the field name or the gene */
        return ["meta", self.conf.clusterField];
    };

    this.fieldNameToIndex = function(fieldName) {
        /* given a meta field name, return its meta table index */
        return cbUtil.findIdxWhereEq(self.conf.metaFields, "name", fieldName);
    };

    this.loadMetaVec = function(fieldIdx, onDone, onProgress) {
    /* get an array of numbers, one per cell, that reflect the meta field contents
     * and an object with some info about the field. call onDone(arr, metaInfo) when done. */
        //var metaInfo = cbUtil.findObjWhereEq(self.conf.metaFields, "name", fieldName);
        var metaInfo = self.conf.metaFields[fieldIdx];
        console.log(metaInfo);
        var fieldName = metaInfo.name;
        var binUrl = cbUtil.joinPaths([self.url, "metaFields", fieldName+".bin.gz"]);
        var arrType = cbUtil.makeType(metaInfo.arrType);
        cbUtil.loadFile(binUrl, arrType, onDone, onProgress, metaInfo);
    };

    this.loadMetaForCell = function(cellIdx, onDone, onProgress) {
    /* for a given cell, call onDone with an array of the metadata values, as strings. */
        // first we need to lookup the offset of the line and its length from the index
        var url = cbUtil.joinPaths([self.url, "meta.index"]);
        var start = (cellIdx*6); // four bytes for the offset + 2 bytes for the line length
        var end   = start+6;

        function lineDone(text) {
            /* called when the line from meta.tsv has been read */
            var fields = text.split("\t");
            fields[fields.length-1] = fields[fields.length-1].trim(); // remove newline
            onDone(fields);
        }

        function offsetDone(arr) {
            /* called when the offset in meta.index has been read */
            var offset = cbUtil.baReadOffset(arr, 0);
            var lineLen = cbUtil.baReadUint16(arr, 4);
            // now get the line from the .tsv file
            var url = cbUtil.joinPaths([self.url, "meta.tsv"]);
            cbUtil.loadFile(url+"?"+cellIdx, null, lineDone, onProgress, null, 
                offset, offset+lineLen);
        }

        // chrome caching sometimes fails with byte range requests, so add cellidx to the URL
        cbUtil.loadFile(url+"?"+cellIdx, Uint8Array, function(byteArr) { offsetDone(byteArr); }, 
            undefined, null, start, end);
    };

    function sortArrOfArr(arr, j) {
        /* sort an array of arrays by the j-th element */
        arr.sort(function(a, b) { 
            return a[j] > b[j] ? 1 : -1;
            });
    }

    function countAndSort(arr) {
        /* count values in array, return an array of [value, count] */
        //var counts = {};
        var counts = new Map();
        for (var i=0; i < arr.length; i++) {
            var num = arr[i];
            //counts[num] = counts[num] ? counts[num] + 1 : 1;
            counts.set(num, (counts.get(num) | 0) + 1);
        }
        var entries = Array.from(counts.entries());
        entries.sort(function(a,b) { return a[0]-b[0]});
        return entries;
    }

    function arrToEnum(arr, counts) {
        /* replace values in array with their enum-index */
        // make a mapping value -> bin
        var valToBin = {};
        sortArrOfArr(counts, 0); // sort by value
        for (var i=0; i<counts.length; i++) {
            var val = counts[i];
            valToBin[val] = i;
        }

        // apply the mapping
        var dArr = new Uint8Array(arr.length);
        for (var i=0; i<arr.length; i++) {
            dArr[i] = valToBin[arr[i]];
        }

        var binInfo = [];
        for (var i=0; i < counts.length; i++) {
            var val = parseFloat(counts[i][0]);
            var count = counts[i][1];
            binInfo.push([val, val, count]);
        }

        var ret = {"dArr":dArr, "binInfo":binInfo};
        return ret;
    }

    // XX is this really faster than a simple iteration?
    function smoolakBS_left(arr, find) {
	// binary search, returns index left of insert point
	// based on https://jsperf.com/binary-search-in-javascript/7
        // insert_left based on https://rosettacode.org/wiki/Binary_search
	var lo = 0;
	var hi = arr.length - 1;
	var i;
	while(lo <= hi) {
	    i = ((lo + hi) >> 1);
	    if(arr[i] >= find) 
		hi = i - 1;
            else
	    //if (arr[i] < find) 
		lo = i + 1;
	    //else 
		//{ return lo; }
	}
	return lo;
    }

    function findBins(numVals, breakVals) {
    /*
    find the bin index for the break defined by breakVals for every value in numVals.
    The comparison uses "<=" ("left"). The first bin is therefore just the value of
    the first break, which makes sense for the most common case, 0, going into
    a special bin.  (for speed, 0 is hardcoded to always go into
    bin0). The caller can decrease break0 for more natural results in cases
    where special treatment of bin0 is not intended, e.g. when 0 does not appear.
    Also returns an array with the count for every bin.
    */
	
        var dArr = new Uint8Array(numVals.length); 
        var binCounts = new Uint32Array(breakVals.length);

        for (var i=0; i<numVals.length; i++) {
            var val = numVals[i];
            var binIdx = 0;
            if (val!==0) // special case for 0, saves some time
                binIdx = smoolakBS_left(breakVals, numVals[i]);
            dArr[i] = binIdx;
            binCounts[binIdx]++;
        }
        return {"dArr":dArr, "binCounts":binCounts};
    }

    function discretizeArray(arr, maxBinCount) {
        /* discretize numeric values to deciles. return an obj with dArr and binInfo */
        /* ported from Python cbAdd:discretizeArray */
        var counts = countAndSort(arr);

        // if we have just a few values, do not do any binning, just count
        if (counts.length < maxBinCount)
            return arrToEnum(arr, counts);

        // make array of count-indices of the breaks
        // e.g. if maxBinCount=10, breakPercs is [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        var breakPercs = [];
        var binPerc = 1.0/maxBinCount;
        for (var i=0; i<maxBinCount; i++)
            breakPercs.push( binPerc*i );

        var countLen = counts.length;
        var breakIndices = [];
        for (var i=0; i < breakPercs.length; i++) {
            var bp = breakPercs[i];
            breakIndices.push( Math.round(bp*countLen) );
        }
        breakIndices.push(countLen-1); // the last break is always a special case

        // make array with values at the break indices
        var breakValues = [];
        for (var i=0; i < breakIndices.length; i++) {
            var breakIdx = breakIndices[i];
            var breakVal = counts[breakIdx][0];
            // for expression vectors without a 0,
            // we don't want special handling of the value 0:
            // force break0 to be 0, so bin0 covers 0-break1
            if (i===0 && breakVal!==0)
                breakVal = 0;
            breakValues.push( breakVal );
        }

        var fb = findBins(arr, breakValues);
        var dArr = fb.dArr;
        var binCounts = fb.binCounts;

        var binInfo = [];
        for (var i=0; i<binCounts.length; i++) {
            var binMax = parseFloat(breakValues[i]);
            var binMin = 0.0;
            if (i===0)
                binMin = binMax;
            else
                binMin = parseFloat(breakValues[i-1]);

            var binCount = binCounts[i];
            binInfo.push( [binMin, binMax, binCount] );
        }

        return {"dArr": dArr, "binInfo": binInfo};
    }

    this.loadExprVec = function(geneSym, onDone, onProgress, binCount) {
    /* given a geneSym (string), retrieve array of deciles and call onDone with the
     * array */
        function onGeneDone(comprData, geneSym) {
            // decompress data and run onDone when ready
            self.exprCache[geneSym] = comprData;

            console.log("Got expression data, size = "+comprData.length+" bytes");
            var buf = pako.inflate(comprData);

            // see python code in cbAdd, function 'exprRowEncode':
            //# The format of a record is:
            //# - 2 bytes: length of descStr, e.g. gene identifier or else
            //# - len(descStr) bytes: the descriptive string descStr
            //# - 132 bytes: 11 deciles, encoded as 11 * 3 floats (=min, max, count)
            //# - array of n bytes, n = number of cells
            
            // read the gene description
            var descLen = cbUtil.baReadUint16(buf, 0);
            var arr = buf.slice(2, 2+descLen);
            var geneDesc = String.fromCharCode.apply(null, arr);

            // read the info about the bins
            //var binInfoLen = 11*3*4;
            //var dv = new DataView(buf.buffer, 2+descLen, binInfoLen);
            //var binInfo = []
            //for (var i=0; i < 11; i++) {
                //var startOfs = i*12;
                //var min = dv.getFloat32(startOfs, true); // true = little-endian
                //var max = dv.getFloat32(startOfs+4, true);
                //var binCount = dv.getFloat32(startOfs+8, true); // number of cells in this bin
                //if (min===0 && max===0 && binCount===0)
                    //break;
                //binInfo.push( [min, max, binCount] );
            //}

            // read the byte array with one bin index per cell
            // var digArr = buf.slice(2+descLen+binInfoLen);
            
            // read the expression array
            var sampleCount = self.conf.sampleCount;
            var matrixType = self.conf.matrixArrType;
            if (matrixType===undefined)
                alert("dataset JSON config file: missing matrixArrType attribute");
            var ArrType = cbUtil.makeType(matrixType);
            var arrData = buf.slice(2+descLen, 2+descLen+(4*sampleCount));
            var exprArr = new ArrType(arrData.buffer);

            console.time("discretize");
            var da = discretizeArray(exprArr, binCount);
            console.timeEnd("discretize");

            onDone(exprArr, da.dArr, geneSym, geneDesc, da.binInfo);
        }

        var offsData = self.geneOffsets[geneSym];
        if (offsData===undefined) {
            alert("cbData.js: "+geneSym+" is not in the expression matrix");
            onDone(null);
        }

        var start = offsData[0];
        var lineLen = offsData[1];

        var end = start + lineLen - 1; // end pos is inclusive

        var url = cbUtil.joinPaths([self.url, "exprMatrix.bin"]);

        if (geneSym in this.exprCache)
            onGeneDone(this.exprCache[geneSym], geneSym);
        else
            cbUtil.loadFile(url+"?"+geneSym, Uint8Array, onGeneDone, onProgress, geneSym, 
                start, end);
    };

    this.loadClusterMarkers = function(markerIndex, clusterName, onDone, onProgress) {
    /* given the name of a cluster, return an array of rows with the cluster-specific genes */
        var url = cbUtil.joinPaths([self.url, "markers", "markers_"+markerIndex, clusterName.replace("/", "_")+".tsv"]);
        cbUtil.loadTsvFile(url, onMarkersDone, {"clusterName":clusterName});

        function onMarkersDone(papaResults, url, otherData) {
            var rows = papaResults.data;
            onDone(rows, otherData);
        }
    };

    this.getMetaFields = function() {
    /* return an array of the meta fields, in the format of the config file:
     * objects with 'name', 'label', 'valCounts', etc */
        return self.conf.metaFields;
    };

    this.getConf = function() {
    /* return an object with a few general settings for the viewer:
     * - alpha: default transparency
     * - radius: circle default radius  */
        return self.conf;
    };

    this.getGenes = function() {
    /* return an object with the geneSymbols */
        return self.geneOffsets;
    };

    this.searchGenes = function(prefix, onDone) {
    /* call onDone with an array of gene symbols that start with prefix (case-ins.)
     * returns an array of objects with .id and .text attributes  */
        var geneList = [];
        for (var geneSym in self.geneOffsets) {
            if (geneSym.toLowerCase().startsWith(prefix))
                geneList.push({"id":geneSym, "text":geneSym});
        }
        onDone(geneList);
    };

    this.getName = function() {
    /* return name of current dataset*/
        if (self.conf!==null)
            return self.conf.name;
        else
            return self.name;
    };

    this.getDefaultClusterFieldIndex = function() {
    /* return field index of default cluster field */
        var idx = cbUtil.findIdxWhereEq(self.conf.metaFields, "name", self.conf.clusterField);
        return idx;
    };

}
