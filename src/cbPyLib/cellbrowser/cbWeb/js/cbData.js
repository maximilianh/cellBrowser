// scDb: a class accessing single cell data from a URL
'use strict';
/*jshint globalstrict: true*/
/* jshint -W104 */  // allow some es6 parts (const) 
/* jshint -W117 */  // ignore undefined classes

/* a module with some helper functions */
var cbUtil = (function () {
    var my = {};

    // the byte range warning message will be shown only once, so we need a global flag
    my.byteRangeWarningShown = false;

    my.joinPaths = function joinPaths(parts, separator) {
    // join paths together, taking care of duplicated /s
      if (parts[0]==="") // ["", "test.txt] should just be test.txt, not /test.txt
          parts.shift();
      return parts.map(function(part) { return part.trim().replace(/(^[\/]*|[\/]*$)/g, ''); }).join(separator || '/');
    };

    my.dumpObj = function (o) {
    /* for debugging */
        console.log(JSON.stringify(o));
    };

    my.keys = function(o, isInt) {
    /* return all keys of object as an array */
        var allKeys = [];
        if (isInt)
            for(var k in o) allKeys.push(parseInt(k));
        else
            for(var j in o) allKeys.push(j);
        return allKeys;
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
                    if (url.search("dataset.json")>-1)
                        alert("Could not find a dataset at "+url+". If you are sure that the link is correct, please "+
                            "contact the administrator of this server, "+
                            "or cells@ucsc.edu if this is running at UCSC. ");
                    else
                        alert("Could not load "+url); 
                onSuccess(null);
            }
        });
    };
    
    my.loadFile = function(url, arrType, onDone, onProgress, otherInfo, start, end) {
        /* load text or binary file with HTTP GET into fileData variable and call function
         * onDone(fileData, otherInfo) when done. 
         * convert the data to arrType. arrType can be either a class like
         * Uint8Array or the value 'string', for normal text or 'comprText'
         * for gzip'ed string.  To switch off type casting, set arrType=null
         *
         * optional: byte range request for start-end (0-based, inclusive).
         * */
        var oReq = new XMLHttpRequest();
        oReq.open("GET", url, true);

        if (arrType==="string")
            oReq.responseType = "text";
        else
            oReq.responseType = "arraybuffer";

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
        /* called when binary file has been loaded. */
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

        // if the user wants only a byte range...
        if (this._start!==undefined) {
            // check if the expected length is OK. Some webservers don't support byte range 
            // requests.
            var expLength = this._end-this._start+1; // byte range is inclusive

            var dataLen = binData.byteLength;
            if (!dataLen)
                dataLen = binData.length;

            if (dataLen < expLength)
                alert("internal error cbData.js: chunk is too small. Does the HTTP server really support byte range requests?");

            if (dataLen > expLength) {
                console.log("Webserver does not support byte range requests, working around it, but this may be slow");
                if (dataLen>30000000 && !my.byteRangeWarningShown) {
                    alert("The webserver of this site does not support byte-range requests. " +
                        "While the cell browser may work to some extent, it will " +
                        " be slower and use more memory than normal. Please contact the administrator who setup this cell browser.");
                    my.byteRangeWarningShown = true;
                    }

                if (this._arrType) {
                    if (this._arrType==="string")
                        // for strings, it's easy to grab a part of it
                        binData = binData.substring(this._start, this._end);
                    else {
                        // it's a bit harder if the data is a buffer
                        binData = new Uint8Array(binData); // buffers don't support slicing, so cast to array first
                        binData = binData.slice(this._start, this._end).buffer; // slice does not include the end
                    }
                }
            }
        }

        //if (url.endsWith(".gz"))
            //binData = pako.ungzip(binData);

        if (this._arrType) {
            if (this._arrType==="comprText") {
                var arr = pako.ungzip(binData);
                // https://stackoverflow.com/questions/6965107/converting-between-strings-and-arraybuffers
                // convert byte array to strin
                // I used the apply() function originally, that's more compatible, but leads to a 
                // stack overflow in bigger datasets
                //binData = String.fromCharCode.apply(null, arr);
                var dec = new TextDecoder("utf-8");
                binData = dec.decode(arr);
            }
            else if (this._arrType!=="string")
                binData = new this._arrType(binData);
        }

        this._onDone(binData, this._otherInfo);
    };

    my.makeType = function(typeStr) {
        /* given a string, return the correct type array for it */
        typeStr = typeStr.toLowerCase();
        if ((typeStr==="double" || typeStr==="float64"))
            return Float64Array;
        if ((typeStr==="float" || typeStr==="float32"))
            return Float32Array;
        else if (typeStr==="int32")
            return Int32Array;
        else if (typeStr==="uint32" || typeStr==="dword")
            return Uint32Array;
        else if (typeStr==="uint16" || typeStr==="word")
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
    self.exprBinCount = 10; 

    self.exprCache = {}; // cached compressed expression arrays
    self.metaCache = {}; // cached compressed meta arrays

    self.quickExpr = {}; // uncompressed expression arrays
    self.allMeta = {};   // uncompressed meta arrays

    // special values representing NaN in data arrays, must match same variables in cellBrowser.py
    var FLOATNAN = Number.NEGATIVE_INFINITY; // NaN and sorting does not work. we want NaN always to be first, so encode as -inf

    this.conf = null;

    this.loadConfig = function(onDone, md5) {
    /* load config and gene offsets from URL and call func when done */

        var doneCount = 0;
        function gotOneFile() {
            doneCount++;
            if (doneCount===2)
                onDone(self.name);
        }

        // load config and call onDone
        var dsUrl = cbUtil.joinPaths([this.url, "dataset.json"]);
        // deactivate the cache - this is a small file that users typically change often
        if (!md5)
           dsUrl = dsUrl+"?"+Math.floor(Math.random()*100000000);
        else
           dsUrl = dsUrl+"?"+md5;
        cbUtil.loadJson(dsUrl, function(data) { self.conf = data; gotOneFile();});

        // start loading gene offsets in background, this takes a while
        var osUrl = cbUtil.joinPaths([this.url, "exprMatrix.json"]);
        cbUtil.loadJson(osUrl, function(data) { self.geneOffsets = data; gotOneFile();}, true);
    };

    this.loadCoords = function(coordIdx, onDone, onProgress) {
    /* load coordinates from URL and call onDone(array of (x,y), coordInfoObj, labelMids) when done */
        //var coordInfo = cbUtil.findObjWhereEq(self.conf.coords, "name", this.coordName);
        var i = 0;
        var binData = null;
        var meta = null;
        var labelMids; // null means: no json file found

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

        if (self.conf.coords.length===0 || self.conf.coords===undefined) {
            alert("There are no coordinates defined in the dataset.json config file. Please add at least one coordinates file to cellbrowser.conf and run cbBuild again with the same html output directory.")
            return;
        }

        var coordInfo = self.conf.coords[coordIdx];
        if (coordInfo===null) {
           alert("Could not find coordinates with name "+this.coordName);
           return;
        }

        var binUrl = cbUtil.joinPaths([self.url, "coords", coordInfo.name, "coords.bin"]);
        if (coordInfo.md5)
            binUrl += "?"+coordInfo.md5
        var arrType = cbUtil.makeType(coordInfo.type || "float64");
        cbUtil.loadFile(binUrl, arrType, binDone, onProgress, coordInfo);

        var jsonUrl = cbUtil.joinPaths([self.url, "coords", coordInfo.name, "clusterLabels.json"]);
        if (coordInfo.labelMd5)
            jsonUrl += "?"+coordInfo.labelMd5
        cbUtil.loadJson(jsonUrl, jsonDone, true);
    };

    this.getDefaultColorField = function() {
        /* return a pair: [0] is "meta" or "gene" and [1] is the field name or the gene */
        if (self.conf.clusterField)
            return ["meta", self.conf.clusterField]; // just for backwards-compat. with old datasets
        else
            return ["meta", self.conf.defColorField];
    };

    this.findMetaInfo= function(findName) {
        /* return meta info field with name, add 'metaIndex' attribute  */
        var metaFieldInfo = self.conf.metaFields;
        for (var i = 0; i < metaFieldInfo.length; i++) {
            var metaInfo = metaFieldInfo[i];
            if (metaInfo.name===findName || metaInfo.label===findName) {
                metaInfo.index = i;
                return metaInfo;
            }
        }
        return null;
    }
    
    this.fieldNameToIndex = function(fieldName) {
        /* given a meta field name, return its meta table index or null */
        var idx = cbUtil.findIdxWhereEq(self.conf.metaFields, "name", fieldName);
        if (idx===null)
            idx = cbUtil.findIdxWhereEq(self.conf.metaFields, "label", fieldName);
        return idx;
    };

    this._startMetaLoad = function(metaInfo, arrType, onMetaDone, onProgress, extraInfo) {
        /* start the loading of a meta data file and call onMetaDone when ready */
        var fieldName = metaInfo.name;
        var binUrl = cbUtil.joinPaths([self.url, "metaFields", fieldName+".bin.gz"]);
        if (metaInfo.md5)
            binUrl += "?"+metaInfo.md5;
        cbUtil.loadFile(binUrl, arrType, 
                onMetaDone,
                onProgress, metaInfo);
    }

    this.loadMetaVec = function(metaInfo, onDone, onProgress, otherInfo) {
    /* get an array of numbers, one per cell, that reflect the meta field contents
     * and an object with some info about the field. call onDone(arr, metaInfo) when done. 
     * Keep all compressed arrays in metaCache;
     * Numerical meta data (int/float) is discretized, binning info is written to metaInfo.binInfo
     * and the original data is added to metaInfo.origVals
     * */

        function onMetaDone(comprBytes, metaInfo) {
            self.metaCache[metaInfo.name] = comprBytes; 
            var ArrType = cbUtil.makeType(metaInfo.arrType);

            var bytes = comprBytes; // some Apache/InternetBrowser combinations silently uncompress
            var comprView = new Uint8Array(comprBytes);
            // magic bytes of gzip are 1f, 8b, and we assume little endianess
            if (comprView[0]===0x1f && comprView[1]===0x8b) {
                try {
                    bytes = pako.ungzip(comprBytes);
                }
                catch(err) {
                    alert("Error when decompressing a file. This has to do with your Apache config or your "+
                             " internet browser. Please contact cells@ucsc.edu, we can help you solve this. "+err);
                }
            }

            var buffer = bytes.buffer;
            var arr = new ArrType(buffer);
            if (metaInfo.arrType==="float32") {
                // numeric arrays have to be binned on the client. They are always floats.
                var discRes = discretizeArray(arr, self.exprBinCount, FLOATNAN);
                metaInfo.origVals = arr; // keep original values, so we can later query for them
                arr = discRes.dArr;
                metaInfo.binInfo = discRes.binInfo;
                //metaInfo.valCounts = binInfoToValCounts(discRes.binInfo);
            }
            onDone(arr, metaInfo, otherInfo);
        }

        //var metaInfo = self.conf.metaFields[fieldIdx];
        //console.log(metaInfo);

        if ((self.allMeta!==undefined) && (metaInfo.name in self.allMeta)) {
            console.log("Found in uncompressed cache");
            onDone(self.allMeta[metaInfo.name], metaInfo, otherInfo);
            return;
        }
        else if ((self.metaCache!==undefined) && (metaInfo.name in self.metaCache)) {
            console.log("Found in compressed cache");
            onMetaDone(self.metaCache[metaInfo.name], metaInfo);
            return;
        }
        else {
            self._startMetaLoad(metaInfo, null, onMetaDone, onProgress, metaInfo);
            return;
        }
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
            cbUtil.loadFile(url+"?"+cellIdx, "string", lineDone, onProgress, null, 
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

    function arrToEnum(arr, counts) {
        /* given an array of numbers, count how often each number appears.
         * return an obj with two keys:
         * - dArr is the new array with the index of each value
         * - binInfo is an array with for each index a tuple of (value, value, count) 
         */
        // replace values in array with their enum-index
        // -> make a mapping value -> bin
        var valToBin = {};
        sortArrOfArr(counts, 0); // sort by value
        for (var i=0; i<counts.length; i++) {
            var val = counts[i][0];
            valToBin[val] = i;
        }

        // apply the mapping
        var dArr = new Uint8Array(arr.length);
        for (let i=0; i<arr.length; i++) {
            dArr[i] = valToBin[arr[i]];
        }

        var binInfo = [];
        for (let i=0; i < counts.length; i++) {
            let val = parseFloat(counts[i][0]);
            var count = counts[i][1];
            binInfo.push([val, val, count]);
        }

        var ret = {"dArr":dArr, "binInfo":binInfo};
        return ret;
    }

    function findBins(numVals, bin0Val, breakVals) {
    /*
    Discretize an array of expression values.

    Find the bin index for the break defined by breakVals for every value in numVals.
    The comparison uses "<=" ("left"). Also returns an array with the count for every bin.
    */
	
        var dArr = new Uint8Array(numVals.length); 
        var binCounts = new Uint32Array(breakVals.length+1); // bin0 is a special bin, so +1

        for (var i=0; i<numVals.length; i++) {
            var val = numVals[i];
            var binIdx = 0;
            if (val!==bin0Val)
                binIdx = smoolakBS_left(breakVals, numVals[i])+1; // bin0 is reserved, so +1
            dArr[i] = binIdx;
            binCounts[binIdx]++;
        }
        return {"dArr":dArr, "binCounts":binCounts};
    }

    function findBins_meta(numVals, breakVals, nanValue) {
    /*
    Discretize an array of any numbers. See findBins_expr
    Put into bins an array of numbers, some of which may be NaN, given the breaks.
    bin0 is reserved for NaN.
    NaN is encoded in numVals as nanValues, which has to smaller than any other value.

    Find the bin index for the break defined by breakVals for every value in numVals.
    The comparison uses "<=" ("left"). 
    Also returns an array with the count for every bin.

    */
	
        var dArr = new Uint8Array(numVals.length); 
        var binCounts = new Uint32Array(breakVals.length-1);
        var breaks = breakVals.slice(2);

        for (var i=0; i<numVals.length; i++) {
            var val = numVals[i];
            var binIdx = 0;
            if (val!==nanValue)
                binIdx = smoolakBS_left(breaks, numVals[i]);
            dArr[i] = binIdx+1;
            binCounts[binIdx]++;
        }
        return {"dArr":dArr, "binCounts":binCounts};
    }

    function discretizeArray(arr, maxBinCount, bin0Val) {
        /* discretize numeric values to deciles. return an obj with dArr and binInfo */
        /* bin0Val is the value that is treated differently, it is kept in its own bin */
        /* is bin0Val is null, switch off special bin0Value handling
        /* ported from Python cbAdd:discretizeArray */
        /* supports NaN special values */
        var counts = countAndSort(arr);

        // if we have just a few values, do not do any binning, just count
        if (counts.length < maxBinCount)
            return arrToEnum(arr, counts);

        // From now on, we treat NAs/0s separately, so remove their counts
        if (counts[0][0]===bin0Val) // NA is always first, as we defined it as -inf
            counts.shift(); // remove first element
        
        // make array of count-indices of the breaks
        // e.g. if maxBinCount=10, breakPercs is [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        var breakPercs = [];
        var binPerc = 1.0/maxBinCount;
        for (var i=0; i<maxBinCount; i++)
            breakPercs.push( binPerc*i );

        // for each percentage, get the index
        var countLen = counts.length;
        var breakIndices = [];
        for (let i=0; i < breakPercs.length; i++) {
            var bp = breakPercs[i];
            breakIndices.push( Math.round(bp*countLen) );
        }
        breakIndices.push(countLen-1); // the last break is always a special case. Here we set it to the last element.
        // now we have 11 breaks for 10 bins

        var minVal = counts[breakIndices[0]][0];

        // make array with values at the break indices except the first one, as comparison is <=
        var breakValues = [];
        for (let i=1; i < breakIndices.length; i++) {
            var breakIdx = breakIndices[i];
            var breakVal = counts[breakIdx][0];
            breakValues.push( breakVal );
        }
        
        var fb = findBins(arr, bin0Val, breakValues);
        var dArr = fb.dArr;
        var binCounts = fb.binCounts;

        // we should have 11 breaks/10bins but 12 values in binCounts
        //assert(len(breakVals)==12));
        //assert(len(binCounts)==11));
        //assert((len(binCounts)+1 == len(breakVals)));

        // convert to format (min, max, count)
        // bin0 has special min/max of "Unknown" or 0
        var binInfo = [];

        var bin0MinMax = "Unknown";
        if (bin0Val===0)
            bin0MinMax = 0;
        binInfo.push( [bin0MinMax, bin0MinMax, binCounts[0]] );

        for (let i=0; i<breakValues.length; i++) {
            var binMin = minVal;
            if (i!==0)
                binMin = parseFloat(breakValues[i-1]);
            var binMax = parseFloat(breakValues[i]);
            var binCount = binCounts[i+1];
            binInfo.push( [binMin, binMax, binCount] );
        }

        return {"dArr": dArr, "binInfo": binInfo};
    }

    this.loadExprAndDiscretize = function(geneSym, onDone, onProgress) {
    /* given a geneSym (string), retrieve array of array put into binCount bins
     * and call onDone with (array, discretizedArray, geneSymbol, geneDesc,
     * binInfo) */

        var binCount = self.exprBinCount;

        function onLoadedVec(exprArr, geneSym, geneDesc) {
            console.time("discretize "+geneSym);
            var specVal = 0;
            var matrixMin = self.getMatrixMin();
            if (matrixMin < 0)
                specVal = null;

            var da = discretizeArray(exprArr, binCount, specVal);
            console.timeEnd("discretize "+geneSym);
            onDone(exprArr, da.dArr, geneSym, geneDesc, da.binInfo);
        }

        this.loadExprVec(geneSym, onLoadedVec, onProgress)
    };

    this.loadExprVec = function(geneSym, onDone, onProgress, otherInfo) {
    /* given a geneSym (string), retrieve array of values and call onDone with
     * (array, geneSym, geneDesc) */
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

            // read the expression array
            var sampleCount = self.conf.sampleCount;
            var matrixType = self.conf.matrixArrType;
            if (matrixType===undefined)
                alert("dataset JSON config file: missing matrixArrType attribute");
            var ArrType = cbUtil.makeType(matrixType);
            var arrData = buf.slice(2+descLen, 2+descLen+(4*sampleCount));
            var exprArr = new ArrType(arrData.buffer);

            onDone(exprArr, geneSym, geneDesc, otherInfo);
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

    this.addCustomMetaField = function(metaInfo) {
        //if (!self.customMeta)
            //self.customMeta = [];
        //self.customMeta.push(metaInfo); 
        metaInfo.isCustom = true;
        self.conf.metaFields.unshift(metaInfo);
    }

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

    this.loadCellIds = function(idxArray, onDone, onProgress) {
        /* Get the cellId strings, the first meta field, for the integer IDs of cells in idxArray.
         * Calls onDone with an array of strings.
         * If idxArray is null, gets all cellIds */

        function mapIdxToId(idxArray) {
            /* Internally cells are referred to by an array of cell indices.
               Translate an array of these numbers to an array of strings */

            if (idxArray===null)
                return self.cellIds;

            var idList = [];
            for (var i=0; i<idxArray.length; i++) {
                var cellIdx = idxArray[i];  
                idList.push(self.cellIds[cellIdx]);
            }
            return idList;
        }

        function onIdsDone (text) {
            /* called when the cellIds are loaded. comprText is the whole gzip'ed text file, one ID per line. */
            self.cellIds = text.split("\n");
            onDone(mapIdxToId(idxArray));
        }

        if (self.cellIds===undefined) {
            // if we haven't loaded them yet, trigger the cellId load
            let cellIdMeta = self.getCellIdMeta();
            self._startMetaLoad(cellIdMeta, "comprText", onIdsDone, onProgress, idxArray)
        }
        else
            onDone(mapIdxToId(idxArray));
    }

    this.loadFindCellIds = function(findIds, onSearchDone, onProgress, hasWildcards) {
        /* load the all cell IDs then convert a given array of IDs to the matching cell indexes.
         * Returns a pair of (matching cell indices, array of not found identifiers) */
        function mapIdsToIdx(findIds) {
            var cellIds = self.cellIds;

            var notFoundIds = [];
            var idArr = [];

            if (hasWildcards) {
                var foundIds = {}; // = set = removes any duplicates
                for (var i=0; i<findIds.length; i++) {
                    var searchId = findIds[i];  
                    var searchRe = new RegExp(searchId);
                    var foundOne = false;
                    for (var cellIdx = 0; cellIdx<cellIds.length; cellIdx++) {
                        var cellId = cellIds[cellIdx];
                        if (searchRe.exec(cellId)!==null) {
                            foundIds[cellIdx] = null;
                            foundOne = true;
                        }
                    }
                    if (!foundOne)
                        notFoundIds.push(searchId);
                }
                idArr = cbUtil.keys(foundIds); // convert to array
                idArr.sort();
            }
            else {
                for (let i=0; i<findIds.length; i++) {
                    let searchId = findIds[i];  
                    var foundIdx = cellIds.indexOf(searchId);
                    if (foundIdx===-1)
                        notFoundIds.push(searchId);
                    else
                        idArr.push(foundIdx);

                }
            }
            return [idArr, notFoundIds];
        }

        function onIdsDone (text) {
            /* called when the cellIds are loaded. comprText is the whole gzip'ed text file, one ID per line. */
            self.cellIds = text.split("\n");
            onSearchDone(mapIdsToIdx(findIds));
        }

        if (self.cellIds===undefined) 
            // trigger the cellId load
            self._startMetaLoad(self.getMetaFields()[0], "comprText", onIdsDone, onProgress)
        else
            onSearchDone(mapIdsToIdx(findIds));
    }

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
            cbUtil.loadFile(url+"?"+cellIdx, "string", lineDone, onProgress, null, 
                offset, offset+lineLen);
        }

        // chrome caching sometimes fails with byte range requests, so add cellidx to the URL
        cbUtil.loadFile(url+"?"+cellIdx, Uint8Array, function(byteArr) { offsetDone(byteArr); }, 
            undefined, null, start, end);
    };

    function countAndSort(arr) {
        /* count values in array, return an array of [value, count]. */
        //var counts = {};
        var counts = new Map(); // this is a relatively recent Javascript feature, but is faster
        for (var i=0; i < arr.length; i++) {
            var num = arr[i];
            //without the Map(), it's a bit slower: counts[num] = counts[num] ? counts[num] + 1 : 1;
            counts.set(num, (counts.get(num) | 0) + 1);
        }
        var entries = Array.from(counts.entries());
        entries.sort(function(a,b) { return a[0]-b[0]});
        return entries;
    }

    // Todo: check one day if this really faster than a simple iteration
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

    this.loadExprVec = function(geneSym, onDone, onProgress, otherInfo) {
    /* given a geneSym (string), retrieve array of values and call onDone with
     * (array, geneSym, geneDesc) */
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

            // read the expression array
            var sampleCount = self.conf.sampleCount;
            var matrixType = self.conf.matrixArrType;
            if (matrixType===undefined)
                alert("dataset JSON config file: missing matrixArrType attribute");
            var ArrType = cbUtil.makeType(matrixType);
            var arrData = buf.slice(2+descLen, 2+descLen+(4*sampleCount));
            var exprArr = new ArrType(arrData.buffer);

            onDone(exprArr, geneSym, geneDesc, otherInfo);
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

    this.getCellIdMeta = function() {
        /* return the cell ID meta field */
        for (let metaInfo of self.conf.metaFields)
            if (!metaInfo.isCustom)
                return metaInfo;
    }

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

    this.preloadAllMeta = function() {
        /* start loading all meta value vectors and add them to db.allMeta. */
        self.allMeta = {};

        function doneMetaVec(arr, metaInfo, otherInfo) {
            self.allMeta[metaInfo.name] = arr; 
            delete self.metaCache[fieldIdx];
        }

        var metaFieldInfo = self.getMetaFields();
        for (var fieldIdx = 0; fieldIdx < metaFieldInfo.length; fieldIdx++) {
           var fieldInfo = metaFieldInfo[fieldIdx];
           if (fieldInfo.type==="uniqueString" || fieldInfo.arr)
               continue;
           self.loadMetaVec(fieldInfo, doneMetaVec);
        }
    }

    this.getMatrixMin = function() {
        /* return the minimum valu in the matrix */
       var validGenes = self.getGenes();
       var matrixMin = 0;
       if ("_range" in validGenes)
           matrixMin = validGenes["_range"][0];
        return matrixMin;
    }

    this.preloadGenes = function(geneSyms, onDone, onProgress) {
       /* start loading the gene expression vectors in the background. call onDone when done. */
       var validGenes = self.getGenes();

       var loadCounter = 0;
       if (geneSyms) {
           for (var i=0; i<geneSyms.length; i++) {
               var sym = geneSyms[i][0];
               if (! (sym in validGenes)) {
                  alert("Error: "+sym+" is in quick genes list but is not a valid gene");
                  continue;
               }

                self.loadExprAndDiscretize(
                   sym, 
                   function(exprVec, discExprVec, geneSym, geneDesc, binInfo) { 
                       self.quickExpr[geneSym] = [discExprVec, geneDesc, binInfo];
                       loadCounter++; 
                       if (loadCounter===geneSyms.length) onDone(); 
                   },
                   onProgress);
           }
       }
    };

    this.loadGeneSetExpr = function(onDone) {
        /* return array of [geneSym, discExprVec, geneDesc, binInfo] */
        var setInfo = [];

        for (var geneInfo of self.conf.quickGenes) {
            var geneSym = geneInfo[0];
            var exprInfo = self.quickExpr[geneSym]; // contains: [discExprVec, geneDesc, binInfo]
            var newInfo = [geneSym, exprInfo[0], exprInfo[1], exprInfo[2]];
            setInfo.push(newInfo);
        }
        onDone(setInfo);
    };

    this.removeAllCustomAnnots = function() {
        /* remove all custom annotation fields */
        var newMetaInfo = [];
        for (let metaField of self.conf.metaFields) {
            if (!metaField.isCustom)
                newMetaInfo.push(metaField);
        }
        self.conf.metaFields = newMetaInfo;
    }

}
