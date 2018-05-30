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
         * onDone(fileData, otherInfo) when done, 
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
        /* called when binary file has been loaded */
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
        cbUtil.loadJson(dsUrl, function(data) { self.conf = data; gotOneFile()});

        // start loading gene offsets, this takes a while
        var osUrl = cbUtil.joinPaths([this.url, "exprMatrixOffsets.json"]);
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
        var fieldName = metaInfo.name;
        var binUrl = cbUtil.joinPaths([self.url, "metaFields", fieldName+".bin"]);
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

    this.loadExprVec = function(geneSym, onDone, onProgress) {
    /* given a geneSym (string), retrieve array of deciles and call onDone with the
     * array */
        function onGeneDone(comprData, geneSym) {
            // decompress data and run onDone when ready
            var buf = pako.inflate(comprData);

            // see python code in cbAdd, function 'exprRowEncode':
            //# The format of a record is:
            //# - 2 bytes: length of descStr, e.g. gene identifier or else
            //# - len(descStr) bytes: the descriptive string descStr
            //# - 132 bytes: 11 deciles, encoded as 11 * 3 floats (=min, max, count)
            //# - array of n bytes, n = number of cells
            
            // read the gene description
            var descLen = cbUtil.baReadUint16(buf, 0);
            var arr = buf.slice(2, 18);
            var geneDesc = String.fromCharCode.apply(null, arr);

            // read the info about the bins
            var binInfoLen = 11*3*4;
            var dv = new DataView(buf.buffer, 2+descLen, binInfoLen);
            var binInfo = []
            for (var i=0; i < 11; i++) {
                var startOfs = i*12;
                var min = dv.getFloat32(startOfs, true); // true = little-endian
                var max = dv.getFloat32(startOfs+4, true);
                var binCount = dv.getFloat32(startOfs+8, true); // number of cells in this bin
                if (min===0 && max===0 && binCount===0)
                    break;
                binInfo.push( [min, max, binCount] );
            }

            // read the byte array with one bin index per cell
            var digArr = buf.slice(2+descLen+binInfoLen);
             
            onDone(digArr, geneSym, geneDesc, binInfo);
        }

        var offsData = self.geneOffsets[geneSym];
        if (offsData===undefined)
            onDone(null);

        var start = offsData[0];
        var lineLen = offsData[1];

        var end = start + lineLen - 1; // end pos is inclusive

        var url = cbUtil.joinPaths([self.url, "exprMatrix.bin"]);
        //cbUtil.loadFile(url+"?"+geneSym, null, onLineDone, onProgress, {'gene':geneSym}, 
         //   start, end);
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
    }

    this.getName = function() {
    /* return name of current dataset*/
        if (self.conf!==null)
            return self.conf.name;
        else
            return self.name;
    }

    this.getDefaultClusterFieldIndex = function() {
    /* return field index of default cluster field */
        var idx = cbUtil.findIdxWhereEq(self.conf.metaFields, "name", self.conf.clusterField);
        return idx;
    };


}
