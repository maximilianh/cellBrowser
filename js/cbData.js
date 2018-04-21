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
        if (this.status !== 200) {
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
            var expLength = this._end-this._start;
            if (binData.byteLength !== expLength) {
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

    self.url = url;
    self._coordsDoneFunc = null;
    self.geneOffsets = null;

    this.conf = null;

    this.init = function(onDone) {
    /* load config from URL and call func when done */
        // load config and call onDone
        var dsUrl = cbUtil.joinPaths([this.url, "dataset.json"]);
        cbUtil.loadJson(dsUrl, function(data) { self.conf = data; });

        // start loading gene offsets, this takes a while
        var osUrl = cbUtil.joinPaths([this.url, "exprMatrixOffsets.json"]);
        cbUtil.loadJson(osUrl, function(data) { self.geneOffsets = data; onDone();});
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
        /* return a pair: [0] is "meta" or "gene" and [1] is the field name */
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
    /* for a given cell, call onDone with a list of the metadata values, as strings. */
        // first we need to lookup the offset of the line and its length from the index
        var url = cbUtil.joinPaths([self.url, "meta.index"]);
        var start = (cellIdx*6); // four bytes for the offset + 2 bytes for the line length
        var end   = start+6;

        function offsetDone(arr) {
            /* called when the offset in meta.index has been read */
            var offset = cbUtil.baReadOffset(arr, 0);
            var lineLen = cbUtil.baReadUint16(arr, 4);
            // now get the line from the .tsv file
            var url = cbUtil.joinPaths([self.url, "meta.tsv"]);
            cbUtil.loadFile(url+"?"+cellIdx, null, onDone, onProgress, null, 
                offset, offset+lineLen);
        }

        // chrome caching sometimes fails with byte range requests, so add cellidx
        cbUtil.loadFile(url+"?"+cellIdx, Uint8Array, function(text) { onDone(text); }, 
            undefined, null, start, end);
    };

    this.loadExprVec = function(geneId, onDone, onProgress) {
    /* given a geneId (string), retrieve array of deciles and call onDone with the
     * array */
        function onGeneDone(comprData, metaData) {
            // decompress data and run onDone on it
            var uncomprData = pako.inflate(comprData);
            onDone(uncomprData, metaData);
        }

        var offsData = self.geneOffsets[geneId];
        if (offsData===undefined)
            onDone(null);

        var start = offsData[0];
        var lineLen = offsData[1];
        var deciles = offsData[2];

        var end = start + lineLen;

        var url = cbUtil.joinPaths([self.url, "exprMatrix.bin"]);
        //cbUtil.loadFile(url+"?"+geneId, null, onLineDone, onProgress, {'gene':geneId}, 
         //   start, end);
        cbUtil.loadFile(url+"?"+geneId, Uint8Array, onGeneDone, onProgress, {'gene':geneId, 'deciles' : deciles}, 
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

    this.getDefaultClusterFieldIndex = function() {
    /* return field index of default cluster field */
        var idx = cbUtil.findIdxWhereEq(self.conf.metaFields, "name", self.conf.clusterField);
        return idx;
    };


}
