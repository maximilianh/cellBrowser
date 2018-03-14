// scDb: a class accessing single cell data from a URL
'use strict';

function scDb(url) {
    this.url = url;

    function joinPaths(parts, separator) {
      return parts.map(function(part) { return part.trim().replace(/(^[\/]*|[\/]*$)/g, ''); }).join(separator || '/');
    }

    this.init = function(successFunc) {
    /* load config from URL and call func when done */
        var dsUrl = joinPaths[url, "dataset.json"];
        var jqxhr = $.getJSON( "example.json", successFunc).fail(function() { alert( "Could not load "+dsUrl ); });
    }

    this.loadCoords = function(coordName) {
    };
}
