Combine results
---------------

You can use `cbTool metaCat` to merge the meta.tsv files from different pipelines (cbScanpy, cbSeurat, etc) into a single one, like this::

    cbTool metaCat myMeta.tsv seuratOut/meta.tsv scanpyOut/meta.tsv ./newMeta.tsv --fixDot

The option --fixDot will work around R's strange habit of replacing special characters in the cell identifiers with ".".
Directories created with ExportToCellbrowser() from R should not have this problem, but others may.

You can start with one of the auto-generated cellbrowser.conf files or start from a fresh one with `cbBuild --init`.
In this cellbrowser.conf, add all the coordinates files from all your pipelines. 
