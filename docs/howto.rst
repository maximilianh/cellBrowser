How To...
____

How to export the tree and data from URD for use in the Cell Browser
^^^^

`URD <https://github.com/farrellja/URD>`_ is an R package that can be used to reconstruct transcriptional
trajectories and then displaying this trajectory as a branching tree. You can export the tree diagram, 
expression data, and metadata from an URD object from within R and then use the resulting files to 
build a cell browser. 

Step 1: Export cell coordinates for the tree
""""

First, we need the coordinates for the cells in relation to the tree:

``write.table(urd_obj@tree$cell.layout, file='urd.coords.tsv', quote=FALSE, sep='\t', col.names = NA)``

Step 2: Export line coordinates for the tree
""""

Next, we need the coordinates for the lines that make up the tree:

``write.table(urd_obj@tree$tree.layout, file='urd.lines.tsv', quote=FALSE, sep='\t', col.names = NA)``

Step 3: Export expression matrix
""""

Export data in MTX format, since it can handle large matrix sizes. MTX consists of three files: 
(1) a sparse matrix, (2) a file of column names, and (3) a file of row names.

Set file names for all three files:

``matrixPath<-file.path('matrix.mtx')``

``genesPath <- file.path("genes.tsv")``

``barcodesPath <- file.path("barcodes.tsv")``

Output to each of the files

(1) MTX sparse matrix

``writeMM(urd_obj@count.data, matrixPath)``

(2) Row names (genes)

``write.table(as.data.frame(cbind(rownames(urd_obj@count.data), 
rownames(urd_obj@count.data))), file=genesPath, sep="\t", 
row.names=F, col.names=F, quote=F)``

(3) Column names (samples/cells)

``write(colnames(urd_obj@count.data), file = barcodesPath)``

Step 4: Convert MTX to tsv.gz
""""

It's easiest to specify a single exprMatrix.tsv.gz file in your cellbrowser.conf later,
so we'll convert our exported MTX to tsv via ``cbTool mtx2tsv``:

``cbTool mtx2tsv matrix.mtx genes.tsv barcodes.tsv exprMatrix.tsv.gz``

Step 5: Export metadata
""""

Metadata annotations are also needed for a cell browser:

``write.table(urd_obj@meta, file='meta.tsv', 
quote=FALSE, sep='\t', col.names = NA)``

Step 6: Export tSNE (optional)
""""

The cell coordinates and lines from steps one and two above satisfy the cell browser's need for a layout, however, 
URD can generate a tSNE layout as part of it's run. You can export these coordinates
for use in the cell browser.

``write.table(urd_obj@tsne.y, file='tsne.coords.tsv', quote=FALSE, sep='\t', col.names = NA)``

Step 7: Create your cellbrowser.conf
""""

Next create the cellbrowser.conf file for your dataset. You can use ``cbBuild --init`` to
place an example cellbrowser.conf (and desc.conf) into your current directory.

You will specifically need to edit these lines to point to the flies that you exported in steps 1-5 above:

::


  exprMatrix="exprMatrix.tsv.gz"

  meta="meta.tsv"

  coords=[
    {
      "file":"urd.coords.tsv",
      "lineFile":"urd.lines.tsv",
      "shortLabel":"URD Trajectory",
      "flipY":True,
      "lineFlipY": True
    },
    {
      "file": "tsne.coords.tsv",
      "shortLabel":"tSNE"
    }
  ]
  
You will still need to set the other `required settings <https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf#L1>`_ in your cellbrowser.conf as well
