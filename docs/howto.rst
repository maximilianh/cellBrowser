How To...
____

How to create a cell browser using a Seurat RDS file
^^^^

You can go from an RDS file to cell browser for a dataset in three easy steps:

Step 1: Export an RDS file of your Seurat object
""""

From within R, run this command to create an RDS file fo your dataset::

  saveRDS(objName, "myDataset.rds")

Step 2: Use cbImportSeurat to export 
""""

Next, you will use cbImportSeurat to create the files needed for a cell browser using the data in the RDS file::

  cbImportSeurat -i myDataset.rds -o myRdsImport -n seurat-import

Note: cbImportSeurat will work with RDS files from Seurat v2 or v3. When importing data, you need to have installed the same version of Seurat that was used to create the RDS file.

Step 3: Build a Cell Browser
""""

Lastly, go into the output directory specified in the cbImportSeurat command and run cbBuild to create the cell browser::

  cd myRdsImport
  cbBuild -o ~/public_html/cb

Or, if you don't have a webserver already, use the built-in one:

  cbBuild -o /myHtmlFiles -p 8888

You should now be able to access your cell browser from the web.

How to use the cell browser export function in Seurat3
^^^^

It is a simple, single-line command to build a web-accessible cell browser from a Seurat object from within R:: 

 ExportToCellbrowser(myObj, dir="myDatasetExport", cb.dir="htdocs", dataset.name="pbmcSmall", port=8080)




How to run a basic Seurat pipeline using cbSeurat
^^^^

Going from an expression matrix to a cell browser by running our basic Seurat pipeline takes two steps:

Step 1: Run cbSeurat on your expression matrix
""""

First, run a Seurat pipeline on your expression matrix using ``cbSeurat``::

  cbSeurat --exprMatrix=myExpressionMatrix.tsv.gz --name=myDataset --outDir=seurat-out

Step 2: Build a Cell Browser
""""

Next, go into the output directory specified in the cbImportSeurat command and run cbBuild to create the cell browser::

  cd seurat-out
  cbBuild -o ~/public_html/cb

Or, if you don't have a webserver already, start the built-in one:

  cbBuild -o /myHtmlFiles -p 8888


How to configue a basic cbSeurat pipeline
^^^^

Running ``cbSeurat`` will run a basic Seurat pipeline with the default settings. ``cbSeurat`` can be configured through a `seurat.conf <https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/seurat.conf>`_.

Step 1: Copy a seurat.conf 
""""

cbSeurat can be used to copy down an example seurat.conf::
 
  cbSeurat --init

Step 2: Edit your seurat.conf
""""

Now that you have a seurat.conf in your current directory, open it up and edit it! If this file is in the same 
directory where you are running ``cbSeurat``, it will be automatically picked up. 



How to create a cell browser using a Scanpy h5ad file
^^^^

Going from an h5ad file to cell browser for a dataset takes two steps:

Step 1: Use cbImportScanpy to export 
""""

First, you will use cbImportScanpy to create the files needed for a cell browser using the data in the RDS file::

  cbImportScanpy -i myDataset.h5ad -o scanpy-import -n my-dataset

Step 2: Build a Cell Browser
""""

Then, go into the output directory specified in the cbImportSeurat command and run cbBuild to create the cell browser::

  cd scanpy-import
  cbBuild -o ~/public_html/cb

Or, if you don't have a webserver already, start the built-in one:

  cbBuild -o /myHtmlFiles -p 8888

You should now be able to access your cell browser from the web or your local computer.


How to convert a Scanpy object within Python
^^^^

It a few simple commands to build a ``cellbrowser.conf`` and all the files you need for a cell
browser. This is particularly useful for Jupyter notebooks. 


Step 1: Export the data needed
""""

Load the cell browser package and export the files from the scanpy object::

 import cellbrowser.cellbrowser as cb
 cb.scanpyToCellbrowser(adata, "scanpyOut", "myScanpyDataset")

Step 2: Build the cell browser
""""

Next, build the dataset::

  cb.build("scanpyOut", "~/public_html/cb")

Step 3: Start (and stop) web server (optional)
""""

This step is only necessary if you don't already have a web server running that is servering up the output of step 2.

Start the web server::

  cb.serve("~/public_html/cb", 8888)

Stop the webserver when you're done::

  cb.stop()


How to run a basic Scanpy pipeline using cbScanpy
^^^^

Going from an expression matrix to a cell browser by running our basic Scanpy pipeline takes two steps:

Step 1: Run cbScanpy on your expression matrix
""""

First, run a Scanpy pipeline on your expression matrix using cbSeurat::

  cbScanpy -e myExpressionMatrix.tsv.gz -n my-scanpy-dataset -o scanpy-out -m cell-annotations.tsv

Step 2: Build a Cell Browser
""""

Next, go into the output directory specified in the ``cbScanpy`` command and build your cell browser::

  cd scanpy-out
  cbBuild -o ~/public_html/cb

How to configue a basic cbScanpy pipeline
^^^^

Running ``cbSeurat`` will run a basic Scanpy pipeline with the default settings. ``cbScanpy`` can be configured through a `scanpy.conf <https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/scanpy.conf>`_.

Step 1: Copy a scanpy.conf 
""""

cbSeurat can be used to copy down an example scanpy.conf::
 
  cbScanpy --init

Step 2: Edit your seurat.conf
""""

Now that you have a scanpy.conf in your current directory, open it up and edit it! If this file is in the same 
directory where you are running ``cbScanpy``, it will be automatically picked up. 



How to export the tree and data from URD for use in the Cell Browser
^^^^

`URD <https://github.com/farrellja/URD>`_ is an R package that can be used to reconstruct transcriptional
trajectories and then displaying this trajectory as a branching tree. You can export the tree diagram, 
expression data, and metadata from an URD object from within R and then use the resulting files to 
build a cell browser. 

Step 1: Export cell coordinates for the tree
""""

First, we need the coordinates for the cells in relation to the tree::

  write.table(urd_obj@tree$cell.layout, file='urd.coords.tsv', quote=FALSE, sep='\t', col.names = NA)

Step 2: Export line coordinates for the tree
""""

Next, we need the coordinates for the lines that make up the tree::

  write.table(urd_obj@tree$tree.layout, file='urd.lines.tsv', quote=FALSE, sep='\t', col.names = NA)

Step 3: Export expression matrix
""""

Export data in MTX format, since it can handle large matrix sizes. MTX consists of three files: 
(1) a sparse matrix, (2) a file of column names, and (3) a file of row names.

Set file names for all three files::

  matrixPath<-file.path('matrix.mtx')
  genesPath <- file.path("genes.tsv")
  barcodesPath <- file.path("barcodes.tsv")

Output to each of the files

(1) MTX sparse matrix:

::

  writeMM(urd_obj@count.data, matrixPath)``

(2) Row names (genes):

::

  write.table(as.data.frame(cbind(rownames(urd_obj@count.data), rownames(urd_obj@count.data))), file=genesPath, sep="\t", row.names=F, col.names=F, quote=F)

(3) Column names (samples/cells):

::

  write(colnames(urd_obj@count.data), file = barcodesPath)

Step 4: Convert MTX to tsv.gz
""""

It's easiest to specify a single exprMatrix.tsv.gz file in your cellbrowser.conf later,
so we'll convert our exported MTX to tsv via ``cbTool mtx2tsv``::

  cbTool mtx2tsv matrix.mtx genes.tsv barcodes.tsv exprMatrix.tsv.gz

Step 5: Export metadata
""""

Metadata annotations are also needed for a cell browser::

  write.table(urd_obj@meta, file='meta.tsv', quote=FALSE, sep='\t', col.names = NA)

Step 6: Export tSNE (optional)
""""

The cell coordinates and lines from steps one and two above satisfy the cell browser's need for a layout, however, 
URD can generate a tSNE layout as part of it's run. You can export these coordinates
for use in the cell browser::

  write.table(urd_obj@tsne.y, file='tsne.coords.tsv', quote=FALSE, sep='\t', col.names = NA)

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

How to start the webserver without building datasets
^^^^

If you have stopped the built-in webserver and want to start it again, without rebuilding the entire dataset, use the cbUpgrade tool:

  cbUpgrade -o /myHtmlFiles -p 8888

