With Scanpy
-----------

If have an expression matrix and would like to do standard
preprocessing, embedding and clustering, you can use our minimal Scanpy pipeline through 
the command line toold ``cbScanpy``.

If you already have a .h5ad file, you can use the program ``cbImportScanpy`` 
to export the data to a directory and then create the Cell Browser html directory with the ``cbBuild``
command.

If you are already using Scanpy, you can convert your anndata Scanpy objects
directly to Cell Browser format and start a webserver, e.g. from Jupyter,
directly with the Python3 function ``cellbrowser.scanpyToCellbrowser(ad, outDir, datasetname)``.

A standard Scanpy pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^

Requirements: Python3 with Scanpy installed, see https://scanpy.readthedocs.io/en/latest/installation.html.
Please make sure that you install the igraph library. It's a requirement for the most basic scanpy features,
but it's not an official requirement of scanpy. The command `pip install scanpy[louvain]` will make sure
that igraph is installed.

We provide a wrapper around Scanpy which runs filtering, PCA,
nearest-neighbors, clustering, t-SNE and UMAP and formats them for cbBuild. An
example file is on our downloads server::

    mkdir ~/cellData
    cd ~/cellData
    rsync -Lavzp genome-test.gi.ucsc.edu::cells/datasets/pbmc3k/ ./pbmc3k/ --progress
    cd pbmc3k

Write an empty **scanpy.conf**::

    cbScanpy --init

Edit the *scanpy.conf* file and adapt it to your needs or just keep the default values. Then run commands like this::
    
    # process the matrix and write results to scanpyout/
    cbScanpy -e filtered_gene_bc_matrices/hg19/matrix.mtx -o scanpyout -n pbmc3k
    # build the html directory from scanpyout/
    cd scanpyout
    cbBuild -o ~/public_html/cb -p 8888

Convert a Scanpy .h5ad
^^^^^^^^^^^^^^^^^^^^^^

The tool for this task is called cbImportScanpy. Specify the input file and the
output directory, then run ``cbBuild`` on the output directory.  There is an
example input file in the cellbrowser github repository::

   cd sample_data/pbmc_small
   cbImportScanpy -i anndata.h5ad -o pbmc3kImportScanpy 
   cd pbmc3kImportScanpy
   cbBuild -o ~/public_html/cb

Convert a Scanpy object
^^^^^^^^^^^^^^^^^^^^^^^

From Jupyter or Python3, create a data directory with the tab-sep files and a basic cellbrowser.conf::

    import cellbrowser.cellbrowser as cb
    cb.scanpyToCellbrowser(adata, "scanpyOut", "myScanpyDataset")

Then, build the cell browser from this output directory into a html directory::

    cb.build("scanpyOut", "~/public_html/cells")

If you don't have a webserver running already, start an http server to serve this directory::

    cb.serve("~/public_html/cells", 8888)

You can later stop this http server::

    cb.stop()

Or from a Unix Shell, build and start the http server::

    cd scanpyOut
    cbBuild -o ~/public_html/cells/ -p 8888

