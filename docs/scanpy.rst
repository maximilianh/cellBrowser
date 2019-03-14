With Scanpy
-----------

If you got a .h5ad file from a collaborator, you can use the program ``cbImportScanpy`` 
to export the data and then create the Cell Browser html directory with the ``cbBuild``
command.

If you are already using Scanpy, you can convert your anndata Scanpy objects
directly to Cell Browser format and start a webserver, e.g. from Jupyter,
without the ``cbBuild`` script or the Unix command line.

If you are not using Scanpy but have an expression matrix and would like to do
preprocessing, embedding and clustering, you can use our minimal Scanpy pipeline.

Convert a Scanpy .h5ad file
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Specify the input file and the output directory, then run ``cbBuild`` on the output directory.
There is an example file in the Cell Browser github repository::

   cd sample_data/pbmc_small
   cbImportScanpy -i anndata.h5ad -o pbmc3kImportScanpy 
   cd pbmc3kImportScanpy
   cbBuild -o ~/public_html/cb

Convert a Scanpy object
^^^^^^^^^^^^^^^^^^^^^^^

From Jupyter or Python3, create a data directory with the tab-sep files and a basic cellbrowser.conf::

    import cellbrowser.cellbrowser as cb
    cb.scanpyToTsv(adata, "scanpyOut", "myScanpyDataset")

Then, build the cell browser into a html directory::

    cb.build("scanpyOut", "~/public_html/cells")

If you don't have a webserver running already, start an http server to serve this directory::

    cb.serve("~/public_html/cells", 8888)

You can later stop this http server::

    cb.stop()

Or from a Unix Shell, build and start the http server::

    cd scanpyOut
    cbBuild -o ~/public_html/cells/ -p 8888

A minimal Scanpy pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^

Requirements: python3 with Scanpy installed, see https://scanpy.readthedocs.io/en/latest/installation.html.
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

Edit the *scanpy.conf* file and adapt it to your needs or just leave the default values. Then::
    
    cbScanpy -e filtered_gene_bc_matrices/hg19/matrix.mtx -o scanpyout -n pbmc3k
    cd scanpyout
    cbBuild -o ~/public_html/cb -p 8888

