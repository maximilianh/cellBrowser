With Scanpy
-----------

There area few different ways to create a cell browser using Scanpy:

* **Run our basic Scanpy pipeline** - with just an expression matrix and ``cbScanpy``, you can the standard preprocessing, embedding, and clustering through Scanpy.
* **Import a Scanpy h5ad file** - create a cell browser from your ``h5ad`` file using the command-line program ``cbImportScanpy``.
* **Use a few Python 3 function** - you can build a cell browser from a Scanpy ``h5ad`` file and start a web server, e.g. from Jupyter, with the Python3 function ``cellbrowser.scanpyToCellbrowser(ad, outDir, datasetname)``.

A standard Scanpy pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^

Requirements: Python3 with Scanpy installed, see their `installation instructions <https://scanpy.readthedocs.io/en/latest/installation.html>`_ for information about setting up Scanpy.
As part of the Scanpy installion process, ensure that the igraph library is also installed.
It's needed for the most basic scanpy features even though it's not an official requirement.
The command ``pip install scanpy[louvain]`` will make sure that igraph is installed.

We provide a wrapper around Scanpy, named ``cbScanpy``, which runs filtering, PCA,
nearest-neighbors, clustering, t-SNE, and UMAP. The individual steps are explained in more detail in 
the `Scanpy PBMC3k tutorial <https://icb-scanpy-tutorials.readthedocs-hosted.com/en/latest/pbmc3k.html>`_.

The output of ``cbScanpy`` is formatted
to be directly usable to build a cell browser with ``cbBuild``. 

You can test ``cbScanpy`` yourself using the following set of steps. 
To process an example dataset, download the 10x pbmc3k expression matrix from our servers::

    mkdir ~/cellData
    cd ~/cellData
    rsync -Lavzp genome-test.gi.ucsc.edu::cells/datasets/pbmc3k/ ./pbmc3k/ --progress
    cd pbmc3k

Next, run the expression matrix ``filtered_gene_bc_matrices/hg19/matrix.mtx`` through Scanpy::
    
    cbScanpy -e filtered_gene_bc_matrices/hg19/matrix.mtx -o scanpyOut -n pbmc3k
    
This will run Scanpy and will fill the directory ``scanpyOut/`` with everything needed
to create a cell browser. After the ``cbScanpy`` script completes, you can build your
cell browser from the output::

    cd scanpyout
    cbBuild -o ~/public_html/cb -p 8888

Changing the defaults using ``scanpy.conf``
"""""

This set of steps will run a basic Scanpy pipeline with the default settings. You can modify the settings
for Scanpy by creating a ``scanpy.conf``::

    cbScanpy --init

You can edit the settings in ``scanpy.conf`` and re-run the ``cbScanpy`` command to generate a new set of
Scanpy output using these new settings.

Convert a Scanpy ``h5ad``
^^^^^^^^^^^^^^^^^^^^^^

If you have run Scanpy and have an output ``h5ad`` file, you can import it 
into a cell browser using the command ``cbImportScanpy``. 

The steps in this section walk you through the process of importing data from a
Scanpy file and then building a cell browser from the output. The steps use an example 
``h5ad`` file available for a small pbmc dataset from our Github repo: 
`anndata.h5ad <https://github.com/maximilianh/cellBrowser/blob/master/sampleData/pbmc_small/anndata.h5ad>`_.

First, use ``cbImportScanpy`` to extract the data from the ``h5ad``::

   cbImportScanpy -i anndata.h5ad -o pbmc3kImportScanpy
   
The ``-i`` option specifies the input ``h5ad`` file and the ``-o`` option specifies
a name for the output directory. You can use the ``-n`` option to change the dataset
name in the cell browser; if it is not specified, it will default to the output
directory name.

The output of ``cbImportScanpy`` will be formatted so that you can immediately
build a cell browser from it. Go into the pbmc3kImportScanpy directory and run
``cbBuild`` to create the cell browser output files::

   cd pbmc3kImportScanpy
   cbBuild -o ~/public_html/cb

Alternatively, you can use the ``--htmlDir`` option for ``cbImportScanpy`` to automatically
run cbBuild for you::

    cbImportScanpy -i anndata.h5ad -o pbmc3kImportScanpy --htmlDir=~/public_html/cb

Convert a Scanpy object
^^^^^^^^^^^^^^^^^^^^^^^

From Jupyter or Python3, you can create a data directory with the necessary
tsv files and a basic ``cellbrowser.conf``::

    import cellbrowser.cellbrowser as cb
    cb.scanpyToCellbrowser(adata, "scanpyOut", "myScanpyDataset")

Here ``adata`` is your Scanpy object, ``scanpyOut`` is your output directory, and
``myScanpyDataset`` is your dataset name.

Then, build the cell browser from this output directory into a html directory::

    cb.build("scanpyOut", "~/public_html/cells")

If you don't have a web server running already, use this function start one to serve up this directory::

    cb.serve("~/public_html/cells", 8888)

You can stop the web server with the function::

    cb.stop()

Or from a Unix shell, you can build and start a web server using ``cbBuild``::

    cd scanpyOut
    cbBuild -o ~/public_html/cells/ -p 8888

