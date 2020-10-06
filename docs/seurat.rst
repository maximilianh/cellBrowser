With Seurat
-----------

There are a number of ways to create a cell browser using Seurat:

* **Import a Seurat rds file** - create a cell browser with the Unix command line tool ``cbImportSeurat``.
* **Using RStudio and a Seurat object** - create a cell browser directly using the ``ExportToCellbrowser()`` R function. 
* **Run our basic Seurat pipeline** - with just an expression matrix, you can run our ``cbSeurat`` pipeline to create a cell browser.

Each of these methods are described in more detail below.

Convert a Seurat ``rds`` or ``.rdata` file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, create an .rds file in R as described in the Seurat tutorial::

    saveRDS(pbmc, "pbmc3k_small.rds")

Next, on the Unix command line, use the ``cbImportSeurat`` script to convert this ``rds``
file into a cell browser::

    cbImportSeurat -i pbmc3k_small.rds -o pbmc3kImport

This works with objects created by versions 2 and 3 of Seurat. ``cbImportSeurat`` can read 
both ``.rds`` and ``.rdata`` files, for `.rdata` it assume the first object is the Seurat object.
Make sure that you have the same major version of Seurat installed that was used to create the
object. You cannot open Seurat2 objects with Seurat3 or vice versa. 
(We often need to switch between Seurat versions and found conda environments very helpful for this.)

The ``-i`` option specifies the input ``rds`` file and the ``-o`` option specifies a name for the output
directory. You can use the ``-n`` option to change the dataset name in the cell browser;
if it is not specified, it will default to the output directory name.

A Seurat object does not contain the marker genes by default, as FindAllMarkers() does not save its output.
You can add it to the object when you save the .rds file with a command like this::

    object@misc$markers <- FindAllMarkers(object)

``cbImportSeurat`` will then use these markers. Otherwise, if ``misc$markers`` is not present in the object, it will
run FindAllMarkers with the default values (Wilcoxon and 0.25 as the cutoff). Alternatively, you can also save the markers
to a tab-separated file yourself and provide this file with the ``--markerFile`` option.

Lastly, go into the ``pbmc3kImport`` directory and run ``cbBuild`` to create the cell browser
output files::

    cd pbmc3kImport
    cbBuild -o ~/public_html/cb
    
Alternatively, you can use the ``--htmlDir`` option for ``cbImportSeurat`` to automatically run ``cbBuild`` for you::

    cbImportSeurat -i pbmc3k_small.rds -o pbmc3kImport --htmlDir=~/public_html/cb

Convert a Seurat object from R
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function ExportToCellbrowser() is already part of Seurat 3. You can install pre-release Seurat3 like this::

    install.packages("devtools")
    devtools::install_github("satijalab/seurat", ref = "release/3.0")

For Seurat 2, you have to load the function with this command::

    source("https://raw.githubusercontent.com/maximilianh/cellBrowser/master/src/cbPyLib/cellbrowser/R/ExportToCellbrowser-seurat2.R")

You can then write a Seurat object to a directory from which you can run cbBuild::

    ExportToCellbrowser(pbmc_small, dir="pbmcSmall", dataset.name="pbmcSmall")

Or, you can build a cell browser from this dataset into the ``htdocs`` directory,
serve the result on port 8080 via http, and open a web browser from within R::

    ExportToCellbrowser(pbmc_small, dir="pbmcSmall", cb.dir="htdocs", dataset.name="pbmcSmall", port=8080)

Writing the expression matrix is somewhat slow. If you have already exported into the same 
output directory before and just updated a part of the cell annotation data
(e.g. clustering), you can use the argument *skip.matrix=TRUE* to save some
time:

    ExportToCellbrowser(pbmc_small, dir="pbmcSmall", dataset.name="pbmcSmall", skip-matrix=TRUE)

Run a basic Seurat pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have never used Seurat before and just want to process an expression matrix
as quickly as possible, this section is for you.

If you do not have R installed yet, we recommend that you install it via conda.
To install miniconda, follow their `installation instructions <https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation>`_.

After setting up conda, install R::

    conda install r

Then, install Seurat::

    conda install -c bioconda r-seurat 

To process an example dataset now, download the 10X pbmc3k expression matrix::

    rsync -Lavzp genome-test.gi.ucsc.edu::cells/datasets/pbmc3k/ ./pbmc3k/ --progress

Now run the expression matrix ``filtered_gene_bc_matrices/hg19/matrix.mtx`` through
Seurat::

    cbSeurat -e filtered_gene_bc_matrices/hg19 --name pbmc3kSeurat -o seuratOut 

This will create a script (``seuratOut/runSeurat.R``), run it through Rscript, and
will fill the directory ``seuratOut/`` with everything needed to create a cell
browser. After the ``cbSeurat`` script completes, you can build your cell browser from the output::

    cd seuratOut
    cbBuild -o ~/public_html/cells

Changing the defaults using ``seurat.conf``
""""""""

This set of steps will run a basic Seurat pipeline with the default settings. You can
modify the settings for Seurat by creating a ``seurat.conf`` file::

    cbSeurat --init

You can edit the settings in ``seurat.conf`` and re-run the ``cbSeurat`` command to
generate a new set of Seurat output using these new settings. 
