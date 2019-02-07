With Seurat
-----------

If you are an RStudio user and have a Seurat object, you can convert
it to html directly without going to the Unix command line. Users of large
servers may prefer to import a Seurat ``rds`` file.  If you have an expression
matrix and no knowledge of Seurat, you can still use our minimal Seurat
pipeline to do some quick inspection of your data with the Cell Browser.

Convert a Seurat2 .rds file
^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use the program ``cbImportSeurat2`` to convert a ``rds`` file to a Cell
Browser. You can create an .rds file in R as described in the Seurat tutorial::

    saveRDS(pbmc, "pbmc3k_small.rds")

Then, on the Unix command line, you specify the input .rds file and the output directory (the name
in the cell browser defaults to the output directory name, but you can change this with -n)::

    cbImportSeurat2 -i pbmc3k_small.rds -o pbmc3kImport

Then go into the directory *pbmc3kImport* and run cbBuild to create the Cell Browser html files::

    cd pbmc3kImport
    cbBuild -o ~/public_html/cb

Convert a Seurat object from R
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function ExportToCellbrowser() is already part of Seurat 3. You can install pre-release Seurat3 like this::

    install.packages("devtools")
    devtools::install_github("satijalab/seurat", ref = "release/3.0")

For Seurat 2, you have to load the function with this command::

    source("https://raw.githubusercontent.com/maximilianh/cellBrowser/master/src/cbPyLib/cellbrowser/R/ExportToCellbrowser-seurat2.R")

You can then write a Seurat object to a directory from which you can run cbBuild::

    ExportToCellbrowser(pbmc_small, dir="pbmcSmall", dataset.name="pbmcSmall")

Or immediately convert the files to html files in the directory ``htdocs`` and
serve the result on port 8080 via http and open a web browser from R::

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
Follow these instructions to install the miniconda installer:
https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation

When conda is installed, install R::

    conda install r

Then, again using conda, install Seurat::

    conda install -c bioconda r-seurat 

To process an example dataset now, download the 10X pbmc3k expression matrix::

    rsync -Lavzp genome-test.gi.ucsc.edu::cells/datasets/pbmc3k/ ./pbmc3k/ --progress

Create a default file ``seurat.conf``::

    cbSeurat --init

You can modify seurat.conf but the default values are good for this dataset.
Now run the expression matrix *filtered_gene_bc_matrices/hg19/matrix.mtx* through
Seurat like this::

    cbSeurat -e filtered_gene_bc_matrices/hg19 --name pbmc3kSeurat -o seuratOut 

This will create a script seuratOut/runSeurat.R, run it through Rscript and
will fill the directory seuratOut/ with everything needed to create a cell
browser. Now you can build your cell browser from the Seurat output::

    cd seuratOut
    cbBuild -o ~/public_html/cells

You can modify the file seurat.conf and rerun the cbSeurat command above.

