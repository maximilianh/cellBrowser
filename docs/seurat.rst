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

Then, on the Unix command line, you specify the input .rds file, the output directory and
the name for the dataset in the cell browser::

    cbImportSeurat2 pbmc3k_small.rds pbmc3kCb pbmc3k-rdsExport

Then go into the directory *pbmc3kCb* and run cbBuild to create the Cell Browser html files::

    cd pbmc3kCb
    cbBuild -o ~/public_html/cb

Convert a Seurat object from R
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function ExportToCellbrowser() is already part of Seurat 3. You can install pre-release Seurat3 like this::

    install.packages("devtools")
    devtools::install_github("satijalab/seurat", ref = "release/3.0")

For Seurat 2, you have to load the function with this command::

    source("https://raw.githubusercontent.com/maximilianh/cellBrowser/master/src/cbPyLib/cellbrowser/R/ExportToCellbrowser-seurat2.R")

You can then write a Seurat object to a directory from which you can run cbBuild::

    ExportToCellbrowser(pbmc_small, dir="pbmcSmall", cb.dir="htdocs", dataset.name="pbmcSmall")

Or immediately convert the files to html and serve the result on port 8080 and open a web browser::

    ExportToCellbrowser(pbmc_small, dir="pbmcSmall", cb.dir="htdocs", dataset.name="pbmcSmall", port=8080)

Writing the expression matrix is somewhat slow. If you have already exported into the same 
output directory before and just updated a part of the cell annotation data
(e.g. clustering), you can use the argument *skip.matrix=TRUE* to save some
time.

Run a basic Seurat pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have never used Seurat before and just want to process an expression matrix
as quickly as possible, this section is for you.

First make sure that you can install the package "hdf5r" in R::

    Rscript -e "install.packages('hdf5r' , dep=TRUE, repos='http://cran.r-project.org/')"

If the above doesn't work, try installing the fake-hdf5r package, which means
that you won't be able to read hdf5 files, but reading .mtx and of course
tab-sep files will still work::

    Rscript -e "install.packages('remotes' , dep=TRUE, repos='http://cran.r-project.org/')"
    Rscript -e "remotes::install_github('UCSF-TI/fake-hdf5r')"

Then install Seurat into your default command line R (not RStudio or another R version you may have)::

    Rscript -e "install.packages(c('Seurat', 'data.table'), dep=TRUE, repos='http://cran.r-project.org/')"

To run an example now, download the 10X pbmc3k expression matrix::

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

