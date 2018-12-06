UCSC Single Cell Browser
========================

Funded by the California Institute of Regenerative Medicine and the
Chan-Zuckerberg Initiative https://www.chanzuckerberg.com/.

The UCSC Cell Browser is a viewer for single cell data. You can click on and
hover over cells to get meta information, search for genes to color on and
click clusters to show cluster-specific marker genes, which in turn are
clickable again.

For a demo of the browser, see http://cells.ucsc.edu

The main script is cbBuild. It is a Python program that takes a gene expression
matrix and related files and converts the output to JSON and binary files to
an output directory which can be put onto a webserver or shown with the built-in
webserver.

The following documentation explains how to create a cellbrowser from a:

* Cellranger directory: command line tool `cbImportCellranger`
* Seurat rds file: command line tool `cbImportSeurat`
* Scanpy h5ad file: command line tool `cbImportScanpy`
* Seurat object: `ExportToCellbrowser()` from R
* Scanpy object: `scanpyToCellbrowser()` from Python
* Expression matrix
  * running a basic Seurat pipeline: `cbSeurat` command line tool
  * running a basic Scanpy pipeline: `cbScanpy` command line tool

This is early research software. You are likely to find bugs. Please open a Github
ticket or email us at cells@ucsc.edu, we can usually fix them quickly.

Installation with pip
---------------------

You need Python2.5+ or Python3+ and pip. On a Mac or any Linux, simply run:

    sudo pip install cellbrowser

On OSX, if this says "command not found", you need to setup pip first:

    sudo easy_install pip

On Linux, if you you're not allowed to run the sudo command, you can install into your user home directory:

    pip install --user cellbrowser
    export PATH=$PATH:~/.local/bin
    
Installation with conda
-----------------------

Alternatively, if you prefer to install through bioconda, since 0.4.23 you can do:

    conda install -c bioconda ucsc-cell-browser
    
Installation: git clone
-----------------------

As an alternative to pip or conda, you can also git clone the repo and
run the command line scripts under cellbrowser/src:

    git clone https://github.com/maximilianh/cellBrowser.git --depth=10
    cd cellBrowser/src
    
Create site with data
---------------------

After installing through one of the methods above, you should be able to run the cbBuild command and see the help messages:

    cbBuild

Here is a small example dataset (Nowakowski et al 2018, fetal brains). The
expression matrix includes only the first 100 genes, otherwise quite a few
features of the browser are used. Download and extract it to the current directory with:

    curl -s https://cells.ucsc.edu/downloads/samples/mini.tgz | tar xvz
    cd mini

You can build a browser consisting of html and other files into the directory
~/public_html/cells/ and serve that directory on port 8888:

    cbBuild -o ~/public_html/cells/ -p 8888

Then point your web browser to http://localhost:8888. To stop the web server, press Ctrl-C. 
You will have to re-run cbBuild again with -p8888 to look at it again.

The file cellbrowser.conf explains all the various settings that are available
in this config file. E.g. you can change the colors, add acronym tables, add
file names, add more marker gene tables, etc.

Putting it onto the internet
============================

To deploy the result onto a real webserver, simply copy all files and directories
under ~/public_html/cells to an empty directory on a webserver and point your
web browser to it. E.g. many universities give their members webspace,
sometimes in a directory called ~/public_html or on a special server. If you
don't have that, contact us or use Cyverse or Amazon S3 to host your files, not
Dropbox. You cannot use Dropbox, iCloud OneDrive or Google Drive, they are not webservers.

To add more datasets, go to the other data directories and run cbBuild
there, with the same output directory. cbBuild will then modify the index.html
in the output directory to show all datasets. Note that the directory that you
provide via -o (or the CBOUT environment variable) is the html directory. The
data for each individual dataset will be copied into subdirectories under this
html directory, one directory per dataset.

Instead of specifying "-o" all the time, you can also add a line like this to
your ~/.bashrc to point to your html directory:
 
    export CBOUT=/var/www

The -p 8888 is optional. A more permanent alternative to the -p option is to
run a webserver on your machine and build directly into its web directory.

On a Mac you can use the Apache that ships with OSX:

    sudo /usr/sbin/apachectl start
    sudo cbBuild -o /Library/WebServer/Documents/cells/

Then you should be able to access your viewer at http://localhost/cells

On Linux, you would install Apache2 (with 'sudo yum install htppd' or 'sudo apt-get install
apache2') and use the directory /var/www/ instead:

    sudo cbBuild -o /var/www/

We hope you do not use this software on Windows. Contact us if you have to.

Convert a Seurat object
=======================

The function ExportToCellbrowser() will be part of Seurat 3. You can install pre-release Seurat3 like this:

    install.packages("devtools")
    devtools::install_github("satijalab/seurat", ref = "release/3.0")

For Seurat 2, you have to load the function with this command:

    source("https://raw.githubusercontent.com/maximilianh/cellBrowser/master/src/R/ExportToCellbrowser-Seurat2.R")

You can then write a Seurat object to a directory from which you can run cbBuild:

    ExportToCellbrowser(pbmc_small, dir="pbmcSmall", cb.dir="htdocs", dataset.name="pbmcSmall")

Or immediately convert the files to html and serve the result on port 8080 and open a web browser:

    ExportToCellbrowser(pbmc_small, dir="pbmcSmall", cb.dir="htdocs", dataset.name="pbmcSmall", port=8080)

A minimal Scanpy pipeline
=========================

Requirements: python3 with Scanpy installed, see https://scanpy.readthedocs.io/en/latest/installation.html.

We provide a wrapper around Scanpy which runs filtering, PCA, nearest-neighbors, clustering, t-SNE and
UMAP and formats them for cbBuild. An example file is on our downloads server:

    mkdir ~/cellData
    cd ~/cellData
    rsync -Lavzp genome-test.gi.ucsc.edu::cells/datasets/pbmc3k/ ./pbmc3k/ --progress
    cd pbmc3k

Write an empty scanpy.conf:

    cbScanpy --init

Edit the scanpy.conf file and adapt it to your needs or just leave the default values. Then:
    
    cbScanpy -e filtered_gene_bc_matrices/hg19/matrix.mtx -o scanpyout -n pbmc3k
    cd scanpyout
    cbBuild -o ~/public_html/cb -p 8888

Currently only the genes are exported that were used by Scanpy and only their
normalized and log'ed value, so the matrix after filtering. This has
advantages, but also disadvantages.  Contact us if you have an opinion on which
expression value should be shown. You can also manually copy your original
expression matrix into the output directory ("scanpyout" in the example) to 
include all genes.

Convert an existing Scanpy object
=================================

The cbScanpy wrapper runs some generic analysis steps with very crude default
values. If you have done the analysis already, you can build a cellbrowser from
your existing Scanpy object.

From Jupyter or Python3, create a data directory with the tab-sep files and a basic cellbrowser.conf:

    import cellbrowser.cellbrowser as cb
    cb.scanpyToTsv(adata, "scanpyOut", "myScanpyDataset")

Then, build the cell browser into a html directory:

    cb.build("scanpyOut", "~/public_html/cells")

If you don't have a webserver running already, start an http server to serve this directory: 

    cb.serve("~/public_html/cells", 8888)

You can later stop this http server:

    cb.stop()

Or from a Unix Shell, build and start the http server:

    cd scanpyOut
    cbBuild -o ~/public_html/cells/ -p 8888

Convert an Seurat object
========================

The function ExportToCellbrowser() will be part of Seurat 3. You can install the pre-release Seurat3 like this:

    install.packages("devtools")
    devtools::install_github("satijalab/seurat", ref = "release/3.0")

For Seurat 2, you have to load the function with this command:

    source("https://raw.githubusercontent.com/maximilianh/cellBrowser/master/src/R/ExportToCellbrowser-seurat2.R")

You can then write a Seurat object to a directory from which you can run cbBuild:

    ExportToCellbrowser(pbmc_small, dir="pbmcSmall", cb.dir="htdocs", dataset.name="pbmcSmall")

Or immediately convert the files to html and serve the result on port 8080 and open a web browser:

    ExportToCellbrowser(pbmc_small, dir="pbmcSmall", cb.dir="htdocs", dataset.name="pbmcSmall", port=8080)

Convert CellRanger results
==========================

Find the cellranger OUT directory, it contains an "analysis" directory and also
a subdirectory "filtered_gene_bc_matrices". This is the directory that is the
input directory for our tool cbImportCellranger. The tool converts the
cellranger files to tab-separated files, then run cbBuild on these.

To import Cellranger .mtx files, we need the scipy package (add --user if you are not admin on your machine):

    pip install scipy

Let's use an example, the pbmc3k cellranger output files from the 10x website:

    rsync -Lavzp genome-test.gi.ucsc.edu::cells/datasets/pbmc3kCellranger/ ./pbmc3kCellranger/ --progress
    cbImportCellranger -i pbmc3kCellranger -o cellrangerOut --name pbmc3k_cellranger
    cd cellrangerOut
    cbBuild -o ~/public_html/cells -p 9999

### Process an expression matrix with a basic Seurat pipeline

First make sure that you can install the package "hdf5r" in R:

    Rscript -e "install.packages('hdf5r' , dep=TRUE, repos='http://cran.r-project.org/')"

If the above doesn't work, try installing the fake-hdf5r package, which means that you won't be able to read 
hdf5 files, but reading .mtx and of course tab-sep files will still work:

    Rscript -e "install.packages('remotes' , dep=TRUE, repos='http://cran.r-project.org/')"
    Rscript -e "remotes::install_github('UCSF-TI/fake-hdf5r')"

Then install Seurat into your default command line R (not RStudio or another R version you may have):

    Rscript -e "install.packages(c('Seurat', 'data.table'), dep=TRUE, repos='http://cran.r-project.org/')"

To run an example now, download the 10X pbmc3k expression matrix:

    rsync -Lavzp genome-test.gi.ucsc.edu::cells/datasets/pbmc3k/ ./pbmc3k/ --progress

Create a default seurat.conf:

    cbSeurat --init

You can modify seurat.conf but the default values are good for this dataset.
Now run the expression matrix filtered_gene_bc_matrices/hg19/matrix.mtx through Seurat like this:

    cbSeurat -e filtered_gene_bc_matrices/hg19 --name pbmc3kSeurat -o seuratOut 

This will create a script seuratOut/runSeurat.R, run it through Rscript and will fill the directory seuratOut/ with everything needed to create a cell browser. Now you can build your cell browser from the Seurat output:

    cd seuratOut
    cbBuild -o 

You can modify the file seurat.conf and rerun the cbSeurat command above.

### Adding a dataset from tab-separated files

Go to the directory with the expression matrix and the cell annotations. Start from a sample cellbrowser.conf:

    cbBuild --init

Then you need at least three but ideally four files, they can be in .tsv or .csv format:

1. The expression matrix, one row per gene, ideally gzipped. The first column
   must be the gene identifier or gene symbol, or ideally
   geneId|symbol. ENSG and ENSMUSG gene identifiers will be translated
   automatically to symbols. The other columns are expression values as
   numbers, one per cell. The number type will be auto-detected (float or int).
   The file must be a header line that describes the columns with the
   identifiers for the cells.

2. The cell annotation meta data table, one row per cell. No need to gzip this
   relatively small file. The first column is the name of the cell and it has
   to match the cell name in the expression matrix. There should be at least
   two columns: one with the name of the cell and one with
   the name of the cluster. Ideally your expression matrix is a tab-separated
   file and has as many cell columns as you have rows in the meta data file
   and they appear in the same order in both files, as cbBuild doesn't have to
   trim the matrix then or reorder the meta file. The meta file has a header
   file, the names of the columns will be refered to later in the cellbrowser.conf file.

3. The coordinates of the cells, often t-SNE coordinates. This file always has three
   columns, (cellName, x, y). The cellName must be the same as in the expression
   matrix and cell annotation meta data file. You can provide multiple files
   in this format, if you have run multiple dimensionality reduction algorithms.

4. The (optional) table with cluster-specific marker genes. The first column is
   the cluster name (from the cell annotation meta file), the second column 
   contains the gene symbol (or gene ID, will be mapped to symbol) and the
   third column is some numeric score (e.g.  p-Value or FDR).  You can add as
   many other columns as you like with additional information about this gene
   or run your table through cbMarkerAnnotate to add information from various
   gene-centric databases to your existing table. Alternatively you can also
   provide the raw Seurat marker gene output. There can be multiple files with
   cluster-specific marker genes, e.g. in case that you are also doing
   differential gene expression analysis or have results from multiple
   algorithms. 

Make sure that all your input files have Unix line endings and fix the line endings if necessary with mac2unix or dos2unix.

    file *.txt *.csv *.tsv *.tab

Edit cellbrowser.conf. Enter the name of the three files for the tags  exprMatrix, meta, coordFiles. If you have
a table with cluster specific genes, put that into clusterFiles.
Enter the value of your cluster name field from the meta annotation file for the tags labelField and clusterField.

From the directory where your cellbrowser.conf is located, run 

    cbBuild -o /tmp/cb -p 8888

Navigate your internet browser to the name of the server (or localhost, if you're running this on your own machine)
followed by :8888, e.g. http://localhost:8888.

This is early testing research software, many things have not been properly tested yet. When you run into problems, just open a ticket or send email to cells@ucsc.edu.

### Combining Seurat, Scanpy and Cellranger result into a single browser

You can use `cbTool metaCat` to merge the meta.tsv files from different pipelines into a single one, like this:

    cbTool metaCat myMeta.tsv seuratOut/meta.tsv scanpyOut/meta.tsv ./newMeta.tsv --fixDot

The option --fixDot will work around R's strange habit of replacing special characters in the cell identifiers with ".".
Directories created with ExportToCellbrowser() should not have this problem, but others may.

You can now take one of the auto-generated cellbrowser.conf files or start from a fresh one with `cbBuild --init`.
In this cellbrowser.conf, add all the coordinates files from all your pipelines. Unfortunately, right now you can
only have a single marker gene list.

### Optional Python modules to install

In cellbrowser.conf you can specify a color file, the format is .tsv or .csv and it has two columns, clusterName<tab>colorCode. If this file contains html color names instead of color codes, you have to install the module webcolors:

    pip install webcolors

To read expression matrices in .mtx format, you have to install scipy:

    pip install scipy

