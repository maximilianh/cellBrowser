UCSC Single Cell Browser
========================

Funded by the California Institute of Regenerative Medicine and the
Chan-Zuckerberg Initiative, see https://www.chanzuckerberg.com/.

This is a viewer for single cell data. It allows you to load an expression
matrix and cell annotation (meta data) file and color the plot by gene or
annotation.

For a demo of the browser, see http://cells.ucsc.edu

This repo contains many different pieces in Python and Javascript, but the main script is
cbBuild. It is a Python script that takes a gene expression matrix and related files and
converts the output to JSON and binary files to the output webserver directory.

Requirements: Python2.6+ or Python3+

Feedback? Open a ticket or email us at cells@ucsc.edu

There is a sample dataset in sampleData/sample1, it's a minimal expression
matrix for a few thousand cells and only the first 100 genes and a bit of meta
data for the cells..

You can build a viewer for it in the directory ~/public_html/cb/ and serve that directory on port 8888:

    cd sampleData/sample1/
    ../../src/cbBuild -o ~/public_html/cb/ -p 8888

The file cellbrowser.conf in sampleData/sample1/ explains all the various settings
that are available in this config file. E.g. you can change the colors, add acronym tables,
add file names, add more marker gene tables, etc.

Then point your web browser to http://localhost:8888. To stop the web server, press Ctrl-C. 
You will have to re-run cbBuild again with -p8888 to look at it again.

To deploy the result onto a real webserver, simply copy all files and directories
under "cells" to an empty directory on a webserver and point your
web browser to it. E.g. many universities give their employees homepage
directories, sometimes in a directory called "~/public_html" or on a special server.

To add more datasets, go to the other data directories and run cbBuild
there, with the same output directory. cbBuild will then modify the index.html
in the output directory to show both datasets (or more). Note that the
directory that you provide via -o or the CBOUT environment variable is the html
directory. The data for each individual dataset will be copied into
subdirectories under this html directory, one directory per dataset.

# Using real webserver

The -p 8888 is optional. A more permanent alternative to the -p option is to
run a webserver on your machine and build directly into its web directory.

On a Mac you can use the Apache that ships with OSX:

    sudo /usr/sbin/apachectl start
    sudo ../../src/cbBuild -o /Library/WebServer/Documents/cells/

Then you should be able to access your viewer at http://localhost/cells

On Linux, you would use the directory /var/www/ instead:

    sudo ../../src/cbBuild -o /var/www/

We hope you do not use this software on Windows. We could make it work, as it's
only Python but we would rather avoid making it work on Windows.

### Process an expression matrix with ScanPy

Requirements: python3 with Scanpy installed, see https://scanpy.readthedocs.io/en/latest/installation.html.

We provide a wrapper around Scanpy which runs filtering, PCA, nearest-neighbors, clustering, t-SNE and
UMAP and formats them for cbBuild. An example file is on our downloads server:

    mkdir ~/cellData
    cd ~/cellData
    rsync -Lavzp genome-test.gi.ucsc.edu::cells/datasets/pbmc3k ./pbmc3k/ --progress
    cd pbmc3k
    ../../cellBrowser/src/cbScanpy -e filtered_gene_bc_matrices/hg19/matrix.mtx -o scanpyout -n pbmc3k
    cd scanpyout
    ../../cellBrowser/src/cbBuild -o ~/public_html/cb

Currently only the genes are exported that were used by Scanpy and only their
normalized and log'ed value. This has advantages, but also disadvantages.  I am
looking forward on your feedback on how this should be done ideally.

### Convert an existing Scanpy object to a cell browser

From Jupyter or Python3, create a data directory with the tab-sep files:

    sys.path.append("cellbrowser/src/cbPyLib")
    import cellbrowser
    # convert to tsv files and create a cellbrowser.conf
    cellbrowser.scanpyToTsv(adata, "scanpyOut", "myDataset")

Then build the cell browser from the Unix shell into a html directory:

    cbBuild -i scanpyOut/cellbrowser.conf -o ~/public_html/cb/

### Convert a CellRanger directory

    cbCellranger -i inputDir -o outputDir
    cd outputDir
    cbBuild

### Convert a Seurat object

    Use src/cbSeurat. More instructions later.

### Optional Python modules to install

In cellbrowser.conf you can specify a color file. If this file contains html color names, you
have to install the module webcolors:

    pip install webcolors

### Adding your dataset

Go to the directory with the expression matrix and the cell annotations. Start from a sample cellbrowser.conf:

    cbBuild --init

Make sure that your files have the correct line endings and fix the line endings if necessary with mac2unix or dos2unix.

    file *.txt *.csv *.tsv *.tab

It's a good idea to gzip your expression matrix at this stage. The expression matrix must have
only one column for the gene identifiers. Ideally your expression matrix is a
tab-separated file and has as many sample columns as you have rows in the meta
data file  and they appear the same order. If this is the case, the conversion of the matrix
is much quicker.

Edit cellbrowser.conf and adapt at least the values for meta, exprMatrix, labelField, clusterField and coordFiles.

From this directory, run 

    cbBuild -o <yourWebserverHtmlDirectory> -p 8888

Navigate your internet browser to the name of the server (or localhost, if you're running this on your own machine)
followed by :8888, e.g. http://localhost:8888.

This is early testing research software, many things have not been properly tested yet. When you run into problems, just open a ticket or send email to cells@ucsc.edu.
