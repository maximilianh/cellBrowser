UCSC Single Cell Browser Demo
=============================

This repo contains two parts:
* a Python script that runs a gene expression matrix through Seurat and
  massages the output data into JSON
* a Javascript viewer that shows the JSON files.

There is a sample dataset in data/quakeBrainGeo1/.

A viewer was created from it with these commands:

    cd data/quakeBrainGeo1/

    ../../cbPrep matrix -n quakeBrainGeo1 -e geneMatrix.tsv --log2 --skip -m meta.tsv -o ~/public_html/cbTest/ -g markerSymbols.txt -l biosample_cell_type

    ../../cbPrep html -o ~/public_html/cbTest

To deploy the result onto a webserver, copy all files ~/public_html/cbTest to an empty directory on a webserver.
Files in ~/public_html/cbTest/build/ are not needed to be copied over, they are
not used by the viewer and will only speed up future "matrix" runs.

Note that the output directory cbTest contains a directory for this dataset (specified via -n) called "quakeBrainGeo1". The
file dataset.json in this directory can be modified to change: 
* long and short labels of this dataset, a description of how the dataset was created
* the default field to show cluster labels for

A subsequent run of "cbPrep html" will then update the index.html with the information in all <subdirectory>/dataset.json files.

Requirements: Seurat 1.4

Installation:

    echo 'export R_LIBS_USER=$HOME/R' >> ~/.bashrc
    source ~/.bashrc
    mkdir -p $R_LIBS_USER
    wget https://github.com/satijalab/seurat/archive/v1.4.0.tar.gz
    R CMD INSTALL v1.4.0.tar.gz -l $R_LIBS_USER
