UCSC Single Cell Browser
========================

For a demo, see http://cells.ucsc.edu

This repo contains different parts:
* cbAdd: a Python script that takes a gene expression matrix and a cell annotation (meta) table
  and converts the output to JSON and binary files
* cbMake: a Python script that finds all results of cbAdd in a directory and creates an index.html
  for it using cellbrowser.js and a few others.

There is a sample dataset in sampleData/sample1.

A viewer was created from it with these commands:

    cd sampleData/sample1/
    ../../cbAdd -o ~/public_html/cbTest
    ../../cbMake -o ~/public_html/cbTest

To deploy the result onto a webserver, copy all files ~/public_html/cbTest to
an empty directory on a webserver.
