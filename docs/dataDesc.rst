Describing datasets
-------------------

A dataset can be described with three HTML files, summary.html, methods.html
and downloads.html.  You can put these in the same directory where
``cellbrowser.conf`` is stored and they will get copied 
along to the webserver and shown in the ``File > Open Dataset...`` dialog.

Howver, when you have many datasets, writing the html files gets repetitive.
This is where datasetDesc.conf is handy, it's a key-value file with the
description of the dataset in a standardized format.

A sample file can be created with the command ``cbBuild --init``.

The following lists all tags that are currently supported.

These tags contain longer text:

- ``title``: title of the dataset, often the paper title
- ``abstract``: a big picture summary of the dataset
- ``unitDesc``: a description of the values / the unit in the expression matrix
  (e.g. 'TPM' or 'log'ed counts')

This tag contains a file name:

- ``image``: picture, usually a 400px-wide thumbnail of the dimensionality reduction

The following tags can contain URLs and optionally, separated with a space, a label for the link:

- ``biorxiv_url``: URL of the pre-print
- ``paper_url``: URL to any website with the fulltext
- ``other_url``: URL to a website that describes the dataset

The following tags contain accession IDs and will be translated to links:

- ``pmid``: Pubmed ID of the publication (CIRM TagsV5)
- ``geo_series``: NCBI GEO series ID (CIRM TagsV5)
- ``dbgap``: NCBI dbGaP accession, starts with phs

The following tags contain text:
- ``submitter``: name and/or email of submitter
- ``lab``: lab and University of submitter
- ``submission_date``: ideally in format year-month-day
- ``version``: version of dataset, a number that is increased over time
