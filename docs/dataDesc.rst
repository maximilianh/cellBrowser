Describing datasets
-------------------

A dataset can be described with three HTML files, summary.html, methods.html
and downloads.html.  You can put these in the same directory where
``cellbrowser.conf`` is stored and they will get copied 
along to the webserver and shown in the ``File > Open Dataset...`` dialog.

Howver, when you have many datasets, writing the html files gets repetitive.
This is where desc.conf is handy, it's a key-value file with the
description of the dataset in a standardized format.

A sample file can be created with the command ``cbBuild --init``.

The following lists all tags that are currently supported.

These tags contain longer text that can include HTML markup:

- ``title``: title of the dataset, often the paper title
- ``abstract``: a big picture summary of the dataset, as a string
- ``methods``: the methods for the dataset, as a string
- ``unitDesc``: a description of the values / the unit in the expression matrix
  (e.g. 'TPM' or 'log'ed counts')

Instead of long strings with HTML content for ``abstract`` and ``methods``, you can also create the
files ``abstract.html`` and ``methods.html``, they will be used instead. Or use the 
statements ``abstractFile`` and ``methodsFile`` to specify other file names. In the HTML, 
you can use text like ``<section>some subtitle</section>`` to split the text into sections.

These tags contains a file name:
- ``image``: usually a 400px-wide thumbnail of the dimensionality reduction
- ``rawMatrixFile``: usually the raw unprocessed matrix. Usually a .zip or .gz file. Also see ``rawMatrixNote``.

The following tags can contain URLs and optionally, separated with a space, a label for the link. If you do 
not specify the label, a default label will be used (e.g. 'Biorxiv Preprint'):

- ``biorxiv_url``: URL of the pre-print
- ``paper_url``: URL to any website with the fulltext
- ``other_url``: URL to a website that describes the dataset

The following tags contain accession IDs and will be translated to links:

- ``pmid``: Pubmed ID of the publication (CIRM TagsV5)
- ``geo_series``: NCBI GEO series ID (CIRM TagsV5)
- ``sra``: NCBI SRA accession
- ``sra_study``: NCBI SRA SRPxxxx accession
- ``doi``: DOI of paper fulltext
- ``dbgap``: NCBI dbGaP accession, starts with phs
- ``bioproject``: NCBI Bioproject accession, a 4-9 digit number, without the PRJNA prefix

The following tags contain just text:

- ``submitter``: name and/or email of submitter
- ``lab``: lab and University of submitter
- ``submission_date``: ideally in format year-month-day
- ``rawMatrixNote``: text to describe the raw matrix, see ``rawMatrixFile``
- ``version``: version of dataset, a simple number (1,2,3,...) that should be increased each time a major change (usually meta data) was received from the lab
