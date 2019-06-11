Annotate Genes
--------------

When you load a basic set of marker genes into your cell browser, they will be
imported as gene symbols along with an associated score, with no extra
annotations. To make the list of markers genes more useful, you can add extra
annotations using the ``cbMarkerAnnotate`` tool. This tool will add information
about in other resources that describe gene expression profiles
(Allen Brain Atlas), diseases they have been linked to (OMIM, HPO, Sfari), and
protein class (HPRD).

Running this script on your marker genes file is very simple::

    cbMarkerAnnotate inFname outFname

The format for ``inFname`` is the same as for standard cell browser marker gene
files, a tsv or csv table with at least three columns, in this order:

1. cluster - needs to match ``labelField`` in ``cellbrowser.conf``.
2. gene - can be a gene symbol or Ensembl gene ID, with or without the version.
3. score - scores are typically "avg_diff" or "p-Value" or similar. Gene 

``cbMarkerAnnotate`` will map Ensembl gene IDs to symbols and then lookup various
gene-related databases to add more columns to ``inFname`` and write the result to
``outFname``. You can then change the marker genes file described by the ``markers``
parameter in your ``cellbrowser.conf`` to point to ``outFname``.

You can also add your own annotations (e.g. number of associated PubMed articles) to
your marker genes files. These will be displayed alongside any other annotations
that you may have added using ``cbMarkerAnnotate``.
