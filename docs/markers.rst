Annotate Genes
--------------

If you load your existing marker genes into the cell browser, they will by
default only be imported as symbols, with no annotations. You may want to 
annotate the marker genes with information about their gene expression profile
(Allan Brain Atlas), the diseases they have been linked to (OMIM, HPO, Sfari),
their protein class (HPRD) or other annotations (number of PubMed publications).

To this end, run the tool cbMarkerAnnotate. The syntax is very simple::

    cbMarkerAnnotate inFname outFname

The format for inFname is the same as for the cellbrowser marker gene files, a
tab-sep or comma-sep table with at least three columns, in this order: cluster,
gene, score. Typical scores are "avg_diff" or "p-Value" or similar. Gene can be
a gene symbol or Ensembl gene ID, with or without the version.

cbMarkerAnnotate will map Ensembl gene IDs to symbols and then lookup various
gene-related databases to add more columns to inFname and write the result to
outFname, in a format that the Cell Browser can easily display.
