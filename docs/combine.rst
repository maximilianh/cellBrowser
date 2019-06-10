``cbTool``: combine and convert your data
=====

The script ``cbTool`` included in the Cell Browser package includes a number of utilties
for combining or converting your data. These different functions and how to use them are
described below. 

Combining results
---------------

Metadata
^^^^^

You can use the ``cbTool metaCat`` utility to merge the metadata files from different
sources or pipelines (e.g. ``cbScanpy`` or ``cbSeurat``) into a single one. A command
to do combine a metadata set of metadata files from Scanpy and Seurat with a separate
one might look like this::

    cbTool metaCat myMeta.tsv seuratOut/meta.tsv scanpyOut/meta.tsv ./newMeta.tsv --fixDot

The resulting file will include the columns from all three of the original files
combined into a new metadata file. Note that ``cbTool metaCat`` assumes that the first
column of each file contains the same cell identifier that it can use to join them.

Matrices 
^^^^^

Similar to ``metaCat`` for combining metadata files, ``cbTool matCat`` can be used to
combine expression matrices from different sources. A command combine the two different
matrixes would look like this::

    cbTool matCat mat1.tsv.gz mat2.tsv.gz exprMatrix.tsv.gz



Converting mtx to tsv
-------

Using ``cbTool mtx2tsv``, you can convert your expression matrix in `matrix market
<https://math.nist.gov/MatrixMarket/formats.html>`_ to tsv format
(with one gene per line and one cell per column)::

    cbTool mtx2tsv matrix.mtx genes.tsv barcodes.tsv exprMatrix.tsv.gz


Fixing R Seurat output
-----

The option ``--fixDot`` for ``cbTool`` will work around R's strange habit of replacing
special characters in the cell identifiers with ".". Directories created with the
``ExportToCellbrowser()`` function from R should not have this problem, but others may. 
