With text files
---------------

The generic way to setup a Cell Browser is from .tsv or .csv text files.

Go to the directory with the expression matrix and the cell annotations. Start from a sample ``cellbrowser.conf``::

    cbBuild --init

The current list of all possible cellbrowser.conf statements can be found in `our Github example <https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf>`.

Then you need at least three but ideally four data files, they can be in .tsv or .csv format::

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
   line, the names of the columns from this line are the ones refered to in the 
   cellbrowser.conf file.

3. The coordinates of the cells, often t-SNE or UMAP coordinates. This file
   always has three columns, (cellName, x, y). The cellName must be the same as in
   the expression matrix and cell annotation meta data file. You can provide
   multiple files in this format, if you have run multiple dimensionality
   reduction algorithms.  You can also specify only a subset of the cells here. In
   this way, you can use a single dimensionality reduction algorithm, but multiple
   subsets of the cells, e.ge. one coord file per tissue. If R has changed your
   cell identifiers (adding dots), you may be able to fix them with ``cbTool metaCat``.

4. The (optional) table with cluster-specific marker genes. The first column is
   the cluster name (from the cell annotation meta file), the second column 
   contains the gene symbol (or Ensembl gene ID, will be mapped to symbol) and
   the third column is some numeric score (e.g.  p-Value or FDR).  You can add
   as many other columns as you like with additional information about this
   gene or run your table through cbMarkerAnnotate to add information from
   various gene-centric databases to your existing table. Alternatively you can
   also provide the raw Seurat marker gene output. There can be multiple files
   with cluster-specific marker genes, e.g. in case that you are also doing
   differential gene expression analysis or have results from multiple
   algorithms. See :ref:`Annotate genes` on how to add link to external
   gene-databases (like Allan Brain Atlas or OMIM) to your marker genes.


Make sure that all your input files have Unix line endings and fix the line
endings if necessary with mac2unix or dos2unix::

    file *.txt *.csv *.tsv *.tab

Edit cellbrowser.conf. Enter the name of the three files with the config
statements exprMatrix, meta, coordFiles. If you have a table with cluster
specific genes, put that into clusterFiles.  Enter the value of your cluster
name field from the meta annotation file for the tags labelField and
clusterField.

From the directory where your cellbrowser.conf is located, run::

    cbBuild -o /tmp/cb -p 8888

Point your internet browser to the name of the server (or localhost, if
you're running this on your own machine) followed by :8888, e.g.
http://localhost:8888.

The output directory (/tmp/cb in the example) can hold multiple datasets. 
If you have a second dataset in another directory that contains cellbrowser.conf,
just make sure that the other cellbrowser.conf specifies a different dataset name 
with `name=xxx`. Then run `cbBuild -o /tmp/cb -p 8888` in the other
directory to add the second dataset to your browser output directory /tmp/cb.
