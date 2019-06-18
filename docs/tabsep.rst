With text files
---------------

The generic way to set up your own cell browser is to start from tab-separated (tsv)
or comma-separated (csv) format text files. The steps on this page assume that you
have already gone through the process of clustering your cells.

The files you will need
^^^^^

You will need the first three files described below in tsv or csv format, the fourth is optional::

1. **Expression matrix**: one row per gene and one column per cell, ideally gzipped. The first column
   must be the gene identifier or gene symbol, or ideally
   geneId|symbol. Ensembl/GENCODE gene identifiers starting with ENSG and ENSMUSG will be translated
   automatically to symbols. The other columns are expression values as
   numbers, one per cell. The number type will be auto-detected (float or int).
   The first line of the file must be a header that includes the cell
   identifiers.

2. **Cell annotation metadata table**, one row per cell. No need to gzip this
   relatively small file. The first column is the name of the cell and it has
   to match the name in the expression matrix. There should be at least
   two columns: one with the name of the cell and one with
   the name of the cluster. To speed up processing of both your expression matrix
   and metadata file, these files should describe the same numbers of cells and be
   in the same order. This allows cbBuild process these files without needing to
   trim the matrix and reorder the metadata file. The metadata file also must have
   a header line.

3. **Cell coordinates**, often t-SNE or UMAP coordinates. This file
   always has three columns, (cellName, x, y). The cellName must be the same as in
   the expression matrix and cell annotation metadata file. If you have run
   multiple dimensionality reduction algorithms, you can specify multiple
   coordinate files in this format. The number rows in these coordinates doesn't
   need to match that of your expression matrix or metadata files, allowing you to
   specify only a subset of the cells. In this way, you can use a single dimensionality
   reduction algorithm, but include multiple subsets and view of the cells,
   e.g. one coordinates file per tissue. Note, if R has changed your
   cell identifiers (e.g. by adding dots), you may be able to fix them by running ``cbTool metaCat``.

4. **Cluster-specific marker genes** (optional). The first column is
   the cluster name (from the cell annotation metadata file), the second column 
   contains the gene symbol (or Ensembl gene ID, which will automatically be mapped
   to the gene symbol), and
   the third column is some numeric score (e.g.  p-Value or FDR). You can add
   as many other columns as you like with additional information about this
   gene. You can also run ``cbMarkerAnnotate`` on this file to add information from
   various gene-centric databases and link-outs to other resources to your existing
   table. See :ref:`Annotate genes` on how to add link to external
   gene-databases (like Allan Brain Atlas or OMIM) to your marker genes.
   If you used Seurat for your clustering, you can just provide the raw Seurat marker gene output.
   You can also specify multiple files of cluster-specific marker genes,
   e.g. in case that you are also doing differential gene expression analysis
   or have results from multiple algorithms. 


Make sure that all your input files have Unix line endings and fix the line
endings if necessary with mac2unix or dos2unix::

    file *.txt *.csv *.tsv *.tab

Setting up your ``cellbrowser.conf``
^^^^^

After you have all of these files in place, go to the directory go to the
directory containing all of these files and run the following command to
copy a sample ``cellbrowser.conf`` into your current directory::

    cbBuild --init

This sample ``cellbrowser.conf`` includes the required settings plus some other useful settings. 
The current list of all possible cellbrowser.conf statements can be found in our `example cellbrowser.conf <https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf>`_.


In your cellbrowser.conf, replace the default values in the config statements:

* ``exprMatrix`` - expression matrix file name
* ``meta`` - cell annotation metadata file name
* ``coords`` - coordinate file names with a layout method label for each
* ``markers`` - cluster-specific marker gene file namew with a label for each
* ``labelField`` and ``clusterField`` - name of cluster field from header line of metadata file

From the directory where your ``cellbrowser.conf`` is located, run::

    cbBuild -o /tmp/cb -p 8888

Point your internet browser to the name of the server (or localhost, if
you're running this on your own machine) followed by :8888, e.g.
http://localhost:8888.

The cell browser output directory (/tmp/cb in this example) can hold multiple datasets. 
If you have a second dataset in another directory that contains cellbrowser.conf,
just make sure that the other cellbrowser.conf specifies a different dataset name 
with `name=xxx`. Then run `cbBuild -o /tmp/cb -p 8888` in the other
directory to add the second dataset to your cell browser output directory /tmp/cb.
