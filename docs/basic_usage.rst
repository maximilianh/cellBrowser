Basic usage
----------

The UCSC Cell Browser tools consist mainly of the Python script ``cbBuild``
that imports existing single cell data from a directory with tab-separated files and generates a
directory of html files. ``cbSeurat`` and ``cbScanpy`` run a very basic single cell pipeline
onto your expression matrix in various formats and create the directory of tab-separated files. 
The ``cbImport-`` tools convert files produced by Cellranger, Seurat and Scanpy
into a directory of tab-separated files (``cbImportCellranger``, ``cbImportScanpy``, etc.). 
There is also a collection of small tools (``cbTool``) to
combine cell annotation files from different pipelines or convert expression matrices.

The main script is ``cbBuild``. It takes a gene expression matrix and related files
and converts the output to JSON and binary files to an output directory which
can be put onto a webserver or used with the built-in webserver. There is no backend
server needed right now for the cell browser, any static webserver at your University
or the ones you can rent from companies will do.

After the installation, you should be able to run the cbBuild command and see
the usage message::

    cbBuild

Here is a small example dataset (Nowakowski et al 2018, fetal brains), this is the
the dataset cortex-dev on `cells.ucsc.edu <http://cells.ucsc.edu/?ds=cortex-dev>`_. The
expression matrix includes only the first 100 genes, otherwise quite a few
features of the browser are used. Download and extract it to the directory
``mini`` with::

    curl -s https://cells.ucsc.edu/downloads/samples/mini.tgz | tar xvz

You can now build a browser consisting of html and other files into the directory
~/public_html/cells/ and serve that directory on port 8888::

    cd mini
    cbBuild -o ~/public_html/cells/ -p 8888

Then point your web browser to ``http://localhost:8888`` or, if you're running
this on a server and not your own computer, replace localhost with the address
of the server. To stop the cbBuild web server, press Ctrl-C.  To keep it running, 
press Ctrl-Z and put it into the background with ``bg``. If you stopped it, you can always run
``cbBuild`` with the same arguments ``-p 8888`` to show it again, it will not re-export 
the whole expression matrix again, if there is one under
``~/public_html/cells/cortex-dev`` already. 

The ``-p PORT`` is option. If you only want to build html files and serve them with your own
webserver, do not specify this option and the cbBuild will only build html files.

Our example `cellbrowser.conf <https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf>`
explains all the various settings that are available in this config file. E.g.
you can change the colors, explain acronyms used in your cluster names,
add file names, add alternative dimensionality reductions, add more marker gene tables, etc. 

The most important setting in cellbrowser.conf is the name of the dataset. The
name of the mini example is 'sample' and this means that the converted expression
matrix and other files will be written to ~/public_html/cells/sample. This means that
you can go to another directory with a cellbrowser.conf file, and run the same cbBuild
command as above, there is no need to change the output directory of cbBuild, it can
contain multiple datasets.
