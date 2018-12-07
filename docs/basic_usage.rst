Basic usage
----------

The UCSC Cell Browser tools consist mainly of the Python script ``cbBuild``
that imports existing single cell data from tab-separated files and generates a
directory of html files. Other tools convert from Cellranger, Seurat and Scanpy
into this format, merge cell annotation files or convert expression matrices.

The main script is ``cbBuild``. It takes a gene expression matrix and related files
and converts the output to JSON and binary files to an output directory which
can be put onto a webserver or used with the built-in webserver.

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
of the server. To stop the web server, press Ctrl-C.  You will have to run
``cbBuild`` with ``-p 8888`` to show it again.

The file ``cellbrowser.conf`` explains all the various settings that are available
in this config file. E.g. you can change the colors, add acronym tables, add
file names, add more marker gene tables, etc.
