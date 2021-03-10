Setup your own
--------------

Overview
^^^^^^^^

The UCSC Cell Browser tool set consists of a number of different scripts to help you set up your own. 
The primary utility being the Python script ``cbBuild`` that will import a set of existing single-cell
data from a directory of tab-separated files and configuration files to generate a
directory of html, json, and css files that can be viewed on the web. The rest of the utilities will
produce output than can be fed directly into ``cbBuild``.

The utilities ``cbSeurat`` and ``cbScanpy`` run a very basic single-cell pipeline on your expression
matrix and will output all the files needed to create a cell browser visualization. 
The ``cbImport*`` (``cbImportCellranger``, ``cbImportScanpy``, etc.) tools convert files
produced by Cellranger, Seurat, and Scanpy into a set of files that you can create a Cell
Browser visualization from. Both the pipeline and import tools are covered in more detail
under their respective sections (With Scanpy, With Seurat, and With Cellranger).
There is also a collection of small tools (``cbTool``) to combine cell annotation files
from different pipelines or convert expression matrices.

Using cbBuild to set up a Cell Browser
^^^^^^^^^^^

The main utility for building your own cell browser is ``cbBuild``. It takes in a gene expression
matrix and a set related files and converts them JSON and binary files outputting them to directory which
can be put onto a web server or used with the built-in webserver. At this time, there is no backend
server needed for a cell browser. You can place the output of  ``cbBuild`` on any static web server at your University
or the ones you can rent from companies will do.

After the installation, you should be able to run the cbBuild command and see
the usage message::

    cbBuild

Example Minimal Cell Browser
""""""

Below are  some instructions to set up a cell browser using a small example dataset based on data from 
`Nowakowski et al. 2017. <https://science.sciencemag.org/content/358/6368/1318.long>`_ and
the cortex-dev dataset on `cells.ucsc.edu <http://cells.ucsc.edu/?ds=cortex-dev>`_. The
expression matrix only includes 100 genes, but it does show off many of the
features of the cell browser. 

First, download and extract it to the directory ``mini`` with::

    curl -ks https://cells.ucsc.edu/downloads/samples/mini.tgz | tar xvz

Next, build a browser consisting of html and other files into the directory
~/public_html/cells/ and serve that directory on port 8888::

    cd mini
    cbBuild -o ~/public_html/cells/ -p 8888

Lastly, point your web browser to ``http://localhost:8888`` to view your minimal cell browser. If you're running
this on a server and not your own computer, replace localhost with the address
of your server. To stop the cbBuild web server, press Ctrl-C. To keep it running in the background, 
press Ctrl-Z and put it into the background with ``bg``. If you have stopped the web server, you
can always run the same ``cbBuild`` command to restart it. Restarting the web server will not re-export 
the entire expression matrix again if there is already one under
``~/public_html/cells/sample``. 

The optional to specify the port, ``-p PORT``, is optional. If you only want to build html files and serve them with your own
web server, do not specify this option and ``cbBuild`` will only build the output files, but won't start a web server.

The example `cellbrowser.conf <https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf>`_
explains all the various settings that are available in this config file. Things 
you can change include the colors for different metadata attributes, explain cluster acronyms used in your cluster names,
add file names, add alternative dimensionality reduction layouts, add more marker gene tables, and more. 

One of the most important settings in cellbrowser.conf is the dataset name. For example, in this 
'mini' example, the dataset name is 'sample'. When you run ``cbBuild``, its output 
files will be written to ~/public_html/cells/sample. You can go to another directory
with a different cellbrowser.conf file and a different dataset name, and if you run the same cbBuild
command as above, the cell browser output files will be copied into a new subdirectory within ~/public_html/cells/. 
A single cbBuild output directory can contain multiple datasets. 
