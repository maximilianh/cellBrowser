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


On your local computer
====

If you are running the cell browser on your local computer, you will likely need to have cbBuild start up a webserver for you. 
You can use the ``-p PORT`` option to specify on which port the webserver will run::

    cbBuild -o ~/public_html/cells/ -p 8888

Pointint your web browser to ``http://localhost:8888`` to view your cell browser. To stop the cbBuild web server, 
press Ctrl-C. To keep it running in the background, press Ctrl-Z and put it into the background with ``bg``. 
If you have stopped the web server, you can always run the same ``cbBuild`` command to restart it. 
Restarting the web server will not re-export the entire expression matrix again if there is already one under
``~/public_html/cells/my-dataset``. 

On a webserver
====

If you are on a webserver, you likely only only want to build the cell browser html files into a web-accessible directory::

  cbBuild -o ~/public_html/cells

Specifying the port is optional is you are running this on a server that is already web-accessible. To view your cell browser, 
navigate to the address for your webserver. 

Customizing your cell browser
^^^^

The example `cellbrowser.conf <https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf>`_
explains all the various settings that are available in this config file. Things 
you can change include the colors for different metadata attributes, explain cluster acronyms used in your cluster names,
add file names, add alternative dimensionality reduction layouts, add more marker gene tables, and more. 

One of the most important settings in cellbrowser.conf is the dataset name. When you run ``cbBuild``,
the output will be written to OUTPUT_DIR/dataset-name. If you specify the same ``-o OUTPUT_DIR`` 
when running cbBuild for a different dataset with a different dataset name and cellbrowser.conf, 
the output will be put into a new subdirectory within ``OUTPUT_DIR``. A single cbBuild output directory 
can contain multiple datasets. 
