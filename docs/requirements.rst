Optional Python modules
-----------------------

There are currently no required Pythom modules for the core Cell Browser script: ``cbBuild``.
However, to take advantage of some of the advanced features or specialized scripts
such as ``cbScanpy`` or ``cbSeurat``, you will need to install some extra packages or tools. 

Custom Colors
^^^^

In your cellbrowser.conf you can specify a color file, e.g. colors.tsv, with custom colors
for your metadata values. The file can be in tsv or csv format and it has two columns,
first metadataValue and then colorCode. If this file contains HTML color names instead
of color codes, you have to install the module webcolors::

    pip install webcolors

Image sizes
^^^^^^

To get the image sizes, cbBuild uses either the "file" command or the "identify" command (for JPEGs). 
You may have to install the ImageMagick package to get the identify command.

Matrices in mtx format
^^^^^^

To read expression matrices in .mtx format, you have to install scipy::

    pip install scipy

``cbScanpy`` and ``cbSeurat``
^^^^^^

``cbScanpy`` requires that Scanpy is installed. See the Scanpy documentation for `installation instructions <https://scanpy.readthedocs.io/en/latest/installation.html>`_. 

``cbSeurat`` requires that both R and Seurat are installed. See the Seurat website for 
`installation instructions <https://satijalab.org/seurat/install.html>`_. 
Note Conda can also be used to install `Seurat <https://anaconda.org/bioconda/r-seurat>`_
and `R <https://anaconda.org/r/r>`_. You can confirm that seurat is installed for R by
typing `Rscript` and looking for Seurat in the output.
