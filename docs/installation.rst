Installation
====

.. contents:: Table of Contents
   :depth: 1
   :backlinks: none

There are currently no required Pythom modules for the core Cell Browser script: ``cbBuild``. 
Installation instructions for a variety of methods are below, with pip being the 
recommended method.

Basic installation
----

With pip (recommended)
^^^^

To install the Cell Browser using pip, you will need Python2.5+ or Python3+ and pip. With these setup, on a Mac or any Linux system, simply run::

    sudo pip install cellbrowser

On Linux, if you are not allowed to run the sudo command, you can install the Cell Browser into your user home directory::

    pip install --user cellbrowser
    export PATH=$PATH:~/.local/bin

You can add the second command to your ~/.profile or ~/.bashrc, this will allow you
to run the Cell Browser commands without having to specify their location.
    
On OSX, if running ``sudo pip`` outputs *command not found*, you will need to setup pip first by running::

    sudo easy_install pip

With conda
^^^^

If you would prefer to install the Cell Browser through bioconda, you can run::

    conda install -c bioconda ucsc-cell-browser
    
There should be conda versions for release 0.4.23 onwards. Please indicate in any bug
reports if you used conda to install.

With git clone
^^^^

Pip is not required to install the Cell Browser. As an alternative to pip or conda, you can also git clone the repo and
run the command line scripts under cellbrowser/src::

    git clone https://github.com/maximilianh/cellBrowser.git --depth=10
    cd cellBrowser/src

With wget or curl
^^^^

You don't use pip, conda or git? You can also download the current master branch::

    wget https://github.com/maximilianh/cellBrowser/archive/master.zip
    unzip master.zip
    cellBrowser-master/src/cbBuild

Notes on Windows Installation
^^^^

First install the Windows Linux subsystem. 
Then open the Windows Linux Subsystem bash terminal and run these commands::

    sudo apt-get update
    sudo apt-get install python-pip
    sudo pip install cellbrowser

Optional modules
----

To take advantage of some of the advanced features or specialized scripts
such as ``cbScanpy`` or ``cbSeurat``, you will need to install some extra packages or tools. 

Image sizes
^^^^^^

To get the image sizes, cbBuild uses either the "file" command or the "identify" command (for JPEGs). 
You may have to install the ImageMagick package to get the identify command.

Custom Colors
^^^^

In your cellbrowser.conf you can specify a file with custom colors
for your metadata values. If this file contains HTML color names instead
of color codes, you have to install the module webcolors::

    pip install webcolors

Matrices in mtx format
^^^^^^

To read expression matrices in .mtx format, you have to install scipy::

    pip install scipy

``cbScanpy`` and ``cbSeurat``
^^^^^^

``cbScanpy`` requires that Scanpy is `installed <https://scanpy.readthedocs.io/en/latest/installation.html>`_. 

``cbSeurat`` requires that both R and `Seurat <https://satijalab.org/seurat/articles/install.html>`_ are installed .

