Installation
------------

Installation with pip
^^^^^^^^^^^^^^^^^^^^^

To install the Cell Browser using pip, you will need Python2.5+ or Python3+ and pip. With these setup, on a Mac or any Linux system, simply run::

    sudo pip install cellbrowser

On Linux, if you are not allowed to run the sudo command, you can install the Cell Browser into your user home directory::

    pip install --user cellbrowser
    export PATH=$PATH:~/.local/bin

You can add the second command to your ~/.profile or ~/.bashrc, this will allow you
to run the Cell Browser commands without having to specify their location.
    
On OSX, if running ``sudo pip`` outputs *command not found*, you will need to setup pip first by running::

    sudo easy_install pip

Installation with conda
^^^^^^^^^^^^^^^^^^^^^^^

If you would prefer to install the Cell Browser through bioconda, you can run::

    conda install -c bioconda ucsc-cell-browser
    
There should be conda versions for release 0.4.23 onwards. The conda version is managed by
Pablo Moreno at the EBI and is often a few releases behind. Please indicate in any bug
reports if you used conda to install.

Installation with git clone
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pip is not required to install the Cell Browser. As an alternative to pip or conda, you can also git clone the repo and
run the command line scripts under cellbrowser/src::

    git clone https://github.com/maximilianh/cellBrowser.git --depth=10
    cd cellBrowser/src

Installation with just wget or curl
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You don't use pip, conda or git? You can also download the current master branch::

    wget https://github.com/maximilianh/cellBrowser/archive/master.zip
    unzip master.zip
    cellBrowser-master/src/cbBuild
