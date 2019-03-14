Installation
------------

Installation with pip
^^^^^^^^^^^^^^^^^^^^^

You need Python2.5+ or Python3+ and pip. On a Mac or any Linux, simply run::

    sudo pip install cellbrowser

On Linux, if you you're not allowed to run the sudo command, you can install into your user home directory::

    pip install --user cellbrowser
    export PATH=$PATH:~/.local/bin

You can add the second command to your ~/.profile or ~/.bashrc.
    
On OSX, if ``sudo pip`` says *command not found*, you need to setup pip first::

    sudo easy_install pip

Installation with conda
^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, if you prefer to install through bioconda, since 0.4.23 you can do::

    conda install -c bioconda ucsc-cell-browser
    
The conda version is managed by Pablo Moreno at the EBI and is often a few releases behind. Please 
indicate in any bug reports if you used conda to install.

Installation with git clone
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pip is not required. As an alternative to pip or conda, you can also git clone the repo and
run the command line scripts under cellbrowser/src::

    git clone https://github.com/maximilianh/cellBrowser.git --depth=10
    cd cellBrowser/src

Installation with just wget or curl
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You don't use pip, conda or git? You can also download the current master branch::

    wget https://github.com/maximilianh/cellBrowser/archive/master.zip
    unzip master.zip
    cellBrowser-master/src/cbBuild
