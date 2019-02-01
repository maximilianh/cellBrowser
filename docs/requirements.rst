Optional Python modules
-----------------------

There are currently no required Pythom modules for the basic Cell Browser script cbBuild.

In cellbrowser.conf you can specify a color file, the format is .tsv or .csv and it has two columns, clusterName<tab>colorCode. If this file contains html color names instead of color codes, you have to install the module webcolors:

    pip install webcolors

To read expression matrices in .mtx format, you have to install scipy:

    pip install scipy


