UCSC Cell Browser
-----------------

The UCSC Cell Browser is a viewer for single cell data. You can click on and
hover over cells to get meta information, search for genes to color on and
click clusters to show cluster-specific marker genes, which in turn are
clickable again.

The main live of the cell browser with various CIRM and HCA datasets is
http://cells.ucsc.edu. We are happy to add your favorite dataset to it, just send us a link or the files to cells@ucsc.edu. 
This documentation describes how you can setup one for
your own data and put the resulting files onto your own webserver.

The UCSC cell browser is funded by grants from the California Institute of Regenerative Medicine
and the `Chan-Zuckerberg Initiative <https://www.chanzuckerberg.com/>`_.

To report problems or look at the source code, see `GitHub <https://github.com/maximilianh/cellBrowser>`_.

This is early research software. You are likely to find bugs. Please open a Github
ticket or email us at cells@ucsc.edu, we can usually fix them quickly.

.. toctree::
   :maxdepth: 1

   installation
   basic_usage
   internet
   tabsep
   seurat
   scanpy
   cellranger
   advanced
   markers
   combine
   cellbrowser_conf.rst
   dataDesc.rst
   requirements
