UCSC Cell Browser
-----------------

The UCSC Cell Browser is a fast, lightweight viewer for single-cell data.
Cells are presented along with metadata and gene expression, with the ability
to color cells by both of these attributes. Additional information, such as
cluster marker genes and selected dataset-relevant genes, can also be displayed
using the Cell Browser.

There is a UCSC Cell Browser website available at http://cells.ucsc.edu, which includes
a handful of datasets from repositories like HCA, CIRM, and GEO as well as user
contributed ones. We are happy to add your favorite dataset to this, you will just need
to send us the files or a link to where we can download them to cells@ucsc.edu. 
If you use the UCSC Cell Browser in your research, please cite
`our Bioinformatics paper <https://dx.doi.org/10.1093/bioinformatics/btab503>`_.
If you are also using data from a specific dataset we host, please also cite
the original authors of that dataset (visible under 'Info & Download' while viewing that dataset).

The documentation on this website describes how you can create a Cell Browser for
your own data and make it available through your own web server.

The UCSC cell browser is funded by grants from the `California Institute for Regenerative Medicine <https://www.cirm.ca.gov/>`_ and the
`Chan-Zuckerberg Initiative <https://www.chanzuckerberg.com/>`_.

To report issues or view the source code, see `GitHub <https://github.com/maximilianh/cellBrowser>`_.

This is early research software. You are likely to run into bugs.
If you do run into any trouble, please open a
`Github issue <https://github.com/maximilianh/cellBrowser/issues/new>`_
or email us at cells@ucsc.edu, we can usually fix them quickly.

.. toctree::
   :maxdepth: 1

   interface.rst
   submission.rst
   installation
   quick_start.rst
   basic_usage
   internet
   howto.rst
   tabsep
   seurat
   scanpy
   cellranger
   advanced
   markers
   combine
   cellbrowser_conf.rst
   dataDesc.rst
   collections.rst
   load.rst
   bulk.rst
   images.rst
