.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
   :target: http://bioconda.github.io/recipes/ucsc-cell-browser/README.html

UCSC Single Cell Browser
========================

The UCSC Cell Browser is a viewer for single cell data. You can click on and
hover over cells to get meta information, search for genes to color on and
click clusters to show cluster-specific marker genes. 

To look at a list of selected single cell datasets, see http://cells.ucsc.edu

To setup your own cell browser, from Cellranger, Seurat, Scanpy or text files 
(tsv/csv), or just a single cell expression matrix, read the documentation
at http://cellbrowser.rtfd.io

If you want us to add a single cell dataset to the website http://cells.ucsc.edu, 
please contact us at cells@ucsc.edu. We are happy to add any dataset.

This is a viewer for a static, precomputed layout. If you're looking for an interative layout, where you can 
move the cells around and run algorithsm interactively, try Chan-Zuckerberg's own cellxgene or Spring.
Another viewer is `Scope <http://scope.aertslab.org/>`_.

Apart from labs at UCSC, UCSF and Gladstone who host their data at
cells.ucsc.edu, other groups have setup their own cell browsers:

* Alexander Misharin Lab, Northwester University, https://www.nupulmonary.org/resources/
* Accelerating Medicine Partnership Consortium, https://immunogenomics.io/cellbrowser/, used in `Zhang et al. 2018 <https://www.biorxiv.org/content/10.1101/351130v1>`_ and `Der et al 2018 <https://www.biorxiv.org/content/10.1101/382846v1>`_
* Sansom Lab, Oxford, https://sansomlab.github.io
* http://caire.ipmc.cnrs.fr/cellbrowser/Differentiation/ Zaragosi group at IPMC CNRS Nice, for the manuscript https://dev.biologists.org/content/early/2019/09/25/dev.177428.abstract

Additional availability
-----------------------

* The preferred installation is via pip https://pypi.org/project/cellbrowser/0.4.12/, see https://cellbrowser.readthedocs.io
* Bioconda: this tool is available to install via `bioconda <https://bioconda.github.io/recipes/ucsc-cell-browser/README.html>`_Note that the conda release is usually a bit outdated relative to the pip release. 
* Biocontainers: there is a biocontainer automatically generated from the bioconda package available `here <https://quay.io/repository/biocontainers/ucsc-cell-browser>`_
* The Seurat3Wizard, demo at http://nasqar.abudhabi.nyu.edu/SeuratV3Wizard, builds a cell browser as its last step
* Galaxy: there is a Galaxy tool for UCSC CellBrowser, which can be installed on any Galaxy instance via its `Galaxy Toolshed entry <https://toolshed.g2.bx.psu.edu/view/ebi-gxa/ucsc_cell_browser>`_ or it can be directly used by users at the `Human Cell Atlas Galaxy instance <https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/ebi-gxa/ucsc_cell_browser/ucsc_cell_browser>`_ or as part of the example workflows, such as the `Human Cell Atlas / Scanpy CellBrowser workflow <https://humancellatlas.usegalaxy.eu/u/pmoreno/w/humancellatlas-scanpy-cellbrowser>`_ or the `EBI Single Cell Expression Atlas / Scanpy / CellBrowser workflow <https://humancellatlas.usegalaxy.eu/u/pmoreno/w/atlas-scanpy-cellbrowser-imported-from-uploaded-file>`_


Funded by the California Institute of Regenerative Medicine and the
Chan-Zuckerberg Initiative https://www.chanzuckerberg.com/.

This is early research software. You are likely to find bugs. Please open a Github
ticket or email us at cells@ucsc.edu, we can usually fix them quickly.
