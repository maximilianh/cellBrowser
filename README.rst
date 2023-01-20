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
at http://cellbrowser.rtfd.io. If you use the UCSC Cell Browser in your research, please cite
`our Bioinformatics paper <https://dx.doi.org/10.1093/bioinformatics/btab503>`_.
If you are also using data from a specific dataset we host, please also cite
the original authors of that dataset (visible under 'Info & Download' while viewing that dataset).

If you want us to add a single cell dataset to the website http://cells.ucsc.edu, 
please contact us at cells@ucsc.edu. We are happy to add any dataset.

This is a viewer for a static, precomputed layout. If you're looking for an interative layout, where you can 
move the cells around and run some algorithms interactively, try Chan-Zuckerberg's own cellxgene or Spring.
A website with both datasets and some analysis is `Scope <http://scope.aertslab.org/>`_.

Many labs host their data at cells.ucsc.edu by sending it to us, but some groups have setup their own cell browsers:

* Alexander Misharin Lab, Northwester University, https://www.nupulmonary.org/resources/
* Accelerating Medicine Partnership Consortium, https://immunogenomics.io/cellbrowser/, used in `Zhang et al. 2018 <https://www.biorxiv.org/content/10.1101/351130v1>`_ and `Der et al 2018 <https://www.biorxiv.org/content/10.1101/382846v1>`_
* Ovarian Cancer Cell Laboratory, Nuffield Department of Women's & Reproductive Health, University of Oxford, https://ovariancancercell.github.io/
* Zemans Lab, Ann Arbor, https://rnabioco.github.io/lung-scrna/
* Unpublished work, unknown lab?, http://covid19ocularsurface.org/
* Chinese University of Hong Kong, Testis cis-element atlas, http://testisatlas.s3-website-us-west-2.amazonaws.com/CB.html
* Sansom Lab, Oxford, https://sansomlab.github.io (they use an older version of the Cell Browser source code, reached out to update with our bugfix but did not get reply), for `Croft et al, Nature 2019 <https://www.nature.com/articles/s41586-019-1263-7>`_ 
* https://www.genomique.eu/cellbrowser/HCA/ Zaragosi group at IPMC CNRS Nice, for the manuscript https://dev.biologists.org/content/early/2019/09/25/dev.177428.abstract
* Same group seems to host a COVID-19 dataset: https://www.genomique.eu/cellbrowser/COVID/
* Bin Ren lab, CAAtlas http://catlas.org/mousebrain/#!/home
* Conrad lab at Charite Berlin: http://singlecell.charite.de/
* `STAB: a spatio-temporal cell atlas of the human brain <https://stab.comp-sysbio.org/tool/cellbrowser/index.html>`_ from  `Song et al NAR 2021 <https://academic.oup.com/nar/article/49/D1/D1029/5911746>`_.
* UCLA: http://mergeomics.research.idre.ucla.edu/PVDSingleCell/
* Lako Lab at Newcastle University, UK: http://retinalstemcellresearch.co.uk/CorneaCellAtlas/ from `Collins et al. 2021. The Ocular Surface. <https://www.sciencedirect.com/science/article/pii/S1542012421000215>`_
* RNA Bioscience Initiative: https://www.pneuroonccellatlas.org/ and https://github.com/rnabioco/lung-scrna
* Paul Gontarz, WUSTL, http://regmedsrv1.wustl.edu/Public_SPACE/pgontarz/Public_html/cellbrower/Exp1/


These papers have cell browsers made at UCSC:

* organoidatlas: https://www.sciencedirect.com/science/article/pii/S221112472030053X
* dros-brain: https://elifesciences.org/articles/50354
* kidney-atlas: https://science.sciencemag.org/content/365/6460/1461.abstract
* allen-celltypes/mouse-cortex: https://www.biorxiv.org/content/10.1101/2020.03.30.015214v1.full
* organoidreportcard: https://www.nature.com/articles/s41586-020-1962-0

Before judging this project by the number of issue tickets or PRs, note that at UCSC we use an internal
ticket system with more features and that a lot of communication with wetlab users is by email at cells@ucsc.edu, as we 
do not require a Github account for feedback. But we do reply to issues here, as you can see from the Github 
account and also use Github for source control.

Additional availability
-----------------------

* The preferred installation is via pip https://pypi.org/project/cellbrowser/, for documentation see https://cellbrowser.readthedocs.io
* Bioconda: this tool is available to install via `bioconda <https://bioconda.github.io/recipes/ucsc-cell-browser/README.html>`_. Note that the conda release is usually a bit outdated relative to the pip release, so use pip if possible. If you cannot use pip, please contact us. 
* Biocontainers: there is a biocontainer automatically generated from the bioconda package available `here <https://quay.io/repository/biocontainers/ucsc-cell-browser>`_
* The Seurat3Wizard, demo at http://nasqar.abudhabi.nyu.edu/SeuratV3Wizard, builds a cell browser as its last step
* Galaxy: there is a Galaxy tool for UCSC CellBrowser, which can be installed on any Galaxy instance via its `Galaxy Toolshed entry <https://toolshed.g2.bx.psu.edu/view/ebi-gxa/ucsc_cell_browser>`_ or it can be directly used by users at the `Human Cell Atlas Galaxy instance <https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/ebi-gxa/ucsc_cell_browser/ucsc_cell_browser>`_ or as part of the example workflows, such as the `Human Cell Atlas / Scanpy CellBrowser workflow <https://humancellatlas.usegalaxy.eu/u/pmoreno/w/humancellatlas-scanpy-cellbrowser>`_ or the `EBI Single Cell Expression Atlas / Scanpy / CellBrowser workflow <https://humancellatlas.usegalaxy.eu/u/pmoreno/w/atlas-scanpy-cellbrowser-imported-from-uploaded-file>`_

This project was funded by the California Institute of Regenerative Medicine and the
Chan-Zuckerberg Initiative https://www.chanzuckerberg.com/. In 2020, it is funded through a supplement to the NHGRI Genome Browser grant.

This is early research software. You are likely to find bugs. Please open a Github
ticket or email us at cells@ucsc.edu, we can usually fix them quickly.

Citation
--------

If you use the UCSC Cell Browser in your work, please cite `Speir et al, Biorxiv 2020 <https://www.biorxiv.org/content/10.1101/2020.10.30.361162v1>`_ 
