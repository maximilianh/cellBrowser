With CellRanger
===============

Find the cellranger OUT directory, it contains an *analysis* directory and also
a subdirectory *filtered_gene_bc_matrices*. The *OUT*
directory is the one for our tool ``cbImportCellranger``. The tool converts the
cellranger files to tab-separated files, then you can run cbBuild on these.

As we are reading Cellranger *mtx* files, we need the scipy package (add --user
if you are not admin on your machine)::

    pip install scipy

Let's use an example, the pbmc3k cellranger output files from the 10x website::

    rsync -Lavzp genome-test.gi.ucsc.edu::cells/datasets/pbmc3kCellranger/ ./pbmc3kCellranger/ --progress
    cbImportCellranger -i pbmc3kCellranger -o cellrangerOut --name pbmc3k_cellranger
    cd cellrangerOut
    cbBuild -o ~/public_html/cells -p 9999
