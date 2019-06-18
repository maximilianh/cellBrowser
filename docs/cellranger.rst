With Cell Ranger
===============

Find the ``cellranger`` **OUT** directory, it should contain an ``analysis`` directory and
a subdirectory ``filtered_gene_bc_matrices``. The **OUT**
directory will be the input for our tool ``cbImportCellranger``. The tool converts the
``cellranger`` files to ones formatted for ``cbBuild``.

As we are reading Cell Ranger ``mtx`` files, we need the scipy package (add ``--user``
if you are not the admin on your machine)::

    pip install scipy

The example below use the pbmc3k ``cellranger`` output files from the 10x website.
First, download the files with the command::

    rsync -Lavzp genome-test.gi.ucsc.edu::cells/datasets/pbmc3kCellranger/ ./pbmc3kCellranger/ --progress

Next, run ``cbImportCellranger`` to convert the ``cellranger`` files into something
that can be used to build a cell browser::

        cbImportCellranger -i pbmc3kCellranger -o cellrangerOut --name pbmc3k_cellranger

The ``-i`` option specifies the input ``cellranger`` directory and the ``-o`` option
specifies a name for the output directory. You can use the ``-n`` option to change the
dataset name in the cell browser; if it is not specified, it will default to the output
directory name.

Lastly, go into the ``cellrangerOut`` directory and run ``cbBuild`` to create a cell browser::
    
    cd cellrangerOut
    cbBuild -o ~/public_html/cells -p 9999
