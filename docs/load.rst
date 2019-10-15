Loading a dataset
-----------------

Once you have found a dataset of interest on https://cells.ucsc.edu, it is
very easy to load it into your favorite analysis environment. (Let us know if 
something if we are missing one below.)

First, download the expression matrix and the meta data, usually in a Unix terminal::

    mkdir adultPancreas
    cd adultPancreas
    wget https://cells.ucsc.edu/adultPancreas/exprMatrix.tsv.gz
    wget https://cells.ucsc.edu/adultPancreas/meta.tsv

Replace "quakePancreas" above with the dataset name of interest, it is shown in
the URL when you open a dataset after "ds=" or in the download instructions or on the dataset
page as the "CellBrowser dataset identifier".

Then open your favorite tool (e.g. RStudio or Jupyter) and follow the instructions below.

Seurat
^^^^^^

Run these commands if you have downloaded the file as above::

    require(Seurat)
    require(data.table)
    setwd("adultPancreas")
    mat <- fread("zcat < exprMatrix.tsv.gz")
    meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
    genes = mat[,1][[1]]
    genes = gsub(".+[|]", "", genes)
    mat = data.frame(mat[,-1], row.names=genes)
    so <- CreateSeuratObject(counts = mat, project = "adultPancreas", meta.data=meta)

Or you can download directly into R, without wget, by replacing the fread commands with these::

    mat <- fread("curl https://cells.ucsc.edu/adultPancreas/exprMatrix.tsv.gz | zcat")
    meta <- data.frame(fread("https://cells.ucsc.edu/adultPancreas/meta.tsv"), row.names=1)

Scanpy
^^^^^^

To create an anndata object in Scanpy::

    import scanpy as sc
    import pandas as pd
    ad = sc.read_text("exprMatrix.tsv.gz")
    meta = pd.read_csv("meta.tsv", sep="\t")
    ad.var = meta

