Loading a dataset
-----------------

Once you have found a dataset of interest on https://cells.ucsc.edu, it is
very easy to load it into your favorite analysis environment. (Let us know if 
something if we are missing one below.)

First, download the expression matrix and the meta data, usually in a Unix terminal::

    wget https://cells.ucsc.edu/quakePancreas/exprMatrix.tsv.gz
    wget https://cells.ucsc.edu/quakePancreas/meta.tsv

Replace "quakePancreas" above with the dataset name of interest, it is shown in
the URL when you open a dataset after "ds=" or in the download instructions.

Then open your favorite tool (e.g. RStudio or Jupyter) and follow the instructions below.

Seurat
^^^^^^

Run these commands if you have downloaded the file as above::

    require(Seurat)
    require(data.tables)
    mat <- fread("zcat < exprMatrix.tsv.gz")
    # or: mat <- read.table(gzfile("exprMatrix.tsv.gz"), header = T)
    meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
    so <- CreateSeuratObject(counts = mat, project = "cellBrowserImport", meta.data=meta)

Or without downloading them first

    require(data.tables)
    mat <- fread("curl https://cells.ucsc.edu/adultPancreas/exprMatrix.tsv.gz | zcat")
    meta <- data.frame(fread("https://cells.ucsc.edu/adultPancreas/meta.tsv"), row.names=1)
    so <- CreateSeuratObject(counts = mat, project = "cellBrowserImport", meta.data=meta)


Scanpy
^^^^^^

    import scanpy as sc
    import pandas as pd
    ad = sc.read_text("exprMatrix.tsv.gz")
    meta = pd.read_csv("meta.tsv", sep="\t")
    ad.var = meta

