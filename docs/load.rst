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
    mat <- fread("exprMatrix.tsv.gz")
    meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
    genes = mat[,1][[1]]
    genes = gsub(".+[|]", "", genes)
    mat = data.frame(mat[,-1], row.names=genes)
    so <- CreateSeuratObject(counts = mat, project = "adultPancreas", meta.data=meta)

Or you can download directly into R, without wget, by replacing the fread and read.table commands above in line 4 and 5 with these::

    mat <- fread("https://cells.ucsc.edu/adultPancreas/exprMatrix.tsv.gz")
    meta <- data.frame(fread("https://cells.ucsc.edu/adultPancreas/meta.tsv"), row.names=1)

If your version of data.tables does not support .gz yet, the fread commands can be changed to this::
 
    # from current directory
    mat <- fread("zcat < exprMatrix.tsv.gz")
    # or direct download:
    mat <- fread("curl https://cells.ucsc.edu/adultPancreas/exprMatrix.tsv.gz | zcat")

If the matrix name is not ``exprMatrix.tsv.gz`` but ``matrix.mtx``, you have to
use Seurat's MTX loader.  In addition to ``matrix.mtx``, make sure to also
download the files ``barcodes.tsv`` and ``genes.tsv`` sometimes
called ``features.tsv``.  If you downloaded these three files and ``meta.tsv`` into a directory ``downloadDir``, 
load them like this::

    require(Seurat)
    setwd("downloadDir")
    mat = Read10X(".")
    meta = read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
    so <- CreateSeuratObject(counts = mat, project = "myProjectName", meta.data=meta)
    
Scanpy
^^^^^^

To create an anndata object in Scanpy if the expression matrix is a .tsv.gz file::

    import scanpy as sc
    import pandas as pd
    ad = sc.read_text("exprMatrix.tsv.gz")
    meta = pd.read_csv("meta.tsv", sep="\t")
    ad.var = meta

If the expression matrix is an MTX file::

    import scanpy as sc
    import pandas as pd
    ad = sc.read_mtx("matrix.mtx.gz")
    meta = pd.read_csv("meta.tsv", sep="\t")
    ad.var = meta

Some datasets use the format identifier|symbol for the ad.obs gene names (e.g. "ENSG0123123.3|HOX3"). To keep only the symbol:

    ad.obs.index = [x.split("|")[1] for x in ad.obs.index.tolist()]
