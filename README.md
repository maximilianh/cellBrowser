UCSC Single Cell Browser
========================

For a demo, see http://cells.ucsc.edu

This repo contains many different pieces in Python and Javascript, but the main script is
cbBuild. It is a Python script that takes a gene expression matrix and related files and
converts the output to JSON and binary files to the output webserver directory.

Requirements: Python2.6+ or Python3+

There is a sample dataset in sampleData/sample1, it's a minimal expression matrix for a few thousand cells and only the first 100 genes and a bit of meta data for the cells..

A viewer was created from it with these commands:

    cd sampleData/sample1/
    ../../src/cbBuild -o ~/public_html/cbTest -p 8888

Then point your web browser to http://localhost:8888

To deploy the result onto a real webserver, copy all files and directories under ~/public_html/cbTest to
an empty directory on a webserver and point your web browser to it.

=== Process an expression matrix ===

Requirements: python3 with Scanpy installed.

We provide a wrapper around Scanpy which runs filtering, PCA, nearest-neighbors, clustering, t-SNE and
UMAP and formats them for cbBuild. An example file is on our downloads server:

    mkdir ~/cellData
    cd ~/cellData
    rsync -Lavzp hgwdev.soe.ucsc.edu::cells/datasets/pbmc3k ./pbmc3k/ --progress
    ../../cellBrowser/src/cbScanpy -e filtered_gene_bc_matrices/hg19/matrix.mtx -o cbScanpyOut/ -n pbmc3k

=== Convert a Scanpy object ===

From Jupyter or Python3:

    sys.path.append("cellbrowser/src/cbLib")
    import cellbrowser
    scanpyToTsv(adata, "scanpyOut")

=== Convert a CellRanger directory ===

=== Convert a Seurat object ===
