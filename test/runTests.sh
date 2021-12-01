#!/bin/bash
set -e
if [ ! -f pbmc_small/exprMatrix.tsv.gz ]; then
        mkdir pbmc_small
        wget https://cells-test.gi.ucsc.edu/downloads/test/pbmc_small/exprMatrix.tsv.gz -O pbmc_small/exprMatrix.tsv.gz
fi
../src/cbScanpy -e pbmc_small/exprMatrix.tsv.gz -o scanpyout -n pbmcSmall 
#../src/cbImportScanpy -i scanpyout/anndata.h5ad -o importScanpyOut
../src/cbSeurat -e pbmc_small/exprMatrix.tsv.gz -o seuratOut -n pbmcSmallSeurat
#../src/cbImportSeurat -i seurat2.rds -o cbImportSeurat
#cd 5kPbmcMtx/ && cbTool mtx2tsv matrix.mtx.gz barcodes.tsv.gz features.tsv.gz exprMatrix.tsv.gz && cd ..
