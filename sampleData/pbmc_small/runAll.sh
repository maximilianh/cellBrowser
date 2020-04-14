#!/bin/bash
set -e
#../../src/cbScanpy -e exprMatrix.tsv.gz -o scanpyout -n pbmcSmall 
#../../src/cbImportScanpy -i scanpyout/anndata.h5ad -o importScanpyOut
#../../src/cbSeurat -e exprMatrix.tsv.gz -o seuratOut -n pbmcSmallSeurat
../../src/cbImportSeurat -i seurat2.rds -o importSeuratOut
cd 5kPbmcMtx/ && cbTool mtx2tsv matrix.mtx.gz barcodes.tsv.gz features.tsv.gz exprMatrix.tsv.gz && cd ..
