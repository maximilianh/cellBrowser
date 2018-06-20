This dataset is from Darmanis et al, PNAS 2015, "A survey of human brain transcriptome diversity at the
single cell level.", https://www.ncbi.nlm.nih.gov/pubmed/26060301

Fastq was downloaded from the SRA and run through kallisto.

Meta-data was downloaded from cirm.ucsc.edu in tagStorm format (fullMeta.tags), converted to
.tsv with tagStormToTab (fullMeta.tsv) and only human-readable columns were written to 
meta.tsv.

A viewer was then created with these commands:

    ../../cbPrep matrix -e geneMatrix.tsv --log2 --skip -m meta.tsv -o ~/public_html/cbTest/ -g markerSymbols.txt -l biosample_cell_type

    ../../cbPrep html -o ~/public_html/cbTest
