aparna:
	cd ../cbData && tar cvzfh ~/public_html/cellBrowser/datasets/aparna.tgz aparna/{dataset.conf,exprMatrix.tsv,meta.tsv,coords.tsv,acronyms.tsv,quickGenes.tsv}

wget-aparna:
	cd ../cbData/ && wget http://hgwdev.soe.ucsc.edu/~max/cellBrowser/datasets/aparna.tgz -O - | tar xvz
