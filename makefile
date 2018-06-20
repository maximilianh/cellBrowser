aparna:
	tar cvzfh ~/public_html/cellBrowser/datasets/aparna.tgz data/aparna/{dataset.conf,exprMatrix.tsv,meta.tsv,coords.tsv,acronyms.tsv,quickGenes.tsv}

wget-aparna:
	cd data && wget http://hgwdev.soe.ucsc.edu/~max/cellBrowser/aparna.tgz -O - | tar xvz
