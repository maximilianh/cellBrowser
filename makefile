aparna:
	mkdir ../cbData -p
	tar cvzfh ~/public_html/cellBrowser/datasets/aparna.tgz ../cbData/aparna/{dataset.conf,exprMatrix.tsv,meta.tsv,coords.tsv,acronyms.tsv,quickGenes.tsv}

get-aparna:
	cd ../cbData/ && wget http://hgwdev.soe.ucsc.edu/~max/cellBrowser/datasets/aparna.tgz -O - | tar xvz
