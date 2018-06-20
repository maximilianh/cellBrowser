aparna:
	tar cvzfh ~/public_html/cellBrowser/datasets/aparna.tgz ../cbData/aparna/{dataset.conf,exprMatrix.tsv,meta.tsv,coords.tsv,acronyms.tsv,quickGenes.tsv}

get-aparna:
	mkdir ../cbData -p
	cd ../cbData/ && wget http://hgwdev.soe.ucsc.edu/~max/cellBrowser/datasets/aparna.tgz -O - | tar xvz
