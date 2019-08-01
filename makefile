aparna:
	mkdir ../cbData -p
	tar cvzfh ~/public_html/cellBrowser/datasets/aparna.tgz ../cbData/aparna/{dataset.conf,exprMatrix.tsv,meta.tsv,coords.tsv,acronyms.tsv,quickGenes.tsv}

get-aparna:
	cd ../cbData/ && wget http://hgwdev.soe.ucsc.edu/~max/cellBrowser/datasets/aparna.tgz -O - | tar xvz

static-data:
	rsync -avzp cellbrowserData/ /hive/data/inside/cells/cellbrowserData/

minisample:
	cd sampleData && tar cvfz /hive/data/inside/cells/samples/mini.tgz --hard-dereference -T miniFiles.txt 
pip:
	rm -rf build/*
	rm -rf dist/*
	rm -rf build/*
	python3 setup.py sdist
	#python3 setup.py sdist 
	twine upload dist/*

