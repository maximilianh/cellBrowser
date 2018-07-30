HPO - human mutant gene phenotypes:
file: hpo_frequent_7Dec17.txt
Downloaded from wget http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_FREQUENT_FEATURES_diseases_to_genes_to_phenotypes.txt -O hpo_frequent_7Dec17.txt
example link: http://compbio.charite.de/hpoweb/showterm?gene=5308 <entrez ID>

OMIM - human rare disease gene:
file: mim2gene.txt
Downloaded from: wget https://omim.org/static/omim/data/mim2gene.txt
Example link: https://omim.org/entry/601542?search=pitx2

Cosmic - cancer genes census:
file: Census_allWed Dec  6 18_35_54 2017.tsv
Downloaded from http://cancer.sanger.ac.uk/census
Example link: http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=ACVR1

HPRD - human protein database:
file: HPRD_molecular_class_081914.txt
Downloaded from http://hprd.org/download

SFARI - autism genes:
file: SFARI-Gene_genes_export06-12-2017.csv
Downloaded from https://gene.sfari.org/tools/
example link: https://gene.sfari.org/database/human-gene/OTX1

HGNC - human gene names:
file: hgnc_complete_set_05Dec17.txt
Downloaded from https://www.genenames.org/cgi-bin/statistics
example link: https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:9005

MGI Mouse/Human homologs:
file: mgi_HGNC_homologene_8Dec17.txt
Download from: http://www.informatics.jax.org/downloads/reports/index.html#homology
Command: wget http://www.informatics.jax.org/downloads/reports/HGNC_homologene.rpt -O mgi_HGNC_homologene_8Dec17.txt
no link, this is only the ortholog table, used for brainspan/Eurexpress human/mouse mapping

Brainspan LMD - human developmental brain laser microdissection microarray:
file: brainspan_genes.csv
Downloaded from: http://www.brainspan.org/static/download.html
Command: wget http://www.brainspan.org/api/v2/well_known_file_download/267666527 -O b.zip
unzip b.zip
mv rows_metadata.csv static/brainspan_genes.csv
example link: http://www.brainspan.org/lcm/search?exact_match=false&search_term=51147&search_type=gene

Brainspan Mouse Development - mouse developmental ISH
file: brainspanMouse_9Dec17.txt
Command: wget 'http://api.brain-map.org/api/v2/data/Gene/query.xml?criteria=products[id$eq3]&num_rows=all' -O - | tr -d ' ' | egrep '<entrez-id|<id>' | awk '/entrez/ { split($1, a, /[<>]/); eId=a[3];} /<id>/ { split($1, a, /[<>]/); id=a[3]; print eId, id}' > brainspanMouse_9Dec17.txt
columns: entrez, brainspanId
example link: http://developingmouse.brain-map.org/gene/show/18505

Eurexpress - whole mount mouse in-situs:
file: eurexpress-7Dec17.txt
Download from: http://www.eurexpress.org/ee/databases/biomart.jsp
example link: http://www.eurexpress.org/ee/databases/assay.jsp?assayID=euxassay_019559

Gencode versions:
gencode22.ens79-80.tab
gencode23.ens81-82.tab
gencode24.ens83-84.tab


