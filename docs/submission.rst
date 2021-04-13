Submitting data to the UCSC Cell Browser
----

At this time, we are happy to host pretty much any single-cell dataset,
regardless of the library prepartion (10x, Smart-seq2, etc), organism 
(human, mouse, zebrafish, etc), or analysis method (Seurat, Scanpy,
Monocle, etc).

A cell browser requires at minimum three things:

* Expression matrix
* Metadata with cell names and cluster field
* 2D Layout coordinates

If you can provide these alongside some description of the dataset,
we'll host it. Contact us at cells@ucsc.edu to get started.

Preparing and sharing your files
^^^^

Before we can make a cell browser for you, you have to share the data
with us. We accept the following file types:

* Seurat RDS, Rdata, or Robj files
* Scanpy h5ad or Loom files
* A collection of tsv or csv files
* The output directory of one of our cbImport* tools
  
After you have your data in one of the formats above, you will have to 
share the data with us.

We prefer data to be shared in a way that is easy for us to download
with something like wget, so the following methods are ideal:

* University or other institutional server space
* AWS
* GEO
* Dropbox (if shared via link, e.g. https://www.dropbox.com/s/NN/my_seurat.rds)

We do accept files via other methods, although they take a little more work for 
us to move to our server. These include methods such as:

* Box
* Google Drive
* Dropbox, other methods of sharing data

Other things we want from you
^^^^

Dataset description
""""

Alongside your submission, it would be great if you filled out a 
`desc.conf <https://cellbrowser.readthedocs.io/dataDesc.html>`_ file. At
the very least, it should have the abstract, methods, and title filled out. 
However, you are welcome to fill out more fields and make it as complete as 
you would like. You can run ``cbBuild --init`` to copy an example desc.conf
into your current directory or you can copy one from 
`our Github repo <https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/desc.conf>`_.

If you have a collection of datasets, please provide basic desc.conf for each of
the datasets in the collection. These don't have to be as comprehensive as the
top-level one for the collection itself. At minimum, they should include two tags:
(1) title and (2) abstract. However, feel free to fill out as many of the desc.conf
fields as you think is helpful for that dataset.

Dataset shortname
""""

It would be great if you could suggest a dataset shortname at the time of
your submission, although we're happy to make one up for you. An ideal short
name fulfills the following requirements

* All lowercase
* Words separated by dashes ("-")
* Four words or less (don't be afraid to abreviate words, e.g. development -> dev)
* Informative

A great example is cortex-dev - it's all lowercase, the two words are separated by 
dashes, it's short at only two words long, and informs you that the dataset is focused on 
cortex development. It fulfills all four points above. 

Other great examples:

* mouse-nervous-system
* skeletal-muscle
* mouse-oligo-het
* covid-hypertension

The short name doesn't have to be perfect, but good enough to communicate something
about your dataset in a few words. 

"Quick Genes"
""""

This is a list of 10-15 genes that you think are important to your dataset(s). In
addition to the list of gene symbols, it is great to have a word or two about why
it's in the list (e.g. "Fst, Paraxial Mesoderm"; "HES1, Fig1D"). If you have a
collection of datasets, you can have one set of genes for every dataset in the
collection or a different set for each. 

This list and descriptions should be in CSV or TSV format. 

Getting your URL
^^^^

After you submit your dataset to us, we will import the data and make a preliminary
version available on our development server. We will work with you to iterate and
make improvements to this version first. Once you give your final approval, we will
push the data to our main site, cells.ucsc.edu. Once there, you will recieve the 
final URL, e.g. cortex-dev.cells.ucsc.edu. This is the URL you should place in your
paper, link to from your lab website, tweet about, etc. Please **do not** put the
url to our development server in your paper, since it is under active development, 
we occasionaly break it.

FAQs
^^^^

Can I share the output of cbBuild with you?
""""

If you are going to share the output of one of our cbImport* tools, we prefer
the directory containing the cellbrowser.conf, desc.conf, etc. The output of 
cbBuild is optimized for web access and display, which makes it difficult if 
not impossible to make changes to the cell browser at a later date (e.g. 
correcting spelling mistakes). If we have access to the desc.conf, cellbrowser.conf, 
and other files, we can easliy make these changes and rebuild the cell browser
if needed. 

Can I keep my dataset private until a later date, but still accessible to reviewers?
""""
Yes, we offer limited methods for keeping datasets private. We can hide datasets from
being listed alongside the others we host. This means that someone would need to know
the URL or dataset name to be able to access your dataset. For example, this means
that someone would need the URL cells.ucsc.edu/?ds=cortex-dev or know the name
(cortex-dev) to access the dataset.
