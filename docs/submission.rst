Submitting data to the UCSC Cell Browser
----

At this time, we are happy to host pretty much any single-cell dataset,
regardlessof the library prepartion (10x, Samrt-seq2, etc), organism 
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

* University or other insitutional server space
* GEO
* Dropbox (if shared via link, e.g. https://www.dropbox.com/s/NN/my_seurat.rds)

We do accept files via other methods, although they take a little more work for 
us to move to our server. These include methods such as:

* Box
* Google Drive
* Dropbox, other methods of sharing data

Other information we want
^^^^
Alongside your submission, it would be great if you filled out a 
`desc.conf <https://cellbrowser.readthedocs.io/dataDesc.html>`_ file. At
the very least, it should have the abstract, methods, and title filled out. 
However, you are welcome to fill out more fields and make it as complete as 
you would like.

FAQs
^^^^

Can I share the output of cbBuild with you?
""""

If you are going to share the output of one of our cbImport* tools, we prefer
the directory containing the cellbrowser.conf, desc.conf, etc. The output of 
cbBuild is optimized for web access and display, which makes it difficult if 
not impossible to make changes to the cell browser at a later date (e.g. 
correcting spelling mistakes). If you have access to the desc.conf, cellbrowser.conf, 
and other files, we can easliy make these changes and rebuild the cell browser
if needed. 

Can I keep my dataset private until a later date, but still accessible to reviewers?
""""
Yes, we offer limited methods for keeping datasets private. We can hide datasets from
being list alongside the others we host. This means that someone would need to know
the URL or dataset name to be able to access your dataset. For example, this means
that someone would need the URL cells.ucsc.edu/?ds=cortex-dev or know the name
(cortex-dev) to access the dataset.

Can you host single
