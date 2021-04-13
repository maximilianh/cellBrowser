# --------- REQUIRED SETTINGS --------------

# example config file with all possible settings

# internal short name, only visible in the URL
# same as the output directory name
# no special chars, no whitespace, please
name = "sample"

# priority determines the order of the datasets
# smallest comes first
priority = 10

# tags are shown in the dataset browser
# current tags:
# smartseq2,10x
tags = ["smartseq2"]

# human-readable name of this dataset
shortLabel="CellBrowser 100-genes demo"

# name of the expression matrix file, genes are rows
exprMatrix="exprMatrix.tsv.gz"

# "gencode-human", "gencode-mouse" or "symbol"
# For "symbol" you can specify which database to use to check
# symbols or, for cbHub, how to map them to the genome.
# 'auto' will automatically detect Ensembl human/mouse  IDs
# and translate to symbols
geneIdType="auto"

# name of the meta data table ("samplesheet). One sample per row. First row is name of sample.
meta="meta.tsv"

# we try to auto-detect the field type of fields in the meta data.
# Sometimes, this doesn't work, e.g. when your cluster ID is a numer
# or your C1 chip ID is a number, but you don't want them binned, you want
# to treat as if they were categories
enumFields = ["c1_cell_id"]

# tsv files with coordinates of every sample in format <sampleId, x, y>
# first the name of the file, then a human readable description
coords=[
    {
            "file":"tsne.coords.tsv", 
            "flipY" : False, # R/Matplotlib files need to be flipped on the Y-axis
            "shortLabel":"t-SNE on WGCNA"
    },
    {
            "file":"subset.coords.tsv", 
            "shortLabel":"neural cells", 
            # you can overlay lines onto the cells, table has to have columns named x1, x2, y1, y2
            "lineFile" : "lines.tsv",
            # you can flip the y-axis of just the lines, relative to the points
            # This was necessary for a user when using the files produced by the URD pseudotime package
            #"lineFlipY" : True,
            # you can automatically switch on coloring on a meta data field whenever a layout is activated
            "colorOnMeta":"neuralCluster"
    },
]

# default field in the meta data table with the name of the cluster
clusterField="WGCNAcluster"

# default field in the meta data table used for the label of the clusters shown by default
labelField="WGCNAcluster"

# --------- OPTIONAL SETTINGS --------------

# genes that are highlighted in your paper can be pre-loaded and are shown as a clickable table on the left
# this is optional but we highly recommend that you define at least 2-3 quick genes, it makes the browser a lot
# more intuitive for users
quickGenesFile = "quickGenes.csv"

# if you want to enforce some order of the values of your enums, e.g. your cluster annotation should be sorted
# in a given order in the display, supply a text file with the values in the right order, one per line.
# You can supply one text file per meta data field
# enumOrder = { "WGCNAcluster" : "clusterorder.txt" }

# tsv files with marker gene lists for the clusters 
# format is (clusterName, geneSymbol, pValue, enrichment) + any additional fields or URLs you want to show
markers=[
    {"file":"markers.tsv", "shortLabel":"Cluster-specific markers"}
]

# do not show this dataset on the dataset list. This can be used for pre-publication data.
# visibility="hide"

# optional: UCSC track hub with the BAM file reads and expression values 
# Alternatively, you can also provide a full link to a UCSC Genome Browser session here
hubUrl="http://cells.ucsc.edu/cortex-dev/hub/hub.txt"

# optional: table with <name><color> for any meta data values
# color is a six-digit hexcode
# name is a any value in the meta data table, e.g. cluster name. Canb be a .tsv or .csv file.
colors="colors.tsv"

# should the cluster labels be shown by default (default: true)
showLabels=True

# the radius of the circles. If not specified, reasonable defaults will be used
#radius = 5
# the alpha/transparency of the circles. If not specified, reasonable defaults will be used.
#alpha = 0.3

# you need short names for your clusters, as there is little space on the plot
# but cell types have complicated and long names
# So you can provide a table with two columns: 1) short cluster name 2) long version
# e.g. EC, endothelial cells
# can be a .tsv or .csv file
acronymFile = "acronyms.tsv"

# the unit of the values in the expression matrix
# any string, shown on genome browser and violin y-Axis
# typical values are: "read count/UMI", "log of read count/UMI", "TPM", "log of TPM", "CPM", "FPKM", "RPKM"
unit = "TPM"   

# format of the numbers in the matrix. 
# 'auto' works in 99% of the cases. Otherwise you can use 'int' for integers and 'float' for  floating point numbers. 
# Use 'forceInt' if your matrix contains only integers but in a format like 3.123e10 
# or the matrix has only integers expressed like 100.000, 200.000, 300.00, ...
matrixType='auto'

# rarely needed: if your expression matrix does not contain genes, but something
# else, like "lipids" or "plankton", you can replace the word "gene" in the
# user interface with another word
# geneLabel = "Lipid"

# the default color palettes for this dataset. By default, we use Paul Tol's
# but you could use other ones, see the URL when you change the palette to see possible values
# defQuantPal = "viridis"
# defCatPal = "rainbow"

# you can optionally show little images for clusters on the tooltip. 
# For now, they have to be PNGs.
# For now, you will have to copy these images to the source destination htdocs directory manually
# right now, only brain-lipids/all-lipids is using this
# clusterPngDir = "clusterImgs"
