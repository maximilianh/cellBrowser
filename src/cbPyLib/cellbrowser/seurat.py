# a wrapper around the R library Seurat

import logging, optparse, sys, glob, os, datetime
from os.path import join, basename, dirname, isfile, isdir, relpath, abspath, getsize, getmtime, expanduser

from .cellbrowser import copyPkgFile, writeCellbrowserConf, pipeLog, makeDir, maybeLoadConfig, errAbort, popen
from .cellbrowser import setDebug, build, isDebugMode, generateHtmls

def parseArgs():
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] -e matrixFile -o outDir - run Seurat and output .tsv files

    If exceptions occur, will automatically start the debugger.
    """)

    parser.add_option("", "--init", dest="init", action="store_true",
        help="copy sample seurat.conf to current directory")

    parser.add_option("-e", "--exprMatrix", dest="exprMatrix", action="store",
            help="gene-cell expression matrix file, possible formats: .mtx, .txt.gz, .csv.gz, .rds. For .mtx, specify the directory where the matrix.mtx file is stored, e.g. filtered_gene_bc_matrices/hg19/")

    parser.add_option("-o", "--outDir", dest="outDir", action="store",
            help="output directory")

    parser.add_option("-c", "--confFname", dest="confFname", action="store", default="seurat.conf", 
            help="config file from which settings are read, default is %default")

    #parser.add_option("-s", "--samplesOnRows", dest="samplesOnRows", action="store_true",
            #help="when reading the expression matrix from a text file, assume that samples are on lines (default behavior is one-gene-per-line, one-sample-per-column)")

    parser.add_option("-n", "--name", dest="name", action="store",
            help="internal name of dataset in cell browser. No spaces or special characters.")

    parser.add_option("-d", "--debug", dest="debug", action="store_true", help="show debug output")

    (options, args) = parser.parse_args()

    if (options.exprMatrix is None and options.outDir is None) and not options.init:
        parser.print_help()
        exit(1)

    setDebug(options.debug)
    return args, options

def runRscript(scriptFname, logFname):
    """ run an R script given as a string through Rscript """
    logging.info("running %s through Rscript" % scriptFname)

    #ofh = open(logFname, "w")

    # removed subprocess. Multiprocessing is notoriously tricky to get right.
    cmd = "time Rscript %s 2>&1 | tee %s" % (scriptFname, logFname)
    #cmd = ["time","Rscript", scriptFname, "&"]
    #proc, stdout = popen(cmd, shell=True)
    #for line in stdout:
        #print(line),
        #ofh.write(line)
    err = os.system(cmd)
    #proc.stdout.close()
    #stat = os.waitpid(proc.pid, 0)
    #err = stat[1]
    assert(err==0)
    logging.info("Wrote logfile of R run to %s" % logFname)

def writeSeurat2Script(conf, inData, tsnePath, clusterPath, markerPath, rdsPath, matrixPath, scriptPath):
    " write the seurat R script to a file "
    cmds = []
    # try to install Seurat if not already installed
    cmds.append("ver = R.Version()")
    cmds.append("if (ver$major<3 || (ver$major==3 && ver$minor < 4)) { error('Sorry, Seurat requires at least R 3.4') } ")
    cmds.append("if (!require('Seurat',character.only = TRUE)) { install.packages(c('Seurat', 'data.table'), dep=TRUE, repos='http://cran.r-project.org/')}")
    cmds.append("library(methods)")
    cmds.append("suppressWarnings(suppressMessages(library(Seurat)))")
    cmds.append('print("Seurat: Reading data")')
    cmds.append('print("Loading input data matrix")')

    writeMatrix = False
    if isfile(inData):
        if inData.endswith(".mtx"):
            logging.error("You specified an .mtx file as the input matrix")
            logging.error("Seurat cannot read that. You need to specify the directory name, with a file matrix.mtx in it.")
            logging.error("The directory also has to contain two other files, genes.tsv and barcodes.tsv")
            logging.error("Please rename your input files if necessary")
            exit(1)

        if inData.endswith(".rds"):
            cmds.append('mat = readRDS(file = "%s")' % inData)
        else:
            if ".csv" in inData:
                cmds.append('mat = read.table("%s",sep=",",header=TRUE,row.names=1)' % inData)
            else:
                cmds.append('mat = read.delim("%s",header=TRUE,row.names=1)' % inData)
    else:
        matrixFname = join(inData, "matrix.mtx")
        if not isfile(matrixFname):
            errAbort("Could not find file %s - sure that you specified the directory of the matrix.mtx file?" % matrixFname)
        cmds.append('mat = Read10X(data.dir="%s")' % inData)

    undoLog = conf.get("undoLog", False)
    if undoLog:
        cmds.append('print("Undoing log2 on matrix, doing 2^mat")')
        cmds.append('mat <- 2^mat')

    minCells = conf.get("minCells", 3)
    minGenes = conf.get("minGenes", 200)

    cmds.append('print("Seurat: Setup")')
    cmds.append('sobj <- CreateSeuratObject(raw.data = mat, min.cells=%d, min.genes=%d)' % (minCells, minGenes))
    #cmds.append('sobj=Setup(nbt,project = "NBT",min.cells = 3,names.field = 2,names.delim = "_",min.genes = 500, do.logNormalize = F, total.expr = 1e4)')
    cmds.append('sobj') # print size of the matrix

    # export the matrix as a proper .tsv.gz
    cmds.append('print("Writing expression matrix to %s")' % matrixPath)
    cmds.append('dataFrame <- as.data.frame(as.matrix(mat))')
    # we MUST USE the raw matrix, so we use the mat object, not sobj, because otherwise
    # some of our markers won't even be in the final matrix. Very strange. Ask Andrew?
    #cmds.append('dataFrame <- as.data.frame(as.matrix(sobj@raw.data))') # raw counts, really not filtered?
    #cmds.append('dataFrame <- as.data.frame(as.matrix(sobj@data))') # log-normalized
    #cmds.append('dataFrame <- as.data.frame(as.matrix(sobj@scale.data))') #  scaled
    cmds.append('z <- gzfile("%s")' % matrixPath)
    cmds.append("write.table(dataFrame, z, quote=FALSE, sep='\t', eol='\n', col.names=NA, row.names=TRUE)")

    # find mito genes, mito-%, and create plots for it
    # XX - ENSG names?
    cmds.append('mito.genes <- grep(pattern = "^MT-", x = rownames(x = sobj@data), value = TRUE)')
    cmds.append('percent.mito <- Matrix::colSums(sobj@raw.data[mito.genes, ])/Matrix::colSums(sobj@raw.data)')
    cmds.append('sobj <- AddMetaData(object = sobj, metadata = percent.mito, col.name = "percent.mito")')
    cmds.append('VlnPlot(object = sobj, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)')
    cmds.append('par(mfrow = c(1, 2))')
    cmds.append('GenePlot(object = sobj, gene1 = "nUMI", gene2 = "percent.mito")')
    cmds.append('GenePlot(object = sobj, gene1 = "nUMI", gene2 = "nGene")')

    # remove cells with > mitoMax and a gene count < minGenes or > maxGenes
    mitoMax=conf.get("mitoMax", 0.05)
    maxGenes=conf.get("maxGenes", 2500)

    # filter based on mito%
    cmds.append('print("Removing cells with more than %f percent of mitochondrial genes")' % mitoMax)
    cmds.append('print("Keeping only cells with a geneCount %d-%d")' % (minGenes, maxGenes))
    cmds.append('sobj <- FilterCells(object = sobj, subset.names = c("nGene", "percent.mito"), low.thresholds = c(%(minGenes)d, -Inf), high.thresholds = c(%(maxGenes)d, %(mitoMax)f))' % locals())
    cmds.append('sobj') # print size of the matrix

    # do the log
    cmds.append('sobj <- NormalizeData(object = sobj, normalization.method = "LogNormalize", scale.factor = 10000)')

    # find most variable genes and gate on it
    minMean = conf.get("varMinMean", 0.0125)
    maxMean = conf.get("varMaxMean", 3)
    minDisp = conf.get("varMinDisp", 0.5)

    cmds.append('sobj <- FindVariableGenes(object = sobj, mean.function = ExpMean, dispersion.function = LogVMR, '
        'x.low.cutoff = %f, x.high.cutoff = %f, y.cutoff = %s)' % (minMean, maxMean, minDisp))
    #cmds.append('length(x = sobj@var.genes)')
    cmds.append('sobj') # print size of the matrix

    # scale
    cmds.append('sobj <- ScaleData(object = sobj, vars.to.regress = c("nUMI", "percent.mito"))')

    # PCA
    cmds.append('print("Running PCA")')
    pcCountConfig = conf.get("pcCount", 30)
    # see: https://github.com/satijalab/seurat/issues/982
    # pca = SeuratObj@dr$pca
    # eigValues = (pca@sdev) ^2 ## 
    # varExplained = eigValues / sum(eigValues)
    #if pcCountConfig == "auto":
        #cmds.append('sobj <- JackStraw(object = sobj, num.replicate = 100, display.progress = FALSE)')
        #cmds.append('JackStrawPlot(object = sobj, PCs = 1:12)')
        #cmds.append('sobj <- RunPCA(object = sobj, pcs.compute = 50, pc.genes = sobj@var.genes, do.print = FALSE, pcs.print = 1:5, genes.print = 5)')
    #else:
    cmds.append('sobj <- RunPCA(object = sobj, pcs.compute = %d, pc.genes = sobj@var.genes, do.print = FALSE, pcs.print = 1:5, genes.print = 5)' % pcCountConfig)
    cmds.append('print("PCA plots")')
    cmds.append('VizPCA(object = sobj, pcs.use = 1:2)')
    cmds.append('PCAPlot(object = sobj, dim.1 = 1, dim.2 = 2)')
    #cmds.append('sobj <- ProjectPCA(object = sobj, do.print = FALSE)')
    cmds.append('PCHeatmap(object = sobj, pc.use = 1:8, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)')

    cmds.append('PCElbowPlot(object = sobj)')

    #if pcCountConfig == "auto":
        #cmds.append('sobj <- JackStraw(object = sobj, num.replicate = 100, display.progress = FALSE)')
        #cmds.append('JackStrawPlot(object = sobj, PCs = 1:12)')
        #cmds.append('pcCount=min(ncol(sobj@dr$pca@pca.rot), 30)')
    #else:
        #cmds.append('pcCount=%d' % pcCountConfig)

    #cmds.append('print("Number of PCs used for Louvain Clustering:")')
    #cmds.append('print(pcCount)')

    louvainRes = conf.get("louvainRes", 0.6)
    cmds.append('print("Finding clusters with resolution %f")' % louvainRes)
    #cmds.append('sobj <- FindClusters(object = sobj, reduction.type = "pca", dims.use = 1:pcCount, resolution = %f, print.output = 0, save.SNN = TRUE)' % (louvainRes))
    cmds.append('sobj <- FindClusters(object = sobj, reduction.type = "pca", resolution = %f, print.output = 0, save.SNN = TRUE)' % (louvainRes))

    cmds.append("PrintFindClustersParams(object = sobj)")

    cmds.append('print("Running t-SNE")')
    #cmds.append("sobj <- RunTSNE(object = sobj, dims.use = 1:pcCount, do.fast = TRUE)")
    # "duplicate" = samples with identicals PC coordinates, more likely with big datasets
    cmds.append("sobj <- RunTSNE(object = sobj, do.fast = TRUE, check_duplicates=FALSE)")
    cmds.append("TSNEPlot(object = sobj, doLabel=T)")

    minMarkerPerc = conf.get("minMarkerPerc", 0.25)
    cmds.append('print("Finding markers")')
    cmds.append('all.markers <- FindAllMarkers(object = sobj, min.pct = %f, only.pos=TRUE, thresh.use=0.25)' % minMarkerPerc)

    cmds.append('print("Saving .rds to %s")' % rdsPath)
    cmds.append('saveRDS(sobj, file = "%s")' % rdsPath)

    # export the data to tsvs
    cmds.append('tsne12 <- FetchData(sobj, c("tSNE_1", "tSNE_2"))')
    cmds.append('write.table(tsne12, "%s", quote=FALSE, sep="\t", col.names=NA)' % tsnePath)
    cmds.append('clusters <- FetchData(sobj,c("ident"))')
    cmds.append('colnames(clusters) <- c("cellId\tCluster")')
    cmds.append('write.table(clusters, "%s", quote=FALSE, sep="\t")' % clusterPath)
    # marker file is the flag file for successful operation
    cmds.append('write.table(all.markers, "%s", quote=FALSE, sep="\t", col.names=NA)' % markerPath)

    cmds.append('PrintCalcParams(sobj, calculation = "RunPCA")')
    cmds.append('PrintCalcParams(sobj, calculation = "RunTSNE")')

    ofh = open(scriptPath, "w")
    for c in cmds:
        ofh.write(c+"\n")
    ofh.close()
    logging.info("Wrote R code to %s" % scriptPath)

def cbSeuratCli():
    global options
    args, options = parseArgs()

    if options.init:
        copyPkgFile("sampleConfig/seurat.conf")
        sys.exit(0)

    if options.name is None:
        errAbort("You need to specify the -n option and provide a name for this dataset for the viewer.")

    if options.exprMatrix is None:
        errAbort("You need to specify the -e option with the file name of the expression matrix.")

    if not isfile(options.confFname):
        logging.warn("%s not found, using default values for all Seurat options. Use --init to create a sample file." % options.confFname)

    inMatrix = options.exprMatrix
    outDir = options.outDir
    inConfFname = options.confFname
    datasetName = options.name

    makeDir(outDir)

    figDir = join(outDir, "figs")

    tsnePath = join(outDir, "tsne.coords.tsv")
    clusterPath = join(outDir, "meta.tsv")
    markerPath = join(outDir, "markers.tsv")
    scriptPath = join(outDir, "runSeurat.R")
    logPath = join(outDir, "analysisLog.txt")
    cbConfPath = join(outDir, "cellbrowser.conf")
    rdsPath = join(outDir, "seurat.rds")
    matrixPath = join(outDir, "exprMatrix.tsv.gz")

    inConf = maybeLoadConfig(inConfFname)
    writeSeurat2Script(inConf, inMatrix, tsnePath, clusterPath, markerPath, rdsPath, matrixPath, scriptPath)
    runRscript(scriptPath, logPath)

    if not isfile(markerPath):
        errAbort("R script did not complete successfully. Check %s and analysisLog.txt." % scriptPath)

    coords = [{'shortLabel':'t-SNE', 'file':'tsne.coords.tsv'}]
    writeCellbrowserConf(datasetName, coords, cbConfPath, args={"clusterField":"Cluster"})

    generateHtmls(datasetName, outDir)

def cbImportSeurat2_parseArgs(showHelp=False):
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] -i input.rds -o outDir [-n datasetName] - convert Seurat2 object to cellbrowser

    Example:
    - %prog -i pbmc3k.rds -o pbmc3kSeurat - convert pbmc3k to directory of tab-separated files
    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages. Also activates debug() in the generated R script.")

    parser.add_option("-i", "--inRds", dest="inRds", action="store",
        help="input .rds file. Required parameter")

    parser.add_option("-o", "--outDir", dest="outDir", action="store",
        help="Output directory. Required parameter")

    parser.add_option("-n", "--name", dest="datasetName", action="store",
        help="Dataset name for generated cellbrowser.conf. If not specified, the last component of -o will be used.")

    parser.add_option("", "--htmlDir", dest="htmlDir", action="store",
        help="do not only convert to tab-sep files but also run cbBuild to"
            "convert the data and add the dataset under htmlDir")

    parser.add_option("-p", "--port", dest="port", action="store", type="int",
            help="only with --htmlDir: start webserver on port to serve htmlDir")

    parser.add_option("-x", "--skipMatrix", dest="skipMatrix", action="store_true",
            default = False,
        help="do not convert the matrix, saves time if the same one has been exported before to the "
        "same outDir directory")

    parser.add_option("-m", "--skipMarkers", dest="skipMarkers", action="store_true",
            default = False,
        help="do not calculate cluster-specific markers with FindAllMarkers(), saves a lot of time")

    parser.add_option("-c", "--clusterField", dest="clusterField", action="store",
        help="Cluster field to color on, by default this is the @ident slot of the Seurat object but it can also be any other meta data field of the @meta.data slot")

    parser.add_option("", "--markerFile", dest="markerFile", action="store",
            help="Instead of calculating cluster markers again, use this file. Format: cluster,gene,pVal + any other fields.")

    (options, args) = parser.parse_args()

    if showHelp:
        parser.print_help()
        exit(1)

    setDebug(options.debug)
    return args, options

def readExportScript(cmds):
    " find the ExportToCellbrowser-seurat2.R script and add the commands between # --- to cmds "
    fname = join(dirname(__file__), "R", "ExportToCellbrowser-seurat2.R")

    blockFound = False
    for line in open(fname):
        if line=="# ---\n":
            blockFound = (not blockFound)
            continue
        if blockFound:
            cmds.append(line.rstrip("\n"))
    return cmds

def writeRScript(cmds, scriptPath, madeBy):
    " write list of R commands to script, prefix with a header "
    ofh = open(scriptPath, "w")
    ofh.write("# generated R code by: cellbrowser %s, %s, username: %s\n" %
            (madeBy, datetime.datetime.now().strftime("%I:%M%p %B %d, %Y"), os.getlogin()))
    for c in cmds:
        ofh.write(c)
        ofh.write("\n")
    ofh.close()

def cbImportSeurat2(rdsPath, outDir, datasetName, options):
    " convert Seurat2 .rds file to tab-sep directory for cellbrowser "
    logging.info("inFname: %s, outDir: %s, datasetName: %s" % (rdsPath, outDir, datasetName))

    makeDir(outDir)

    skipMatrix = options.skipMatrix
    skipMarkers = options.skipMarkers
    clusterField = options.clusterField
    markerFile = options.markerFile
    if markerFile is None:
        markerFileStr = "NULL"
    else:
        markerFileStr = '"%s"' % markerFile

    tsnePath = join(outDir, "tsne.coords.tsv")
    metaPath = join(outDir, "meta.tsv")
    markerPath = join(outDir, "markers.tsv")
    scriptPath = join(outDir, "runSeurat.R")
    matrixPath = join(outDir, "exprMatrix.tsv.gz")
    logPath = join(outDir, "analysisLog.txt")

    cmds = ["require(methods)"] # for the 'slots()' function
    cmds.append("require(Seurat)")

    cmds = readExportScript(cmds)

    cmds.append("message('Reading %s')" % rdsPath)
    cmds.append("sobj <- readRDS(file='%s')" % rdsPath)
    skipStr = str(skipMatrix).upper()
    skipMarkerStr = str(skipMarkers).upper()
    if isDebugMode():
        cmds.append("debug(ExportToCellbrowser)")
    if clusterField is None:
        clusterStr = 'NULL'
    else:
        clusterStr = "'%s'" % clusterField

    cmds.append("message('Exporting Seurat data to %s')" % outDir)
    cmds.append("ExportToCellbrowser(sobj, '%s', '%s', markers.file = %s, cluster.field=%s, skip.expr.matrix = %s, skip.markers = %s, all.meta=TRUE)" %
            (outDir, datasetName, markerFileStr, clusterStr, skipStr, skipMarkerStr))

    writeRScript(cmds, scriptPath, "cbImportSeurat2")
    runRscript(scriptPath, logPath)
    if not isfile(metaPath):
        errAbort("R script did not complete successfully. Check %s and analysisLog.txt." % scriptPath)

    #cbConfPath = join(outDir, "cellbrowser.conf")
    #coords = [{'shortLabel':'t-SNE', 'file':'tsne.coords.tsv'}]
    #writeCellbrowserConf(datasetName, coords, cbConfPath, args={"clusterField":"Cluster"})

    generateHtmls(datasetName, outDir)

def cbImportSeurat2Cli():
    " convert .rds to directory "
    args, options = cbImportSeurat2_parseArgs()

    inFname = options.inRds
    outDir = options.outDir
    if None in [inFname, outDir]:
        cbImportSeurat2_parseArgs(showHelp=True)
        errAbort("You need to specify at least an input rds file and an output directory")

    datasetName = options.datasetName
    if datasetName is None:
        datasetName = basename(outDir.rstrip("/"))

    cbImportSeurat2(inFname, outDir, datasetName, options)

    if options.port and not options.htmlDir:
        errAbort("--port requires --htmlDir")

    if options.htmlDir:
        build(outDir, options.htmlDir, port=options.port)
