# a wrapper around the R library Seurat

import logging, optparse, sys, glob, os, datetime, shutil
from os.path import join, basename, dirname, isfile, isdir, relpath, abspath, getsize, getmtime, expanduser, splitext

from .cellbrowser import copyPkgFile, writeCellbrowserConf, pipeLog, makeDir, maybeLoadConfig, errAbort, popen
from .cellbrowser import setDebug, build, isDebugMode, generateHtmls, runCommand

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

    parser.add_option("", "--threadCount", dest="threadCount", action="store", type="int",
            help="Number of threads to use via the future library. Default is not use multithreading, so there is no requirement for future library")

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

def checkRVersion():
    " make sure that we have > 3.4"
    logging.debug("Checking R version")
    cmd = ["Rscript", "-e",  "ver = R.Version(); if (ver$major<3 || (ver$major==3 && ver$minor < 4)) { stop('Sorry, Seurat requires at least R 3.4') } "]
    runCommand(cmd)

def getSeuratVersion():
    " return current Seurat version "
    logging.debug("Checking Seurat version")
    cmd = ["Rscript", "-e",  "message(packageVersion('Seurat'))"]
    proc, stdout = popen(cmd, useStderr=True)
    if proc is None:
        return None
    else:
        verStr = stdout.readline()
        logging.debug("version is %s", verStr)
        return verStr

def writeCbSeuratScript(conf, inData, tsnePath, clusterPath, markerPath, rdsPath, matrixPath, scriptPath,
        datasetName, outDir, threadCount):
    " write the seurat R script to a file "
    checkRVersion()

    version = getSeuratVersion()
    if version is None:
        logging.warn("Seurat is not installed, trying to install it now")
        cmd = ["Rscript", "-e", "install.packages(c('Seurat', 'data.table'), dep=TRUE, repos='http://cran.r-project.org/')"]
        runCommand(cmd)
        version = getSeuratVersion()
        if version is None:
            errAbort("Could not find or install Seurat")

    isSeurat2 = False
    isSeurat3 = False
    isSeurat4 = False

    verStr = version.split(".")[0]
    if verStr=="2":
        isSeurat2 = True
    elif verStr=="3":
        isSeurat3 = True
    elif verStr=="4":
        isSeurat4 = True
    else:
        assert(False) # unsupported Seurat version

    cmds = []
    # try to install Seurat if not already installed
    cmds.append("ver = R.Version()")
    cmds.append("if (ver$major<3 || (ver$major==3 && ver$minor < 4)) { error('Sorry, Seurat requires at least R 3.4') } ")
    cmds.append("message('Seurat loaded, version ', packageVersion('Seurat'))")
    cmds.append('print("Seurat: Reading data")')
    cmds.append("library(methods)")
    cmds.append("suppressWarnings(suppressMessages(library(Seurat)))")
    if threadCount!=None and threadCount > 0:
        cmds.append('print("Using %d cores")")' % threadCount)
        cmds.append("library(future)")
        cmds.append('plan(strategy = "multicore", workers = %d)' % threadCount)
    cmds.append('print("Seurat: Reading data")')
    cmds.append('print("Loading input data matrix")')
    readExportScript(cmds) # add the ExportToCellbrowser function to cmds

    writeMatrix = False
    if isfile(inData):
        if inData.endswith(".mtx"):
            logging.error("You specified an .mtx file as the input matrix")
            logging.error("Seurat cannot read that. You need to specify the directory name, with a file matrix.mtx in it.")
            logging.error("The directory also has to contain two other files, genes.tsv and barcodes.tsv")
            logging.error("They have to be called like this. Please rename your input files if necessary.")
            exit(1)

        if inData.endswith(".rds"):
            cmds.append('mat = readRDS(file = "%s")' % inData)
        else:
            if ".csv" in inData:
                cmds.append('print("Loading csv file with read.table")')
                cmds.append('mat = read.table("%s",sep=",",header=TRUE,row.names=1)' % inData)
            else:
                cmds.append('print("Loading tsv file with read.delim")')
                cmds.append('mat = read.delim("%s",header=TRUE,row.names=1)' % inData)
    else:
        matrixFname = join(inData, "matrix.mtx")
        if not isfile(matrixFname):
            errAbort("Could not find file %s - sure that you specified the directory of the matrix.mtx file?" % matrixFname)
        cmds.append('mat = Read10X(data.dir="%s")' % inData)

    undoLog = conf.get("undoLog", False)
    if undoLog:
        cmds.append('print("Undoing log2 on matrix, doing 2^mat")')
        cmds.append('mat <- 2^(mat-1)')

    minCells = conf.get("minCells", 3)
    minGenes = conf.get("minGenes", 200)

    cmds.append('print("Seurat: Setup")')
    if isSeurat3:
        cmds.append('sobj <- CreateSeuratObject(mat)')
    elif isSeurat4:
        cmds.append('sobj <- CreateSeuratObject(counts = mat, min.cells=%d, min.genes=%d)' % (minCells, minGenes))
    else:
        cmds.append('sobj <- CreateSeuratObject(raw.data = mat, min.cells=%d, min.genes=%d)' % (minCells, minGenes))
    cmds.append('sobj') # print size of the matrix

    # find mito genes, mito-%, and create plots for it
    # XX - ENSG names?
    if isSeurat3 or isSeurat4:
        cmds.append('sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")')
    else:
        cmds.append('mito.genes <- grep(pattern = "^MT-", x = rownames(x = sobj@data), value = TRUE)')
        cmds.append('percent.mito <- Matrix::colSums(sobj@raw.data[mito.genes, ])/Matrix::colSums(sobj@raw.data)')
        cmds.append('sobj <- AddMetaData(object = sobj, metadata = percent.mito, col.name = "percent.mito")')

    if isSeurat2:
        cmds.append('VlnPlot(object = sobj, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)')
        cmds.append('par(mfrow = c(1, 2))')
        cmds.append('GenePlot(object = sobj, gene1 = "nUMI", gene2 = "percent.mito")')
        cmds.append('GenePlot(object = sobj, gene1 = "nUMI", gene2 = "nGene")')


    # remove cells with > mitoMax and a gene count < minGenes or > maxGenes
    mitoMax=conf.get("mitoMax", 0.05)
    maxGenes=conf.get("maxGenes", 2500)

    # filter based on mito%
    cmds.append('print("Removing cells with more than %f percent of mitochondrial genes")' % (100*mitoMax))
    cmds.append('print("Keeping only cells with a geneCount %d-%d")' % (minGenes, maxGenes))

    if isSeurat3 or isSeurat4:
        cmds.append('sobj <- subset(sobj, subset = nFeature_RNA > %(minGenes)d & nFeature_RNA < %(maxGenes)d & percent.mt < %(mitoMax)f)' % locals())
    else:
        cmds.append('sobj <- FilterCells(object = sobj, subset.names = c("nGene", "percent.mito"), low.thresholds = c(%(minGenes)d, -Inf), high.thresholds = c(%(maxGenes)d, %(mitoMax)f))' % locals())
    cmds.append('sobj') # print size of the matrix

    # do the log
    cmds.append('sobj <- NormalizeData(object = sobj, normalization.method = "LogNormalize", scale.factor = 10000)')

    if isSeurat3 or isSeurat4:
        cmds.append('sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)')
    else:
        # find most variable genes and gate on it
        minMean = conf.get("varMinMean", 0.0125)
        maxMean = conf.get("varMaxMean", 3)
        minDisp = conf.get("varMinDisp", 0.5)
        cmds.append('sobj <- FindVariableGenes(object = sobj, mean.function = ExpMean, dispersion.function = LogVMR, '
            'x.low.cutoff = %f, x.high.cutoff = %f, y.cutoff = %s)' % (minMean, maxMean, minDisp))
    cmds.append('sobj') # print size of the matrix

    # scale
    if isSeurat3 or isSeurat4:
        cmds.append('sobj <- ScaleData(object = sobj)')
    else:
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
    if isSeurat3 or isSeurat4:
        cmds.append('sobj <- RunPCA(sobj, npcs = %d)' % pcCountConfig)
        cmds.append('VizDimLoadings(sobj, dims = 1:2, reduction = "pca")')
        cmds.append('VizDimLoadings(sobj, dims = 2:3, reduction = "pca")')
    else:
        cmds.append('sobj <- RunPCA(object = sobj, pcs.compute = %d, pc.genes = sobj@var.genes, do.print = FALSE, pcs.print = 1:5, genes.print = 5)' % pcCountConfig)
        cmds.append('print("PCA plots")')
        cmds.append('VizPCA(object = sobj, pcs.use = 1:2)')
        cmds.append('PCAPlot(object = sobj, dim.1 = 1, dim.2 = 2)')
        cmds.append('PCHeatmap(object = sobj, pc.use = 1:8, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)')
        cmds.append('PCElbowPlot(object = sobj)')

    louvainRes = conf.get("louvainRes", 0.6)
    cmds.append('print("Finding clusters with resolution %f")' % louvainRes)

    if isSeurat3 or isSeurat4:
        cmds.append('sobj <- FindNeighbors(sobj, dims=%d)' % (min(10, pcCountConfig)))
        cmds.append('sobj <- FindClusters(sobj, resolution = %f)' % louvainRes)
    else:
        cmds.append('sobj <- FindClusters(object = sobj, reduction.type = "pca", resolution = %f, print.output = 0, save.SNN = TRUE)' % (louvainRes))
        cmds.append("PrintFindClustersParams(object = sobj)")

    cmds.append('print("Running t-SNE")')
    perplexity = str(conf.get("perplexity", 30))

    doUmap = conf.get("doUmap", False)
    if isSeurat3 or isSeurat4:
        cmds.append("sobj <- RunTSNE(sobj, perplexity=%s)" % perplexity)
        if doUmap:
            cmds.append("sobj <- RunUMAP(sobj, dims=1:%s)" % str(pcCountConfig))
    else:
        cmds.append("sobj <- RunTSNE(object = sobj, do.fast = TRUE, check_duplicates=FALSE, perplexity=%s)" % perplexity)
        cmds.append("TSNEPlot(object = sobj, doLabel=T)")

    minMarkerPerc = conf.get("minMarkerPerc", 0.25)
    cmds.append('print("Finding markers")')

    if isSeurat3 or isSeurat4:
        cmds.append('all.markers <- FindAllMarkers(object = sobj)')
        cmds.append('sobj@misc[["markers"]] <- all.markers')

    cmds.append('print("Saving .rds to %s")' % rdsPath)
    cmds.append('saveRDS(sobj, file = "%s")' % rdsPath)

    cmds.append("message('Exporting Seurat data object to cbBuild directory %s')" % outDir)
    cmds.append("ExportToCellbrowser(sobj, '%s', '%s', use.mtx=T, matrix.slot='%s')" %
            (outDir, datasetName, "counts"))

    writeRScript(cmds, scriptPath, "cbSeurat")

def cbSeurat2Cli():
    " stay backwards compatible "
    cbSeuratCli()

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
    threadCount = options.threadCount

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

    confArgs = {"clusterField":"Cluster"}
    if inConf.get("skipMatrixExport", False):
        confArgs["exprMatrix"] = inMatrix
        matrixPath = None

    writeCbSeuratScript(inConf, inMatrix, tsnePath, clusterPath, markerPath, rdsPath, matrixPath, scriptPath,
            datasetName, outDir, threadCount)
    runRscript(scriptPath, logPath)

    if not isfile(markerPath):
        errAbort("R script did not complete successfully. Check %s and analysisLog.txt." % scriptPath)

    coords = [{'shortLabel':'t-SNE', 'file':'tsne.coords.tsv'}]


    writeCellbrowserConf(datasetName, coords, cbConfPath, args=confArgs)

    generateHtmls(datasetName, outDir, desc = {"supplFiles": [{"label":"Seurat RDS", "file":basename(rdsPath)}]})
    copyPkgFile("sampleConfig/desc.conf", outDir)

def cbImportSeurat_parseArgs(showHelp=False):
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] -i input.rds -o outDir [-n datasetName] - convert Seurat object to cellbrowser

    Example:
    - %prog -i pbmc3k.rds -o pbmc3kSeurat - convert pbmc3k to directory of tab-separated files
    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages. Also activates debug() in the generated R script.")

    parser.add_option("-i", "--in", dest="inFname", action="store",
        help="input file. Required parameter")

    parser.add_option("-o", "--outDir", dest="outDir", action="store",
        help="Output directory. Required parameter")

    parser.add_option("-n", "--name", dest="datasetName", action="store",
        help="Dataset name for generated cellbrowser.conf. If not specified, the last component of -o will be used.")

    parser.add_option("-f", "--inFormat", dest="inFormat", action="store",
            help="the format of the input file. Either 'rds' or 'rdata'. If nothing is specified and the input filename ends with .rdata or .robj, the option defaults to rdata, otherwise rds is used.")

    parser.add_option("", "--htmlDir", dest="htmlDir", action="store",
        help="do not only convert to tab-sep files but also run cbBuild to"
            "convert the data and add the dataset under htmlDir")

    parser.add_option("-p", "--port", dest="port", action="store", type="int",
            help="only with --htmlDir: start webserver on port to serve htmlDir")

    parser.add_option("-x", "--skipMatrix", dest="skipMatrix", action="store_true",
            default = False,
        help="do not convert the matrix, saves time if the same one has been exported before to the "
        "same outDir directory")

    parser.add_option("", "--threads", dest="threadCount", action="store", type="int", default=0,
            help="activate multiprocess strategy, default thread count is 0, which uses no multithreading")

    parser.add_option("-m", "--skipMarkers", dest="skipMarkers", action="store_true",
            default = False,
        help="do not calculate cluster-specific markers with FindAllMarkers(), saves a lot of time")

    parser.add_option("-c", "--clusterField", dest="clusterField", action="store",
        help="Cluster field to color on, by default this is the @ident slot of the Seurat object but it can also be any other meta data field of the @meta.data slot")

    parser.add_option("", "--markerFile", dest="markerFile", action="store",
            help="Instead of calculating cluster markers again, use this file. Format: cluster,gene,pVal + any other fields. Or alternatively the native Seurat cluster markers format, as created by write.table")

    parser.add_option("", "--useMtx", dest="useMtx", action="store_true",
            help="Write a .mtx.gz file, instead of a tsv.gz file. Necessary for big datasets.")

    parser.add_option("-s", "--matrixSlot", dest="matrixSlot", action="store", default="counts,data,scale.data",
            help="Export this slot of the matrix. Can be 'counts', 'scale.data' or 'data'. Default is %default. The first value is the default matrix for coloring the plot.")

    parser.add_option("", "--assay", dest="assay", action="store",
            help="Select the Seurat3 assay to export.")

    (options, args) = parser.parse_args()

    if showHelp:
        parser.print_help()
        exit(1)

    setDebug(options.debug)
    return args, options

def readExportScript(cmds):
    " find the cellbrowser.R script and add the commands between # --- to cmds "
    fname = join(dirname(__file__), "R", "cellbrowser.R")

    blockFound = False
    for line in open(fname):
        cmds.append(line.rstrip("\n"))

    # the exactly same R export source file is also part of seurat-wrappers now.
    # we want to have only a single source code file, and seurat-wrappers code 
    # cannot use require, so we add the require commands here
    cmds.insert(0, "message(R.version$version.string)")
    cmds.insert(0, "library(Matrix, warn.conflicts=FALSE)")
    cmds.insert(0, "library(R.utils, warn.conflicts=FALSE)")
    cmds.insert(0, "library(R.oo, warn.conflicts=FALSE)")

    assert(len(cmds)!=0)
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
    logging.info("Wrote R script to %s" % scriptPath)

def findObjectVersion(outDir):
    " look for analysisLog.txt and search for the SeuratObject line in there "
    logFname = join(outDir, "analysisLog.txt")
    logging.debug("Checking %s to find the Seurat version used to write the object file" % logFname)
    if not isfile(logFname):
        logging.info("Cannot find %, no Seurat version available" % logFname)
        return ""

    for l in open(logFname):
        if l.startswith("ObjectVersion"):
            ver = l.strip().split(" = ")[1].strip()
            return ver

    return ""

def cbImportSeurat(inFname, outDir, datasetName, options):
    " convert Seurat 2 or 3 .rds/.rdata file to tab-sep directory for cellbrowser "
    logging.info("inFname: %s, outDir: %s, datasetName: %s" % (inFname, outDir, datasetName))

    makeDir(outDir)

    skipMatrix = options.skipMatrix
    skipMarkers = options.skipMarkers
    clusterField = options.clusterField
    inFormat = options.inFormat
    assay = options.assay
    markerFile = options.markerFile
    if markerFile is None:
        markerFileStr = "NULL"
    else:
        markerFileStr = '"%s"' % markerFile

    tsnePath = join(outDir, "tsne.coords.tsv")
    metaPath = join(outDir, "meta.tsv")
    markerPath = join(outDir, "markers.tsv")
    scriptPath = join(outDir, "runSeurat.R")
    matrixPath = join(outDir, "matrix.mtx")
    logPath = join(outDir, "analysisLog.txt")

    cmds = ["require(methods)"] # for the 'slots()' function
    cmds.append("require(Seurat)")

    cmds = readExportScript(cmds)

    if inFormat is None:
        inExt = splitext(inFname.lower())[1]
        if inExt in [".robj", ".rdata"]:
            inFormat="rdata"
        else:
            inFormat = "rds"

    if inFormat=="rds":
        cmds.append("message('Reading %s as .rds file')" % inFname)
        cmds.append("sobj <- readRDS(file='%s')" % inFname)
    else:
        cmds.append("message('Reading %s as .RData file, assuming that the first object in the file is the Seurat object')" % inFname)
        cmds.append("names <- load('%s')" % inFname)
        cmds.append("sobj <- get(names[1])")

    cmds.append("if (class(sobj)!='seurat' && class(sobj)[1]!='Seurat') { stop('The input .rds file does not seem to contain a Seurat object') }")
    skipStr = str(skipMatrix).upper()
    skipMarkerStr = str(skipMarkers).upper()

    if isDebugMode():
        cmds.append("debugonce(ExportToCellbrowser)")

    if clusterField is None:
        clusterStr = 'NULL'
    else:
        clusterStr = "'%s'" % clusterField

    useMtx = "FALSE"
    if options.useMtx:
        useMtx = "TRUE"

    matrixSlot = "counts,dataScale"
    if options.matrixSlot:
        matrixSlot = options.matrixSlot

    if assay:
        cmds.append("message('Setting Seurat assay to %s')" % assay)
        cmds.append("DefaultAssay(sobj) = '%s'" % assay)

    cmds.append("message('Exporting Seurat data to %s')" % outDir)
    cmds.append("ExportToCellbrowser(sobj, '%s', '%s', markers.file = %s, cluster.field=%s, skip.expr.matrix = %s, skip.markers = %s, use.mtx=%s, matrix.slot='%s')" %
            (outDir, datasetName, markerFileStr, clusterStr, skipStr, skipMarkerStr, useMtx, matrixSlot))

    writeRScript(cmds, scriptPath, "cbImportSeurat")
    runRscript(scriptPath, logPath)
    if not isfile(metaPath):
        errAbort("R script did not complete successfully. Check %s and analysisLog.txt." % scriptPath)

    descDict = None
    if inFormat=="rds":
        rdsOutPath = join(outDir, basename(inFname))
        logging.info("Copying %s to %s" % (inFname, rdsOutPath))
        shutil.copyfile(inFname, rdsOutPath)
        objectVersion = findObjectVersion(outDir)
        descDict = {"supplFiles": [{"label":"Seurat %s RDS" % objectVersion, "file":basename(inFname)}]}

    cbConfPath = join(outDir, "cellbrowser.conf")

    generateHtmls(datasetName, outDir, desc = descDict)

def cbImportSeuratCli():
    " convert .rds to directory "
    args, options = cbImportSeurat_parseArgs()

    inFname = options.inFname
    outDir = options.outDir
    if None in [inFname, outDir]:
        cbImportSeurat_parseArgs(showHelp=True)
        errAbort("You need to specify at least an input rds file and an output directory")

    datasetName = options.datasetName
    if datasetName is None:
        datasetName = basename(abspath(outDir).rstrip("/"))

    cbImportSeurat(inFname, outDir, datasetName, options)

    if options.port and not options.htmlDir:
        errAbort("--port requires --htmlDir")

    if options.htmlDir:
        build(outDir, options.htmlDir, port=options.port)
