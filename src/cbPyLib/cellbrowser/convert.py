# various format converters for single cell data:
# - cellranger, mtx to tsv, matConcat, etc

import logging, optparse, io, sys, os, shutil, operator, glob, re
from collections import defaultdict

from .cellbrowser import runGzip, openFile, errAbort, setDebug, moveOrGzip, makeDir, iterItems
from .cellbrowser import mtxToTsvGz, writeCellbrowserConf, getAllFields, readMatrixAnndata
from .cellbrowser import anndataToTsv, loadConfig, sanitizeName, lineFileNextRow, scanpyToCellbrowser, build

from os.path import join, basename, dirname, isfile, isdir, relpath, abspath, getsize, getmtime, expanduser

def cbToolCli_parseArgs(showHelp=False):
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] mtx2tsv|matCat|metaCat - convert various single-cell related files

    mtx2tsv   - convert matrix market to .tsv.gz
    matCat - merge expression matrices with one line per gene into a big matrix.
        Matrices must have identical genes in the same order and the same number of
        lines. Handles .csv files, otherwise defaults to tab-sep input. gzip OK.
    metaCat - concat meta tables

    Examples:
    - %prog mtx2tsv matrix.mtx genes.tsv barcodes.tsv exprMatrix.tsv.gz - convert .mtx to .tsv.gz file
    - %prog matCat mat1.tsv.gz mat2.tsv.gz exprMatrix.tsv.gz - concatenate expression matrices
    - %prog metaCat meta.tsv seurat/meta.tsv scanpy/meta.tsv newMeta.tsv - merge meta matrices
    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages")
    parser.add_option("", "--fixDots", dest="fixDots", action="store_true",
        help="try to fix R's mangling of various special chars to '.' in the cell IDs")


    (options, args) = parser.parse_args()

    if showHelp:
        parser.print_help()
        exit(1)

    setDebug(options.debug)
    return args, options

def cbToolCli():
    " run various tools from the command line "
    args, options = cbToolCli_parseArgs()

    if len(args)<=1:
        cbToolCli_parseArgs(showHelp=True)
        sys.exit(1)

    cmd = args[0]

    cmds = ["mtx2tsv", "matCat", "metaCat"]

    if cmd=="mtx2tsv":
        mtxFname = args[1]
        geneFname = args[2]
        barcodeFname = args[3]
        outFname = args[4]
        mtxToTsvGz(mtxFname, geneFname, barcodeFname, outFname)
    elif cmd=="matCat":
        inFnames = args[1:-1]
        outFname = args[-1]
        matCat(inFnames, outFname)
    elif cmd=="metaCat":
        inFnames = args[1:-1]
        outFname = args[-1]
        metaCat(inFnames, outFname, options)
    else:
        errAbort("Command %s is not a valid command. Valid commands are: %s" % (cmd, ", ".join(cmds)))

def matCat(inFnames, outFname):
    tmpFname = outFname+".tmp"
    ofh = openFile(tmpFname, "w")

    ifhs = []

    sep = "\t"
    if ".csv" in inFnames[0]:
        sep = ","

    for fn in inFnames:
        ifhs.append( openFile(fn) )

    headers, colCounts = getAllFields(ifhs, sep)
    ofh.write("\t".join(headers))
    ofh.write("\n")

    for i, ifh in enumerate(ifhs):
        logging.info("File %s: %d columns with values" % (ifhs[i].name, colCounts[i]-1)) 

    fileCount = len(inFnames)
    progressStep = 1000

    doProcess = True
    lineCount = 0

    while doProcess:
        geneId = None

        lineCount += 1
        if lineCount % progressStep == 0:
            logging.info("Processed %d rows" % (lineCount))

        for i, ifh in enumerate(ifhs):
            lineStr = ifh.readline()
            # check if we're at EOF
            if lineStr=='':
                doProcess = False
                break

            fields = lineStr.rstrip('\r\n').split(sep)
            fields = [x.strip('"') for x in fields]
            if (len(fields)!=colCounts[i]): # check number of columns against header count
                raise Exception("Illegal number of fields: file %s, line %d, field count: %d, expected %d" % 
                    (ifh.name, lineCount, len(fields), colCounts[i]))

            # check the gene ID
            if i==0: # get the gene ID from the first file
                geneId = fields[0]
                allVals = [geneId]
            else:
                #assert(geneId == fields[0]) # if this fails, then the rows are not in the same order
                if geneId != fields[0]:
                    print("Error: File %s has gene IDs that is out of order." % ifh.name)
                    print("Expected geneID: %s, got geneID: %s" % (geneId, fields[0]))
                    sys.exit(1)

            allVals.extend(fields[1:])

        ofh.write("\t".join(allVals))
        ofh.write("\n")

    # make sure that we've read through all files
    for ifh in ifhs:
        assert(ifh.readline()=='') # a file has still lines left to read?

    ofh.close()
    moveOrGzip(tmpFname, outFname)
    logging.info("Wrote %d lines (not counting header)" % lineCount)

def metaCat(inFnames, outFname, options):
    " merge all tsv/csv columns in inFnames into a new file, outFname. Column 1 is ID to join on. "
    allHeaders = ["cellId"]
    allRows = defaultdict(dict) # cellId -> fileIdx -> list of non-ID fields
    fieldCounts = {} # fileIdx -> number of non-ID fields
    allIds = set() # set with all cellIds

    for fileIdx, fname in enumerate(inFnames):
        headers = None
        for row in lineFileNextRow(fname):
            if headers is None:
                headers = row._fields[1:]
                logging.info("Reading %s, %d columns:  %s" % (fname, len(headers), repr(headers)))
                allHeaders.extend(headers)
                fieldCounts[fileIdx] = len(headers)

            cellId = row[0]
            allIds.add(cellId)
            allRows[cellId][fileIdx] = row[1:]

    nonCharRe = re.compile(r'[^a-zA-Z0-9]')

    if options.fixDots:
        # first create a map realId -> rId
        rToReal = {}
        for cellId, rowData in iterItems(allRows):
            rId = nonCharRe.sub(".", cellId)
            if rId!=cellId and rId in allRows:
                assert(rId not in rToReal) # Uh-oh! Mapping R <-> realId is not simple. May need some manual work.
                rToReal[rId]=cellId
        logging.info("Found %d cell IDs that look like they've been mangled by R" % (len(rToReal)))

        # now merge the R-id entries into the real ID entries
        newRows = {}
        for rId, readlId in iterItems(rToReal):
            allRows[cellId].update(allRows[rId])
        for rId in rToReal.keys():
            del allRows[rId]
        logging.info("Merged back rIds into normal data")

    tmpFname = outFname+".tmp"
    ofh = openFile(tmpFname, "w")
    ofh.write("\t".join(allHeaders))
    ofh.write("\n")

    for cellId, rowData in iterItems(allRows):
        row = []
        for fileIdx in range(0, len(inFnames)):
            if fileIdx in rowData:
                row.extend(rowData[fileIdx])
            else:
                row.extend([""]*fieldCounts[fileIdx])

        ofh.write(cellId)
        ofh.write("\t")
        ofh.write("\t".join(row))
        ofh.write("\n")

    ofh.close()
    os.rename(tmpFname, outFname)
    logging.info("Wrote %d lines (not counting header)" % len(allRows))

def crangerToCellbrowser(datasetName, inDir, outDir):
    " convert cellranger output to a cellbrowser directory "
    makeDir(outDir)
    # copy over the clusters
    clustFname = join(inDir, "analysis/clustering/graphclust/clusters.csv")
    if not isfile(clustFname):
        logging.warn("Cannot find %s" % clustFname)
        clustFname = join(inDir, "analysis/kmeans/10_clusters/clusters.csv")
        logging.warn("Using%s instead" % clustFname)
    metaFname = join(outDir, "meta.csv")
    shutil.copy(clustFname, metaFname)

    # copy over the t-SNE coords
    tsneFname = join(inDir, "analysis/tsne/2_components/projection.csv")
    if not isfile(tsneFname):
        tsneFname = join(inDir, "analysis/tsne/projection.csv")
    coordFname = join(outDir, "tsne.coords.csv")
    shutil.copy(tsneFname, coordFname)

    # copy over the markers
    dgeFname = join(inDir, "analysis/diffexp/graphclust/differential_expression.csv")
    if not isfile(dgeFname):
        logging.warn("Not found: %s" % dgeFname)
        dgeFname = join(inDir, "analysis/kmeans/10_clusters/differential_expression.csv")
        logging.warn("Using instead: %s" % dgeFname)

    markerFname = join(outDir, "markers.tsv")
    geneFname = join(outDir, "gene2sym.tsv")
    crangerSignMarkers(dgeFname, markerFname, geneFname, 0.01, 100)

    # convert the matrix
    outExprFname = join(outDir, "exprMatrix.tsv.gz")
    mask1 = join(inDir, "filtered_gene_bc_matrices/*/matrix.mtx")
    logging.info("Looking for %s" % mask1)
    matFnames = glob.glob(mask1)
    if len(matFnames)!=0:
        assert(len(matFnames)==1)
        matFname = matFnames[0]
        barcodeFname = matFname.replace("matrix.mtx", "barcodes.tsv")
        geneFname = matFname.replace("matrix.mtx", "genes.tsv")
        mtxToTsvGz(matFname, geneFname, barcodeFname, outExprFname)
    else:
        mask2 = join(inDir, "*_filtered_gene_bc_matrices_h5.h5")
        logging.info("Looking for %s" % mask2)
        matFnames = glob.glob(mask2)
        if len(matFnames)==0:
            errAbort("Could not find matrix, neither as %s nor as %s" % (mask1, mask2))
        import scanpy.api as sc
        logging.info("Reading matrix %s" % matFname)
        adata = readMatrixAnndata(matFname)
        anndataToTsv(adata, outExprFname)

    confName = join(outDir, "cellbrowser.conf")
    coordDescs = [{"file":"tsne.coords.csv", "shortLabel":"CellRanger t-SNE"}]
    confArgs =  {"meta" : "meta.csv", "clusterField" : "Cluster", "tags" : ["10x"]}
    writeCellbrowserConf(datasetName, coordDescs, confName, confArgs)

    crangerWriteMethods(inDir, outDir, matFname)
    generateDownloads(datasetName, outDir)

def crangerWriteMethods(inDir, outDir, matFname):
    htmlFname = join(outDir, "methods.html")
    if isfile(htmlFname):
        logging.info("%s exists, not overwriting" % htmlFname)
        return

    import csv
    csvMask = join(inDir, "*_metrics_summary.csv")
    csvFnames = glob.glob(csvMask)
    if len(csvFnames)==0:
        logging.warn("Cannot find %s, not writing %s" % (csvMask, htmlFname))
        return

    qcVals = list(csv.DictReader(open(csvFnames[0])))[0]

    ofh = open(htmlFname, "w")
    ofh.write("This dataset was imported from a CellRanger analysis directory with cbCellranger.<p><p>")
    ofh.write("<p><b>QC metrics reported by CellRanger:</b></p>\n")

    for key, value in iterItems(qcVals):
        ofh.write("%s: %s<br>\n" % (key, value))
    ofh.close()
    logging.info("Wrote %s" % ofh.name)

def generateDownloads(datasetName, outDir):
    htmlFname = join(outDir, "downloads.html")
    if isfile(htmlFname):
        logging.info("%s exists, not overwriting" % htmlFname)
        return

    ofh = open(htmlFname, "w")
    ofh.write("<b>Expression matrix:</b> <a href='%s/exprMatrix.tsv.gz'>exprMatrix.tsv.gz</a><p>\n" % datasetName)

    cFname = join(outDir, "cellbrowser.conf")
    if isfile(cFname):
        conf = loadConfig(cFname)
        if "unit" in conf:
            ofh.write("Unit of expression matrix: %s<p>\n" % conf["unit"])

    ofh.write("<b>Cell meta annotations:</b> <a href='%s/meta.tsv'>meta.tsv</a><p>" % datasetName)

    coordDescs = conf["coords"]
    for coordDesc in coordDescs:
        coordLabel = coordDesc["shortLabel"]
        cleanName = sanitizeName(coordLabel.replace(" ", "_"))
        coordFname = cleanName+".coords.tsv.gz"
        ofh.write("<b>%s coordinates:</b> <a href='%s/%s'>%s</a><br>" % (coordLabel, datasetName, coordFname, coordFname))

    rdsFname = join(datasetName, "seurat.rds")
    if isfile(rdsFname):
        ofh.write("<b>Seurat R data file:</b> <a href='%s'>seurat.rds</a><p>" % rdsFname)

    scriptFname = join(datasetName, "runSeurat.R")
    if isfile(scriptFname):
        ofh.write("<b>Seurat R analysis script:</b> <a href='%s'>runSeurat.R</a><p>" % scriptFname)

    logFname = join(datasetName, "analysisLog.txt")
    if isfile(logFname):
        ofh.write("<b>Analysis Log File:</b> <a href='%s'>analysisLog.txt</a><p>" % logFname)

def crangerSignMarkers(dgeFname, markerFname, geneFname, maxPval, maxGenes):
    " convert cellranger diff exp file to markers.tsv file "
    # Old Cellranger 1 format:
    # Gene ID,Gene Name,Cluster 1 Weight,Cluster 1 UMI counts/cell,Cluster 2 Weight,Cluster 2 UMI counts/cell,Cluster 3 Weight,Cluster 3 UMI counts/cell,Cluster 4 Weight,Cluster 4 UMI counts/cell,Cluster 5 Weight,Cluster 5 UMI counts/cell,Cluster 6 Weight,Cluster 6 UMI counts/cell,Cluster 7 Weight,Cluster 7 UMI counts/cell,Cluster 8 Weight,Cluster 8 UMI counts/cell,Cluster 9 Weight,Cluster 9 UMI counts/cell,Cluster 10 Weight,Cluster 10 UMI counts/cell

    ofh = open(markerFname, "w")
    ofh.write("cluster\tgene\tAdj. P-Value\tLog2 fold change\tMean UMI Counts\n")

    clusterToGenes = defaultdict(list)

    # read the significant markers and their p-Values
    for line in open(dgeFname):
        if line.startswith("Gene ID"):
            if "Cluster 1 Weight" in line:
                fieldsProCluster = 2
                fileVersion = 1
            else:
                fieldsProCluster = 3
                fileVersion = 2
            continue
        row = line.rstrip("\n\r").split(",")
        geneId = row[0]
        sym = row[1]
        clusterCount = int((len(row)-2) / fieldsProCluster)
        for clusterIdx in range(0, clusterCount):
            startField = (clusterIdx*fieldsProCluster)+2
            if fileVersion==2:
                mean = float(row[startField])
                fc = float(row[startField+1])
                pVal = float(row[startField+2])
            else:
                fc = float(row[startField])
                mean = float(row[startField+1])
                pVal = 0.0001

            if pVal < maxPval:
                clusterToGenes[clusterIdx+1].append((fc, mean, pVal, sym))

    # write out the markers
    for clusterId, clusterGenes in iterItems(clusterToGenes):
        clusterGenes.sort(key=operator.itemgetter(2)) # sort by fold change
        maxIdx = min(maxGenes, len(clusterGenes))
        for i in range(0, maxIdx):
            fc, mean, pVal, sym = clusterGenes[i]
            ofh.write("%d\t%s\t%g\t%f\t%f\n" % (clusterId, sym, pVal, fc, mean))

    ofh.close()
    logging.info("Wrote %s" % ofh.name)

def cbCellrangerCli():
    args, options = cbCellrangerCli_parseArgs()

    if options.outDir is None or options.inDir is None or options.datasetName is None:
        logging.error("You have to specify at least an input and an output directory and a dataset name.")
        cbCellrangerCli_parseArgs(showHelp=True)

    crangerToCellbrowser(options.datasetName, options.inDir, options.outDir)

def cbCellrangerCli_parseArgs(showHelp=False):
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] -i cellRangerDir -o outputDir - convert the cellranger output to cellbrowser format and create a cellranger.conf file

    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages")

    parser.add_option("-i", "--inDir", dest="inDir", action="store", help="input folder with the cellranger analysis output. This is the directory with the two directories 'analysis' and 'filtered_gene_bc_matrices'")
    parser.add_option("-o", "--outDir", dest="outDir", action="store", help="output directory")
    #parser.add_option("-g", "--geneSet", dest="geneSet", action="store", help="geneset, e.g. gencode28 or gencode-m13 or similar. Default: %default", default="gencode24")
    parser.add_option("-n", "--name", dest="datasetName", action="store", help="name of the dataset")

    (options, args) = parser.parse_args()

    if showHelp:
        parser.print_help()
        exit(1)

    setDebug(options.debug)

    return args, options

def cbImportScanpy_parseArgs(showHelp=False):
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] input.h5ad outDir datasetName - convert Scanpy AnnData object to cellbrowser

    Example:
    - %prog pbmc3k.h5ad pbmc3kScanpy pbmc3kScanpy - convert pbmc3k to directory with tab-separated files
    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages")

    parser.add_option("", "--htmlDir", dest="htmlDir", action="store",
        help="do not only convert to tab-sep files but also run cbBuild to"
            "convert the data and add the dataset under htmlDir")

    parser.add_option("-p", "--port", dest="port", action="store", type="int",
            help="only with --htmlDir: start webserver on port to serve htmlDir")

    parser.add_option("-m", "--skipMatrix", dest="skipMatrix", action="store_true",
        help="do not convert the matrix, saves time if the same one has been exported before to the "
        "same directory")

    (options, args) = parser.parse_args()

    if showHelp:
        parser.print_help()
        exit(1)

    setDebug(options.debug)
    return args, options

def cbImportScanpyCli():
    " convert h5ad to directory "
    args, options = cbImportScanpy_parseArgs()

    if len(args)<3:
        cbImportScanpy_parseArgs(showHelp=True)
        sys.exit(1)

    inFname, outDir, datasetName = args

    import anndata
    ad = anndata.read_h5ad(inFname)
    scanpyToCellbrowser(ad, outDir, datasetName, skipMatrix=options.skipMatrix)

    if options.port and not options.htmlDir:
        errAbort("--port requires --htmlDir")

    if options.htmlDir:
        build(outDir, options.htmlDir, port=options.port)
