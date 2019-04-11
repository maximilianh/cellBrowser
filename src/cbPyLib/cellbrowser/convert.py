# various format converters for single cell data:
# - cellranger, mtx to tsv, matcat, metaCat etc

import logging, optparse, io, sys, os, shutil, operator, glob, re, json
from collections import defaultdict

from .cellbrowser import runGzip, openFile, errAbort, setDebug, moveOrGzip, makeDir, iterItems
from .cellbrowser import mtxToTsvGz, writeCellbrowserConf, getAllFields, readMatrixAnndata
from .cellbrowser import anndataToTsv, loadConfig, sanitizeName, lineFileNextRow, scanpyToCellbrowser, build
from .cellbrowser import generateHtmls

from os.path import join, basename, dirname, isfile, isdir, relpath, abspath, getsize, getmtime, expanduser

def cbToolCli_parseArgs(showHelp=False):
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] mtx2tsv|matCat|metaCat - convert various single-cell related files

    mtx2tsv   - convert matrix market to .tsv.gz
    matCat - merge expression matrices with one line per gene into a big matrix.
        Matrices must have identical genes in the same order and the same number of
        lines. Handles .csv files, otherwise defaults to tab-sep input. gzip OK.
    metaCat - concat/join meta tables on the first (cell ID) field

    Examples:
    - %prog mtx2tsv matrix.mtx genes.tsv barcodes.tsv exprMatrix.tsv.gz - convert .mtx to .tsv.gz file
    - %prog matCat mat1.tsv.gz mat2.tsv.gz exprMatrix.tsv.gz - concatenate expression matrices
    - %prog metaCat meta.tsv seurat/meta.tsv scanpy/meta.tsv newMeta.tsv - merge meta matrices
    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages")
    parser.add_option("", "--fixDots", dest="fixDots", action="store_true",
        help="try to fix R's mangling of various special chars to '.' in the cell IDs")
    parser.add_option("", "--first", dest="first", action="store",
        help="only for metaCat: names of fields to order first, comma-sep, e.g. disease,age. Not cellId, that's always the first field")


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

def reorderFields(row, firstFields):
    " reorder the row to have firstFields first "
    #logging.debug("Reordering row to have %s fields first" % firstFields)
    newRow = [row[0]]
    for idx in firstFields:
        newRow.append(row[idx])

    for i in range(1, len(row)):
        if i not in firstFields:
            newRow.append(row[i])

    return newRow

def metaCat(inFnames, outFname, options):
    " merge all tsv/csv columns in inFnames into a new file, outFname. Column 1 is ID to join on. "
    allHeaders = ["cellId"]
    allRows = defaultdict(dict) # cellId -> fileIdx -> list of non-ID fields
    fieldCounts = {} # fileIdx -> number of non-ID fields
    allIds = set() # set with all cellIds

    firstFields = []
    if options.first!="" and options.first is not None:
        firstFields = options.first.split(",")

    for fileIdx, fname in enumerate(inFnames):
        logging.info("Reading %s" % fname)
        headers = None
        for row in lineFileNextRow(fname, headerIsRow=True):
            if headers is None:
                headers = row[1:]
                logging.info("Reading %s, %d columns:  %s" % (fname, len(headers), repr(headers)))
                allHeaders.extend(headers)
                fieldCounts[fileIdx] = len(headers)
                continue

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

    # find the field indices of the fields we're putting first
    firstFieldIdx = []
    for h in firstFields:
        try:
            idx = allHeaders.index(h)
        except ValueError:
            errAbort("Field %s specified on command line is not in the meta data file" % repr(h))
        firstFieldIdx.append(idx)

    tmpFname = outFname+".tmp"
    ofh = openFile(tmpFname, "w")
    allHeaders = reorderFields(allHeaders, firstFieldIdx)
    ofh.write("\t".join(allHeaders))
    ofh.write("\n")

    for cellId, rowData in iterItems(allRows):
        row = [cellId]
        for fileIdx in range(0, len(inFnames)):
            if fileIdx in rowData:
                row.extend(rowData[fileIdx])
            else:
                row.extend([""]*fieldCounts[fileIdx])

        row = reorderFields(row, firstFieldIdx)
        ofh.write("\t".join(row))
        ofh.write("\n")

    ofh.close()
    os.rename(tmpFname, outFname)
    logging.info("Output field order is: %s" % allHeaders)
    logging.info("Wrote %d lines (not counting header)" % len(allRows))

def importCellrangerMatrix(inDir, outDir):
    " convert the cellranger 2 or 3 matrix, use the .mtx or .h5 file, whatever is available "
    outExprFname = join(outDir, "exprMatrix.tsv.gz")
    # cellranger 2: filtered_feature_bc_matrix/<db>/matrix.mtx and not gzipped
    mask1 = join(inDir, "filtered_gene_bc_matrices/*/matrix.mtx")
    # cellranger 3: filtered_gene_bc_matrices and gzipped
    mask2 = join(inDir, "filtered_feature_bc_matrix/matrix.mtx.gz")
    # h5 file as fallback, requires scanpy
    # cellranger3: filtered_feature_bc_matrix.h5
    mask3 = join(inDir, "*filtered_*matri*.h5") # should work for cr 1, 2 and 3
    matFnames1 = glob.glob(mask1)
    matFnames2 = glob.glob(mask2)
    matFnames3 = glob.glob(mask3)
    logging.info("Looking for %s (cellranger 1/2) or %s (cellranger3) or %s" % (mask1, mask2, mask3))

    # ugly code warning
    if len(matFnames1)!=0:
        # cellranger 1/2 files are not gziped
        assert(len(matFnames1)==1)
        matFname = matFnames1[0]
        logging.info("Found %s" % matFname)
        barcodeFname = matFname.replace("matrix.mtx", "barcodes.tsv")
        geneFname = matFname.replace("matrix.mtx", "genes.tsv")
        mtxToTsvGz(matFname, geneFname, barcodeFname, outExprFname)
    elif len(matFnames2)!=0:
        # cellranger 3 files are gzipped
        assert(len(matFnames2)==1)
        matFname = matFnames2[0]
        logging.info("Found %s" % matFname)
        barcodeFname = matFname.replace("matrix.mtx.gz", "barcodes.tsv.gz")
        geneFname = matFname.replace("matrix.mtx.gz", "features.tsv.gz")
        mtxToTsvGz(matFname, geneFname, barcodeFname, outExprFname)
    elif len(matFnames3)!=0:
        matFnames3 = glob.glob(mask3)
        logging.info("Found %s" % matFnames3)
        assert(len(matFnames3)==1)
        matFname = matFnames3[0]

        import scanpy.api as sc
        logging.info("Reading matrix %s" % matFname)
        adata = readMatrixAnndata(matFname)
        anndataToTsv(adata, outExprFname)
    else:
        errAbort("Could not find matrix, neither as %s, %s nor %s" % (mask1, mask2, mask3))
        logging.info("Looking for %s" % mask3)

    return matFname

def crangerToCellbrowser(datasetName, inDir, outDir, noMat):
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

    if not noMat:
        matFname = importCellrangerMatrix(inDir, outDir)
    else:
        matFname = "unknown"

    confName = join(outDir, "cellbrowser.conf")
    coordDescs = [{"file":"tsne.coords.csv", "shortLabel":"CellRanger t-SNE"}]
    confArgs =  {"meta" : "meta.csv", "clusterField" : "Cluster", "tags" : ["10x"]}
    writeCellbrowserConf(datasetName, coordDescs, confName, confArgs)

    crangerWriteMethods(inDir, outDir, matFname)
    generateHtmls(datasetName, outDir)

def crangerWriteMethods(inDir, outDir, matFname):
    htmlFname = join(outDir, "methods.html")
    if isfile(htmlFname):
        logging.info("%s exists, not overwriting" % htmlFname)
        return

    import csv
    csvMask = join(inDir, "*_metrics_summary.csv")
    csvFnames = glob.glob(csvMask)

    jsonFname = join(inDir, "summary.json")
    if len(csvFnames)==0:
        if isfile(jsonFname):
            qcVals = json.load(open(jsonFname))
        else:
            logging.warn("Cannot find %s nor %s, not writing %s" % (csvMask, jsonFname, htmlFname))
            return
    else:
        qcVals = list(csv.DictReader(open(csvFnames[0])))[0]

    ofh = open(htmlFname, "w")
    ofh.write("<p>This dataset was imported from a CellRanger analysis directory with cbCellranger.</p>\n")
    ofh.write("<p><b>QC metrics reported by CellRanger:</b></p>\n")

    for key, value in iterItems(qcVals):
        ofh.write("<i>%s:</i> %s<br>\n" % (key, str(value)))
    ofh.close()
    logging.info("Wrote %s" % ofh.name)

def crangerSignMarkers(dgeFname, markerFname, geneFname, maxPval, maxGenes):
    " convert cellranger diff exp file to markers.tsv file "
    # Old Cellranger 1 format:
    # Gene ID,Gene Name,Cluster 1 Weight,Cluster 1 UMI counts/cell,Cluster 2 Weight,Cluster 2 UMI counts/cell,Cluster 3 Weight,Cluster 3 UMI counts/cell,Cluster 4 Weight,Cluster 4 UMI counts/cell,Cluster 5 Weight,Cluster 5 UMI counts/cell,Cluster 6 Weight,Cluster 6 UMI counts/cell,Cluster 7 Weight,Cluster 7 UMI counts/cell,Cluster 8 Weight,Cluster 8 UMI counts/cell,Cluster 9 Weight,Cluster 9 UMI counts/cell,Cluster 10 Weight,Cluster 10 UMI counts/cell

    ofh = open(markerFname, "w")
    ofh.write("cluster\tgene\tAdj. P-Value\tLog2 fold change\tMean UMI Counts\n")

    clusterToGenes = defaultdict(list)

    # read the significant markers and their p-Values
    logging.info("Reading %s" % dgeFname)
    ifh = open(dgeFname)

    # determine file format
    line1 = ifh.readline()
    # cellranger 3
    fileVersion = None
    if line1.startswith("Feature ID"):
        fieldsProCluster = 3
        fileVersion = 3
    # cellranger 1 or 2
    elif line1.startswith("Gene ID"):
        if "Cluster 1 Weight" in line1:
            fieldsProCluster = 2
            fileVersion = 1
        else:
            fieldsProCluster = 3
            fileVersion = 2
    assert(fileVersion is not None) # unknown cellranger version?

    for line in ifh:
        row = line.rstrip("\n\r").split(",")
        geneId = row[0]
        sym = row[1]
        clusterCount = int((len(row)-2) / fieldsProCluster)
        for clusterIdx in range(0, clusterCount):
            startField = (clusterIdx*fieldsProCluster)+2
            if fileVersion>=2:
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

    crangerToCellbrowser(options.datasetName, options.inDir, options.outDir, options.noMat)

def cbCellrangerCli_parseArgs(showHelp=False):
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] -i cellRangerDir -o outputDir - convert the cellranger output to cellbrowser format and create a cellranger.conf file

    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages")

    parser.add_option("-i", "--inDir", dest="inDir", action="store", help="input folder with the cellranger analysis output. This is the directory with the two directories 'analysis' and 'filtered_gene_bc_matrices'")
    parser.add_option("-o", "--outDir", dest="outDir", action="store", help="output directory")
    #parser.add_option("-g", "--geneSet", dest="geneSet", action="store", help="geneset, e.g. gencode28 or gencode-m13 or similar. Default: %default", default="gencode24")
    parser.add_option("-n", "--name", dest="datasetName", action="store", help="name of the dataset. No spaces or special characters.")
    parser.add_option("-m", "--noMat", dest="noMat", action="store_true", help="do not export the matrix again, saves some time if you changed something small since the last run")

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
    - %prog -i pbmc3k.h5ad -o pbmc3kScanpy - convert pbmc3k to directory with tab-separated files
    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages")

    parser.add_option("-i", "--inFile", dest="inFile", action="store",
        help="input .h5ad file. Required parameter")

    parser.add_option("", "--proc", dest="useProc", action="store_true",
        help="when exporting, do not use the raw input data, instead use the normalized and corrected matrix scanpy. This has no effect if the anndata.raw attribute is not used in the anndata object")

    parser.add_option("-o", "--outDir", dest="outDir", action="store",
        help="Output directory. Required parameter")

    parser.add_option("-n", "--name", dest="datasetName", action="store",
        help="Dataset name for generated cellbrowser.conf. If not specified, the last component of -o will be used.")

    parser.add_option("", "--htmlDir", dest="htmlDir", action="store",
        help="do not only convert to tab-sep files but also run cbBuild to"
            "convert the data and add the dataset under htmlDir")

    parser.add_option("-p", "--port", dest="port", action="store", type="int",
            help="only with --htmlDir: start webserver on port to serve htmlDir")

    parser.add_option("", "--markerField", dest="markerField", action="store",
            help="name of the marker genes field, default: %default", default="rank_genes_groups")

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

    if len(args)==0 and not options.inFile and not options.outDir:
        cbImportScanpy_parseArgs(showHelp=True)
        sys.exit(1)

    inFname = options.inFile
    outDir = options.outDir
    if None in [inFname, outDir]:
        errAbort("You need to specify at least an input .h5ad file and an output directory")

    datasetName = options.datasetName
    if datasetName is None:
        datasetName = basename(outDir.rstrip("/"))

    markerField = options.markerField

    import anndata
    ad = anndata.read_h5ad(inFname)
    scanpyToCellbrowser(ad, outDir, datasetName, skipMatrix=options.skipMatrix, useRaw=(not options.useProc),
            markerField=markerField)
    generateHtmls(datasetName, outDir)

    if options.port and not options.htmlDir:
        errAbort("--port requires --htmlDir")

    if options.htmlDir:
        build(outDir, options.htmlDir, port=options.port)

