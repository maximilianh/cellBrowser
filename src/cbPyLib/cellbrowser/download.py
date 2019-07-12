# download data from single cell repos

import logging, optparse, io, sys, os, shutil, operator
from collections import defaultdict, OrderedDict, Counter

from .cellbrowser import runGzip, openFile, errAbort, setDebug, moveOrGzip, makeDir, iterItems
from .cellbrowser import mtxToTsvGz, writeCellbrowserConf, parseIntoColumns, readMatrixAnndata, splitOnce
from .cellbrowser import generateHtmls, writePyConf

from os.path import join, basename, dirname, isfile, isdir, relpath, abspath, getsize, getmtime, expanduser

def cbGet_parseArgs(showHelp=False):
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] <repo> <accession> - import single cell data from repositories

    The only valid repo right now is: ebi

    Examples:
    - %prog ebi -a E-GEOD-100058 -o myEbiDataset - import E-GEOD-10058 to myEbiDataset/
    - %prog ebi -i E-GEOD-100058-mirror/ -o myEbiDataset - convert E-GEOD-10058 FTP mirror to myEbiDataset/
    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages")
    parser.add_option("-a", "--accession", dest="acc", action="store",
        help="Accession in source database")
    parser.add_option("-i", "--inDir", dest="inDir", action="store",
        help="Input directory")
    parser.add_option("-o", "--outDir", dest="outDir", action="store",
        help="Output directory")

    (options, args) = parser.parse_args()

    if showHelp:
        parser.print_help()
        exit(1)

    setDebug(options.debug)
    return args, options

def mirrorFtp(hostName, ftpDir, localDir):
    " mirror an ftp directory to local disk "
    makeDir(localDir)
    from ftplib import FTP
    logging.info("FTP: connecting to %s, dir %s" % (hostName, ftpDir))
    ftp = FTP(hostName)
    ftp.login()
    ftp.cwd(ftpDir)

    ls = ftp.nlst()
    count = len(ls)
    curr = 0
    logging.info("Found {} files".format(count))
    for fn in ls:
        if fn=="complete":
            continue
        curr += 1
        logging.info('Downloading file {} ... {} of {} ...'.format(fn, curr, count))
        outPath = join(localDir, fn)
        ftp.retrbinary('RETR ' + fn, open(outPath, 'wb').write)

    ftp.quit()

    complFn = join(localDir, "complete")
    open(complFn, "w").close() # = create 0-size file
    logging.info("FTP download complete.")

#def downloadFromEbi(acc, outDir):
    #" download files to outDir/orig and then create cell browser files in outDir/ "
    #for fExt in [".clusters.tsv"]:
        #fname = acc+fExt
        #fileUrl = baseUrl+"/"+fname
        #localPath = join(outDir, fname)

def parseClusters(fname):
    """ find the optimal K number in the clusters file and the cluster assignment
    return a the value of K and a dictionary with sampleName -> clusterNumber
    """
    logging.info("Parsing %s" % fname)
    bestK = None
    for line in open(fname):
        row = line.rstrip("\n\r").split("\t")
        if row[0]=="sel.K":
            sampleNames = row[2:]
            continue
        selK = row[0]
        lastRow = row
        if selK=="TRUE":
            bestK = row[1]
            clusters = row[2:]
            break

    if bestK is None:
        logging.warn("No best K found, using last K")
        bestK = lastRow[1]
        clusters = lastK[2:]

    return bestK, dict(zip(sampleNames, clusters))

def parseCondSdrf(cellIds, fname):
    """ parse the condensed sdrf into cellId -> dict of key-val """
    cellIds = set(cellIds)
    cellMeta = defaultdict(dict)
    for line in open(fname):
        #E-MTAB-7303             ERR2847884      characteristic  organism        Homo sapiens    http://purl.obolibrary.org/obo/NCBITaxon_9606
        #E-MTAB-4850		SAMEA50486668	characteristic	age	42 year
        row = line.rstrip("\n\r").split("\t")
        dsId, _, cellId, _, fieldName, val = row[:6]
        if cellId not in cellIds:
            continue
        cellMeta[cellId][fieldName] = val
    return cellMeta

def parseSdrfFieldNames(sdrfFname):
    " return list of field names in sdrf file "
    ifh = open(sdrfFname)

    #Source Name	Comment[ENA_SAMPLE]	Comment[BioSD_SAMPLE]	Characteristics[organism]	Characteristics[disease]	Characteristics[disease staging]	Characteristics[age]	Unit[time unit]	Characteristics[sex]	Characteristics[organism part]	Characteristics[cell type]	Characteristics[cell line]	Characteristics[single cell well quality]	Comment[single cell quality]	Material Type	Protocol REF	Protocol REF	Extract Name	Comment[LIBRARY_LAYOUT]	Comment[LIBRARY_SELECTION]	Comment[LIBRARY_SOURCE]	Comment[LIBRARY_STRAND]	Comment[LIBRARY_STRATEGY]	Comment[NOMINAL_LENGTH]	Comment[NOMINAL_SDEV]	Comment[ORIENTATION]	Comment[single cell isolation]	Comment[library construction]	Comment[input molecule]	Comment[primer]	Comment[end bias]	Protocol REF	Performer	Assay Name	Comment[technical replicate group]	Technology Type	Comment[ENA_EXPERIMENT]	Scan Name	Comment[SUBMITTED_FILE_NAME]	Comment[ENA_RUN]	Comment[FASTQ_URI]	Comment[SPOT_LENGTH]	Comment[READ_INDEX_1_BASE_COORD]	Factor value[cell line]	Factor Value[stimulus]	Factor Value[single cell identifier]

    fieldNames = ifh.readline().strip("\n\r").split("\t")

    newNames = []
    for fieldName in fieldNames:
        if "[" in fieldName:
            fieldName = fieldName.split("[")[1].rstrip("]")
        newNames.append(fieldName)
    return newNames

#def parseSdrf(srdfFname):
    #""" parse EBI srdf file to list of fieldNames. Remove fields with only a single value.
    #return a dict with the fields and values where all lines contain the same value.
    #Also remove filename and URI fields.
    #"""
    #logging.debug("Parsing %s" % srdfFname)
    #columns = parseIntoColumns(srdfFname)
    ##cleanList = []
    #fieldList = []
    #genInfo = {}
    #for fieldName, values in columns:
        #if "[" in fieldName:
            #fieldName = fieldName.split("[")[1].rstrip("]")

        #valSet = set(values)

        #logging.debug(fieldName)
        #logging.debug(valSet)
        #if len(valSet)==1:
            #logging.debug("Field has only a single value, adding to info dict: %s=%s" % (fieldName, values[0]))
            #genInfo[fieldName] = values[0]
        #else:
            #cleanList.append( (fieldName, values) )
            #fieldList.append(fieldName)

    #return fieldList, genInfo

def parseIdf(fname):
    " parse idf and return as dict key = value "
    logging.debug("Parsing %s" % fname)
    data = OrderedDict()
    for line in open(fname, "Ut"):
        if line.startswith("\n"):
            continue
        line = line.rstrip("\n")
        key, val = splitOnce(line, "\t")
        data[key] = val
    return data

def toStr(s):
    " prettify raw IDF tab-sep list to a normal string "
    return s.strip("\t").replace("\t", ", ")

def translateIdf(idf):
    " convert dict read from IDF to a simple key-value format for the cell browser desc.conf "
    desc = OrderedDict()
    desc["title"] = toStr(idf["Investigation Title"])
    desc["abstract"] = toStr(idf["Experiment Description"])
    desc["design"] = toStr(idf["Experimental Design"])
    desc["expType"] = toStr(idf["Comment[EAExperimentType]"])
    desc["methods"] = toStr(idf["Protocol Description"])
    desc["curator"] = toStr(idf["Comment[EACurator]"])+", EBI Single Cell Expression Atlas"
    desc["submitter"] = toStr(idf["Person First Name"])+" "+toStr(idf["Person Last Name"])+" <"+toStr(idf["Person Email"])+">"
    desc["submission_date"] = toStr(idf["Public Release Date"])
    if "PubMed ID" in idf:
        desc["pmid"] = toStr(idf["PubMed ID"])
    if "Publication DOI" in idf:
        desc["doi"] = toStr(idf["Publication DOI"])

    if "Comment[ArrayExpressAccession]" in idf:
        desc["arrayexpress"] = toStr(idf["Comment[ArrayExpressAccession]"])

    if "Comment [SecondaryAccession]" in idf:
        desc["ena_project"] = toStr(idf["Comment [SecondaryAccession]"])

    return desc

def convertCoords(coordFname, newCoordFname):
    " simply reorder columns"
    ofh = open(newCoordFname, "w")
    for line in open(coordFname):
        row = line.rstrip("\n\r").split("\t")
        newRow = (row[2], row[0], row[1])
        ofh.write("\t".join(newRow))
        ofh.write("\n")
    ofh.close()
    logging.info("Wrote %s" % newCoordFname)

def convertMarkers(inFname, outFname):
    " reorder columns and reheader and return list of one top marker per cluster "
    #names	groups	scores	logfc	pvals	pvals_adj
    #ENSG00000112096	0	9.120421	1.4149776	3.5174541399729827e-15	5.52275474517158e-11
    ofh = open(outFname, "w")
    clusterToMarkers = defaultdict(list)
    for line in open(inFname):
        row = line.rstrip("\n\r").split("\t")
        #cluster, pval, auroc, gene = row
        geneId, cluster, score, logfc, pval, pvalAdj = row
        if row[0]=="names":
            newRow = ["Cluster", "Gene", "p-Value Adj.", "p-Value", "logFC", "score"]
        else:
            newRow = (row[1], row[0], row[5], row[4], row[3], row[2])
            clusterToMarkers[cluster].append( (float(pvalAdj), geneId) )

        ofh.write("\t".join(newRow))
        ofh.write("\n")

    ofh.close()
    logging.info("Wrote %s" % outFname)

    topMarkers = {}
    for cluster, pValGenes in iterItems(clusterToMarkers):
        pValGenes.sort()
        topGene = pValGenes[0][1]
        topMarkers[cluster] = topGene

    return topMarkers

def writeQuickGenes(topGenes, outFname):
    " given cluster -> gene, create a quickGenes file "
    ofh = open(outFname, "w")
    for cluster, gene in iterItems(topGenes):
        ofh.write("%s\tCluster %s\n"% (gene, cluster))
    ofh.close()
    logging.info("Wrote %s" % outFname)

def makeBasicCbConf(name, shortLabel, hasInferred, matrixFname="exprMatrix.tsv.gz"):
    " just fill a basic conf dict and return it "
    c = OrderedDict()
    c["name"] = name
    c["shortLabel"] = shortLabel
    c["exprMatrix"] = matrixFname
    c["meta"] = "meta.tsv"
    c["#priority"] = "10"
    c["tags"] = ["ebi"]
    c["enumFields"] = "cluster"

    if hasInferred:
        clusterField = "inferred cell type"
        c["acroFname"] = "acronyms.tsv"
    else:
        clusterField = "cluster"

    c["clusterField"] = clusterField
    c["labelField"] = clusterField

    c["#unit"] = "TPM"
    c["coords"] = [{"file":"tsne_perp25.coords.tsv", "shortLabel" : "t-SNE Perp=25"}]
    c["markers"] = [{"file":"markers.tsv", "shortLabel":"Cluster-specific genes"}]
    c["quickGenesFile"] = "quickGenes.tsv"
    return c

def writeMeta(clusters, cellIds, sdrfFields, boringFields, cellMeta, metaFname):
    " write a meta.tsv file, add in the clusters "
    #idIdx = None
    #fieldNames = []
    #for fieldIdx, (fieldName, vals) in enumerate(metaColumns):
        #if fieldName=="ENA_RUN":
            #idIdx = fieldIdx
            #break
        #fieldNames.append(fieldName)

    #if idIdx is None:
        #errAbort("Couldn't find the main ID field in the sdrf file. Possible field names: %s" % ",".join(fieldNames))

    # basic check of the ID field
    #idName = metaColumns[idIdx][0]
    #idList = metaColumns[idIdx][1]
    #if len(idList)!=len(set(idList)):
        #idCounts = Counter(idList)
        #errAbort("field %s was identified as the ID field, but it is not unique" % idCounts)

    # put the id field first
    #newOrder = []
    #newOrder.append(metaColumns[idIdx])
    #for i in range(0, idList):
        #if i==idIdx or i==0: # remove the first source name field
            #continue
        #newOrder.append(metaColumns[i])

    # add the cluster field
    #clustVals = []
    #for cellId in newOrder[0][1]:
        #clustId = clusters[cellId]
        #clustVals.append(clustId)

    #newOrder.insert(1, ("Cluster", clustVals))

    #headers = [f[0] for f in newOrder]

    #fieldValues = [x[1] for x in newOrder]
    #for i in range(0, len(idList)):
        #row = []
        #for field in fieldValues:
            #row.append(field[i])
        #ofh.write("\t".join(row))
        #ofh.write("\n")
    #ofh.close()
    #logging.info("Wrote %s" % newCoordFname)

    logging.debug("Fields in the SDRF: %s" % sdrfFields)
    logging.debug("Boring fields: %s" % boringFields.keys())

    allFields = set()
    for cellId in cellIds:
        meta = cellMeta[cellId]
        allFields.update(meta.keys())
    logging.debug("All fields from the condensed: %s" % allFields)

    fieldOrder = []
    doneFields = set()
    for fn in sdrfFields:
        if fn not in boringFields and fn in allFields and fn not in doneFields:
            fieldOrder.append(fn)
            doneFields.add(fn)

    logging.debug("final field order: %s" % fieldOrder)
    assert(len(set(fieldOrder))==len(fieldOrder)) # fields must not appear twice

    ofh = open(metaFname, "w")
    ofh.write("cell ID\tcluster\t")
    ofh.write("\t".join(fieldOrder))
    ofh.write("\n")

    for cellId, meta in iterItems(cellMeta):
        row = [cellId]
        clustId = clusters.get(cellId, "Unknown")
        row.append(clustId)
        for fieldName in fieldOrder:
            val = meta.get(fieldName, "")
            row.append(val)
        ofh.write("\t".join(row))
        ofh.write("\n")
    ofh.close()
    logging.info("Wrote %s" % ofh.name)

def findBoringFields(cellMetaDict):
    """ return a dict of all fieldNames that have only a single value. The keys of this dict are the field names, 
    the values of this dict are the single values."""
    fieldValues = defaultdict(set)
    for cellId, cellMeta in iterItems(cellMetaDict):
        for fieldName, val in iterItems(cellMeta):
            fieldValues[fieldName].add(val)

    boringFields = {}
    for fieldName, vals in iterItems(fieldValues):
        if fieldName.endswith("_FILE_NAME") or fieldName.endswith("_URI") or fieldName=="single cell identifier":
            logging.debug("Skipping field %s" % fieldName)
            boringFields[fieldName] = None
            continue
        if len(vals)==1:
            boringFields[fieldName] = list(vals)[0]
    return boringFields

def writeAcronyms(cellMeta, acroFname):
    " extract acronyms from inferred cell name field "
    ofh = open(acroFname, "w")
    for cellId, metaDict in iterItems(cellMeta):
        cellName = metaDict.get("inferred cell type")
        if cellName and "(" in cellName:
            longName, shortName = cellName.split("(")
            longName = longName.strip()
            shortName = shortName.strip(")")
            metaDict["inferred cell type"] = shortName
            ofh.write(shortName+"\t"+longName+"\n")
    ofh.close()
    logging.info("Wrote %s" % ofh.name)
    return cellMeta

def convertFromEbi(origDir, acc, outDir):
    " convert single cell expression files to cell browser format "
    # cluster is a special meta data attribute, parse it into a dict
    clusterFname = join(origDir, acc+".clusters.tsv")
    bestK, clusters = parseClusters(clusterFname)

    # best K determines marker file name
    markerFname = join(origDir, acc+".marker_genes_%s.tsv" % bestK)
    newMarkerFname = join(outDir, "markers.tsv")
    topGenes = convertMarkers(markerFname, newMarkerFname)

    quickGeneFname = join(outDir, "quickGenes.tsv")
    writeQuickGenes(topGenes, quickGeneFname)

    # make sure we use the same cellID order in meta and matrix
    sampleNameFname = join(origDir, acc+".aggregated_filtered_counts.mtx_cols")
    sampleNames = open(sampleNameFname).read().splitlines()

    # meta data
    sdrfFname = join(origDir, acc+".sdrf.txt")
    #fieldNames, sdrfInfo = parseSdrf(sdrfFname)
    fieldNames = parseSdrfFieldNames(sdrfFname) # we're doing this to get them in the proper order
    condSdrfFname = join(origDir, acc+".condensed-sdrf.tsv")
    cellMeta = parseCondSdrf(sampleNames, condSdrfFname)
    boringFields = findBoringFields(cellMeta)

    # default perp is 25
    coordFname = join(origDir, acc+".tsne_perp_25.tsv")
    newCoordFname = join(outDir, "tsne_perp25.coords.tsv")
    convertCoords(coordFname, newCoordFname)

    # dataset descriptor
    idfFname = join(origDir, acc+".idf.txt")
    idfInfo = parseIdf(idfFname)
    datasetDesc = translateIdf(idfInfo)
    datasetDesc["custom"] = boringFields
    descFname = join(outDir, "desc.conf")
    writePyConf(datasetDesc, descFname)

    # write acronyms and add the short names to the meta data
    hasInferred = ("inferred cell type" in fieldNames)
    if hasInferred:
        acroFname = join(outDir, "acronyms.tsv")
        cellMeta = writeAcronyms(cellMeta, acroFname)

    # write the meta 
    metaFname = join(outDir, "meta.tsv")
    writeMeta(clusters, sampleNames, fieldNames, boringFields, cellMeta, metaFname)

    # cellbrowser.conf file
    confFname = join(outDir, "cellbrowser.conf")
    conf = makeBasicCbConf(acc, datasetDesc["title"], hasInferred)
    writePyConf(conf, confFname)

    mtxName = join(origDir, acc+".aggregated_filtered_counts.mtx.gz")
    sampleNamesFname = join(origDir, acc+".aggregated_filtered_counts.mtx_cols")
    geneNamesFname = join(origDir, acc+".aggregated_filtered_counts.decorated.mtx_rows")
    outMatName = join(outDir, "exprMatrix.tsv.gz")
    #runGzip(inMatName, outMatName)
    if not isfile(outMatName):
        mtxToTsvGz(mtxName, geneNamesFname, sampleNamesFname, outMatName)

def cbGetCli():
    " run downloaders "
    args, options = cbGet_parseArgs()

    if len(args)<1:
        cbGet_parseArgs(showHelp=True)
        sys.exit(1)

    outDir = options.outDir
    if outDir is None:
        outDir = acc

    cmd = args[0]

    cmds = ["ebi"]

    if cmd=="ebi":
        acc = options.acc
        if options.inDir:
            origDir = options.inDir
        else:
            origDir = join(outDir, "orig")
            flagFname = join(origDir, "complete")
            logging.debug("Checking %s" % flagFname)
            if not isfile(flagFname):
                hostName = "ftp.ebi.ac.uk"
                ftpDir = "/pub/databases/microarray/data/atlas/sc_experiments/%s" % acc
                mirrorFtp(hostName, ftpDir, origDir)
        convertFromEbi(origDir, acc, outDir)
    else:
        errAbort("Repository %s is not valid. Valid repos are: %s" % (cmd, ", ".join(cmds)))
