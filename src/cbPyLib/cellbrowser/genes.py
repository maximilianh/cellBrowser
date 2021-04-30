# functions to guess the gene model release given a list of gene IDs
# tested on python3 and python2
import logging, sys, optparse, string, glob, gzip, json
from io import StringIO
#from urllib.request import urlopen
from urllib.request import Request, urlopen
from collections import defaultdict
from os.path import join, basename, dirname, isfile

from .cellbrowser import sepForFile, getStaticFile, openFile, splitOnce, setDebug, getStaticPath
from .cellbrowser import getGeneSymPath, downloadUrlLines, getSymToGene, getGeneBedPath, errAbort, iterItems
from .cellbrowser import findCbData, readGeneSymbols, getGeneJsonPath, getDownloadsUrl

# ==== functions =====
def cbGenes_parseArgs():
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] command - download gene model files and auto-detect the version.

    Commands:
    fetch <geneType> - download pre-built geneId -> symbol table from UCSC
    fetch <assembly>.<geneType> - download pre-built gene models and symbol table from UCSC
    build <assembly>.<geneType> - Download a gene model file from UCSC, pick one transcript per gene and save to ~/cellbrowserData/genes/<db>.<geneType>.bed.gz and <geneType>.symbols.tsv.gz

    Run "fetch" or "build" without arguments to list the available files at UCSC.

    ls - list all available (built or downloaded)  gene models on this machine

    guess <inFile> <organism> - Guess best gene type. Reads the first tab-sep field from inFile and prints genetypes sorted by % of matching unique IDs to inFile.

    Examples:
    %prog fetch                   # show the files that are available
    %prog fetch gencode-34        # geneId -> symbol mapping for human gencode relase 34
    %prog fetch hg38.gencode-34   # gene -> chrom mapping for human gencode relase 34
    %prog build                   # show the files that are available
    %prog build mm10 gencode-M25
    %prog ls
    %prog guess genes.txt mouse
    %prog index # only used at UCSC to prepare the files for 'guess'
    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true", help="show debug messages")
    (options, args) = parser.parse_args()

    if args==[] or (args[0]=="guess" and len(args)==1):
        parser.print_help()
        exit(1)

    setDebug(options.debug)
    return args, options

# ----------- main --------------
def parseSignatures(org, geneIdType):
    " return dict with gene release -> list of unique signature genes "
    ret = {}
    logging.info("Parsing gencode release signature genes")
    fname = getStaticFile("genes/%s.%s.unique.tsv.gz" % (org, geneIdType))
    logging.info("Parsing %s" % fname)
    genes = set()
    verToGenes = {}
    for line in openFile(fname):
        if line.startswith("#"):
            continue
        version, geneIds = line.rstrip("\n").split('\t')
        geneIds = set(geneIds.split("|"))
        verToGenes[version] = geneIds

    return verToGenes
        
def guessGeneIdType(genes):
    " return tuple organism / identifier type "
    logging.debug("Trying to guess organism and identifier type (syms or ids)")
    gene1 = list(genes)[0]
    if gene1.startswith("ENSG"):
        return "human", "ids"
    if gene1.startswith("ENSMUS"):
        return "mouse", "ids"

    upCount = 0
    for g in genes:
        if g.isupper():
            upCount += 1

    logging.debug("%d of %d genes are uppercase" % (upCount, len(genes)))
    if upCount/float(len(genes)) > 0.8:
        return "human", "syms"
    else:
        return "mouse", "syms"

def parseGenes(fname):
    " return gene IDs in column 1 of file "
    fileGenes = set()
    headDone = False
    logging.info("Parsing first column from %s" % fname)
    sep = sepForFile(fname)
    for line in openFile(fname):
        if not headDone:
            headDone = True
            continue
        geneId = splitOnce(line[:50], sep)[0]
        geneId = geneId.strip("\n").strip("\r").strip()
        #fileGenes.add(geneId.split('.')[0].split("|")[0])
        fileGenes.add(geneId.split("|")[0])
    logging.info("Read %d genes" % len(fileGenes))
    return fileGenes

def guessGencodeVersion(fileGenes, signGenes):
    logging.info("Number of genes that are only in a gene model release:")
    #diffs = []
    infos = []
    for version, uniqGenes in signGenes.items():
        intersection = list(fileGenes.intersection(uniqGenes))
        share = 100.0 * (float(len(intersection)) / len(uniqGenes))
        intLen = len(intersection)
        geneCount = len(uniqGenes)
        infoStr = "release "+version+": %0.2f%%, %d out of %d" % (share, len(intersection), len(uniqGenes))
        if len(intersection)!=0:
            expStr = ", ".join(intersection[:5])
            infoStr += (" e.g. "+ expStr)
        #logging.info(infoStr)
        infos.append((share, version, intLen, geneCount, infoStr))

        #diffs.append((len(intersection), version))

    infos.sort(reverse=True)
    bestVersion = infos[0][1]

    for info in infos:
        share, version, intLen, geneCount, infoStr = info
        print(infoStr)

    return bestVersion

def guessGencode(fname, org):
    inGenes = set(parseGenes(fname))
    guessOrg, geneType = guessGeneIdType(inGenes)
    if org is None:
        org = guessOrg
    logging.info("Looks like input gene list is from organism %s, IDs are %s" % (org, geneType))
    signGenes = parseSignatures(org, geneType)
    bestVersion = guessGencodeVersion(inGenes, signGenes)
    print("Best %s Gencode release\t%s" % (org, bestVersion))

    allIds = readGeneSymbols(bestVersion)
    if geneType=="syms":
        allIds = allIds.values()
    notFoundIds = inGenes - set(allIds)
    print("%d of the genes in the input are not part of %s" % (len(notFoundIds), bestVersion))
    print("Examples: %s" % " ".join(list(notFoundIds)[:50]))

def buildSymbolTable(geneType):
    if geneType.startswith("gencode"):
        release = geneType.split("-")[1]
        rows = iterGencodePairs(release)
    else:
        errAbort("unrecognized gene type '%s'" % geneType)

    outFname = getStaticPath(getGeneSymPath(geneType))
    writeRows(rows, outFname)

def iterGencodePairs(release, doTransGene=False):
    " generator, yields geneId,symbol or transId,geneId pairs for a given gencode release"
    # e.g. trackName = "wgEncodeGencodeBasicV34"
    #attrFname = trackName.replace("Basic", "Attrs").replace("Comp", "Attrs")
    #assert(release[1:].isdigit())
    db = "hg38"
    if release[0]=="M":
        db = "mm10"
    if release in ["7", "14", "17", "19"] or "lift" in release:
        db = "hg19"
    url = "https://hgdownload.cse.ucsc.edu/goldenPath/%s/database/wgEncodeGencodeAttrsV%s.txt.gz" %  (db, release)
    logging.info("Downloading %s" % url)
    doneIds = set()
    for line in downloadUrlLines(url):
        row = line.rstrip("\n").split("\t")

        if doTransGene:
            # key = transcript ID, val is geneId
            key = row[4]
            val = row[0]
            val = val
        else:
            # key = geneId, val is symbol
            key = row[0]
            key = key
            val = row[1]

        if key not in doneIds:
            yield key, val
            doneIds.add(key)

def iterGencodeBed(db, release):
    " generator, yields a BED12+1 with a 'canonical' transcript for every gencode comprehensive gene "
    transToGene = dict(iterGencodePairs(release, doTransGene=True))

    url = "http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/wgEncodeGencodeCompV%s.txt.gz" % (db, release)
    logging.info("Downloading %s" % url)
    geneToTransList = defaultdict(list)
    for line in downloadUrlLines(url):
        row = tuple(line.split('\t'))
        transId = row[1]
        geneId = transToGene[transId]
        score = int(''.join(c for c in geneId if c.isdigit())) # extract only the xxx part of the ENSGxxx ID
        geneToTransList[geneId].append( (score, row) )

    logging.info("Picking one transcript per gene")
    for geneId, transList in iterItems(geneToTransList):
        transList.sort() # prefer older transcripts
        canonTransRow = transList[0][1]
        binIdx, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, cdsStartStat, cdsEndStat, exonFrames = canonTransRow
        blockStarts = []
        blockLens = []
        for exonStart, exonEnd in zip(exonStarts.split(","), exonEnds.split(",")):
            if exonStart=="":
                continue
            blockSize = int(exonEnd)-int(exonStart)
            blockStarts.append(exonStart)
            blockLens.append(str(blockSize))
        newRow = [chrom, txStart, txEnd, geneId, score, strand, cdsStart, cdsEnd, exonCount, ",".join(blockLens), ",".join(blockStarts), name2]
        yield newRow

def writeRows(rows, outFname):
    with openFile(outFname, "wt") as ofh:
        for row in rows:
            ofh.write("\t".join(row))
            ofh.write("\n")
    logging.info("Wrote %s" % outFname)

def buildLocusBed(db, geneType):
    " build a BED file with a 'canonical' transcript for every gene and a json file for it "
    if geneType.startswith("gencode"):
        release = geneType.split("-")[1]
        rows = iterGencodeBed(db, release)
    else:
        errAbort("Unknown gene model type: %s" % geneType)

    outFname = getStaticPath(getGeneBedPath(db, geneType))
    writeRows(rows, outFname)

    jsonFname = getStaticPath(getGeneJsonPath(db, geneType))
    bedToJson(db, geneType, jsonFname)

def listModelsLocal():
    " print all gene models on local machine "

    dataDir = join(findCbData(), "genes")
    print("Local cell browser genes data directory: %s" % dataDir)
    fnames = glob.glob(join(dataDir, "*.symbols.tsv.gz"))
    names = [basename(x).split(".")[0] for x in fnames]
    print("Installed gene/symbol mappings:")
    print("\n".join(names))
    print()

    fnames = glob.glob(join(dataDir, "*.bed.gz"))
    names = [basename(x).replace(".bed.gz","") for x in fnames]
    print("Installed gene/chrom-location mappings:")
    print("\n".join(names))

def iterBedRows(db, geneIdType):
    " yield BED rows of gene models of given type "
    fname = getStaticPath(getGeneBedPath(db, geneIdType))
    logging.info("Reading BED file %s" % fname)
    with openFile(fname) as ofh:
        for line in ofh:
            row = line.rstrip("\n\r").split("\t")
            yield row

def parseApacheDir(lines):
    fnames = []
    for l in lines:
        hrefCount = l.count("<a href=")
        if hrefCount==1:
            if "Parent Directory<" in l: 
                continue
            fname = l.split('<a href="')[1].split('"')[0]
            fnames.append(fname)
    return fnames

def listModelRemoteFetch():
    " print all gene models that can be downloaded "
    url = join(getDownloadsUrl(), "genes")
    lines = downloadUrlLines(url)
    fnames = parseApacheDir(lines)
    geneFnames = [f.replace(".bed.gz","") for f in fnames if f.endswith(".bed.gz")]
    symFnames = [f.replace(".symbols.tsv.gz", "") for f in fnames if f.endswith(".symbols.tsv.gz")]

    sep = "\n"

    print("Pre-built gene model mapping files available for 'fetch' at %s" % url)
    print(sep.join(geneFnames))
    #for g in geneFnames:
        #print(g.replace(".bed.gz",""))

    print()
    print("Pre-built geneId/symbol tables available for 'fetch' at %s" % url)
    print(sep.join(symFnames))
    #for g in symFnames:
        #print(g.replace(".symbols.tsv.gz", ""))

def listModelRemoteBuild():
    sep = "\n"
    urls = [("hg38", "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/"),
            ("mm10", "https://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/"),
            ("hg19", "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/")
            ]

    allNames = defaultdict(list)
    for db, url in urls:
        print()
        print("Files available for 'build' for assembly %s (%s)" % (db, url))
        lines = downloadUrlLines(url)
        fnames = parseApacheDir(lines)
        geneFnames = [x for x in fnames if x.startswith("wgEncodeGencodeAttrs") and x.endswith(".txt.gz")]
        relNames = [x.replace("wgEncodeGencodeAttrsV", "gencode-").replace(".txt.gz", "") for x in geneFnames]
        allNames[db].extend(relNames)
        print(sep.join(relNames))

    #for db, names in allNames.items():
        #for name in names:
            ##print("%s\t%s" % (db, name))
            #print(name)

def keepOnlyUnique(dictSet):
    """ give a dict with key -> set, return a dict with key -> set, but only with elements in the set that
    that don't appear in any other set
    """
    uniqVals = {}
    for key1, origVals in dictSet.items():
        vals = set(list(origVals))

        for key2 in dictSet.keys():
            if key1==key2:
                continue
            vals = vals - dictSet[key2]
        uniqVals[key1] = vals

    setList = list(dictSet.values())
    allCommon = set.intersection(*setList)
    return uniqVals, len(allCommon)

def writeUniqs(dictSet, outFname):
    " wrote to output file in format <key>tab<comma-sep-list of vals> "
    logging.info("Writing to %s" % outFname)
    with openFile(outFname, "wt") as ofh:
        for key, vals in dictSet.items():
            ofh.write("%s\t%s\n" % (key, "|".join(vals)))

def uniqueIds(org):
    """ find unique identifiers in all symbols and geneIds of infileMask and write to
    outBase.{syms,ids}.unique.syms.tsv.gz
    """
    logging.info("Processing: %s" % org)
    infileMask = "gencode*.symbols.tsv.gz"
    dataDir = join(findCbData(), "genes")
    fnames = glob.glob(join(dataDir, infileMask))
    allSyms = {}
    allIds = {}
    for fname in fnames:
        baseName = basename(fname)
        if "lift" in baseName or "mouse" in baseName or "human" in baseName:
            continue
        if org=="human" and "M" in baseName:
            continue
        if org=="mouse" and not "M" in baseName:
            continue
        geneType = basename(fname).split(".")[0]
        logging.info("Reading %s" % fname)

        syms = set()
        ids = set()
        for line in openFile(fname):
            row = line.rstrip("\n").split("\t")
            geneId, sym = row[:2]
            syms.add(sym)
            ids.add(geneId)
        allSyms[geneType] = syms
        allIds[geneType] = ids

    # force refseq into this
    syms = []
    fname = getStaticFile("genes/entrez-%s.symbols.tsv.gz" % (org))
    for line in openFile(fname):
        if line.startswith("#"):
            continue
        geneId, sym = line.rstrip("\n").split("\t")
        syms.append(sym)
    #verToGenes["refseq"] = set(syms)
    allSyms["entrez"] = set(syms)

    logging.info("Finding unique values")
    uniqSyms, commonSyms = keepOnlyUnique(allSyms)
    uniqIds, commonIds = keepOnlyUnique(allIds)
    logging.info("%d symbols and %d geneIds are shared among all releases" % (commonSyms, commonIds))

    writeUniqs(uniqSyms, join(dataDir, org+".syms.unique.tsv.gz"))
    writeUniqs(uniqIds, join(dataDir, org+".ids.unique.tsv.gz"))

def bedToJson(db, geneIdType, jsonFname):
    " convert BED file to more compact json file: chrom -> list of (start, end, strand, gene) "
    geneToSym = readGeneSymbols(geneIdType)

    # index transcripts by gene
    bySym = defaultdict(dict)
    for row in iterBedRows(db, geneIdType):
        chrom, start, end, geneId, score, strand = row[:6]
        sym = geneToSym[geneId]
        start = int(start)
        end = int(end)
        transLen = end-start
        bySym[sym].setdefault(chrom, []).append( (transLen, start, end, strand, geneId) )

    symLocs = defaultdict(list)
    for sym, chromDict in bySym.items():
        for chrom, transList in chromDict.items():
            transList.sort(reverse=True) # take longest transcript per chrom
            _, start, end, strand, transId = transList[0]
            symLocs[chrom].append( (start, end, strand, sym) )

    sortedLocs = {}
    for chrom, geneList in symLocs.items():
        geneList.sort()
        sortedLocs[chrom] = geneList

    ofh = open(jsonFname, "wt")
    outs = json.dumps(sortedLocs)
    #md5 = hashlib.md5(outs.encode("utf8")).hexdigest()[:10]
    ofh.write(outs)
    ofh.close()
    logging.info("Wrote %s" % jsonFname)

    #fileInfo[code] = {"label":label, "file" : jsonFname, "md5" :md5}

def buildGuessIndex():
    " read all gene model symbol files from the data dir, and output <organism>.unique.tsv.gz "
    dataDir = join(findCbData(), "genes")
    uniqueIds("human")
    uniqueIds("mouse")

def fetch(fileDesc):
    " download symbol or gene files to local dir "
    if "." in fileDesc:
        # user wants a gene model file
        ext = "bed.gz"
    else:
        ext = "symbols.tsv.gz"
    fname = getStaticFile("genes/%s.%s" % (fileDesc, ext), verbose=True)
    return

def cbGenesCli():
    args, options = cbGenes_parseArgs()

    command = args[0]
    if command=="guess":
        fname = args[1]
        org = None
        if len(args)==3:
            org = args[2]
        guessGencode(fname, org)

    elif command=="fetch":
        if len(args)==1:
            listModelRemoteFetch()
        else:
            arg = args[1]
            fetch(arg)

    elif command=="syms": # undocumented
        geneType = args[1]
        buildSymbolTable(geneType)

    elif command=="build":
        if len(args)==1:
            listModelRemoteBuild()
        else:
            db, geneType = args[1:]
            buildSymbolTable(geneType)
            buildLocusBed(db, geneType)

    elif command=="ls":
        listModelsLocal()

    elif command=="index":
        buildGuessIndex()

    elif command=="json": # undocumented
        db, geneType, outFname = args[1:]
        bedToJson(db, geneType, outFname)
    else:
        errAbort("Unrecognized command: %s" % command)

