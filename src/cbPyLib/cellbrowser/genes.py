# functions to guess the gene model release given a list of gene IDs
# tested on python3 and python2
import logging, sys, optparse, string, glob, gzip
from io import StringIO
#from urllib.request import urlopen
from urllib.request import Request, urlopen
from collections import defaultdict
from os.path import join, basename, dirname, isfile

from .cellbrowser import sepForFile, getStaticFile, openFile, splitOnce, setDebug, getStaticPath
from .cellbrowser import getGeneSymPath, downloadUrlLines, getSymToGene, getGeneBedPath, errAbort, iterItems

# ==== functions =====
def cbGenes_parseArgs():
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] command - download gene model files and auto-detect the version.
    syms <geneType> - Download a table with geneId <-> symbol from a database.
    locs <assembly> <geneType> - Download a gene model file from UCSC, pick one transcript per gene and save to ~/cellbrowserData/genes/<db>.<geneType>.bed.
    guess - Guess Ensembl/Gencode version given a file with list of gene IDs. Reads the first tab-sep field from filename and reports best gencode version.
    
    Examples:
    %prog syms gencode-34
    %prog locs hg38 gencode-34
    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true", help="show debug messages")
    (options, args) = parser.parse_args()

    if args==[]:
        parser.print_help()
        exit(1)

    setDebug(options.debug)
    return args, options

# ----------- main --------------
def parseSignatures(org, geneIdType):
    " return dict with gene release -> list of unique signature genes "
    ret = {}
    logging.info("Parsing gencode release signature genes")
    fname = getStaticFile("genes/gencode-%s.guessVersion.%s.tsv.gz" % (org, geneIdType))
    logging.info("Parsing %s" % fname)
    genes = set()
    verToGenes = {}
    for line in openFile(fname):
        if line.startswith("#"):
            continue
        version, geneIds = line.rstrip("\n").split('\t')
        geneIds = geneIds.split(",")
        verToGenes[version] = geneIds
    return verToGenes
        
def guessGeneIdType(genes):
    " return tuple organism / identifier type "
    gene1 = list(genes)[0]
    if gene1.startswith("ENSG"):
        return "human", "acc"
    if gene1.startswith("ENSMUS"):
        return "mouse", "acc"
    if gene1.upper()==gene1:
        return "human", "sym"
    else:
        return "mouse", "sym"

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
        fileGenes.add(geneId.split('.')[0].split("|")[0])
    logging.info("Read %d genes" % len(fileGenes))
    return fileGenes

def guessGencodeVersion(fileGenes, signGenes):
    logging.info("Number of genes that are specific for gene model release:")
    diffs = []
    for version, uniqGenes in signGenes.items():
        intersection = list(fileGenes.intersection(uniqGenes))
        infoStr = "release "+version+": %d out of %d" % (len(intersection), len(uniqGenes))
        if len(intersection)!=0:
            expStr = ", ".join(intersection[:5])
            infoStr += (" e.g. "+ expStr)
        logging.info(infoStr)
        diffs.append((len(intersection), version))

    diffs.sort(reverse=True)
    bestVersion = diffs[0][1]
    return bestVersion

def guessGencode(fname):
    inGenes = set(parseGenes(fname))
    org, geneType = guessGeneIdType(inGenes)
    logging.info("Looks like input gene list is from organism %s, IDs are %s" % (org, geneType))
    signGenes = parseSignatures(org, geneType)
    bestVersion = guessGencodeVersion(inGenes, signGenes)
    print("Best %s Gencode release\t%s" % (org, bestVersion))

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
    assert(release.isdigit())
    db = "hg38"
    url = "https://hgdownload.cse.ucsc.edu/goldenPath/%s/database/wgEncodeGencodeAttrsV%s.txt.gz" %  (db, release)
    logging.info("Downloading %s" % url)
    doneIds = set()
    for line in downloadUrlLines(url):
        row = line.rstrip("\n").split("\t")

        # strip the .x version off from the geneId
        if doTransGene:
            # key = transcript ID, val is geneId
            key = row[4]
            val = row[0]
            val = val.split('.')[0]
        else:
            # key = geneId, val is symbol
            key = row[0]
            key = key.split('.')[0]
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
    " build a BED file with a 'canonical' transcript for every gene "
    if geneType.startswith("gencode"):
        release = geneType.split("-")[1]
        rows = iterGencodeBed(db, release)
    else:
        errAbort("Unknown gene model type: %s" % geneType)

    outFname = getStaticPath(getGeneBedPath(db, geneType))
    writeRows(rows, outFname)

def cbGenesCli():
    args, options = cbGenes_parseArgs()

    command = args[0]
    if command=="guess":
        fname = args[1]
        guessGencode(fname)
    elif command=="syms":
        geneType = args[1]
        buildSymbolTable(geneType)
    elif command=="locs":
        options = args[1:]
        buildLocusBed(*options)
    else:
        errAbort("Unrecognized command: %s" % command)

