# functions to guess the gene model release given a list of gene IDs
# tested on python3 and python2
import logging, sys, optparse, string, glob
from collections import defaultdict
from os.path import join, basename, dirname, isfile

import cellbrowser

# ==== functions =====
def cbGuessGencode_parseArgs():
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("usage: %prog [options] filename - guess Gencode version of a list of gene IDs. Reads the first tab-sep field from filename and reports best gencode version.")

    parser.add_option("-d", "--debug", dest="debug", action="store_true", help="show debug messages")
    (options, args) = parser.parse_args()

    if args==[]:
        parser.print_help()
        exit(1)

    if options.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    return args, options

# ----------- main --------------
def parseSignatures(org, geneIdType):
    " return dict with gene release -> list of unique signature genes "
    ret = {}
    logging.info("Parsing gencode release signature genes")
    fname = cellbrowser.getStaticFile("genes/gencode-%s.guessVersion.%s.tsv.gz" % (org, geneIdType))
    logging.info("Parsing %s" % fname)
    genes = set()
    verToGenes = {}
    for line in cellbrowser.openFile(fname):
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
    sep = cellbrowser.sepForFile(fname)
    for line in cellbrowser.openFile(fname):
        if not headDone:
            headDone = True
            continue
        geneId = cellbrowser.splitOnce(line[:50], sep)[0]
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

def cbGuessGencodeCli():
    args, options = cbGuessGencode_parseArgs()

    fname = args[0]
    inGenes = set(parseGenes(fname))
    org, geneType = guessGeneIdType(inGenes)
    logging.info("Looks like input gene list is from organism %s, IDs are %s" % (org, geneType))
    signGenes = parseSignatures(org, geneType)
    bestVersion = guessGencodeVersion(inGenes, signGenes)
    print("Best %s Gencode release\t%s" % (org, bestVersion))
