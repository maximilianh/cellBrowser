#!/usr/bin/env python

# this library mostly contains functions that convert tab-sep files
# (=single cell expression matrix and meta data) into the binary format that is read by the
# javascript viewer cbWeb/js/cellbrowser.js and cbData.js.
# Helper functions here allow importing data from other tools, e.g. cellranger or scanpy.

# requires at least python2.6, version tested was 2.6.6
# should work with python2.5, not tested
# works on python3, version tested was 3.6.5
# all functions related to cbScanpy() require python3, as scanpy requires at least python3

import logging, sys, optparse, struct, json, os, string, shutil, gzip, re, unicodedata
import zlib, math, operator, doctest, copy, bisect, array, glob, io, time, subprocess
import hashlib, timeit, datetime, keyword
from distutils import spawn
from collections import namedtuple, OrderedDict
from os.path import join, basename, dirname, isfile, isdir, relpath, abspath, getsize, getmtime, expanduser
from time import gmtime, strftime

try:
    # python2
    from urlparse import urljoin
    from urllib2 import urlopen
except:
    # python3
    from urllib.parse import urljoin
    from urllib.request import urlopen

try:
    # python2.7+
    from collections import defaultdict
    from collections import Counter
except:
    # python2.6 has no defaultdict or Counter yet
    from backport_collections import defaultdict # error? -> pip2 install backport-collections
    from backport_collections import Counter # error? -> pip2 install backport-collections

# We do not require numpy but numpy is around 30-40% faster in serializing arrays
# So use it if it's present
numpyLoaded = False
try:
    import numpy as np
except:
    numpyLoaded = False
    logging.debug("Numpy could not be loaded. This is fine but matrix export may be 1/3 slower.")

# older numpy versions don't have tobytes()
if numpyLoaded:
    try:
        np.ndarray.tobytes
    except:
        numpyLoaded = False
        logging.error("Numpy version too old. Falling back to normal Python array handling.")

isPy3 = False
if sys.version_info >= (3, 0):
    isPy3 = True

# directory to static data files, e.g. gencode tables
# By default, this is ~/cbData, or alternatively /usr/local/share/cellbrowser 
# or the directory in the environment variable CBDATA, see findCbData()
#dataDir = join(dirname(__file__), "..", "cbData")
dataDir = None

# the default html dir, used if the --htmlDir option is set but empty
# this variable is initialized below (no forward declaration in Python)
# just before cbBuild_parseArgs
defOutDir = None

CBHOMEURL = "https://cells.ucsc.edu/downloads/cellbrowserData/"
#CBHOMEURL = "http://localhost/downloads/cellbrowserData/"

# a special value that is used for both x and y to indicate that the cell should not be shown
# must match the same value in maxPlot.js
HIDDENCOORD = 12345

# special value representing NaN in floating point arrays
# must match the same value in cellBrowser.js
FLOATNAN = float('-inf') # NaN and sorting does not work. we want NaN always to be first, so encode as -inf
# special value representing NaN in integer arrays, again, we want this to be first after sorting
# must match the same value in cellBrowser.js
INTNAN = -2**16

# how many md5 characters to keep in version identifiers. We load all files using their md5 to address
# internet browser caching
MD5LEN = 10

coordLabels = {
    #  generic igraph neighbor-based layouts
    "fa": "ForceAtlas2",
    "fr": "Fruchterman Reingold",
    "grid_fr": "Grid Fruchterman Reingold", # deactivated for now due to https://github.com/igraph/python-igraph/issues/152
    "kk": "Kamadi Kawai",
    "lgl": "Large Graph Layout", # looks useless
    "drl": "DrL Distributed Recursive Layout",
    "rt": "Reingold Tilford tree", # doesn't look useful

    # special scanpy layouts
    "tsne" : "t-SNE",
    "umap" : "UMAP",
    "pagaFa" : "PAGA/ForceAtlas2",
    "pagaFr" : "PAGA/Fruchterman-Reingold",
    "phate" : "PHATE"
}

recommendedLayouts = ["fa", "fr", "kk", "drl", "tsne", "umap", "pagaFa", "phate"]

metaLabels = {
    "louvain" : "Louvain Cluster",
    "percent_mito" : "Percent Mitochond.",
    "n_genes" : "Expressed Genes",
    "n_counts" : "UMI Count"
}

# ==== functions =====

debugDone = False
debugMode = False

def setDebug(doDebug):
    " activate debugging if needed "
    global debugDone
    global debugMode
    if debugDone:
        return

    if doDebug:
        logging.basicConfig(level=logging.DEBUG)
        logging.getLogger().setLevel(logging.DEBUG)
        debugMode = True
    else:
        logging.basicConfig(level=logging.INFO)
        logging.getLogger().setLevel(logging.INFO)
    debugDone = True

def isDebugMode():
    return debugMode

def makeDir(outDir):
    if not isdir(outDir):
        logging.info("Creating %s" % outDir)
        os.makedirs(outDir)

def errAbort(msg):
        logging.error(msg)
        sys.exit(1)

def iterItems(d):
    " wrapper for iteritems for all python versions "
    if isPy3:
        return d.items()
    else:
        return d.iteritems()

def nextEl(d):
    " wrapper for next() for all python versions "
    if isPy3:
        return next(d)
    else:
        return d.next()

def splitOnce(s, sep, splitCount=1):
    " split s only once on sep, for all python versions "
    if isPy3:
        tup = s.split(sep, maxsplit=splitCount)
    else:
        tup = string.split(s, sep, maxsplit=splitCount)
    return tup

def which(prog):
    " return path of program in PATH "
    try:
        import distutils.spawn
        return distutils.spawn.find_executable(prog)
    except:
        return shutil.which(prog)
    assert(False)

def findCbData():
    """ return the name of the dataDir directory:
    This is the directory in the environment variable CBDATA. If this variable is not set,
    we use <directory-of-this-library/cellbrowserData> or /usr/local/share/cellbrowser, if any of these exist.
    If all of that fails, we use ~/cellbrowserData.
    """
    global dataDir
    if dataDir is not None:
        return dataDir
    
    envCbData = os.environ.get("CBDATA")
    if envCbData is not None:
        logging.debug("CBDATA variable found, points to %s" % envCbData)
        dataDir = envCbData
    else:
        baseDir = dirname(__file__) # = directory of this library
        libRelDir = abspath(join(baseDir, "../../../cellbrowserData"))
        for altDir in [libRelDir, "/usr/local/share/cellbrowser", "/usr/local/cellbrowser"]:
            if isdir(altDir):
                logging.debug("Found data directory in %s" % altDir)
                dataDir = altDir
                break

        if dataDir is None:
            dataDir = expanduser("~/cellbrowserData")
            logging.debug("Using %s as data directory" % dataDir)
            if not isdir(dataDir):
                makeDir(dataDir)

    return dataDir

def downloadStaticFile(remotePath, localPath):
    " download a file from CBHOMEURL/<remotePath> to localPath "
    localDir = dirname(localPath)
    makeDir(localDir)

    remoteUrl = urljoin(CBHOMEURL, remotePath)
    logging.info("Downloading %s to %s..." % (remoteUrl, localPath))
    data = urlopen(remoteUrl).read()

    localTmp = localPath+".download"
    ofh = open(localTmp, "wb")
    ofh.write(data)
    ofh.close()
    os.rename(localTmp, localPath)

def getStaticFile(relPath):
    """ get the full path to a static file in the dataDir directory (~/cbData or $CBDATA, by default, see above).
    If the file is not present, it will be downloaded from https://cells.ucsc.edu/downloads/cbData/<pathParts>
    and copied onto the local disk under dataDir
    """
    dataDir = findCbData()
    absPath = join(dataDir, relPath)
    if not isfile(absPath):
        logging.info("%s not found" % absPath)
        downloadStaticFile(relPath, absPath)

    return absPath

def openStaticFile(relPath, mode="r"):
    " download static file and return an open file handle to it "
    absPath = getStaticFile(relPath)
    fh = openFile(absPath, mode)
    return fh

def staticFileNextRow(relPath):
    " yield rows from a static file, possibly downloading it first "
    fh = openStaticFile(relPath)
    for row in lineFileNextRow(fh):
        yield row

def copyPkgFile(relPath, outDir=None, values=None):
    """ copy file from directory under the current package directory to outDir or current directory
    Don't overwrite if the file is already there.
    """
    if outDir is None:
        outDir = os.getcwd()
    baseDir = dirname(__file__) # = directory of this script
    srcPath = join(baseDir, relPath)
    destPath = join(outDir, basename(relPath))
    if isfile(destPath):
        logging.info("%s already exists, not overwriting" % destPath)
    else:
        if values is None:
            shutil.copyfile(srcPath, destPath)
            logging.info("Wrote %s" % destPath)
            # egg-support commented out for now
            #s = pkg_resources.resource_string(__name__, srcPath)
            #ofh = open(destPath, "wb")
            #ofh.write(s)
            #ofh.close()
        else:
            logging.debug("Using %s as template for %s, vars %s" % (srcPath, destPath, values))
            data = open(srcPath).read()
            dataNew = data % values
            ofh = open(destPath, "w")
            ofh.write(dataNew)
            ofh.close()
            logging.info("Wrote %s" % destPath)

def execfile(filepath, globals=None, locals=None):
    " version of execfile for both py2 and py3 "
    logging.debug("Executing %s" % filepath)
    if globals is None:
        globals = {}
    globals.update({
        "__file__": filepath,
        "__name__": "__main__",
    })
    with open(filepath, 'rb') as file:
        exec(compile(file.read(), filepath, 'exec'), globals, locals)

def loadConfig(fname, requireTags=['name', 'coords', 'meta', 'exprMatrix']):
    """ parse python in fname and return variables as dictionary.
    add the directory of fname to the dict as 'inDir'.
    """
    logging.debug("Loading settings from %s" % fname)
    g = {}
    l = OrderedDict()
    execfile(fname, g, l)

    conf = l

    for rt in requireTags:
        if not rt in conf:
            errAbort("The input configuration has to define the %s statement" % rt)
        if rt=="tags":
            if type(conf["tags"])!=type([]):
                errAbort("'tags' in input config file must be a list")

    conf["inDir"] = dirname(fname)

    return conf

def maybeLoadConfig(confFname):
    if isfile(confFname):
        conf = loadConfig(confFname, requireTags=[])
    else:
        logging.debug("Could not find %s, not loading config file" % confFname)
        conf = OrderedDict()
    return conf

cbConf = None
def getConfig(tag, defValue=None):
    " get a global cellbrowser config value from ~/.cellbrowser.conf "
    global cbConf
    if cbConf is None:
        confPath = expanduser("~/.cellbrowser.conf")
        cbConf = maybeLoadConfig(confPath)

    ret = cbConf.get(tag, defValue)
    return ret

def main_parseArgs():
    " arg parser for __main__, only used internally "
    parser = optparse.OptionParser("""usage: %prog serve outDir port
            serve outDir via http on port
    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages")

    (options, args) = parser.parse_args()

    if len(args)==0:
        parser.print_help()
        exit(1)

    setDebug(options.debug)

    return args, options

# ---- GLOBAL ----
defOutDir = getConfig("htmlDir")
if defOutDir is None:
    defOutDir = os.environ.get("CBOUT")

if defOutDir is not None:
    defOutDir = expanduser(defOutDir)

# ---- GLOBAL END ----

def cbBuild_parseArgs(showHelp=False):
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] -i cellbrowser.conf -o outputDir - add a dataset to the single cell viewer directory

    If you have previously built into the same output directory with the same dataset and the
    expression matrix has not changed its filesize, this will be detected and the expression
    matrix will not be copied again. This means that an update of a few meta data attributes
    is quite quick.

    """)

    parser.add_option("", "--init", dest="init", action="store_true",
        help="copy sample cellbrowser.conf and desc.conf to current directory")

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages")

    parser.add_option("-i", "--inConf", dest="inConf", action="append",
        help="a cellbrowser.conf file that specifies labels and all input files, default %default, can be specified multiple times")

    parser.add_option("-o", "--outDir", dest="outDir", action="store", help="output directory, default can be set through the env. variable CBOUT or ~/.cellbrowser.conf, current value: %default", default=defOutDir)

    parser.add_option("-p", "--port", dest="port", action="store",
        help="if build is successful, start an http server on this port and serve the result via http://localhost:port", type="int")

    (options, args) = parser.parse_args()

    if showHelp:
        parser.print_help()
        exit(1)

    setDebug(options.debug)

    return args, options

def cbUpgrade_parseArgs():
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("usage: %prog [options] outDir - copy all relevant js/css files into outDir, look for datasets in it and create index.html")

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages")
    parser.add_option("-o", "--outDir", dest="outDir", action="store",
        help="output directory, default can be set through the env. variable CBOUT, current value: %default",
        default=defOutDir)
    parser.add_option("", "--dev", dest="devMode", action="store_true",
        help="only for developers: do not add version to js/css links")

    (options, args) = parser.parse_args()

    if options.outDir==None:
        parser.print_help()
        exit(1)

    setDebug(options.debug)
    return args, options

def cbScanpy_parseArgs():
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] -e matrixFile -o outDir -n datasetName - run scanpy and output .tsv files
    """)

    parser.add_option("-e", "--exprMatrix", dest="exprMatrix", action="store",
            help="gene-cell expression matrix file, possible formats: .csv, .h5, .mtx, .txt, .tab, .loom, .h5ad. Existing meta data from .loom and .h5ad will be kept and exported.")

    parser.add_option("-m", "--meta", dest="meta", action="store",
            help="name of cell meta data table. A table like tsv or csv format, first row has cellId and the cellId must match a sample name in the expression matrix. Optional but required if you use --inCluster. 'inMeta' in scanpy.conf")

    parser.add_option("-o", "--outDir", dest="outDir", action="store",
            help="output directory")

    parser.add_option("-n", "--name", dest="name", action="store",
            help="internal name of dataset in cell browser. No spaces or special characters.")

    parser.add_option("", "--init", dest="init", action="store_true",
            help="copy sample scanpy.conf to current directory")

    parser.add_option("-s", "--samplesOnRows", dest="samplesOnRows", action="store_true",
            help="when reading the expression matrix from a text file, assume that samples are on lines (default behavior is one-gene-per-line, one-sample-per-column). Also in scanpy.conf.")

    parser.add_option("-c", "--confFname", dest="confFname", action="store", default="scanpy.conf",
            help="config file from which settings are read, default is %default")

    #parser.add_option("", "--inMeta", dest="inMeta", action="store",
            #help="Existing meta data to read into the scanpy anndata.obs system. A .tsv or .csv file. Use in combination with --inCluster to get marker genes for existing cell type clusters. Can also be set in scanpy.conf.")

    parser.add_option("", "--inCluster", dest="inCluster", action="store",
            help="Do not louvain-cluster, but use this meta field (=obs) when calculating marker genes. The default is to use the louvain clustering results. Also in scanpy.conf.")

    parser.add_option("", "--copyMatrix", dest="copyMatrix", action="store_true",
            help="Instead of reading the input matrix into scanpy and then writing it back out, just copy the input matrix. Only works if the input matrix is gzipped and in the right format and a tsv or csv file, not mtx or h5-based files.")

    #parser.add_option("", "--skipMatrix", dest="skipMatrix", action="store_true",
            #help="Do not write the matrix. You will have to copy it manually.")

    parser.add_option("-g", "--genome", dest="genome", action="store",
            help="when reading 10X HDF5 files, the genome to read. Default is %default. Use h5ls <h5file> to show possible genomes", default="GRCh38")

    #parser.add_option("-m", "--metaFields", dest="metaFields", action="store",
            #help="optional list of comma-separated meta-fields to export from the annData object. All fields are exported by default.")

    parser.add_option("", "--test",
        dest="test",
        action="store_true", help="run doctests")
    parser.add_option("-d", "--debug", dest="debug", action="store_true", help="open an iPython shell when an exception occurs. also output debug messages")

    (options, args) = parser.parse_args()

    if options.test:
        import doctest
        doctest.testmod()
        sys.exit(0)

    if (options.exprMatrix is None or options.outDir is None or options.name is None) and not options.init:
        print("Please specify at least the expression matrix (-e), the output directory (-o) and a name (-n)")
        parser.print_help()
        exit(1)

    setDebug(options.debug)
    return args, options

kwSet = set(keyword.kwlist)

def lineFileNextRow(inFile, utfHacks=False, headerIsRow=False):
    """
    parses tab-sep file with headers in first line
    yields collection.namedtuples
    strips "#"-prefix from header line
    utfHacks forces all chars to latin1 and removes anything that doesn't fit into latin1
    """

    if isinstance(inFile, str):
        # input file is a string = file name
        fh = openFile(inFile, mode="rtU")
        sep = sepForFile(inFile)
    else:
        fh = inFile
        sep = "\t"

    line1 = fh.readline()
    line1 = line1.strip("\n\r").lstrip("#")
    if utfHacks:
        line1 = line1.decode("latin1")
        # skip special chars in meta data and keep only ASCII
        line1 = unicodedata.normalize('NFKD', line1).encode('ascii','ignore')
    headers = line1.split(sep)

    if len(headers)>=255:
        errAbort("Cannot read more than 255 columns. Are you sure that this file is in the correct format?"
                " It may have the wrong line endings and may require treatment with dos2unix or mac2unix. "
                " Or it may be the wrong file type for this input, e.g. an expression matrix instead of a "
                " coordinate file.")

    # python does unfortunately not accept reserved names as named tuple names
    # We append a useless string to avoid errors
    if len(kwSet.intersection(headers))!=0:
        newHeaders = []
        for h in headers:
            if h in kwSet:
                h = h+"_Value"
            newHeaders.append(h)
        headers = newHeaders

    #headers = [x if x!="" else "noName" for x in headers]
    if headers[0]=="": # R does not name the first column by default
        headers[0]="rowName"

    if "" in headers:
        logging.error("Found empty cells in header line of %s" % inFile)
        logging.error("This often happens with Excel files. Make sure that the conversion from Excel was done correctly. Use cut -f-lastColumn to remove empty trailing columns.")
        assert(False)

    # Python does not accept headers that start with a digit
    filtHeads = []
    for h in headers:
        if h[0].isdigit():
            filtHeads.append("x"+h)
        else:
            filtHeads.append(h)
    headers = filtHeads

    origHeaders = headers
    headers = [sanitizeHeader(h) for h in headers]

    if headerIsRow:
        yield origHeaders

    Record = namedtuple('tsvRec', headers)

    for line in fh:
        if line.startswith("#"):
            continue
        if utfHacks:
            line = line.decode("latin1")
            # skip special chars in meta data and keep only ASCII
            line = unicodedata.normalize('NFKD', line).encode('ascii','ignore')
        #line = line.decode("latin1")
        # skip special chars in meta data and keep only ASCII
        #line = unicodedata.normalize('NFKD', line).encode('ascii','ignore')
        line = line.rstrip("\r\n")
        fields = line.split(sep)

        if sep==",":
            fields = [x.lstrip('"').rstrip('"') for x in fields]

        try:
            rec = Record(*fields)
        except Exception as msg:
            logging.error("Exception occurred while parsing line, %s" % msg)
            logging.error("Filename %s" % fh.name)
            logging.error("Line was: %s" % line)
            logging.error("Does number of fields match headers?")
            logging.error("Headers are: %s" % headers)
            raise Exception("header count: %d != field count: %d wrong field count in line %s" % (len(headers), len(fields), line))
        yield rec

def parseOneColumn(fname, colName):
    " return a single column from a tsv as a list "
    ifh = open(fname)
    sep = sepForFile(fname)
    headers = ifh.readline().rstrip("\r\n").split(sep)
    colIdx = headers.index(colName)
    vals = []
    for line in ifh:
        row = line.rstrip("\r\n").split(sep)
        vals.append(row[colIdx])
    return vals

def parseIntoColumns(fname):
    " parse tab sep file vertically, return as a list of (headerName, list of values) "
    ifh = open(fname)
    sep = "\t"
    headers = ifh.readline().rstrip("\r\n").split(sep)
    if headers[0]=="":
        headers[0]="cell_id" # some tolerance, for R
    for i, h in enumerate(headers):
        if h=="":
            errAbort("Header '%s' of column %d is empty. Please fix the meta data file and give every column a name" %
                    (h, i))
    colsToGet = range(len(headers))

    columns = []
    for h in headers:
        columns.append([])

    for line in ifh:
        row = line.rstrip("\r\n").split(sep)
        for colIdx in colsToGet:
            columns[colIdx].append(row[colIdx])
    return zip(headers, columns)

def openFile(fname, mode="rt"):
    if fname.endswith(".gz"):
        mode = mode.replace("U", "")
        if isPy3:
            fh = gzip.open(fname, mode, encoding="latin1")
        else:
            fh = gzip.open(fname, mode)
    else:
        if isPy3:
            fh = io.open(fname, mode)
        else:
            fh = open(fname, mode)
    return fh

def parseDict(fname):
    """ parse text file in format key<tab>value and return as dict key->val """
    d = {}

    fh = openFile(fname)

    sep = "\t"
    if fname.endswith(".csv"):
        sep = ","

    for line in fh:
        key, val = line.rstrip("\r\n").split(sep)
        d[key] = val
    return d

def readGeneToSym(fname):
    " given a file with geneId,symbol return a dict geneId -> symbol. Strips anything after . in the geneId "
    if fname.lower()=="none":
        return None

    logging.info("Reading gene,symbol mapping from %s" % fname)

    # Jim's files and CellRanger files have no headers, they are just key-value
    line1 = openFile(fname).readline().rstrip("\r\n")
    fieldCount = len(line1.split('\t'))
    if "geneId" not in line1:
        logging.debug("geneID,symbol file does not start with 'geneId', parsing as key,value")
        d = parseDict(fname)
    # my new files are smaller and have headers
    elif line1=="#geneId\tsymbol" or fieldCount==2:
        d = {}
        for row in lineFileNextRow(fname):
            if row.symbol=="":
                continue
            geneId = row.geneId
            if geneId.startswith("EN") and "." in geneId:
                geneId = geneId.split(".")[0]
            d[geneId]=row.symbol
    else:
        assert(False) # symbols file does not have a header like #geneId<tab>symbol
    logging.debug("Found symbols for %d genes" % len(d))
    return d

def getDecilesList_np(values):
    deciles = np.percentile( values, [0,10,20,30,40,50,60,70,80,90,100] )
    return deciles

def bytesAndFmt(x):
    """ how many bytes do we need to store x values and what is the sprintf
    format string for it?
    """

    if x > 65535:
        assert(False) # field with more than 65k elements or high numbers? Weird meta data.

    if x > 255:
        return "Uint16", "<H" # see javascript typed array names, https://developer.mozilla.org/en-US/docs/Web/JavaScript/Typed_arrays
    else:
        return "Uint8", "<B"

#def getDecilesWithZeros(numVals):
#    """ return a pair of the deciles and their counts.
#    Counts is 11 elements long, the first element holds the number of zeros,
#    which are treated separately
#
#    >>> l = [0,0,0,0,0,0,0,0,0,0,0,0,1,2,3,4,5,6,7,8,9,10]
#    >>> getDecilesWithZeros(l)
#     ([1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [12, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
#    """
#    nonZeros  = [x for x in numVals if x!=0.0]
#
#    zeroCount = len(numVals) - len(nonZeros)
#    deciles   = getDecilesList_np(nonZeros)
#
#    decArr = np.searchsorted(deciles, numVals)
#    decCounts(deciles, nonZeros)
#
#    decCounts.insert(0, zeroCount)
#    return deciles, decCounts, newVals

def findBins(numVals, breakVals, hasNan):
    """
    find the right bin index defined by breakVals for every value in numVals.
    Special handling for the last value. The comparison uses "<=". The first
    break is assumed to be the minimum of numVals and is therefore ignored.
    Also returns an array with the count for every bin. hasNan triggers a special
    mode where the first bin is reserved for NaNs (encoded as -inf)
    >>> findBins([1,1,1,2,2,2,3,3,4,4,5,5,6,6], [1, 2,3,5,6])
    ([0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 3, 3], [6, 2, 4, 2])
    """
    bArr = []
    binCounts = [0]*(len(breakVals)-1)

    # NaNs mean that the non-NaN bins are all +1
    if hasNan:
        breaks = breakVals[2:]
    else:
        breaks = breakVals[1:]

    if hasNan:
        for x in numVals:
            # we use -inf for the NaN value everywhere, sorting is undefined in lists that contain NaN
            if math.isinf(x):
                binIdx = 0
            else:
                binIdx = bisect.bisect_left(breaks, x)+1
            bArr.append(binIdx)
            binCounts[binIdx]+=1
    else:
        for x in numVals:
            binIdx = bisect.bisect_left(breaks, x) # comparison operator is <=
            bArr.append(binIdx)
            binCounts[binIdx]+=1

    return bArr, binCounts

def countBinsBetweenBreaks(numVals, breakVals):
    """ count how many numVals fall into the bins defined by breakVals.
    Special handling for the last value. Comparison uses "<=". The first
    break is assumed to be the minimum of numVals.
    Also returns an array with the bin for every element in numVals
    >>> countBinsBetweenBreaks([1,1,1,2,2,2,3,3,4,4,5,5,6,6], [1,2,3,5,6])
    ([6, 2, 4, 2], [0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 3, 3])
    """

    binCounts = []
    binCount = 0
    i = 1
    dArr = []
    for x in numVals:
        if x <= breakVals[i]:
            binCount+=1
        else:
            binCounts.append(binCount)
            binCount = 1
            i += 1
        dArr.append(i-1)

    binCounts.append(binCount)

    assert(len(dArr)==len(numVals))
    assert(len(binCounts)==len(breakVals)-1)
    return binCounts, dArr

def discretizeArray(numVals, fieldMeta):
    """
    discretize numeric values based on quantiles.
    """
    maxBinCount = 10
    counts = Counter(numVals).most_common()
    counts.sort() # sort by value, not count

    if len(counts) < maxBinCount:
        # if we have just a few values, do not do any binning
        binCounts = [y for x,y in counts]
        values = [x for x,y in counts]

        valToBin = {}
        for i, x in enumerate(values):
            valToBin[x] = i

        dArr = [valToBin[x] for x in numVals]

        fieldMeta["binMethod"] = "raw"
        fieldMeta["values"] = values
        fieldMeta["binCounts"] = binCounts
        return dArr, fieldMeta

    # ten breaks
    breakPercs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    countLen = len(counts)
    breakIndices = [int(round(bp*countLen)) for bp in breakPercs]
    # as with all histograms, the last break is always a special case (0-based array)
    breakIndices.append(countLen-1)
    breakVals = [counts[idx][0] for idx in breakIndices]

    # NaNs are encoded as -inf so they always are the first break
    # The first non-NaN value is at index 1
    # If we have NaNs, we need one more bin, with the first non-Nan value
    hasNan = False
    if math.isinf(breakVals[0]):
        hasNan = True
        breakVals.insert(1, counts[1][0])
        logging.info("Field has NaN/Unknown values")
    logging.debug("Breaks are: %s" % breakVals)

    dArr, binCounts = findBins(numVals, breakVals, hasNan)
    logging.info("Number of values per decile-bin: %s" % binCounts)

    # we should have 11 breaks/10 bins, or 12 breaks/11bins if we have NaN elements
    assert((not hasNan and len(breakVals)==11) or (hasNan and len(breakVals)==12))
    assert((not hasNan and len(binCounts)==10) or (hasNan and len(binCounts)==11))
    assert((len(binCounts)+1 == len(breakVals)))

    fieldMeta["binMethod"] = "quantiles"
    fieldMeta["binCounts"] = binCounts
    if math.isinf(breakVals[0]): # -infinity is not valid in JSON
        breakVals[0] = "Unknown"
    fieldMeta["breaks"] = breakVals

    return dArr, fieldMeta

def discretizeNumField(numVals, fieldMeta, numType):
    " given a list of numbers, add attributes to fieldMeta that describe the binning scheme "
    digArr, fieldMeta = discretizeArray(numVals, fieldMeta)

    #deciles, binCounts, newVals = getDecilesWithZeros(numVals)

    fieldMeta["arrType"] = "uint8"
    fieldMeta["_fmt"] = "<B"
    return digArr, fieldMeta

def typeForStrings(strings):
    """ given a list of strings, determine if they're all ints or floats or strings
    """
    floatCount = 0
    intCount = 0
    for val in strings:
        try:
            newVal = int(val)
            intCount += 1
        except:
            try:
                newVal = float(val)
                floatCount += 1
            except:
                return "string"

    if floatCount!=0:
        return "float"
    return "int"

emptyVals = ["", "null", "none", "None", "unknown", "nd", "n.d.", "Unknown", "NaN", "NA", "undefined", "Na"]

def likeEmptyString(val):
    " returns true if string is a well-known synonym of 'unknown' or 'NaN'. ported from cellbrowser.js "
    return val.strip() in emptyVals

def floatToIntList(vals):
    " convert a list of floats to a integers, take care of -inf values "
    newVals = []
    for x in vals:
        if x==FLOATNAN:
            newVals.append(INTNAN)
        else:
            newVals.append(int(x))
    return newVals


def guessFieldMeta(valList, fieldMeta, colors, forceEnum):
    """ given a list of strings, determine if they're all int, float or
    strings. Return fieldMeta, as dict, and a new valList, with the correct python type
    - 'type' can be: 'int', 'float', 'enum' or 'uniqueString'
    - if int or float: 'deciles' is a list of the deciles
    - if uniqueString: 'maxLen' is the length of the longest string
    - if enum: 'values' is a list of all possible values
    - if colors is not None: 'colors' is a list of the default colors
    """
    unknownCount = 0
    intCount = 0
    floatCount = 0
    valCounts = defaultdict(int)
    newVals = []
    for val in valList:
        fieldType = "string"
        newVal = val

        if likeEmptyString(val):
            unknownCount += 1
            newVal = FLOATNAN
        else:
            try:
                newVal = int(val)
                intCount += 1
                floatCount += 1
            except:
                try:
                    newVal = float(val)
                    floatCount += 1
                except:
                    pass

        newVals.append(newVal)
        valCounts[val] += 1

    valToInt = None
    assert(len(newVals)==len(valList))

    if len(valCounts)==1:
        logging.warn("Field %s contains only a single value" % fieldMeta["name"])


    if intCount+unknownCount==len(valList) and not forceEnum:
        # field is an integer
        #newVals = floatToIntList(newVals)
        #newVals, fieldMeta = discretizeNumField(numVals, fieldMeta, "int")
        #assert(min(newVals) > -2**32) # please contact us if you need very big numbers
        #assert(max(newVals) < 2**32)  # please contact us if you need very big numbers
        #minVal = min(newVals)
        #maxVal = max(newVals)
        #if minVal > -2**16 and maxVal < 2**16:
            #fieldMeta["arrType"] = "int32"
            #fieldMeta["_fmt"] = "<l" # signed long, 4 bytes
        #elif minVal >= 1 and maxVal < 2**32:
            #fieldMeta["arrType"] = "uint32" # unsigned long, 4 bytes
            #fieldMeta["_fmt"] = "<L"

        # JS supports only 32bit signed ints so we store integers as floats
        newVals = [float(x) for x in newVals]
        fieldMeta["arrType"] = "float32"
        fieldMeta["_fmt"] = "<f"
        fieldMeta["type"] = "int"

    elif floatCount+unknownCount==len(valList) and not forceEnum:
        # field is a floating point number: convert to decile index
        newVals = [float(x) for x in newVals]
        #newVals, fieldMeta = discretizeNumField(numVals, fieldMeta, "float")
        fieldMeta["arrType"] = "float32"
        fieldMeta["_fmt"] = "<f"
        fieldMeta["type"] = "float"

    elif len(valCounts)==len(valList) and not forceEnum:
        # field is a unique string
        fieldMeta["type"] = "uniqueString"
        maxLen = max([len(x) for x in valList])
        fieldMeta["maxSize"] = maxLen
        fieldMeta["_fmt"] = "%ds" % (maxLen+1)
        newVals = valList

    else:
        # field is an enum - convert to enum index
        fieldMeta["type"] = "enum"
        valArr = list(valCounts.keys())

        valCounts = list(sorted(valCounts.items(), key=operator.itemgetter(1), reverse=True)) # = (label, count)
        #valCounts = valCounts.items()
        if colors!=None:
            colArr = []
            foundColors = 0
            notFound = set()
            for val, _ in valCounts:
                if val in colors:
                    colArr.append(colors[val])
                    foundColors +=1
                else:
                    notFound.add(val)
                    colArr.append(None) # maybe I should fail hard here?

            if foundColors > 0:
                fieldMeta["colors"] = colArr
                if len(notFound)!=0:
                    logging.warn("No default color found for field values %s. Set these to defaults." % notFound)

        fieldMeta["valCounts"] = valCounts
        fieldMeta["arrType"], fieldMeta["_fmt"] = bytesAndFmt(len(valArr))
        valToInt = dict([(y[0],x) for (x,y) in enumerate(valCounts)]) # dict with value -> index in valCounts
        newVals = [valToInt[x] for x in valList] #

    fieldMeta["diffValCount"] = len(valCounts)

    return fieldMeta, newVals

def writeNum(col, packFmt, ofh):
    " write a list of numbers to a binary file "

def cleanString(s):
    " returns only alphanum characters in string s "
    newS = []
    for c in s:
        if c.isalnum():
            newS.append(c)
    return "".join(newS)

def moveOrGzip(inFname, outFname):
    " if outFname has .gz, runGzip, otherwise just move file over "
    if outFname.endswith(".gz"):
        runGzip(inFname, outFname)
    else:
        logging.debug("Renaming %s to %s" % (inFname, outFname))
        os.rename(inFname, outFname)

def runGzip(fname, finalFname=None):
    " compress fname and move to finalFname when done "
    logging.debug("Compressing %s" % fname)
    cmd = "gzip -f %s" % fname
    runCommand(cmd)
    gzipFname = fname+".gz"

    if finalFname==None:
        return gzipFname

    if isfile(finalFname):
        os.remove(finalFname)
    logging.debug("Renaming %s to %s" % (gzipFname, finalFname))
    os.rename(gzipFname, finalFname)
    if isfile(fname):
        os.remove(fname)
    return finalFname

def addLongLabels(acronyms, fieldMeta):
    """ add a 'longLabel' meta info entry, a list of longer strings for enum fields, if any shortLabels have a longLabel 
    in the acronyms dict of shortLabel = longLabel """
    if acronyms is None:
        return fieldMeta

    foundCount = 0
    for shortLabel, count in fieldMeta["valCounts"]:
        if shortLabel in acronyms:
            foundCount+=1

    if foundCount > 0:
        longLabels = []
        for shortLabel, count in fieldMeta["valCounts"]:
            longLabel = acronyms.get(shortLabel)
            if longLabel is None:
                logging.warning("Field %s: value %s has no long label through the acronyms" %
                        (fieldMeta["label"], shortLabel))
                longLabel = shortLabel
            longLabels.append(longLabel)
        fieldMeta["longLabels"] = longLabels

    return fieldMeta

def metaToBin(inConf, outConf, fname, colorFname, outDir, enumFields):
    """ convert meta table to binary files. outputs fields.json and one binary file per field.
    adds names of metadata fields to outConf and returns outConf
    """
    logging.info("Converting to numbers and compressing meta data fields")
    makeDir(outDir)

    colData = parseIntoColumns(fname)

    colors = parseColors(colorFname)
    acronyms = readAcronyms(inConf, outConf)

    # the user inputs the enum fields cellbrowser.conf as their real names, but internally, unfortunately
    # we have to strip special chars so fix the user's field names to our format
    sanEnumFields = []
    if enumFields is not None:
        sanEnumFields = [sanitizeHeader(n) for n in enumFields]

    fieldInfo = []
    validFieldNames = set()
    for colIdx, (fieldName, col) in enumerate(colData):
        logging.debug("Meta data field index %d: '%s'" % (colIdx, fieldName))
        validFieldNames.add(fieldName)

        forceEnum = (fieldName in sanEnumFields)
        # very dumb heuristic to recognize fields that should not be treated as numbers but as enums
        # res.0.6 is the default field name for Seurat clustering. Field header sanitizing changes it to
        # res_0_6 which is not optimal, but namedtuple doesn't allow dots in names
        if "luster" in fieldName or "ouvain" in fieldName or (fieldName.startswith("res_") and "_" in fieldName):
            forceEnum=True

        cleanFieldName = cleanString(fieldName)
        binName = join(outDir, cleanFieldName+".bin")

        fieldMeta = OrderedDict()
        fieldMeta["name"] = cleanFieldName
        fieldMeta["label"] = fieldName

        fieldMeta, binVals = guessFieldMeta(col, fieldMeta, colors, forceEnum)

        fieldType = fieldMeta["type"]

        if fieldType=="enum":
            fieldMeta = addLongLabels(acronyms, fieldMeta)

        if "metaOpt" in inConf and fieldName in inConf["metaOpt"]:
            fieldMeta["opt"] = inConf["metaOpt"][fieldName]

        packFmt = fieldMeta["_fmt"]

        # write the binary file
        binFh = open(binName, "wb")
        if fieldMeta["type"]!="uniqueString":
            for x in binVals:
                binFh.write(struct.pack(packFmt, x))
        else:
            for x in col:
                if isPy3:
                    binFh.write(bytes("%s\n" % x, encoding="ascii"))
                else:
                    binFh.write("%s\n" % x)
        binFh.close()

        runGzip(binName)
        zippedName = binName+".gz"

        fieldMeta["md5"] = md5WithPython(zippedName)[:MD5LEN]

        del fieldMeta["_fmt"]
        fieldInfo.append(fieldMeta)
        if "type" in fieldMeta:
            logging.info(("Field %(name)s: type %(type)s, %(diffValCount)d different values" % fieldMeta))
        else:
            logging.info(("Field %(name)s: type %(type)s, %(diffValCount)d different values, max size %(maxSize)d " % fieldMeta))

    return fieldInfo, validFieldNames

def iterLineOffsets(ifh):
    """ parse a text file and yield tuples of (line, startOffset, endOffset).
    endOffset does not include the newline, but the newline is not stripped from line.
    """
    line = True
    start = 0
    while line!='':
       line = ifh.readline()
       end = ifh.tell()-1
       if line!="":
           yield line, start, end
       start = ifh.tell()

class MatrixTsvReader:
    " open a .tsv or .csv file and yield rows via iterRows. gz and csv OK."

    def __init__(self, geneToSym=None):
        " can automatically translate to symbols, if dict geneId -> sym is provided "
        self.geneToSym = geneToSym

    def open(self, fname, matType=None, usePyGzip=False):
        " open file and guess field sep and format of numbers (float or int) "

        if which("gunzip")==None:
            logging.warn("Gunzip not in PATH, falling back to Python's built-in")
            usePyGzip = True
        # Note: The encoding for popen below is utf8. Is this a good choice?
        # Does performance change if we don't do unicode strings but
        # byte strings instead? Does unicode eat up the performance gain
        # of gunzip -c ?
        logging.debug("Opening %s" % fname)
        self.fname = fname
        if fname.endswith(".gz"):
            # python's gzip is slower, but does not return an error if we read just 
            # a single line
            if usePyGzip:
                self.ifh = openFile(fname)
            else:
                cmd = ['gunzip', '-c', fname]
                proc, stdout = popen(cmd, doWait=False)
                self.ifh = stdout # faster implementation and in addition uses 2 CPUs
        else:
            self.ifh = io.open(fname, "r", encoding="utf8") # utf8 performance? necessary for python3?

        self.sep = sepForFile(fname)
        logging.debug("Field separator is %s" % repr(self.sep))

        headLine = self.ifh.readline()
        headLine = headLine.rstrip("\r\n")
        self.sampleNames = headLine.split(self.sep)[1:]
        self.sampleNames = [x.strip('"') for x in self.sampleNames]
        assert(len(self.sampleNames)!=0)
        logging.debug("Read %d sampleNames, e.g. %s" % (len(self.sampleNames), self.sampleNames[0]))

        self.oldRows = []
        if matType is None:
            self.matType = self.autoDetectMatType(10)
            logging.info("Auto-detect: Numbers in matrix are of type '%s'", self.matType)
        else:
            logging.debug("Pre-determined: Numbers in matrix are '%s'" % matType)
            self.matType = matType

    def close(self):
        self.ifh.close()

    def getMatType(self):
        return self.matType

    def getSampleNames(self):
        return self.sampleNames

    def autoDetectMatType(self, n):
        " check if matrix has 'int' or 'float' data type by looking at the first n genes"
        # auto-detect the type of the matrix: int vs float
        logging.info("Auto-detecting number type of %s" % self.fname)
        geneCount = 0

        self.matType = "float" # iterRows needs this attribute

        matType = "int"
        for geneId, sym, a in self.iterRows():
            self.oldRows.append( (geneId, sym, a) )
            geneCount+=1
            if numpyLoaded:
                a_int = a.astype(int)
                hasOnlyInts = np.array_equal(a, a_int)
                if not hasOnlyInts:
                    matType = "float"
                    break
            else:
                for x in a:
                    frac, whole = math.modf(x)
                    if frac != 0.0:
                        matType = "float"
                        break
                if matType=="float":
                    break
            if geneCount==n:
                break

        if geneCount==0:
            errAbort("empty expression matrix?")
        logging.debug("Matrix type is: %s" % matType)
        return matType

    def iterRows(self):
        " yield (geneId, symbol, array) tuples from gene expression file. "
        if self.matType == "float":
            npType = "float32"
        else:
            npType = "int32"

        # during auto-detection, we've already read a few lines
        logging.debug("spooling back %d saved rows" % len(self.oldRows))
        for (geneId, sym, arr) in self.oldRows:
            # for the integer case though we have to fix up the type now
            if self.matType=="int":
                if numpyLoaded:
                    arr = arr.astype(int)
                else:
                    arr = [int(x) for x in arr]
            yield (geneId, sym, arr)

        skipIds = 0
        doneGenes = set()
        lineNo = 0
        sep = self.sep
        sampleCount = len(self.sampleNames)
        geneToSym = self.geneToSym
        for line in self.ifh:
            self.lineLen = len(line)

            gene, rest = splitOnce(line.rstrip("\r\n"), sep)
            gene = gene.strip('"')

            if numpyLoaded:
                arr = np.fromstring(rest, dtype=npType, sep=sep, count=sampleCount)
            else:
                if self.matType=="int":
                    #try:
                        arr = [int(x) for x in rest.split(sep)]
                    #except ValueError as ex:
                    #logging.warn("Cannot parse expression matrix. This may be due to the numbers incorrectly auto-detected as integers even if they are floating point numbers. Set matrixType='float' in cellbrowser.conf to fix this and re-run cbBuild. The exact error message was: %s" % ex)
                    #arr = map(int, rest.split(sep)) # this doesn't work in python3, requires list(), so slower
                else:
                    arr = [float(x) for x in rest.split(sep)]
                    #arr = map(float, rest.split(sep))

            if "|" in gene:
                gene, symbol = gene.split("|")
            else:
                if geneToSym is None:
                    symbol = gene
                else:
                    symbol = geneToSym.get(gene.split(".")[0])
                    logging.debug("%s -> %s" % (gene, symbol))
                    #if symbol is None:
                        #symbol = geneToSym.get(gene)

                    if symbol is None:
                        skipIds += 1
                        logging.warn("line %d: %s is not a valid gene ID, check geneIdType setting in cellbrowser.conf" % (lineNo, gene))
                        symbol = gene

                    if symbol.isdigit():
                        logging.warn("line %d in gene matrix: gene identifier %s is a number. If this is indeed a gene identifier, you can ignore this warning. Otherwise, your matrix may have no gene ID in the first column and you will have to fix the matrix." % (lineNo, symbol))

            if symbol in doneGenes:
                logging.warn("line %d: Gene %s/%s is duplicated in matrix, using only first occurrence for symbol, kept second occurrence with original geneId" % (lineNo, gene, symbol))
                symbol = gene

            doneGenes.add(gene)

            lineNo += 1

            logging.debug("Yielding gene %s, sym %s, %d fields" % (gene, symbol, len(arr)))
            yield gene, symbol, arr

        if skipIds!=0:
            logging.warn("Kept %d genes as original IDs, due to duplication or unknown ID" % skipIds)

    def iterRowsWithOffsets(self):
        " like iterRows, but also return offset and line length "
        offset = self.ifh.tell()
        for gene, sym, row in self.iterRows():
            yield gene, sym, row, offset, self.lineLen
            offset = self.ifh.tell()

def getDecilesList(values):
    """ given a list of values, return the 10 values that define the 10 ranges for the deciles
    """
    if len(values)==0:
        return None

    valCount = len(values)
    binSize = float(valCount-1) / 10.0; # width of each bin, in number of elements, fractions allowed

    values = list(sorted(values))

    # get deciles from the list of sorted values
    deciles = []
    pos = 0
    for i in range(10): # 10 bins means that we need 10 limits, the last limit is at 90%
        pos = int (binSize * i)
        if pos > valCount: # this should not happen, but maybe it can, due to floating point issues?
            logging.warn("decile exceeds 10, binSize %d, i %d, len(values) %d" % (binSize, i, len(values)))
            pos = len(values)
        deciles.append ( values[pos] )
    return deciles

def findBin(ranges, val):
    """ given an array of values, find the index i where ranges[i] < val <= ranges[i+1]
    ranges have to be sorted.
    This is a dumb brute force implementation - maybe binary search is faster, if ever revisit this again
    Someone said up to 10 binary search is not faster.
    """
    if val==0: # speedup
        return 0
    for i in range(0, len(ranges)):
        if (val < ranges[i]):
            return i
    # if doesn't fit in anywhere, return beyond last possible index
    return i+1

def discretizeArr_uniform(arr, fieldMeta):
    """ given an array of numbers, get min/max, create 10 bins between min and max then
    translate the array to bins and return the list of bins
    """
    arrMin = min(arr)
    arrMax = max(arr)
    stepSize = (arrMax-arrMin)/10.0

    dArr = [0]*len(arr)
    binCounts = [0]*10
    for i, x in enumerate(arr):
        binIdx = int(round((x - arrMin)/stepSize))
        if x == arrMax:
            binIdx = 9
        assert(binIdx <= 9)
        dArr[i] = binIdx
        binCounts[binIdx]+=1

    fieldMeta["binMethod"] = "uniform"
    fieldMeta["minVal"] = arrMin
    fieldMeta["maxVal"] = arrMax
    fieldMeta["stepSize"] = stepSize
    fieldMeta["binCounts"] = binCounts
    return dArr, fieldMeta

def digitize_py(arr, matType):
    """ calculate deciles ignoring 0s from arr, use these deciles to digitize the whole arr,
    return (digArr, zeroCount, bins).
    bins is an array of (min, max, count)
    There are at most 11 bins and bin0 is just for the value zero.
    For bin0, min and max are both 0.0

    matType can be "int" or "float".
    If it is 'int' and arr has only <= 11 values, will not calculate deciles, but rather just
    count the numbers and use them to create bins, one per number.
    #>>> digitize_py([1,1,1,1,1,2,3,4,5,6,4,5,5,5,5], "float")
    """
    if matType=="int":
        valCount = len(set(arr))
        if valCount <= 11: # 10 deciles + 0s
            counts = Counter(arr).most_common()
            counts.sort()

            valToIdx = {}
            for i, (val, count) in enumerate(counts):
                valToIdx[val] = i

            digArr = [valToIdx[x] for x in arr]
            bins = []
            for val, count in counts:
                bins.append( (val, val, count) )
            return digArr, bins

    noZeroArr = [x for x in arr if x!=0]
    zeroCount = len(arr) - len(noZeroArr)
    deciles = getDecilesList(noZeroArr) # there are 10 limits for the 10 deciles, 0% - 90%
    deciles.insert(0, 0) # bin0 is always for the zeros
    # we now have 11 limits
    assert(len(deciles)<=11)

    # digitize and count bins
    digArr = []
    binCounts = len(deciles)*[0]
    for x in arr:
        binIdx = findBin(deciles, x)
        # bin1 is always empty, so move down all other indices
        if binIdx>0:
            binIdx-=1
        digArr.append(binIdx)
        binCounts[binIdx]+=1

    # create the bin info
    bins = []
    if zeroCount!=0:
        bins.append( [float(0), float(0), float(zeroCount)])

    for i in range(1, len(deciles)):
        minVal = deciles[i-1]
        maxVal = deciles[i]
        count = binCounts[i]
        # skip empty bins
        #if count!=0:
        bins.append( [float(minVal), float(maxVal), float(count)] )

    # add the maximum value explicitly, more meaningful
    bins[-1][1] = np.amax(arr)
    return digArr, bins

def digitizeArr(arr, numType):
    if numpyLoaded:
        return digitize_np(arr, numType)
    else:
        return digitize_py(arr, numType)

def binEncode(bins):
    " encode a list of at 11 three-tuples into a string of 33 floats (little endian)"
    # add (0,0,0) elements to bins until it has 11 elements "
    padBins = copy.copy(bins)
    for i in range(len(bins), 11):
        padBins.append( (0.0, 0.0, 0.0) )
    #print len(padBins), padBins, len(padBins)
    assert(len(padBins)==11)

    strList = []
    for xMin, xMax, count in padBins:
        strList.append( struct.pack("<f", xMin) )
        strList.append( struct.pack("<f", xMax) )
        strList.append( struct.pack("<f", count) )
    ret = "".join(strList)
    assert(len(ret)==11*3*4)
    return ret

def digitize_np(arr, matType):
    """ hopefully the same as digitize(), but using numpy
    #>>> digitize_np([1,2,3,4,5,6,4,1,1,1], "int")
    #>>> digitize_np([0,0,0,1,1,1,1,1,2,3,4,5,6,4,5,5,5,5], "float")
    #>>> digitize_np([1,1,1,1,1,2,3,4,5,6,4,5,5,5,5], "float")
    """

    # meta data comes in as a list
    if not type(arr) is np.ndarray:
        arr = np.array(arr)

    if matType=="int":
        # raw counts mode:
        # first try if there are enough unique values in the array
        # if there are <= 10 values, deciles make no sense,
        # so simply enumerate the values and map to bins 0-10
        binCounts = np.bincount(arr)
        nonZeroCounts = binCounts[np.nonzero(binCounts)] # remove the 0s
        if nonZeroCounts.size <= 11:
            logging.debug("we have read counts and <11 values: not using quantiles, just enumerating")
            posWithValue = np.where(binCounts != 0)[0]
            valToBin = {}
            bins = []
            binIdx = 0
            #for val, count in enumerate(binCounts):
                #if count!=0:
            for val in posWithValue:
                count = binCounts[val]
                bins.append( (val, val, count) )
                valToBin[val] = binIdx
                binIdx += 1
            # map values to bin indices, from stackoverflow
            digArr = np.vectorize(valToBin.__getitem__)(arr)
            return digArr, bins

    logging.debug("calculating deciles")
    # calculate the deciles without the zeros, otherwise
    # the 0s completely distort the deciles
    #noZero = np.copy(arr)
    #nonZeroIndices = np.nonzero(arr)
    noZero = arr[np.nonzero(arr)]

    # gene not expressed -> do nothing
    if noZero.size==0:
        logging.debug("expression vector is all zeroes")
        return np.zeros(arr.size, dtype=np.int8), [(0.0, 0.0, arr.size)]

    deciles = np.percentile( noZero, [0,10,20,30,40,50,60,70,80,90] , interpolation="lower")
    # make sure that we always have a bin for the zeros
    deciles = np.insert(deciles, 0, 0)
    logging.debug("deciles are: %s" % str(deciles))

    # now we have 10 limits, defining 11 bins
    # but bin1 will always be empty, as there is nothing between the value 0 and the lowest limit
    digArr = np.searchsorted(deciles, arr, side="right")
    # so we decrease all bin indices that are not 0
    np.putmask(digArr, digArr>0, digArr-1)
    binCounts = np.bincount(digArr)

    bins = []
    zeroCount = binCounts[0]

    # bin0 is a bit special
    if zeroCount!=0:
        bins.append( [float(0), float(0), zeroCount] )

    for i in range(1, len(deciles)):
        binCount = binCounts[i]
        #if binCount==0:
            #continue
        minVal = deciles[i-1]
        maxVal = deciles[i]
        bins.append( [minVal, maxVal, binCount] )

    bins[-1][1] = np.amax(arr)
    #print bins, len(digArr), digArr
    return digArr, bins

def maxVal(a):
    if numpyLoaded:
        return np.amax(a)
    else:
        return max(a)

def discretExprRowEncode(geneDesc, binInfo, digArr):
    " encode geneDesc, deciles and array of decile indixes into a single string that can be read by the .js code "
    # The format of a record is:
    # - 2 bytes: length of descStr, e.g. gene identifier or else
    # - len(descStr) bytes: the descriptive string descStr
    # - 132 bytes: 11 deciles, encoded as 11 * 3 floats (=min, max, count)
    # - array of n bytes, n = number of cells
    decChrList = [chr(x) for x in digArr]
    decStr = "".join(decChrList)
    geneIdLen = struct.pack("<H", len(geneDesc))

    binStr = binEncode(binInfo)
    geneStr = geneIdLen+geneDesc+binStr+decStr

    geneCompr = zlib.compress(geneStr)
    logging.debug("compression factor of %s: %f, before %d, after %d"% (geneDesc, float(len(geneCompr)) / len(geneStr), len(geneStr), len(geneCompr)))

    return geneCompr

def exprEncode(geneDesc, exprArr, matType):
    """ convert an array of numbers of type matType (int or float) to a compressed string of
    floats
    The format of a record is:
    - 2 bytes: length of descStr, e.g. gene identifier or else
    - len(descStr) bytes: the descriptive string descStr
    - array of n 4-byte floats (n = number of cells)
    """
    geneDesc = str(geneDesc) # make sure no unicode
    geneIdLen = struct.pack("<H", len(geneDesc))

    # on cortex-dev, numpy was around 30% faster. Not a huge difference.
    if numpyLoaded:
        exprStr = exprArr.tobytes()
    else:
        if matType=="float":
            arrType = "f"
        elif matType=="int":
            arrType = "I"
        else:
            assert(False) # internal error
        exprStr = array.array(arrType, exprArr).tostring()

    if isPy3:
        geneStr = geneIdLen+bytes(geneDesc, encoding="ascii")+exprStr
    else:
        geneStr = geneIdLen+geneDesc+exprStr

    geneCompr = zlib.compress(geneStr)

    fact = float(len(geneCompr)) / len(geneStr)
    logging.debug("raw - compression factor of %s: %f, before %d, after %d"% (geneDesc, fact, len(geneStr), len(geneCompr)))
    return geneCompr

def matrixToBin(fname, geneToSym, binFname, jsonFname, discretBinFname, discretJsonFname, matType=None):
    """ convert gene expression vectors to vectors of deciles
        and make json gene symbol -> (file offset, line length)
    """
    logging.info("converting %s to %s and writing index to %s, type %s" % (fname, binFname, jsonFname, matType))
    #logging.info("Shall expression values be log-transformed when transforming to deciles? -> %s" % (not skipLog))
    logging.info("Compressing gene expression vectors...")

    tmpFname = binFname + ".tmp"
    ofh = open(tmpFname, "wb")

    discretTmp = discretBinFname + ".tmp"
    discretOfh = open(discretTmp, "w")

    discretIndex = {}
    exprIndex = {}

    skipIds = 0
    highCount = 0

    matReader = MatrixTsvReader(geneToSym)
    matReader.open(fname, matType=matType)

    if matType is None:
        matType = matReader.getMatType()

    sampleNames = matReader.getSampleNames()

    symCounts = defaultdict(int)
    geneCount = 0
    for geneId, sym, exprArr in matReader.iterRows():
        geneCount += 1

        symCounts[sym]+=1
        if symCounts[sym] > 800:
            errAbort("The gene symbol %s appears more than 800 times in the expression matrix. "
                    "Are you sure that the matrix is in the right format? Each gene should be on a row. "
                    "The gene ID must be in the first column and "
                    "can optionally include the gene symbol, e.g. 'ENSG00000142168|SOD1'. " % sym)

        #if maxVal(exprArr) > 200:
            #highCount += 1

        logging.debug("Processing %s, symbol %s" % (geneId, sym))
        exprStr = exprEncode(geneId, exprArr, matType)
        exprIndex[sym] = (ofh.tell(), len(exprStr))
        ofh.write(exprStr)

        if geneCount % 1000 == 0:
            logging.info("Wrote compressed expression values for %d genes" % geneCount)

    discretOfh.close()
    ofh.close()

    #if highCount==0:
        #logging.warn("No single value in the matrix is > 200. It looks like this matrix has been log'ed before. Our recommendation for visual inspection is to not transform matrices, but that is of course up to you.")
        #logging.error("Rerun with --skipLog.")
        #sys.exit(1)

    if len(exprIndex)==0:
        errAbort("No genes from the expression matrix could be mapped to symbols."
            "Are you sure these are Ensembl IDs? Adapt geneIdType in cellbrowser.conf.")

    jsonOfh = open(jsonFname, "w")
    json.dump(exprIndex, jsonOfh)
    jsonOfh.close()

    jsonOfh = open(discretJsonFname, "w")
    json.dump(discretIndex, jsonOfh)
    jsonOfh.close()

    os.rename(tmpFname, binFname)
    os.rename(discretTmp, discretBinFname)

    return matType

def sepForFile(fname):
    if fname.endswith(".csv") or fname.endswith(".csv.gz") or fname.endswith(".csv.Z"):
        sep = ","
    else:
        sep = "\t"
    logging.debug("Separator for %s is %s" %  (fname, repr(sep)))
    return sep

def indexMeta(fname, outFname):
    """ index a tsv by its first field. Writes binary data to outFname.
        binary data is (offset/4 bytes, line length/2 bytes)
    """
    ofh = open(outFname, "wb")
    logging.info("Indexing meta file %s to %s" % (fname, outFname))
    ifh = open(fname)
    sep = sepForFile(fname)
    headerDone = False
    for line, start, end in iterLineOffsets(ifh):
        if not headerDone:
            headerDone = True
            continue

        row = splitOnce(line, sep)
        field1 = row[0]

        lineLen = end - start
        assert(lineLen!=0)
        assert(lineLen<65535) # meta data line cannot be longer than 2 bytes
        ofh.write(struct.pack("<L", start))
        ofh.write(struct.pack("<H", lineLen))
    ofh.close()

# ----------- main --------------

def parseColors(fname):
    " parse color table and return as dict value -> color "
    if fname==None:
        return {}

    if not isfile(fname):
        logging.warn("File %s does not exist" % fname)
        return None

    colDict = parseDict(fname)
    newDict = {}
    for metaVal, color in iterItems(colDict):
        if color.lower()=="color":
            continue

        color = color.strip().strip("#") # hbeale had a file with trailing spaces

        isHex = True
        if len(color)!=6: # colors can be no more than six hex digits
            isHex = False
        else:
            for c in color:
                if (c not in "0123456789ABCDEFabcdef"):
                    isHex = False
                    break

        if not isHex:
            logging.debug("Not a six-digit hex color code. Trying to map '%s' to a hex color" % color)
            import webcolors # error? -> pip install webcolors
            try:
                color = webcolors.name_to_hex(color, spec='css3').lstrip("#")
            except ValueError:
                # R knows more colors, like deeppink4. We simply map to deeppink for now
                # there does not seem to be a good table with R colors in Python yet
                color = "".join([c for c in color if not c.isdigit()])
                color = webcolors.name_to_hex(color, spec='css3').lstrip("#")

        newDict[metaVal] = color
    return newDict

def parseScaleCoordsAsDict(fname, useTwoBytes, flipY):
    """ parse tsv file in format cellId, x, y and return as dict (cellId, x, y)
    Optionally flip the y coordinates to make it more look like plots in R, for people transitioning from R.
    """
    logging.debug("Parsing coordinates from %s. FlipY=%s, useTwoBytes=%s" % (fname, flipY, useTwoBytes))
    coords = []
    maxY = 0
    minX = 2^32
    minY = 2^32
    maxX = -2^32
    maxY = -2^32
    skipCount = 0

    # parse and find the max values
    warn1Done = False
    warn2Done = False
    for row in lineFileNextRow(fname):
        if (len(row)<3):
            if not warn1Done:
                errAbort("file %s needs to have at least three columns" % fname)
                warn1Done = True
        if (len(row)>3): # coord file has to have three rows (cellId, x, y), we just ignore the headers
            if not warn2Done:
                logging.warn("file %s has more than three columns. Everything beyond column 3 will be ignored" % fname)
                warn2Done = True
        cellId = row[0]
        x = float(row[1])
        y = float(row[2])
        coords.append( (cellId, x, y) )

        # special values (12345,12345) mean "unknown cellId"
        if x!=HIDDENCOORD and y!=HIDDENCOORD:
            minX = min(x, minX)
            minY = min(y, minY)
            maxX = max(x, maxX)
            maxY = max(y, maxY)

    if useTwoBytes is None:
        if len(coords)>100000:
            useTwoBytes = True
            logging.info("More than 100k cells, so automatically enabling two-byte encoding for coordinates")
        else:
            useTwoBytes = False

    if useTwoBytes:
        scaleX = 65535/(maxX-minX)
        scaleY = 65535/(maxY-minY)

    newCoords = {}
    for cellId, x, y in coords:
        if useTwoBytes:
            x = int(scaleX * (x - minX))
            y = int(scaleY * (y - minY))
            if flipY:
                y = 65535 - y
        else:
            if flipY:
                y = maxY - y

        newCoords[cellId] = (x, y)

    return newCoords

def sliceRow(row, skipFields):
    " yield all fields, except the ones with an index in skipFields "
    for i, val in enumerate(row):
        if i not in skipFields:
            yield val

def metaReorder(matrixFname, metaFname, fixedMetaFname):
    """ check and reorder the meta data, has to be in the same order as the
    expression matrix, write to fixedMetaFname. Remove single-value fields. """

    logging.info("Checking and reordering meta data to %s" % fixedMetaFname)
    metaSampleNames = readSampleNames(metaFname)

    if matrixFname is not None:
        matrixSampleNames = readHeaders(matrixFname)[1:]
    else:
        matrixSampleNames=metaSampleNames

    # check that there is a 1:1 sampleName relationship
    mat = set(matrixSampleNames)
    meta = set(metaSampleNames)
    if len(meta)!=len(metaSampleNames):
        logging.error("sample names in the meta data differ in length from the sample names in the matrix: %d sample names in the meta data, %d sample names in the matrix" % (len(meta), len(metaSampleNames)))
        sys.exit(1)

    if len(mat.intersection(meta))==0:
        logging.error("Meta data and expression matrix have no single sample name in common. Sure that the expression matrix has one gene per row?")
        sys.exit(1)

    metaNotMatrix = meta - mat
    matrixNotMeta = mat - meta
    stop = False
    mustFilterMatrix = False
    if len(matrixNotMeta)!=0:
        logging.warn("%d sample names are in the meta data, but not in the expression matrix. Examples: %s" % (len(matrixNotMeta), list(matrixNotMeta)[:10]))
        logging.warn("These samples will be removed from the meta data")
        matrixSampleNames = [x for x in matrixSampleNames if x in meta]
        mustFilterMatrix = True

    if len(metaNotMatrix)!=0:
        logging.warn("%d sample names are in the expression matrix, but not in the meta data. Examples: %s" % (len(metaNotMatrix), list(metaNotMatrix)[:10]))
        logging.warn("These samples will be removed from the expression matrix")

    # filter the meta data file
    logging.info("Data contains %d samples/cells" % len(matrixSampleNames))

    # slurp in the whole meta data, keep track of which fields contain only a single value
    tmpFname = fixedMetaFname+".tmp"
    ofh = open(tmpFname, "w")
    metaToRow = {}
    sep = sepForFile(metaFname)
    fieldValues = defaultdict(set)
    for lNo, line in enumerate(open(metaFname, "rtU")):
        row = line.rstrip("\r\n").split(sep)
        if lNo==0:
            headers = row
            continue
        row = line.rstrip("\r\n").split(sep)
        metaToRow[row[0]] = row

        for fieldIdx, val in enumerate(row):
            fieldValues[fieldIdx].add(val)

    # find fields that contain only a single value
    skipFields = set()
    for fieldIdx, values in iterItems(fieldValues):
        #logging.debug("fieldIdx %d, values %s" % (fieldIdx, values))
        if len(values)==1:
            logging.info("Field %d, '%s', has only a single value. Removing this field from meta data." %
                    (fieldIdx, headers[fieldIdx] ))
            skipFields.add(fieldIdx)

    # write the header line, removing unused fields
    ofh.write("\t".join(sliceRow(headers, skipFields)))
    ofh.write("\n")

    # and write the rows in the right order, also removing unused fields
    for matrixName in matrixSampleNames:
        ofh.write("\t".join(sliceRow(metaToRow[matrixName], skipFields)))
        ofh.write("\n")
    ofh.close()

    os.rename(tmpFname, fixedMetaFname)

    return matrixSampleNames, mustFilterMatrix

def writeCoords(coordName, coords, sampleNames, coordBinFname, coordJson, useTwoBytes, coordInfo, textOutName):
    """ write coordinates given as a dictionary to coordBin and coordJson, in the order of sampleNames
    Also return as a list.
    """
    tmpFname = coordBinFname+".tmp"
    logging.info("Writing coordinates for %s" % (coordName))
    logging.debug("Writing coordinates to %s and %s" % (coordBinFname, coordJson))
    binFh = open(tmpFname, "wb")

    minX = 2^32
    minY = 2^32
    maxX = -2^32
    maxY = -2^32
    xVals = []
    yVals = []

    textOutTmp = textOutName+".tmp"
    textOfh = open(textOutTmp, "w")

    missNames = []

    for sampleName in sampleNames:
        coordTuple = coords.get(sampleName)
        if coordTuple is None:
            #logging.warn("sample name %s is in meta file but not in coordinate file %s, setting to (12345,12345)" % (sampleName, coordName))
            missNames.append(sampleName)
            x = HIDDENCOORD
            y = HIDDENCOORD
        else:
            x, y = coordTuple
            textOfh.write("%s\t%f\t%f\n" % (sampleName, x, y))

        # special values (12345,12345) mean "unknown cellId"
        if x!=HIDDENCOORD and y!=HIDDENCOORD:
            minX = min(x, minX)
            minY = min(y, minY)
            maxX = max(x, maxX)
            maxY = max(y, maxY)

        # all little endian
        if useTwoBytes:
            binFh.write(struct.pack("<H", x))
            binFh.write(struct.pack("<H", y))
        else:
            binFh.write(struct.pack("<f", x))
            binFh.write(struct.pack("<f", y))

        xVals.append(x)
        yVals.append(y)

    if len(missNames)!=0:
        logging.info("%s: %d cells have meta but no coord. E.g. %s" % \
            (coordName, len(missNames), missNames[:3]))

    binFh.close()
    os.rename(tmpFname, coordBinFname)

    md5 = md5WithPython(coordBinFname)
    coordInfo["md5"] = md5[:MD5LEN]

    coordInfo["minX"] = minX
    coordInfo["maxX"] = maxX
    coordInfo["minY"] = minY
    coordInfo["maxY"] = maxY
    if useTwoBytes:
        coordInfo["type"] = "Uint16"
    else:
        coordInfo["type"] = "Float32"

    textOfh.close()
    runGzip(textOutTmp, textOutName)

    logging.debug("Wrote %d coordinates to %s and %s" % (len(sampleNames), coordBinFname, textOutName))
    return coordInfo, xVals, yVals

def runCommand(cmd, verbose=False):
    " run command "
    if verbose:
        logging.info("Running %s" % cmd)
    else:
        logging.debug("Running %s" % cmd)

    if type(cmd)==type([]):
        err = subprocess.call(cmd)
    else:
        err = os.system(cmd)

    if err!=0:
        errAbort("Could not run: %s" % cmd)
    return 0

def copyMatrixTrim(inFname, outFname, filtSampleNames, doFilter, geneToSym, matType):
    """ copy matrix and compress it. If doFilter is true: keep only the samples in filtSampleNames
    Returns the format of the matrix, "float" or "int", or None if not known
    """
    if not doFilter and not ".csv" in inFname.lower():
        logging.info("Copying/compressing %s to %s" % (inFname, outFname))

        # XX stupid .gz heuristics... 
        if inFname.endswith(".gz"):
            cmd = "cp \"%s\" \"%s\"" % (inFname, outFname)
        else:
            cmd = "cat \"%s\" | gzip -c > %s" % (inFname, outFname)
        ret = runCommand(cmd)

        if ret!=0 and isfile(outFname):
            os.remove(outFname)
            sys.exit(1)
        return None

    sep = "\t"

    logging.info("Copying+reordering+trimming %s to %s, keeping %d columns with sample ID in meta" % (inFname, outFname, len(filtSampleNames)))
    logging.debug("matrix type: %s" % matType)

    matIter = MatrixTsvReader(geneToSym)
    matIter.open(inFname, matType=matType)

    sampleNames = matIter.getSampleNames()

    keepFields = set(filtSampleNames)
    keepIdx = []
    for i, name in enumerate(sampleNames):
        if name in keepFields:
            keepIdx.append(i)

    assert(len(keepIdx)!=0)
    logging.debug("Keeping %d fields" % len(keepIdx))

    tmpFname = outFname+".tmp"

    ofh = openFile(tmpFname, "w")
    ofh.write("gene\t")
    ofh.write("\t".join(filtSampleNames))
    ofh.write("\n")

    count = 0
    for geneId, sym, exprArr in matIter.iterRows():
        newRow = [geneId+"|"+sym]
        for idx in keepIdx:
            newRow.append(str(exprArr[idx]))
        ofh.write("\t".join(newRow))
        ofh.write("\n")
        count += 1
        if count%1000==0:
            logging.info("Wrote %d text rows" % count)
    ofh.close()
    matIter.close()

    #tmpFnameGz = outFname+".tmp.gz"
    #runCommand("gzip -c %s > %s " % (tmpFname, tmpFnameGz))
    #os.remove(tmpFname)
    #os.rename(tmpFnameGz, outFname)
    runGzip(tmpFname, outFname)

    return matIter.getMatType()

def convIdToSym(geneToSym, geneId, printWarning=True):
    if geneToSym is None:
        return geneId
    else:
        if geneId not in geneToSym:
            if printWarning:
                logging.warn("Could not convert geneId %s to symbol" % geneId)
            return geneId
        return geneToSym[geneId]

def to_camel_case(snake_str):
    components = snake_str.split('_')
    # We capitalize the first letter of each component except the first one                                     # with the 'title' method and join them together.
    return components[0] + ''.join(x.title() for x in components[1:])

def sanitizeName(name):
    " remove all nonalpha chars, allow underscores "
    assert(name!=None)
    #newName = to_camel_case(name.replace(" ", "_"))
    newName = ''.join([ch for ch in name if (ch.isalnum() or ch=="_")])
    if newName!=name:
        logging.debug("Sanitizing %s -> %s" % (repr(name), newName))
    assert(len(newName)!=0)
    return newName

def sanitizeHeader(name):
    " for tab-sep tables: replace nonalpha chars with  underscores "
    assert(name!=None)
    #newName = to_camel_case(name.replace(" ", "_"))
    newName = re.sub("[^a-zA-Z0-9_]","_", name)
    newName = re.sub("^_","", newName)  # remove _ prefix
    logging.debug("Sanitizing %s -> %s" % (repr(name), newName))
    return newName

def parseMarkerTable(filename, geneToSym):
    " parse marker gene table and return dict clusterName -> list of rows (geneId, geneSym, otherFields...)"
    logging.debug("Reading cluster markers from %s" % (filename))
    ifh = openFile(filename)

    seuratLine = '\tp_val\tavg_logFC\tpct.1\tpct.2\tp_val_adj\tcluster\tgene'
    seuratLine2 = '"","p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene"'
    seuratLine3 = ",p_val,avg_logFC,pct.1,pct.2,p_val_adj,cluster,gene"
    headerLine = ifh.readline().rstrip("\r\n")

    sep = sepForFile(filename)
    if headerLine == seuratLine or headerLine == seuratLine2 or headerLine == seuratLine3:
        logging.debug("Cluster marker file %s was recognized to be in Seurat format" % filename)
        # field 0 is not the gene ID, it has some weird suffix appended.
        headers = ["rowNameFromR", "pVal", "avg. logFC", "PCT1", "PCT2", "pVal adj.", "Cluster", "Gene"]
        geneIdx = 7
        scoreIdx = 1
        clusterIdx = 6
        otherStart = 2
        otherEnd = 6
    else:
        logging.info("Reading %s: assuming marker file format (cluster, gene, score) + any other fields" % filename)
        headers = headerLine.split(sep)
        clusterIdx = 0
        geneIdx = 1
        scoreIdx = 2
        otherStart = 3
        otherEnd = 9999

    otherHeaders = headers[otherStart:otherEnd]
    logging.debug("Other headers: %s" % otherHeaders)

    data = defaultdict(list)
    otherColumns = defaultdict(list)
    for line in ifh:
        row = line.rstrip("\r\n").split(sep)
        clusterName = row[clusterIdx]
        geneId = row[geneIdx]
        scoreVal = float(row[scoreIdx])
        otherFields = row[otherStart:otherEnd]

        for colIdx, val in enumerate(otherFields):
            otherColumns[colIdx].append(val)

        geneSym = convIdToSym(geneToSym, geneId, printWarning=False)

        newRow = []
        newRow.append(geneId)
        newRow.append(geneSym)
        newRow.append(scoreVal)
        newRow.extend(otherFields)
        data[clusterName].append(newRow)

    # annotate otherColumns with their data type, separated by |. This is optional and only used for sorting
    # in the user interface
    otherColType = {}
    for colIdx, vals in iterItems(otherColumns):
        otherColType[colIdx] = typeForStrings(vals)

    otherHeadersWithType = []
    for colIdx, header in enumerate(otherHeaders):
        if otherColType[colIdx]!="string":
            header = header+"|"+otherColType[colIdx]
        otherHeadersWithType.append(header)
    logging.debug("Other headers with type: %s" % otherHeadersWithType)

    newHeaders = ["id", "symbol", headers[scoreIdx]+"|float"]
    newHeaders.extend(otherHeadersWithType)

    if len(data) > 200:
        errAbort("Your marker file has more than 200 clusters. Are you sure that this is correct? The input format is (clusterName, geneSymName, Score), is it possible that you have inversed the order of cluster and gene?")

    # determine if the score field is most likely a p-value, needed for sorting
    revSort = True
    scoreHeader = headers[scoreIdx].lower()
    logging.debug("score field has name '%s'" % scoreHeader)
    if scoreHeader in ["p_val", "p-val", "p.val", "pval", "fdr"]:
        logging.debug("score field name '%s' looks like a p-Value, sorting normally" % scoreHeader)
        revSort = False

    for clusterName, rows in iterItems(data):
        rows.sort(key=operator.itemgetter(2), reverse=revSort)

    return data, newHeaders

def splitMarkerTable(filename, geneToSym, outDir):
    """ split .tsv on first field and create many files in outDir with columns 2-end.
        Returns the names of the clusters and a dict topMarkers with clusterName -> list of three top marker genes.
    """
    topMarkerCount = 5

    if filename is None:
        return

    data, newHeaders = parseMarkerTable(filename, geneToSym)

    logging.debug("Splitting cluster markers into directory %s" % (outDir))
    fileCount = 0
    sanNames = set()
    topMarkers = {}
    for clusterName, rows in iterItems(data):
        logging.debug("Cluster: %s" % clusterName)
        sanName = sanitizeName(clusterName)
        assert(sanName not in sanNames) # after sanitation, cluster names must be unique
        sanNames.add(sanName)

        outFname = join(outDir, sanName+".tsv")
        logging.debug("Writing %s" % outFname)
        ofh = open(outFname, "w")
        ofh.write("\t".join(newHeaders))
        ofh.write("\n")
        for row in rows:
            row[2] = "%0.5E" % row[2] # limit score to 5 digits
            ofh.write("\t".join(row))
            ofh.write("\n")

        topSyms = [row[1] for row in rows[:topMarkerCount]]
        topMarkers[clusterName] = topSyms

        ofh.close()

        runGzip(outFname)

        fileCount += 1
    logging.info("Wrote %d .tsv.gz files into directory %s" % (fileCount, outDir))
    return data.keys(), topMarkers

#def guessConfig(options):
    #" guess reasonable config options from arguments "
    #conf = {}
    #conf.name = dirname(options.matrix)

    #if options.inDir:
        #inDir = options.inDir
        #metaFname = join(inDir, "meta.tsv")
        #matrixFname = join(inDir, "exprMatrix.tsv")
        #coordFnames = [join(inDir, "tsne.coords.tsv")]
        #markerFname = join(inDir, "markers.tsv")
        #if isfile(markerFname):
            #markerFnames = [markerFname]
        #else:
            #markerFnames = None
#
        #acronymFname = join(inDir, "acronyms.tsv")
        #if isfile(acronymFname):
            #otherFiles["acronyms"] = [acronymFname]
#
        #markerFname = join(inDir, "markers.tsv")
        #if isfile(acronymFname):
            #otherFiles["markerLists"] = [markerFname]
    #return conf

def copyDatasetHtmls(inDir, outConf, datasetDir):
    " copy dataset description html files to output directory. Add md5s to outConf. "
    filesToCopy = []

    outConf["descMd5s"] = {}

    for fileBase in ["summary.html", "methods.html", "downloads.html", "thumb.png", "protocol.pdf", "desc.conf"]:
        inFname = makeAbs(inDir, fileBase)
        if not isfile(inFname):
            logging.info("%s does not exist" % inFname)
        else:
            #copyFiles.append( (fname, "summary.html") )
            outPath = join(datasetDir, fileBase)
            logging.debug("Copying %s -> %s" % (inFname, outPath))
            shutil.copyfile(inFname, outPath)
            outConf["descMd5s"][fileBase.split(".")[0]] = md5ForFile(inFname)[:MD5LEN]

def copyImage(inDir, summInfo, datasetDir):
    """ copy image to datasetDir and write size to summInfo["imageWidth"] and "imageHeight" """
    inFname = join(inDir, summInfo["image"])
    logging.debug("Copying %s to %s" % (inFname, datasetDir))
    #shutil.copy(inFname, datasetDir)
    shutil.copyfile(inFname, join(datasetDir, basename(inFname)))

    cmd = ["file", inFname]
    proc, stdout = popen(cmd)
    imgLine = stdout.readline()
    # thumb.png: PNG image data, 400 x 267, 8-bit/color RGBA, non-interlaced
    sizeStr = imgLine.split(":")[1].split(", ")[1]
    width, height = sizeStr.split(" x ")
    summInfo["image"] = (summInfo["image"], width, height)
    return summInfo

def readFile(fname, encoding="utf8"):
    " read entire file with encoding "
    if isPy3:
        fh = open(fname, "r", encoding=encoding)
        text = fh.read()
    else:
        fh = open(fname, "r")
        text = fh.read().decode(encoding)
    fh.close()
    return text

def readFileIntoDict(summInfo, key, inDir, fname, mustExist=False, encoding="utf8"):
    " return file with encoding as string into dictionary "
    fname = join(inDir, fname)
    if not isfile(fname):
        if mustExist:
            errAbort("%s does not exist" % fname)
        else:
            return

    text = readFile(fname)
    summInfo[key] = text

def writeDatasetDesc(inDir, outConf, datasetDir, coordFiles=None):
    " write a json file that describes the dataset abstract/methods/downloads, easier than summary/methods/downloads.html "
    confFname = join(inDir, "datasetDesc.conf")
    if not isfile(confFname):
        confFname = join(inDir, "desc.conf")

    if not isfile(confFname):
        logging.debug("Could not find %s" % confFname)
        return False

    if "fileVersions" not in outConf:
        outConf["fileVersions"] = {}

    outConf["fileVersions"]["desc"] = getFileVersion(confFname)
    outPath = join(datasetDir, "desc.json")

    summInfo = loadConfig(confFname, requireTags=[])

    if coordFiles:
        summInfo["coordFiles"] = coordFiles

    # try various ways to get the abstract and methods html text
    readFileIntoDict(summInfo, "abstract", inDir, "abstract.html")
    readFileIntoDict(summInfo, "abstract", inDir, "summary.html")
    readFileIntoDict(summInfo, "methods", inDir, "methods.html")
    if "abstractFile" in summInfo:
        readFileIntoDict(summInfo, "abstract", inDir, summInfo["abstractFile"], mustExist=True)
        del summInfo["abstractFile"]
    if "methodsFile" in summInfo:
        readFileIntoDict(summInfo, "methods", inDir, summInfo["methodsFile"], mustExist=True)
        del summInfo["methodsFile"]

    # import the unit description from cellbrowser.conf
    if "unit" in outConf and not "unitDesc" in summInfo:
        summInfo["unitDesc"] = outConf["unit"]

    # copy over the raw matrix file, usually this is a zip or gzip file
    if "rawMatrixFile" in summInfo:
        rawInPath = join(inDir, summInfo["rawMatrixFile"])
        rawOutPath = join(datasetDir, summInfo["rawMatrixFile"])
        if not isfile(rawOutPath) or getsize(rawInPath)!=getsize(rawOutPath):
            logging.info("Copying %s to %s" % (rawInPath, rawOutPath))
            shutil.copyfile(rawInPath, rawOutPath)

    # need the collection info, too
    if "collections" in outConf and not "collections" in summInfo:
        summInfo["collections"] = outConf["collections"]

    if "image" in summInfo:
        summInfo = copyImage(inDir, summInfo, datasetDir)

    writeJson(summInfo, outPath)

    if "descMd5s" not in outConf:
        outConf["descMd5s"] = {}

    outConf["descMd5s"]["datasetDesc"] = md5ForFile(confFname)[:MD5LEN]
    logging.debug("Wrote %s" % outPath)
    return True

def makeAbs(inDir, fname):
    " return absolute path of fname under inDir "
    if fname is None:
        return None
    return abspath(join(inDir, fname))

def makeAbsDict(conf, key):
    " given list of dicts with key 'file', assume they are relative to inDir and make their paths absolute "
    inDir = conf["inDir"]
    dicts = conf[key]
    for d in dicts:
        d["file"] = makeAbs(inDir, d["file"])
    return dicts

def parseTsvColumn(fname, colName):
    " parse a tsv file and return column as a pair (values, assignment row -> index in values) "
    logging.info("Parsing column %s from %s" % (colName, fname))
    vals = parseOneColumn(fname, colName)

    newVals = []
    valToInt = {}
    maxIdx = -1
    for v in vals:
        if v not in valToInt:
            maxIdx+=1
            valToInt[v] = maxIdx
        idx = valToInt[v]
        newVals.append(idx)


    # inverse key/val dict
    intToVal = {}
    for k, v in iterItems(valToInt):
        intToVal[v] = k

    valArr = []
    for i in range(0, maxIdx+1):
        valArr.append(intToVal[i])

    return newVals, valArr

def makeMids(xVals, yVals, labelVec, labelVals, coordInfo):
    """
    calculate the positions (centers) for the cluster labels
    given a coord list and a vector of the same size with the label indices, return a list of [x, y, coordLabel]
    """
    logging.debug("Making cluster labels for %s" % labelVals)
    assert(len(xVals)==len(labelVec)==len(yVals))

    # prep the arrays
    clusterXVals = []
    clusterYVals = []
    for i in range(len(labelVals)):
        clusterXVals.append([])
        clusterYVals.append([])
    assert(len(clusterXVals)==len(labelVals))

    # sort the coords into separate arrays, one per cluster
    for i in range(len(labelVec)):
        #for (x, y), clusterIdx in zip(coords, labelVec):
        clusterIdx = labelVec[i]
        x = xVals[i]
        y = yVals[i]
        if x==HIDDENCOORD and y==HIDDENCOORD:
            continue
        clusterXVals[clusterIdx].append(x)
        clusterYVals[clusterIdx].append(y)

    midInfo = []
    for clustIdx, xList in enumerate(clusterXVals):
        clusterName = labelVals[clustIdx]
        if len(xList)==0:
            midInfo.append([HIDDENCOORD, HIDDENCOORD, clusterName])
            continue

        yList = clusterYVals[clustIdx]
        # get the midpoint of this cluster
        midX = sum(xList) / float(len(xList))
        midY = sum(yList) / float(len(yList))

        if len(xList)<3:
            midInfo.append([midX, midY, clusterName])
            continue

        # take only the best 70% of the points closest to the midpoints
        xyDist = []
        for x, y in zip(xList, yList):
            dist = math.sqrt((x-midX)**2+(y-midY)**2)
            xyDist.append( (dist, x, y) )
        xyDist.sort()
        xyDistBest = xyDist[:int(0.7*len(xyDist))]

        # now recalc the midpoint
        xSum = sum([x for dist, x, y in xyDistBest])
        ySum = sum([y for dist, x, y in xyDistBest])
        fixMidX = xSum / float(len(xyDistBest))
        fixMidY = ySum / float(len(xyDistBest))

        midInfo.append([fixMidX, fixMidY, clusterName])

    # make some minimal effort to reduce overlaps
    #spanX = coordInfo['maxX'] - coordInfo['minX']
    #spanY = coordInfo['maxY'] - coordInfo['minY']
    #tickX = spanX / 1000 # rough guess how much one pixel could be on
    #tickY = spanY / 1000 # the final screen
    #for i, (midX1, midY1, clusterName1) in enumerate(midInfo):
        #print "first", i, midX1, midY1, clusterName1
        #for j, (midX2, midY2, clusterName2) in enumerate(midInfo[i+1:]):
            #print "second", j, midX2, midY2, clusterName1, clusterName2
            #distX = abs(midX2-midX1)
            #distY = abs(midY2-midY1)
            #print distX, distY
            ## if distance between two labels too short:
            #dist = math.sqrt((((midX2-midX1)/tickX)**2+((midY2-midY1)/tickY)**2))
            #print "dist in pixels", dist
            #if dist< 30:
                #print "moving"
                #print "before", midInfo[j]
                ## move the first label slightly downwards
                #midInfo[j][1] = midY1 + 5 * tickY
                #print "after", midInfo[j]

    return midInfo

def readHeaders(fname):
    " return headers of a file "
    logging.info("Reading headers of file %s" % fname)
    ifh = openFile(fname, "rtU")
    line1 = ifh.readline().rstrip("\r\n")
    sep = sepForFile(fname)
    row = line1.split(sep)
    row = [x.rstrip('"').lstrip('"') for x in row]
    logging.debug("Found %d fields, e.g. %s" % (len(row), row[:3]))
    return row

def parseGeneInfo(geneToSym, fname):
    """ parse a file with three columns: symbol, desc (optional), pmid (optional).
    Return as a dict symbol -> [description, pmid] """
    if fname is None:
        return {}
    logging.info("Parsing %s" % fname)
    validSyms = None
    if geneToSym is not None:
        validSyms = set()
        for gene, sym in iterItems(geneToSym):
            validSyms.add(sym)

    sep = sepForFile(fname)
    geneInfo = []
    for line in openFile(fname):
        if line.startswith("#"):
            continue
        hasDesc = False
        hasPmid = False
        if line.startswith("symbol"):
            continue
        row = line.rstrip("\r\n").split(sep)
        sym = row[0]
        if validSyms is not None and sym not in validSyms:
            sym = geneToSym.get(sym)
            if sym is None:
                logging.error("'%s' is not a valid gene gene symbol, skipping it" % sym)
                continue

        info = [sym]
        if len(row)>1:
            info.append(row[1])
        if len(row)>2:
            info.append(row[2])
        geneInfo.append(info)
    return geneInfo

def readSampleNames(fname):
    " read only the first column of fname, strip the headers "
    logging.info("Reading sample names from %s" % fname)
    sampleNames = []
    i = 1
    doneNames = set()
    for row in lineFileNextRow(fname):
        metaName = row[0]
        if metaName=="":
            errAbort("invalid sample name - line %d in %s: sample name (first field) is empty" % (i, fname))
        if metaName in doneNames:
            errAbort("sample name duplicated - line %d in %s: sample name %s (first field) has been seen before" % (i, fname, metaName))

        doneNames.add(metaName)
        sampleNames.append(row[0])
        i+=1
    logging.debug("Found %d sample names, e.g. %s" % (len(sampleNames), sampleNames[:3]))
    return sampleNames

def guessGeneIdType(inputObj):
    """ Accepts a list of gene IDs or a matrix file name.
    returns 'gencode-human', 'gencode-mouse' or 'symbols' depending on the first gene
    """
    if type(inputObj)==type(list()):
        geneIds = inputObj
    else:
        matrixFname = inputObj
        matIter = MatrixTsvReader()
        matIter.open(matrixFname, usePyGzip=True)

        geneId, sym, exprArr = nextEl(matIter.iterRows())
        matIter.close()
        geneIds = [geneId]

    geneId = geneIds[0]

    if geneId.startswith("ENSG"):
        geneType = "gencode-human"
    elif geneId.startswith("ENSMUSG"):
        geneType =  "gencode-mouse"
    elif geneId.startswith("KH2013:"):
        geneType = "ciona-kh2013"
    elif geneId.isdigit():
        geneType = "entrez-human"
    else:
        geneType = "symbols"

    logging.info("Auto-detected gene IDs type: %s" % (geneType))
    return geneType

def convertExprMatrix(inConf, outMatrixFname, outConf, metaSampleNames, geneToSym, outDir, needFilterMatrix):
    """ trim a copy of the expression matrix for downloads, also create an indexed
    and compressed version
    """
    matType = inConf.get("matrixType")

    # step1: copy expression matrix, so people can download it, potentially
    # removing those sample names that are not in the meta data
    matrixFname = getAbsPath(inConf, "exprMatrix")
    outConf["fileVersions"]["inMatrix"] = getFileVersion(matrixFname)
    try:
        matType = copyMatrixTrim(matrixFname, outMatrixFname, metaSampleNames, needFilterMatrix, geneToSym, matType)
    except ValueError:
        logging.warn("This is rare: mis-guessed the matrix data type, trying again and using floating point numbers. To avoid this message in the future, you can set matrixType='float' in cellbrowser.conf.")
        matType = copyMatrixTrim(matrixFname, outMatrixFname, metaSampleNames, needFilterMatrix, geneToSym, "float")

    # step2: compress matrix and index to file
    binMat = join(outDir, "exprMatrix.bin")
    binMatIndex = join(outDir, "exprMatrix.json")
    discretBinMat = join(outDir, "discretMat.bin")
    discretMatrixIndex = join(outDir, "discretMat.json")

    matType = matrixToBin(outMatrixFname, geneToSym, binMat, binMatIndex, discretBinMat, discretMatrixIndex, matType=matType)

    if matType=="int":
        outConf["matrixArrType"] = "Uint32"
    elif matType=="float":
        outConf["matrixArrType"] = "Float32"
    else:
        assert(False)

    outConf["fileVersions"]["outMatrix"] = getFileVersion(outMatrixFname)

def copyConf(inConf, outConf, keyName):
    " copy value of keyName from inConf dict to outConf dict "
    if keyName in inConf:
        outConf[keyName] = inConf[keyName]

def convertCoords(inConf, outConf, sampleNames, outMeta, outDir):
    " convert the coordinates "
    coordFnames = makeAbsDict(inConf, "coords")

    flipY = inConf.get("flipY", False)
    useTwoBytes = inConf.get("useTwoBytes", None)

    hasLabels = False
    if "labelField" in inConf and inConf["labelField"] is not None:
        hasLabels = True
        clusterLabelField = inConf["labelField"]
        labelVec, labelVals = parseTsvColumn(outMeta, clusterLabelField)
        outConf["labelField"] = clusterLabelField

    outFnames = []
    newCoords = []
    for coordIdx, inCoordInfo in enumerate(coordFnames):
        coordFname = inCoordInfo["file"]
        coordLabel = inCoordInfo["shortLabel"]
        logging.info("Parsing coordinates for "+coordLabel)
        coords = parseScaleCoordsAsDict(coordFname, useTwoBytes, flipY)
        coordName = "coords_%d" % coordIdx
        coordDir = join(outDir, "coords", coordName)
        makeDir(coordDir)
        coordBin = join(coordDir, "coords.bin")
        coordJson = join(coordDir, "coords.json")

        coordInfo = OrderedDict()
        coordInfo["name"] = coordName
        coordInfo["shortLabel"] = coordLabel
        if "radius" in inCoordInfo:
            coordInfo["radius"] = inCoordInfo["radius"]
        if "colorOnMeta" in inCoordInfo:
            coordInfo["colorOnMeta"] = inCoordInfo["colorOnMeta"]

        cleanName = sanitizeName(coordLabel.replace(" ", "_"))
        textOutBase = cleanName+".coords.tsv.gz"
        textOutName = join(outDir, textOutBase)
        outFnames.append(textOutBase)
        coordInfo, xVals, yVals = writeCoords(coordLabel, coords, sampleNames, coordBin, coordJson, useTwoBytes, coordInfo, textOutName)

        if hasLabels:
            logging.info("Calculating cluster midpoints for "+coordLabel)
            clusterMids= makeMids(xVals, yVals, labelVec, labelVals, coordInfo)
            clusterOrder = orderClusters(clusterMids, outConf)

            clusterInfo = {}
            clusterInfo["labels"] = clusterMids
            clusterInfo["order"] = clusterOrder

            clusterLabelFname = join(coordDir, "clusterLabels.json")
            midFh = open(clusterLabelFname, "w")
            json.dump(clusterInfo, midFh, indent=2)
            logging.debug("Wrote cluster labels, midpoints and order to %s" % clusterLabelFname)
            addMd5(coordInfo, clusterLabelFname, keyName="labelMd5")

        newCoords.append( coordInfo )

    outConf["coords"] = newCoords
    copyConf(inConf, outConf, "labelField")
    copyConf(inConf, outConf, "useTwoBytes")
    return outFnames, labelVals

def readAcronyms(inConf, outConf):
    " read the acronyms and save them into the config "
    inDir = inConf["inDir"]
    fname = inConf.get("acroFname")
    if fname is not None:
        fname = makeAbs(inDir, fname)
        if not isfile(fname):
            logging.warn("%s specified in config file, but does not exist, skipping" % fname)
            acronyms = None
        else:
            acronyms = parseDict(fname)
            logging.info("Read %d acronyms from %s" % (len(acronyms), fname))
            #outConf["acronyms"] = acronyms
        return acronyms

def checkClusterNames(markerFname, clusterNames, clusterLabels, doAbort):
    " make sure that the cluster names from the meta.tsv match the ones from markers.tsv "
    markerClusters = set(clusterNames)
    labelClusters = set(clusterLabels)
    notInLabels = markerClusters - labelClusters
    notInMarkers = labelClusters - markerClusters
    if len(notInMarkers)!=0:
        msg = ("%s: the following cluster names are in the meta file but not in the marker file: %s. "+
        "Please fix one of the files, clicks onto a label will otherwise not work.") % (markerFname, notInMarkers)
        if doAbort:
            errAbort(msg)
        else:
            logging.error(msg)

    if len(notInLabels)!=0:
        logging.warn(("%s: the following cluster names are in the marker file but not in the meta file: %s. "+
                "Users may not notice the problem, but it may indicate an erronous meta data file.") % \
                (markerFname, notInLabels))

def convertMarkers(inConf, outConf, geneToSym, clusterLabels, outDir):
    """ split the marker tables into one file per cluster and add filenames as 'markers' in outConf
    also add the 'topMarkers' to outConf, the top five markers for every cluster.
    """
    markerFnames = []
    if "markers" in inConf:
        markerFnames = makeAbsDict(inConf, "markers")

    newMarkers = []
    #doAbort = True # only the first marker file leads to abort, we're more tolerant for the others
    doAbort = False # temp hack
    topMarkersDone = False
    for markerIdx, markerInfo in enumerate(markerFnames):
        markerFname = markerInfo["file"]
        markerLabel = markerInfo["shortLabel"]

        clusterName = "markers_%d" % markerIdx # use sha1 of input file ?
        markerDir = join(outDir, "markers", clusterName)
        makeDir(markerDir)

        clusterNames, topMarkers = splitMarkerTable(markerFname, geneToSym, markerDir)
        
        # only use the top markers of the first marker file
        if not topMarkersDone:
            outConf["topMarkers"] = topMarkers
            topMarkersDone = True

        checkClusterNames(markerFname, clusterNames, clusterLabels, doAbort)
        doAbort = False

        newDict = {"name" : sanitizeName(clusterName), "shortLabel" : markerLabel}
        if "selectOnClick" in markerInfo:
            newDict["selectOnClick"] = markerInfo["selectOnClick"]
        newMarkers.append( newDict )

    outConf["markers"] = newMarkers

def readQuickGenes(inConf, geneToSym, outConf):
    quickGeneFname = inConf.get("quickGenesFile")
    if quickGeneFname:
        fname = getAbsPath(inConf, "quickGenesFile")
        quickGenes = parseGeneInfo(geneToSym, fname)
        outConf["quickGenes"] = quickGenes
        logging.info("Read %d quick genes from %s" % (len(quickGenes), fname))

def getFileVersion(fname):
    data = {}
    data["fname"] = fname
    addMd5(data, fname, shortMd5=False)
    data["size"] = getsize(fname)
    data["mtime"] = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(getmtime(fname)))
    return data

def checkFieldNames(outConf, fieldNames, validFieldNames):
    " make sure that all fieldNames in outConf are valid field names. errAbort is not. "
    for fn in fieldNames:
        if fn not in outConf:
            continue

        if not outConf[fn] in validFieldNames:
            errAbort("Config statement '%s' contains an invalid field name, '%s'. Valid meta field names are: %s" % \
                (fn, outConf[fn], ", ".join(validFieldNames)))

def convertMeta(inConf, outConf, outDir):
    """ convert the meta data to binary files. The new meta is re-ordered, so it's in the same
    order as the samples in the expression matrix.
    """
    if not "fileVersions" in outConf:
        outConf["fileVersions"] = {}

    metaFname = getAbsPath(inConf, "meta")
    outConf["fileVersions"]["inMeta"] = getFileVersion(metaFname)

    metaDir = join(outDir, "metaFields")
    makeDir(metaDir)
    metaIdxFname = join(outDir, "meta.index")

    finalMetaFname = join(outDir, "meta.tsv")

    matrixFname = getAbsPath(inConf, "exprMatrix")
    sampleNames, needFilterMatrix = metaReorder(matrixFname, metaFname, finalMetaFname)

    outConf["sampleCount"] = len(sampleNames)
    outConf["matrixWasFiltered"] = needFilterMatrix

    colorFname = inConf.get("colors")
    enumFields = inConf.get("enumFields")
    fieldConf, validFieldNames = metaToBin(inConf, outConf, finalMetaFname, colorFname, metaDir, enumFields)
    outConf["metaFields"] = fieldConf

    checkFieldNames(outConf, ["violinField", "clusterField", "labelField"], validFieldNames)

    indexMeta(finalMetaFname, metaIdxFname)

    logging.info("Kept %d cells present in both meta data file and expression matrix" % len(sampleNames))

    outConf["fileVersions"]["outMeta"] = getFileVersion(finalMetaFname)

    return sampleNames, needFilterMatrix, finalMetaFname

def readGeneSymbols(geneIdType, matrixFnameOrGeneIds):
    " return geneToSym, based on gene tables "
    if geneIdType==None or geneIdType=="auto":
        geneIdType = guessGeneIdType(matrixFnameOrGeneIds)

    if geneIdType.startswith('symbol'):
        return None

    geneIdTable = getStaticFile(join("genes", geneIdType+".symbols.tsv.gz"))
    geneToSym = readGeneToSym(geneIdTable)
    return geneToSym

def readMitos(idType):
    ' return the gene IDs of all mitochondrial genes. 37 for human for all gencode versions '
    logging.debug("Creating mito list for gene ID type %s" % idType)

    if idType.startswith("symbol"):
        human = readGeneSymbols("gencode-human", None)
        mouse = readGeneSymbols("gencode-human", None)
        mitos = []
        for sym in human.values():
            if sym.startswith("MT-"):
                mitos.append(sym)
        for sym in mouse.values():
            if sym.startswith("mt-"):
                mitos.append(sym)
        logging.info("Built list of human and mouse mitochondrial symbols, %d symbols found" % len(mitos))
        if len(mitos)>100:
            errAbort("Too many mitos found. %s" % mitos)
        return mitos

    geneToSym = readGeneSymbols(idType, None)
    mitos = []
    allGeneIds = []
    for geneId, sym in iterItems(geneToSym):
        allGeneIds.append(geneId)
        if sym.lower().startswith("mt-"):
            mitos.append(geneId)
    if len(mitos)==0:
        errAbort("Could not find any mitochondrial genes in cell browser gene lists for gene ID type %s. Example gene IDs from input file: %s" % (idType, ",".join(allGeneIds[:10])))
    logging.debug("Found %d mitochondrial genes for %s, e.g. %s" % (len(mitos), idType, mitos[0]))
    return mitos

def getAbsPath(conf, key):
    " get assume that value of key in conf is a filename and use the inDir value to make it absolute "
    return abspath(join(conf["inDir"], conf[key]))

def popen(cmd, shell=False, useStderr=False, doWait=True):
    " run command and return proc object with its stdout attribute  "
    if useStderr:
        if isPy3:
            proc = subprocess.Popen(cmd, stderr=subprocess.PIPE, encoding="utf8", shell=shell)
        else:
            proc = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=shell)
    else:
        if isPy3:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, encoding="utf8", shell=shell)
        else:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=shell)

    if doWait:
        proc.wait()

    if doWait and proc.returncode!=0:
        logging.debug("Cmd %s return non-zero exit code %s" % (cmd, str(proc.returncode)))
        return None, None
    else:
        if useStderr:
            return proc, proc.stderr
        else:
            return proc, proc.stdout

def getMd5Using(md5Cmd, fname):
    " posix command line tool is much faster than python "
    logging.debug("Getting md5 of %s using %s command line tool" % (fname, md5Cmd))
    cmd = [md5Cmd, fname]
    logging.debug("Cmd: %s" % cmd)
    proc, stdout = popen(cmd)
    md5 = stdout.readline()
    #proc.stdout.close()
    #stat = os.waitpid(proc.pid, 0)
    #err = stat[1]
    #assert(err==0)
    return md5

def md5ForList(l):
    " given a list of strings, return their md5 "
    hash_md5 = hashlib.md5()
    for s in l:
        hash_md5.update(s.encode("utf8"))
    md5 = hash_md5.hexdigest()
    return md5

def md5WithPython(fname):
    " get md5 using python lib. OK for small files, very slow for big files. "
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            hash_md5.update(chunk)
    md5 = hash_md5.hexdigest()
    return md5

def md5ForFile(fname, isSmall=False):
    " return the md5sum of a file. Use a command line tool, if possible. "
    logging.debug("Getting md5 of %s" % fname)
    if spawn.find_executable("md5sum")!=None and not isSmall:
        md5 = getMd5Using("md5sum", fname).split()[0]
    elif spawn.find_executable("md5")!=None and not isSmall:
        md5 = getMd5Using("md5", fname).split()[-1]
    else:
        md5 = md5WithPython(fname)
    return md5

def matrixOrSamplesHaveChanged(datasetDir, inMatrixFname, outMatrixFname, outConf):
    """ compare filesize stored in datasetDir/cellbrowser.json.bak with file
    size of inMatrixFname and also compare the sample names with the sample names in
    outMatrixFname
    """
    logging.info("Determining if %s needs to be created" % outMatrixFname)
    if not isfile(outMatrixFname):
        logging.info("%s does not exist." % outMatrixFname)
        return True

    confName = join(datasetDir, "dataset.json")
    if not isfile(confName):
        logging.info("%s does not exist. This looks like the first run with this output directory" % confName)
        return True

    try:
        lastConf = json.load(open(confName))
    except json.decoder.JSONDecodeError:
        errAbort("Is the file %s broken? Please remove the file and run this command again." % confName)

    if not "fileVersions" in lastConf or not "inMatrix" in lastConf["fileVersions"] \
        or not "outMatrix" in lastConf["fileVersions"]:
            logging.warn("Internal error? Missing 'fileVersions' tag in %s" % confName)
            return True

    oldMatrixInfo = lastConf["fileVersions"]["inMatrix"]
    origSize = oldMatrixInfo ["size"]
    nowSize = getsize(inMatrixFname)
    matrixIsSame = (origSize==nowSize)

    if not matrixIsSame:
        logging.info("input matrix has input file size that is different from prevously processed matrix, have to reindex the expression matrix. Old file: %s, current file: %d" % (oldMatrixInfo, nowSize))
        return True
    outConf["fileVersions"]["inMatrix"] = oldMatrixInfo
    outConf["fileVersions"]["outMatrix"] = lastConf["fileVersions"]["outMatrix"]
    outConf["matrixArrType"] = lastConf["matrixArrType"]

    # this obscure command gets file with the the cell identifiers in the dataset directory
    sampleNameFname = join(datasetDir, "metaFields", outConf["metaFields"][0]["name"]+".bin.gz")
    logging.debug("Reading meta sample names from %s" % sampleNameFname)

    # python3 has 'text mode' but python2 doesn't have that so decode explicitely
    metaSampleNames = []
    for line in gzip.open(sampleNameFname, "r"):
        metaSampleNames.append(line.decode("utf8").rstrip("\n\r"))

    outMatrixFname = join(datasetDir, "exprMatrix.tsv.gz")
    matrixSampleNames = readHeaders(outMatrixFname)[1:]
    assert(matrixSampleNames!=0)

    if metaSampleNames!=matrixSampleNames:
        logging.info("meta sample samples from previous run are different from sample names in current matrix, have to reindex the matrix. Counts: %d vs. %d" % (len(metaSampleNames), len(matrixSampleNames)))
        return True

    logging.info("current input matrix looks identical to previously processed matrix, same file size, same sample names")
    return False

def readJson(fname, keepOrder=False):
    " read .json and return as a dict "
    if keepOrder:
        customdecoder = json.JSONDecoder(object_pairs_hook=OrderedDict)
        inStr = readFile(fname)
        data = customdecoder.decode(inStr)
    else:
        data = json.load(open(fname))
    return data

def orderClusters(labelCoords, outConf):
    " given the cluster label coordinates, order them by similarity "
    #labelCoords = readJson(clusterLabelFname)

    # create dict with label1 -> list of (dist, label2)
    dists = defaultdict(list)
    for i, (x1, y1, label1) in enumerate(labelCoords):
        for x2, y2, label2 in labelCoords[i:]:
            dist = math.sqrt((x2-x1)**2+(y2-y1)**2)
            dists[label1].append((dist, label2))
            dists[label2].append((dist, label1))

    labelOrder = []
    doneLabels = set()

    currLabel = labelCoords[0][2]
    for i in range(0, len(labelCoords)):
        otherLabels = dists[currLabel]
        otherLabels.sort()
        for dist, label in otherLabels:
            if not label in doneLabels:
                currLabel = label
                labelOrder.append(currLabel)
                doneLabels.add(currLabel)
                break

    return labelOrder

def convertDataset(inConf, outConf, datasetDir):
    """ convert everything needed for a dataset to datasetDir, write config to outConf.
    If the expression matrix has not changed since the last run, and the sampleNames are the same,
    it won't be converted again.
    """
    # some config settings are passed through unmodified to the javascript
    for tag in ["name", "shortLabel", "radius", "alpha", "priority", "tags",
        "clusterField", "xenaPhenoId", "xenaId", "hubUrl", "showLabels", "ucscDb",
        "unit", "violinField", "visibility", "collections"]:
        copyConf(inConf, outConf, tag)
        if tag in ["collections"]:
            if tag in inConf and not type(inConf[tag])==type([]):
                errAbort("Error in cellbrowser.conf: '%s' must be a list" % (tag))
        if tag=="visibility":
            if tag in inConf and inConf[tag] not in ["hide", "show"]:
                errAbort("Error in cellbrowser.conf: '%s' can only have values: 'hide' or 'show'" % (tag))

    if " " in inConf["name"]:
        errAbort("Sorry, please no whitespace in the dataset 'name' in the .conf file")

    # convertMeta also compares the sample IDs between meta and matrix
    # outMeta is a reordered & trimmed tsv version of the meta table
    sampleNames, needFilterMatrix, outMeta = convertMeta(inConf, outConf, datasetDir)

    geneToSym = None

    outMatrixFname = join(datasetDir, "exprMatrix.tsv.gz")
    geneToSym = -1 # None would mean "there are no gene symbols to map to"
    inMatrixFname = getAbsPath(inConf, "exprMatrix")
    doMatrix = matrixOrSamplesHaveChanged(datasetDir, inMatrixFname, outMatrixFname, outConf)

    if doMatrix:
        geneToSym = readGeneSymbols(inConf.get("geneIdType"), inMatrixFname)
        convertExprMatrix(inConf, outMatrixFname, outConf, sampleNames, geneToSym, datasetDir, needFilterMatrix)
        # in case script crashes after this, keep the current state of the config
        writeConfig(inConf, outConf, datasetDir)
    else:
        logging.info("Matrix and meta sample names have not changed, not indexing matrix again")

    coordFiles, clusterLabels = convertCoords(inConf, outConf, sampleNames, outMeta, datasetDir)

    foundConf = writeDatasetDesc(inConf["inDir"], outConf, datasetDir, coordFiles)
    if not foundConf:
        copyDatasetHtmls(inConf["inDir"], outConf, datasetDir)

    if geneToSym==-1:
        geneToSym = readGeneSymbols(inConf.get("geneIdType"), inMatrixFname)

    convertMarkers(inConf, outConf, geneToSym, clusterLabels, datasetDir)

    readQuickGenes(inConf, geneToSym, outConf)

    # need to able to see quickly how a dataset is described
    if ("descMd5s" in outConf) and ("datasetDesc" in outConf["descMd5s"]):
        outConf["hasFiles"] = ["datasetDesc"]

def writeAnndataCoords(anndata, fieldName, outDir, filePrefix, fullName, desc):
    " write embedding coordinates from anndata object to outDir, the new filename is <prefix>_coords.tsv "
    import pandas as pd
    fileBase = filePrefix+"_coords.tsv"
    fname = join(outDir, fileBase)

    existNames = getObsmKeys(anndata)
    altName1 = "X_"+fieldName
    altName2 = "X_draw_graph_"+fieldName
    if fieldName not in existNames:
        if altName1 in existNames:
            fieldName = altName1
        elif altName2 in existNames:
            fieldName = altName2
        else:
            logging.debug('Couldnt find coordinates for %s, tried keys %s, %s and %s' % (fullName, fieldName, altName1, altName2))
            return

    logging.info("Writing %s coords to %s" % (fullName, fname))
    fa2_coord=pd.DataFrame(anndata.obsm[fieldName],index=anndata.obs.index)
    fa2_coord.columns=['x','y']
    fa2_coord.to_csv(fname,sep='\t')
    desc.append( {'file':fileBase, 'shortLabel': fullName} )

def writeCellbrowserConf(name, coordsList, fname, addMarkers=True, args={}):
    for c in name:
        assert(c.isalnum() or c in ["-", "_"]) # only digits and letters are allowed in dataset names

    metaFname = args.get("meta", "meta.tsv")
    clusterField = args.get("clusterField", "Louvain Cluster")
    coordStr = json.dumps(coordsList, indent=4)
    matrixFname = args.get("exprMatrix", "exprMatrix.tsv.gz")

    conf = """
# This is a bare-bones, auto-generated cellbrowser config file.
# Look at https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf
# for a list of possible options
name='%(name)s'
shortLabel='%(name)s'
exprMatrix='%(matrixFname)s'
#tags = ["10x", 'smartseq2']
meta='%(metaFname)s'
geneIdType='symbols'
clusterField='%(clusterField)s'
labelField='%(clusterField)s'
enumFields=['%(clusterField)s']
coords=%(coordStr)s
#quickGenesFile='quickGenes.csv'
#alpha=0.3
#radius=2
""" % locals()

    if addMarkers:
        conf += 'markers = [{"file": "markers.tsv", "shortLabel":"Cluster Markers"}]\n'

    if "geneToSym" in args:
        conf += "geneToSym='%s'\n" % args["geneToSym"]

    if isfile(fname):
        logging.info("Not overwriting %s, file already exists." % fname)
        return

    ofh = open(fname, "w")
    ofh.write(conf)
    ofh.close()
    logging.info("Wrote %s" % ofh.name)

def anndataMatrixToTsv(ad, matFname, usePandas=False, useRaw=False):
    " write ad expression matrix to .tsv file and gzip it "
    import pandas as pd
    import scipy.sparse

    if useRaw and ad.raw is None:
        logging.warn("The option to export raw expression data is set, but the scanpy object has no 'raw' attribute. Exporting the processed scanpy matrix. Some genes may be missing.")

    if useRaw and ad.raw is not None:
        mat = ad.raw.X
        var = ad.raw.var
        logging.info("Processed matrix has size (%d cells, %d genes)" % (mat.shape[0], mat.shape[1]))
        logging.info("Using raw expression matrix")
    else:
        mat = ad.X
        var = ad.var

    rowCount, colCount = mat.shape
    logging.info("Writing scanpy matrix (%d cells, %d genes) to %s" % (rowCount, colCount, matFname))
    tmpFname = matFname+".tmp"

    logging.info("Transposing matrix") # necessary, as scanpy has the samples on the rows
    mat = mat.T

    if scipy.sparse.issparse(mat):
        logging.info("Converting csc matrix to row-sparse matrix")
        mat = mat.tocsr() # makes writing to a file ten times faster, thanks Alex Wolf!

    if usePandas:
        logging.info("Converting anndata to pandas dataframe")
        data_matrix=pd.DataFrame(mat, index=var.index.tolist(), columns=ad.obs.index.tolist())
        logging.info("Writing pandas dataframe to file (slow?)")
        data_matrix.to_csv(tmpFname, sep='\t', index=True)
    else:
        # manual writing row-by-row should save quite a bit of memory
        logging.info("Writing gene-by-gene, without using pandas")
        ofh = open(tmpFname, "w")
        sampleNames = ad.obs.index.tolist()
        ofh.write("gene\t")
        ofh.write("\t".join(sampleNames))
        ofh.write("\n")

        # when reading 10X files, read_h5 puts the geneIds into a separate field
        # and uses only the symbol. We prefer ENSGxxxx|<symbol> as the gene ID string
        if "gene_ids" in var:
            geneIdObj = var["gene_ids"]
            geneIdAndSyms = zip(geneIdObj.values, geneIdObj.index)
            genes = [x+"|"+y for (x,y) in geneIdAndSyms]
        elif "gene_symbols" in var:
            geneIdObj = var['gene_symbols']
            geneIdAndSyms = zip(geneIdObj.index, geneIdObj.values)
            genes = [x+"|"+y for (x,y) in geneIdAndSyms]
        else:
            genes = var.index.tolist()

        logging.info("Writing %d genes in total" % len(genes))
        for i, geneName in enumerate(genes):
            if i % 2000==0:
                logging.info("Wrote %d genes" % i)
            ofh.write(geneName)
            ofh.write("\t")
            if scipy.sparse.issparse(mat):
                logging.debug("Converting csr row to dense")
                row = mat.getrow(i).todense()
            else:
                row = mat[i,:]

            # Note: Faster is not possible. savetxt() is super slow, as it is not C code
            # We used to have this 'np.savetxt(ofh, row, fmt="%g", delimiter="\t", newline="\t")'
            # But tofile() is several times faster:
            row.tofile(ofh, sep="\t", format="%.7g")
            ofh.write("\n")

        ofh.close()

    if matFname.endswith(".gz"):
        runGzip(tmpFname, matFname)
    else:
        os.rename(tmpFname, matFname)

def makeDictDefaults(inVar, defaults):
    " convert inVar to dict if necessary, defaulting to our default labels "
    if type(inVar) is dict:
        return inVar

    d = {}
    for val in inVar:
        d[val] = defaults.get(val, val)
    return d

def scanpyToCellbrowser(adata, path, datasetName, metaFields=None, clusterField=None,
        nb_marker=50, doDebug=False, coordFields=None, skipMatrix=False, useRaw=False,
        markerField='rank_genes_groups'):
    """
    Mostly written by Lucas Seninge, lucas.seninge@etu.unistra.fr

    Given a scanpy object, write dataset to a dataset directory under path.

    This function export files needed for the ucsc cells viewer from the Scanpy Anndata object
    :param anndata: Scanpy AnnData object where information are stored
    :param path : Path to folder where to save data (tsv tables)
    :param clusterField: name of cluster field, used for labeling/default coloring, default is 'louvain'
    :param metaFields: list of metadata names (string) to export
    from the AnnData object (other than 'louvain' to also save (eg: batches, ...)).
    This can also be a dict of name -> label, if you want to have more human-readable names.
    :param nb_marker: number of cluster markers to store. Default: 50

    """
    if doDebug is not None:
        setDebug(doDebug)

    if not isdir(path):
        makeDir(path)

    confName = join(path, "cellbrowser.conf")
    if isfile(confName):
        logging.warn("%s already exists. Overwriting existing files." % confName)

    import numpy as np
    import pandas as pd
    import anndata

    if not skipMatrix:
        matFname = join(path, 'exprMatrix.tsv.gz')
        anndataMatrixToTsv(adata, matFname, useRaw=useRaw)

    if coordFields=="all" or coordFields is None:
        coordFields = coordLabels
    coordsFields = makeDictDefaults(coordFields, coordLabels)

    coordDescs = []

    for layoutCode, layoutName in coordFields.items():
        writeAnndataCoords(adata, layoutCode, path, layoutCode, layoutName, coordDescs)

    if len(coordDescs)==0:
        raise ValueError("No valid embeddings were found in anndata.obsm but at least one array of coordinates is required. Keys that were tried: %s" % (coordFields))

    ##Check for cluster markers
    if markerField not in adata.uns:
        logging.warn("Couldnt find list of cluster marker genes in the h5ad file in adata.uns with the key '%s'. "
        "As a result, your cell browser will not include marker genes. "
        "From Python, try running sc.tl.rank_genes_groups(adata) to "
        "create the cluster annotation." % markerField)
        addMarkers = False
    else:
        top_score=pd.DataFrame(adata.uns[markerField]['scores']).loc[:nb_marker]
        top_gene=pd.DataFrame(adata.uns[markerField]['names']).loc[:nb_marker]
        marker_df= pd.DataFrame()
        #for i in range(len(top_score.columns)):
        for clustName in top_score.columns:
            topScoreCol = top_score[[clustName]]
            topGeneCol = top_gene[[clustName]]
            concat=pd.concat([topScoreCol,topGeneCol],axis=1,ignore_index=True)
            concat['cluster_name']=clustName
            col=list(concat.columns)
            col[0],col[-2]='z_score','gene'
            concat.columns=col
            marker_df=marker_df.append(concat)

        #Rearranging columns -> Cluster, gene, score
        cols=marker_df.columns.tolist()
        cols=cols[::-1]
        marker_df=marker_df[cols]
        #Export
        fname = join(path, "markers.tsv")
        logging.info("Writing %s" % fname)
        pd.DataFrame.to_csv(marker_df,fname,sep='\t',index=False)
        addMarkers = True

    ##Save metadata
    if metaFields is None:
        metaFields = list(adata.obs.columns.values)
    else:
        # check that field names exist
        for name in metaFields:
            if name not in adata.obs.keys():
                logging.warn('There is no annotation field with the name `%s`.' % name)
                if name not in ["percent_mito", "n_genes", "n_counts"]:
                    # tolerate missing default fields
                    raise ValueError()

    metaFields = makeDictDefaults(metaFields, metaLabels)

    meta_df=pd.DataFrame()

    for metaKey, metaLabel in iterItems(metaFields):
        logging.debug("getting meta field: %s -> %s" % (metaKey, metaLabel))
        if metaKey not in adata.obs:
            logging.warn(str(metaKey) + ' field is not present in the AnnData.obs object')
        else:
            temp=adata.obs[[metaKey]]
            #if metaKey in metaLabels:
            #logging.debug("Using new name %s -> %s" % (metaKey, metaLabels[metaKey]))
            #temp.name = metaLabel
            meta_df=pd.concat([meta_df,temp],axis=1)

    meta_df.rename(metaFields, axis=1, inplace=True)
    fname = join(path, "meta.tsv")
    meta_df.to_csv(fname,sep='\t')

    if clusterField is None:
        clusterField = 'louvain'

    argDict = {}
    if clusterField:
        clusterFieldLabel = metaFields.get(clusterField, clusterField)
        argDict['clusterField'] = clusterFieldLabel

    if addMarkers:
        generateQuickGenes(path)
        argDict['quickGenes'] = "quickGenes.tsv"

    if isfile(confName):
        logging.info("%s already exists, not overwriting. Remove and re-run command to recreate." % confName)
    else:
        writeCellbrowserConf(datasetName, coordDescs, confName, addMarkers=addMarkers, args=argDict)

def writeJson(data, outFname):
    """ https://stackoverflow.com/a/37795053/233871 """
    # Write JSON file
    tmpName = outFname+".tmp"
    with io.open(tmpName, 'w', encoding='utf8') as outfile:
        #str_ = json.dumps(data, indent=2, sort_keys=True,separators=(',', ': '), ensure_ascii=False)
        if isPy3:
            str_ = json.dumps(data, indent=2, separators=(',', ': '), ensure_ascii=False)
        else:
            str_ = json.dumps(data, indent=2, separators=(',', ': '), ensure_ascii=False, encoding="utf8")
        outfile.write(str_)
    os.rename(tmpName, outFname)
    logging.info("Wrote %s" % outFname)

def writePyConf(idfInfo, descFname):
    " write python-style conf file "
    ofh = open(descFname, "w")
    for key, val in iterItems(idfInfo):
        ofh.write("%s = %s\n" % (key, repr(val)))
    ofh.close()
    logging.info("Wrote %s" % descFname)

def writeConfig(inConf, outConf, datasetDir):
    " write dataset summary info to json file. Also keep a copy of the input config. "
    # keep a copy of the original config in the output directory for debugging later
    confName = join(datasetDir, "cellbrowser.json.bak")
    writeJson(inConf, confName)
    logging.info("Wrote %s" % confName)

    outConfFname = join(datasetDir, "dataset.json")
    writeJson(outConf, outConfFname)
    logging.info("Wrote %s" % outConfFname)

def startHttpServer(outDir, port):
    " start an http server on localhost serving outDir on a given port "
    port = int(port)

    import RangeHTTPServer
    try:
        # py3
        import http.server as SimpleHTTPServer
        from http.server import HTTPServer
    except:
        # py2
        import SimpleHTTPServer
        from BaseHTTPServer import HTTPServer

    #server_address = ('localhost', port) # use this line to allow only access from localhost
    server_address = ('', port) # by default, we allow access from anywhere
    HandlerClass = RangeHTTPServer.RangeRequestHandler
    HandlerClass.protocol_version = "HTTP/1.0"
    httpd = HTTPServer(server_address, HandlerClass)

    sa = httpd.socket.getsockname()
    ipAddr = sa[0]
    if ipAddr=="0.0.0.0":
        ipAddr = "127.0.0.1"

    outDir = os.path.expanduser(outDir)
    os.chdir(outDir)
    print("Serving "+outDir+" on port %s" % str(port))
    print("Point your internet browser to http://"+ipAddr+":"+str(sa[1])+" (or the IP address of this server)")
    sys.stderr = open("/dev/null", "w") # don't show http status message on console
    httpd.serve_forever()

def cbBuild(confFnames, outDir, port=None):
    " stay compatible with old name "
    build(confFnames, outDir, port)

def build(confFnames, outDir, port=None, doDebug=False, devMode=False):
    " build browser from config files confFnames into directory outDir and serve on port "
    if outDir=="" or outDir==None:
        outDir = defOutDir

    setDebug(doDebug)
    if type(confFnames)==type(""):
        # it's very easy to forget that the input should be a list so we tolerate a single string
        confFnames = [confFnames]

    outDir = expanduser(outDir)
    datasets = []
    for inConfFname in confFnames:
        if isdir(inConfFname):
            inConfFname = join(inConfFname, "cellbrowser.conf")
        inConf = loadConfig(inConfFname)

        dsName = inConf["name"]
        datasets.append(dsName)

        datasetDir = join(outDir, dsName)
        makeDir(datasetDir)

        outConfFname = join(outDir, "dataset.conf")
        #if onlyMeta:
            #outConf = json.parse(open(outConfFname)) # reuse the old config
        #else:
        outConf = OrderedDict()

        convertDataset(inConf, outConf, datasetDir)

        outConf["fileVersions"]["conf"] = getFileVersion(abspath(inConfFname))
        outConf["md5"] = calcMd5ForDataset(outConf)

        writeConfig(inConf, outConf, datasetDir)

    cbUpgrade(outDir, datasets)

    if port:
        print("Interrupt this process, e.g. with Ctrl-C, to stop the webserver")
        startHttpServer(outDir, int(port))

pidFname = "cellbrowser.pid"

def savePid():
    " save PID to file "
    import tempfile
    tempDir = tempfile.gettempdir()
    tempName = join(tempDir, pidFname)
    tempFh = open(tempName, "w")
    pid = os.getpid()
    tempFh.write("%d" % pid)
    tempFh.close()
    return pid

def loadPid():
    " return previously saved PID and delete the pid file "
    import tempfile
    tempDir = tempfile.gettempdir()
    tempName = join(tempDir, pidFname)
    if not isfile(tempName):
        return None
    pid = int(open(tempName).read())
    os.remove(tempName)
    return pid

def serve(outDir, port):
    " forking from R/ipython/jupyter is not a good idea at all. Instead, we start a shell that runs Python. "
    if outDir is None:
        raise Exception("html outDir must be set if a port is set")

    port = int(port)

    stop()
    cmd = [sys.executable, __file__, "cbServe", "'"+outDir+"'", str(port)]
    cmdStr = " ".join(cmd) + "&"
    logging.debug("Running command through shell: '%s'" % cmdStr)
    ret = os.system(cmdStr)
    if ret!=0:
        print("ERROR: Could not run command '%s'" % cmdStr)

def serveDirect(outDir, port):
    """ run a webserver on localhost with a port, fork/detach as a daemon, and
    write PID into <tempdir>/cellbrowser.pid. NOT USED FOR NOW - Jupyther/ipython/Rstudio don't like this  """
    # mostly copied from
    # http://code.activestate.com/recipes/66012-fork-a-daemon-process-on-unix/
    # do the UNIX double-fork magic, see Stevens' "Advanced 
    # Programming in the UNIX Environment" for details (ISBN 0201563177)
    try:
        pid = os.fork()
        if pid > 0:
            # exit first parent
            #sys.exit(0)
            return
    except OSError as e:
        logging.error("fork #1 failed: %d (%s)" % (e.errno, e.strerror))
        sys.exit(1)

    # decouple from parent environment
    os.chdir("/")
    os.setsid()
    os.umask(0)

    # do second fork
    try:
        pid = os.fork()
        if pid > 0:
            # exit from second parent
            #logging.error( "Daemon PID %d" % pid )
            sys.exit(0)
    except OSError as e:
        logging.error("fork #2 failed: %d (%s)" % (e.errno, e.strerror) )
        sys.exit(1)

    # start the daemon main loop
    pid = savePid()
    logging.info("Starting webserver, process ID: %d" % pid)
    startHttpServer(outDir, port)

def stop():
    " kill the process with PID in <tempdir>/cellbrowser.pid "
    import signal
    pid = loadPid()
    if pid is None:
        logging.debug("No saved process ID found, nothing to kill")
        return

    try:
        os.kill(pid, signal.SIGTERM) #or signal.SIGKILL 
    except OSError as ex:
        print("Unable to kill cellbrowser process with PID %d" % pid)
        return

    print("Successfully killed cellbrowser process with PID %d" % pid)

def cbBuildCli():
    " command line interface for dataset converter, also copies the html/js/etc files "
    args, options = cbBuild_parseArgs()

    confFnames = options.inConf
    if confFnames==None:
        confFnames = ["cellbrowser.conf"]

    if options.init:
        copyPkgFile("sampleConfig/cellbrowser.conf")
        copyPkgFile("sampleConfig/desc.conf")
        sys.exit(1)

    for fname in confFnames:
        if not isfile(fname):
            logging.error("File %s does not exist." % fname)
            cbBuild_parseArgs(showHelp=True)
    if options.outDir is None:
        logging.error("You have to specify at least the output directory via -o or set the env. variable CBOUT or set htmlDir in ~/.cellbrowser.conf.")
        cbBuild_parseArgs(showHelp=True)

    outDir = options.outDir
    #onlyMeta = options.onlyMeta
    port = options.port

    build(confFnames, outDir, port)

def readMatrixAnndata(matrixFname, samplesOnRows=False, genome="hg38"):
    " read an expression matrix and return an adata object. Supports .mtx, .h5 and .tsv (not .tsv.gz) "
    import scanpy as sc
    if matrixFname.endswith(".mtx"):
        import pandas as pd
        logging.info("Loading expression matrix: mtx format")
        adata = sc.read(matrixFname, cache=False).T

        mtxDir = dirname(matrixFname)
        adata.var_names = pd.read_csv(join(mtxDir, 'genes.tsv'), header=None, sep='\t')[1]
        adata.obs_names = pd.read_csv(join(mtxDir, 'barcodes.tsv'), header=None)[0]

    elif matrixFname.endswith(".h5"):
        import h5py
        ifh = h5py.File(matrixFname,'r')
        groups = list(ifh.keys())
        ifh.close()

        if genome not in groups:
            errAbort("The file %s does not have expression info for the genome %s. Possible genomes are: %s. "
                     "Choose one of these and and specify it with the option -g" %
                     (matrixFname, genome, groups))

        logging.info("Loading expression matrix: 10X h5 format")
        adata = sc.read_10x_h5(matrixFname, genome=genome)

    else:
        logging.info("Loading expression matrix: scanpy-supported format, like h5ad, loom, tab-separated, etc.")
        adata = sc.read(matrixFname, first_column_names=True)
        if not samplesOnRows:
            logging.debug("Scanpy defaults to samples on lines, so transposing the expression matrix")
            adata = adata.T

    return adata

def addMd5(d, fname, keyName="md5", isSmall=False, shortMd5=True):
    " add a key 'md5' to dict d with first MD5LEN letters of fname "
    logging.debug("Getting md5 of %s" % fname)
    md5 = md5ForFile(fname, isSmall=isSmall)[:MD5LEN]
    if shortMd5:
        md5 = md5[:MD5LEN]
    d[keyName] = md5

def addCollections(collDir, datasets):
    """ get all names of collections, find all their cellbrowser.conf files under collDir,
    append to 'datasets', mark datasets with collections as hidden and return new list of datasets
    """
    # make map from dataset name to dataset info dict
    nameToDs = {}
    for ds in datasets:
        assert(ds["name"] not in nameToDs)
        nameToDs[ds["name"]] = ds

    # make a map from collName -> dataset info dicts
    # and dict with collection name -> list of dataset md5s
    collToDatasets = defaultdict(list)
    dsMd5s = defaultdict(list)
    for ds in datasets:
        for collName in ds.get("collections", []):
            collToDatasets[collName].append(nameToDs[ds["name"]])
            dsMd5s[collName].append(ds["md5"])
            # a dataset that is part of a collection is hidden from the main list by default
            if not "visibility" in ds:
                logging.debug("Setting dataset %s to hide as it's part of collection %s" % (ds["name"], collName))
                ds["visibility"] = "hide"

    if collDir is None:
        if len(collToDatasets)>0:
            dsAssignment = {}
            for collName, dsList in iterItems(collToDatasets):
                dsAssignment[collName] = []
                for ds in dsList:
                    dsAssignment[collName].append(ds["name"])
            errAbort("Some datasets defined collections but collDir is not defined in ~/.cellbrowser.conf. Found collection-dataset assignment: %s" % dsAssignment)
        else:
            logging.debug("collDir not defined in ~/.cellbrowser.conf, but also no collections used")
            return datasets, {}
    else:
        logging.debug("collDir is set to %s" % collDir)

    logging.debug("Scanning %s for collections" % collDir)
    collDir = expanduser(collDir)
    collLabels = {}
    for collName, collDatasets in iterItems(collToDatasets):
        datasetNames = [ds["name"] for ds in collDatasets]
        subPath = join(collDir, collName)
        if not isdir(subPath):
            errAbort("%s is not a directory but datasets %s refer to a collection of this name" % (subPath, datasetNames))
        fname = join(subPath, "cellbrowser.conf")
        if not isfile(fname):
            #logging.debug("Cannot find %s" % fname)
            errAbort("%s does not exist but datasets %s refer to a collection of this name" % (fname, datasetNames))

        conf = loadConfig(fname, requireTags=["name", "shortLabel"])

        conf["name"] = collName
        conf["isCollection"] = True
        conf["datasetCount"] = len(collToDatasets[collName])
        if isfile(join(subPath, "desc.conf")) or isfile(join(subPath, "datasetDesc.conf")):
            conf["hasFiles"] = ["datasetDesc"]

        if "collections" in conf:
            errAbort("File %s defines a collection but is itself a collection. This is not allowed." % subPath)

        fileVersions = {}
        fileVersions["config"] = getFileVersion(fname)
        conf["fileVersions"] = fileVersions

        collLabels[collName] = conf["shortLabel"]
        #conf["md5"] = collMd5s[collName]

        datasets.append(conf)

    # make map from collName -> list of only primary datasets (not others)
    primCollDatasets = defaultdict(list)
    for ds in datasets:
        if "collections" in ds:
            primCollName = ds["collections"][0]
            primCollDatasets[primCollName].append(ds)

    # add info about siblings to all datasets with the same primary collection
    # it's easier to have all linked dataset labels in the dataset than having to load these all the time
    #dsSiblingLabels = defaultdict(list)
    #for collName, dsList in iterItems(primCollDatasets):
            #sibInfo = []
            #for ds in dsList:
                #sibInfo.append({"name":ds["name"], "shortLabel":ds["shortLabel"]})
            #for ds in dsList:
                #ds["sameCollection"] = sibInfo

    return datasets, primCollDatasets

def calcMd5ForDataset(datasetDesc):
    " make a combined md5 for the dataset, based on the config file, meta, expression and the html or dataset descs "
    md5s = []
    if "fileVersions" in datasetDesc:
        fileVers = datasetDesc["fileVersions"]
        for key, fileInfo in iterItems(fileVers):
            md5s.append(fileInfo["md5"])

    if "descMd5s" in datasetDesc:
        descs = datasetDesc["descMd5s"]
        md5s.append (descs.get("summary", ""))
        md5s.append (descs.get("methods", ""))
        md5s.append (descs.get("downloads", ""))

    return md5ForList(md5s)[:MD5LEN]

def findDatasets(outDir):
    """ search all subdirs of outDir for dataset.json files and return their
    contents as a list A dataset description is a list with three members: A
    label, the base URL and a longer description that can contain html.
    The attribute "priority" can be used to enforce an order on the datasets
    """
    datasets = []
    dsNames = defaultdict(list)
    for subDir in os.listdir(outDir):
        if not isdir(join(outDir, subDir)):
            continue
        if subDir.endswith(".skip"):
            continue
        fname = join(outDir, subDir, "dataset.json")
        if not isfile(fname):
            continue

        datasetDesc = readJson(fname, keepOrder=True)

        if "isCollection" in datasetDesc:
            logging.debug("Dataset %s is a collection, so not parsing now" % datasetDesc["name"])
            continue
        #datasetDesc = json.load(open(fname))
        customdecoder = json.JSONDecoder(object_pairs_hook=OrderedDict)
        inStr = open(fname).read()
        if not isPy3:
            inStr = inStr.decode("utf8")
        inStr = readFile(fname)

        datasetDesc = customdecoder.decode(inStr)

        if not "md5" in datasetDesc:
            datasetDesc["md5"] = calcMd5ForDataset(datasetDesc)

        if not "name" in datasetDesc: # every dataset has to have a name
            errAbort("Dataset %s must have a 'name' field." % datasetDesc)

        dsName = datasetDesc["name"]
        if dsName in dsNames:
            errAbort("Duplicate name: %s appears in these directories: %s and %s" % \
                  (dsName, dsNames[dsName], subDir))
        dsNames[dsName].append(subDir)

        #assert("shortLabel" in datasetDesc)
        if not "shortLabel" in datasetDesc:
            datasetDesc["shortLabel"] = datasetDesc["name"]

        datasetDesc["baseUrl"] = subDir+"/"
        datasets.append(datasetDesc)
    datasets = list(sorted(datasets, key=lambda k: k.get('priority', 10)))
    logging.info("Found %d datasets" % len(datasets))
    return datasets

def copyAllFiles(fromDir, subDir, toDir):
    " copy all files in fromDir/subDir to toDir/subDir "
    outDir = join(toDir, subDir)
    makeDir(outDir)
    for filename in glob.glob(join(fromDir, subDir, '*')):
    # egg-support commented out for now
    #for filename in pkg_resources.resource_listdir(__name__, join(fromDir, subDir)):
        if isdir(filename):
            continue
        #fullPath = join(fromDir, subDir, filename)
        fullPath = filename
        logging.debug("Copying %s to %s" % (fullPath, outDir))
        # copy uses chmod() which we don't want
        #shutil.copy(filename, outDir)
        dstPath = join(outDir, basename(filename))
        shutil.copyfile(fullPath, dstPath)
        #s = pkg_resources.resource_string(__name__, filename)
        #outFname = join(outDir, filename)
        #ofh = open(outFname, "wb")
        #ofh.write(s)
        #ofh.close()

def writeCollectionFiles(collDir, datasets, collToDatasets, onlyDatasets, outDir):
    """ copy over the desc.conf files for all collections and generate md5s for them.
    re-write all collection children, too.
    """
    if len(collToDatasets.keys())==0:
        logging.debug("No collections used")
        return

    if collDir is None:
        logging.warn("Collections are mentioned in at least one cellbrowser.conf file but 'collDir' is not defined in ~/.cellbrowser.conf. Not adding collections.")
        logging.warn("Collections found: %s" % collToDatasets)
        return

    collDir = expanduser(collDir)

    if not isdir(collDir):
        logging.info("directory %s not found. To use collections, please read the documentation on how to set them up." % collDir)
        return

    # make map from dataset name to dataset info dict
    nameToDs = {}
    for ds in datasets:
        #if (ds["name"] in nameToDs):
            #errAbort("dataset name %s already used by %s" % (ds["name"], nameToDs[ds["name"]]))
        nameToDs[ds["name"]] = ds

    for collName, datasets in iterItems(collToDatasets):
        inDir = join(collDir, collName)

        collOutDir = join(outDir, collName)
        if not isdir(collOutDir):
            makeDir(collOutDir)

        # make md5 for a collection: based on all datasets and the description files
        collMd5s = []
        collInfo = OrderedDict() # description of collection, goes into .json file at the end

        for key, val in iterItems(nameToDs[collName]):
            collInfo[key] = val

        foundConf = writeDatasetDesc(inDir, collInfo, collOutDir)
        if foundConf:
            collMd5s.append(collInfo["descMd5s"]["datasetDesc"])

        summFname = join(inDir, "summary.html")
        if isfile(summFname):
            logging.debug("Copying %s to %s" % (summFname, collOutDir))
            #shutil.copy(summFname, collOutDir) # copy() uses chmod
            shutil.copyfile(summFname, join(collOutDir, basename(summFname)))
            collMd5s.append(md5ForFile(summFname))

        for ds in datasets:
            collMd5s.append( ds["md5"] )

        collMd5 = md5ForList(collMd5s)[:MD5LEN]

        summDsList = summarizeDatasets(datasets)
        collInfo["datasets"] = summDsList
        if collName not in nameToDs:
            errAbort("collection %s is used by datasets %s, but there was no cellbrowser.conf found under %s" %
                    (collName, [ds["name"] for ds in datasets], inDir))

        # copy the full list of linked datasets into every child dataset
        # and write these JSON files again
        linkedDsList = summarizeDatasets(summDsList, minimal=True)

        for dsInfo in summDsList:
            dsName = dsInfo["name"]
            childInfo = nameToDs[dsName]
            childInfo["datasets"] = linkedDsList
            childFname = join(outDir, dsName, "dataset.json")
            writeJson(childInfo, childFname)

        collInfo["md5"] = collMd5
        nameToDs[collName]["md5"] = collMd5 # update the md5 in the big list, for the main summary

        jsonFname = join(collOutDir, "dataset.json")
        writeJson(collInfo, jsonFname)

def copyStatic(baseDir, outDir):
    " copy all js, css and img files to outDir "
    logging.info("Copying js, css and img files to %s" % outDir)
    imgDir = join(outDir, "img")

    copyAllFiles(baseDir, "ext/images", outDir)
    copyAllFiles(baseDir, "img", outDir)
    copyAllFiles(baseDir, "ext", outDir)
    copyAllFiles(baseDir, "js", outDir)
    copyAllFiles(baseDir, "css", outDir)

def writeVersionedLink(ofh, mask, webDir, relFname, addVersion=True):
    " write sprintf-formatted mask to ofh, but add ?md5 to jsFname first. Goal is to force cache reload in browser. "
    # hack for jquery - avoid jquery button overriding any other button function
    if relFname.endswith("jquery-ui.min.js"):
        ofh.write("""<script>
  $.fn.bsButton = $.fn.button.noConflict();
  $.fn.bsTooltip = $.fn.tooltip.noConflict();
</script>
""")

    filePath = join(webDir, relFname)
    md5 = md5WithPython(filePath)
    if addVersion:
        verFname = relFname+"?"+md5[:MD5LEN]
    else:
        verFname = relFname
    outLine = mask % verFname
    ofh.write(outLine+"\n")

def writeGaScript(ofh, gaTag):
    " write Google analytics script to ofh "
    ofh.write("""
<!-- Global site tag (gtag.js) - Google Analytics -->
<script async src="https://www.googletagmanager.com/gtag/js?id=%s"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', '%s');
</script>
<!-- END - Google Analytics -->

""" % (gaTag, gaTag))

def summarizeDatasets(datasets, minimal=False):
    " keep only the most important fields of a list of datasets and return them as a list of dicts "
    dsList = []
    for ds in datasets:
        summDs = {
                "shortLabel" : ds["shortLabel"],
                "name" : ds["name"],
                "md5" : ds["md5"],
                }

        if not minimal:
            if "sampleCount" in ds:
                summDs["sampleCount"] = ds["sampleCount"]
            else:
                summDs["isCollection"] = True
                summDs["datasetCount"] = ds["datasetCount"]

            for optListTag in ["tags", "collections", "hasFiles"]:
                if optListTag in ds:
                    assert(type(ds[optListTag])==type([])) # has to be a list
                    summDs[optListTag] = ds[optListTag]

            #for optTag in ["visibility"]:
                #if optTag in ds:
                    #summDs[optTag] = ds[optTag]

        dsList.append(summDs)
    return dsList

def makeIndexHtml(baseDir, datasets, outDir, devMode=False):
    " make the index.html, copy over all .js and related files and add their md5s "
    logging.debug("Development mode is: %s" % repr(devMode))

    dsList = summarizeDatasets(datasets)
    indexFname = join(baseDir, "html", "index.html")
    datasetListJs = "var datasets = "+json.dumps(dsList, sort_keys=True, indent=4, separators=(',', ': '))+";"

    newFname = join(outDir, "index.html")
    tmpFname = newFname+".tmp"
    ofh = open(tmpFname, "w")

    ofh.write("<!doctype html>\n")
    ofh.write('<html>\n')
    ofh.write('<head>\n')
    ofh.write('<meta charset="utf-8">\n')
    ofh.write('<title>UCSC Cell Browser</title>\n')

    cssFnames = ["ext/jquery-ui-1.12.1.css", "ext/spectrum-1.8.0.css", "ext/jquery.contextMenu.css",
        "ext/slick.grid.css", "ext/slick.examples.css",
        #"ext/jquery-ui-1.11.3.custom.css",
        "ext/jquery.tipsy.1.0.3.min.css", "ext/bootstrap.min.css",
        "ext/introjs.2.4.0.min.css", "ext/bootstrap-submenu.min.css",
        "ext/bootstrap-dropmenu.min.css", "ext/font-awesome.css",
        "ext/googleMaterialIcons.css", "ext/chosen.1.8.2.min.css",
        "ext/select2.4.0.4.min.css",
        "ext/selectize.bootstrap3.css",
        #"ext/selectize.0.12.4.min.css",
        "ext/OverlayScrollbars.min.css", # 1.6.2, from https://cdnjs.com/libraries/overlayscrollbars
        "css/cellBrowser.css"
        ]

    addVersion = not devMode

    for cssFname in cssFnames:
        writeVersionedLink(ofh, '<link rel="stylesheet" href="%s">', baseDir, cssFname, addVersion=addVersion)

    jsFnames = ["ext/FileSaver.1.1.20151003.min.js", "ext/jquery.3.1.1.min.js",
        "ext/palette.js", "ext/spectrum.min.js", "ext/jsurl2.js",
        "ext/chosen.jquery.min.js", "ext/mousetrap.min.js",
        "ext/jquery.contextMenu.js", "ext/jquery.ui.position.min.js",
        "ext/jquery.tipsy.min.js", "ext/intro.min.js", "ext/papaparse.min.js",
        "ext/bootstrap.min.js", "ext/bootstrap-submenu.js", "ext/pako_inflate.min.js",
        "ext/FastBitSet.js", "ext/hamster.js", "ext/split.js", "ext/normalizeWheel.js",
        "ext/tablesort.js", "ext/tablesort.number.min.js", "ext/Chart.bundle.min.js",
        "ext/chartjs-chart-box-and-violin-plot.js", "ext/jquery-ui.min.js",
        "ext/select2.min.js",
        "ext/selectize.js", # 0.12.16
        "ext/jquery.sparkline.min.js",
        "ext/jquery.overlayScrollbars.min.js", # 1.6.2 from https://cdnjs.com/libraries/overlayscrollbars
        "ext/jquery.event.drag-2.3.0.js", # for slickgrid 2.4.5
        "ext/lz-string.js",  # 1.4.4, https://raw.githubusercontent.com/pieroxy/lz-string/master/libs/lz-string.js
        "ext/slick.core.js",
        "ext/slick.cellrangedecorator.js", "ext/slick.cellrangeselector.js", "ext/slick.cellselectionmodel.js",
        "ext/slick.editors.js", "ext/slick.formatters.js", "ext/slick.grid.js",
        "ext/tiny-queue.js", "ext/science.v1.js", "ext/reorder.v1.js",  # commit d51dda9ad5cfb987b9e7f2d7bd81bb9bbea82dfe
        "js/cellBrowser.js", "js/cbData.js", "js/maxPlot.js", "js/maxHeat.js",
        ]

    # at UCSC, for grant reports, we need to get some idea how many people are using the cell browser
    ofh.write('<script async defer src="https://cells.ucsc.edu/js/cbTrackUsage.js"></script>\n')

    for jsFname in jsFnames:
        writeVersionedLink(ofh, '<script src="%s"></script>', baseDir, jsFname, addVersion=addVersion)

    if getConfig("gaTag") is not None:
        gaTag = getConfig("gaTag")
        writeGaScript(ofh, gaTag)

    ofh.write('</head>\n')
    ofh.write('<body>\n')
    ofh.write('<script>\n')
    ofh.write(datasetListJs)
    #ofh.write(colListJs)
    ofh.write('\n')
    ofh.write('cellbrowser.loadData(datasets);\n')
    ofh.write('</script>\n');
    ofh.write('</body>\n')
    ofh.write('</html>\n')

    ofh.close()
    os.rename(tmpFname, newFname)

    datasetLabels = [x["name"] for x in dsList]
    logging.info("Wrote %s, added datasets: %s" % (newFname, " - ".join(datasetLabels)))

def removeHiddenDatasets(datasets):
    """ visibility="hidden" removes datasets from the list, e.g. during paper review """
    newDsList = []
    for ds in datasets:
        if ds.get("visibility")=="hide":
            logging.debug("Dataset %s is set to hide, skipping" % ds["name"])
            continue
        newDsList.append(ds)
    return newDsList

def cbMake(outDir, onlyDatasets=None, devMode=False):
    cbUpgrade(outDir, onlyDatasets, devMode)

def cbUpgrade(outDir, onlyDatasets=None, devMode=False):
    " create index.html in outDir and copy over all other static files "
    baseDir = dirname(__file__) # = directory of this script
    webDir = join(baseDir, "cbWeb")
    copyStatic(webDir, outDir)
    datasets = findDatasets(outDir)

    collDir = getConfig("collDir")
    for ds in datasets:
        print(ds["name"])
    datasets, primCollDatasets = addCollections(collDir, datasets)

    writeCollectionFiles(collDir, datasets, primCollDatasets, onlyDatasets, outDir)

    datasets = removeHiddenDatasets(datasets)

    makeIndexHtml(webDir, datasets, outDir, devMode=devMode)


def cbUpgradeCli():
    " command line interface for copying over the html and js files and recreate index.html "
    args, options = cbUpgrade_parseArgs()
    outDir = options.outDir

    if outDir is None:
        errAbort("You have to specify at least the output directory or set the environment variable CBOUT.")
    if len(args)!=0:
        errAbort("This command does not accept arguments without options. Did you mean: -o <outDir> ? ")

    cbUpgrade(outDir, devMode=options.devMode)

def parseGeneLocs(geneType):
    """
    return dict with geneId -> list of bedRows
    bedRows have (chrom, start, end, geneId, score, strand)
    """
    fname = getStaticFile(join("genes", geneType+".genes.bed"))
    ret = defaultdict(list)
    for line in open(fname):
        row = line.rstrip("\r\n").split('\t')
        name = row[3].split(".")[0]
        ret[name].append(row)
    return ret

def extractMatrix(inMatrixFname, hubMatrixFname):
    if isfile(hubMatrixFname):
        logging.info("Not extracting to %s, file already exists" % hubMatrixFname)
    else:
        logging.info("Extracting matrix to %s" % hubMatrixFname)
        cmd = "gunzip -c %s > %s" % (inMatrixFname, hubMatrixFname)
        runCommand(cmd)

def getSizesFname(genome):
    " return chrom.sizes filename for db "
    fname = getStaticFile(join("genomes", genome+".sizes"))
    assert(isfile(fname))
    return fname

def pipeLog(msg):
    logging.info(msg)

# for iphython debug command line
def excepthook(type, value, traceback):
    from IPython import embed
    embed()

def checkLayouts(conf):
    """ it's very easy to get the layout names wrong: check them and handle the special value 'all' """
    if "doLayouts" not in conf:
        return ["tsne", "umap"]

    doLayouts = conf["doLayouts"]
    if doLayouts=="all":
        doLayouts = recommendedLayouts

    for l in doLayouts:
        if l not in coordLabels:
            errAbort("layout name %s is not valid. Valid layouts: %s" % (l, str(coordLabels)))

    return doLayouts

class Tee(object):
    " file like object that writes to stdout and also to a file "
    def __init__(self, fname):
        self.terminal = sys.stdout
        self.log = open(fname, "a")
        self.fname = fname

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass

def addMetaToAnnData(adata, fname):
    " parse .csv or .tsv and add the meta data as fields to adata.obs "
    logging.info("Adding meta data from %s to anndata object" % fname)
    import pandas as pd
    df1 = adata.obs
    df2 = pd.read_csv(fname, sep="\t", index_col=0)

    ids1 = set(df1.index)
    ids2 = set(df2.index)
    commonIds = ids1.intersection(ids2)
    logging.info("Meta data from %s has %d cell identifiers in common with anndata" % (fname, len(commonIds)))
    if len(commonIds)==0:
        errAbort("Vales in first column in file %s does not seem to match the cell IDs from the expression matrix" % fname)

    df3 = df1.join(df2, how="left")
    logging.debug("list of column names in merged meta data: %s"% ",".join(list(df3.columns)))
    adata.obs = df3

    #if len(adata.obs)!=len(df):
        #errAbort("The number of cells in the expression matrix does not match the number of lines in '%s'" % fname)

    # re-order meta data to be in the same order as the matrix
    #df = df.reindex(adata.obs.index)

    #for colName in list(df.columns):
        #logging.debug("Adding column %s to adata.obs" % colName)
        #adata.obs[colName] = df[colName].astype("category")
    #adata.obs["Cluster"] = df["Cluster"].astype("category")
    #adata.obs["Cluster"] = pd.Categorical(df["Cluster"], ordered=True)
    return adata

def getObsmKeys(adata):
    "get the keys of the obsm object. Has this changed in newer versions of anndata? "
    try:
        obsmKeys = adata.obsm.dtype.names # this used to work
    except:
        obsmKeys = list(adata.obsm.keys()) # this seems to work with newer versions
    return obsmKeys

def cbScanpy(matrixFname, inMeta, inCluster, confFname, figDir, logFname):
    """ run expr matrix through scanpy, output a cellbrowser.conf, a matrix and the meta data.
    Return an adata object. Optionally keeps a copy of the raw matrix in adata.raw """
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import warnings
    warnings.filterwarnings("ignore")

    conf = maybeLoadConfig(confFname)

    doLayouts = checkLayouts(conf)

    sc.settings.set_figure_params(dpi=200)
    sc.settings.file_format_figs = 'png'
    sc.settings.plot_suffix=''
    sc.settings.autosave=True
    sc.settings.autoshow=False
    sc.settings.figdir=figDir
    #sc.settings.logfile=logFname
    #sc.settings.logDir=join(outDir, "cbScanpy.log")

    sys.stdout = Tee(logFname) # we want our log messages and also scanpy messages into one file

    pipeLog("cbScanpy $Id$")
    pipeLog("Input file: %s" % matrixFname)
    #printLog("Output directory: %s" % outDir)
    pipeLog("Start time: %s" % datetime.datetime.now())
    sc.logging.print_versions()

    start = timeit.default_timer()
    adata = readMatrixAnndata(matrixFname, samplesOnRows=options.samplesOnRows, genome=options.genome)

    # set default values
    conf["useRaw"] = conf.get("useRaw", True)
    conf["doExp"] = conf.get("doExp", False)
    conf["doTrimCells"] = conf.get("doTrimCells", True)
    conf["minGenes"] = conf.get("minGenes", 200)
    conf["doTrimGenes"] = conf.get("doTrimGenes", True)
    conf["minCells"] = conf.get("minCells", 3)
    conf["doFilterMito"] = conf.get("doFilterMito", True)
    conf["geneIdType"] = conf.get("geneIdType", "auto")
    conf["mitoMax"] = conf.get("mitoMax", 0.05)
    conf["doFilterGenes"] = conf.get("doFilterGenes", True)
    conf["filterMaxGenes"] = conf.get("filterMaxGenes", 15000)
    conf["filterMinGenes"] = conf.get("filterMinGenes", 10)
    conf["doNormalize"] = conf.get("doNormalize", True)
    conf["countsPerCell"] = conf.get("countsPerCell", 10000)
    conf["doLog"] = conf.get("doLog", True)
    conf["doTrimVarGenes"] = conf.get("doTrimVarGenes", True)
    conf["varMinMean"] = conf.get("varMinMean", 0.0125)
    conf["varMaxMean"] = conf.get("varMaxMean", 3)
    conf["varMinDisp"] = conf.get("varMinDisp", 0.5)
    conf["doRegress"] = conf.get("doRegress", True)
    conf["regressMax"] = conf.get("regressMax", 10)
    conf["pcCount"] = conf.get("pcCount", "auto")
    conf["doLayouts"] = doLayouts
    conf["doLouvain"] = conf.get("doLouvain", True)
    conf["louvainNeighbors"] = int(conf.get("louvainNeighbors", 6))
    conf["louvainRes"] = float(conf.get("louvainRes", 1.0))
    conf["doMarkers"] = conf.get("doMarkers", True)
    conf["markerCount"] = int(conf.get("markerCount", 20))
    conf["inMeta"] = conf.get("inMeta", inMeta)
    conf["inCluster"] = conf.get("inCluster", inCluster)

    if conf["inMeta"]:
        adata = addMetaToAnnData(adata, inMeta)

    sampleCount = len(adata.obs)

    bigDataset = False
    if sampleCount < 140000:
        bigDataset = True

    #useRaw = conf["useRaw"]
    #if useRaw and not bigDataset:
        #adata.raw = adata # this is doing much more than assigning, it calls implicitely a function that copies
        ## a few things around. See the anndata source code under basic.py
    #else:
        #bigDataset = True
        #logging.info("Big dataset: not keeping a .raw copy of the data to save memory during analysis")

    if conf["doExp"]:
        pipeLog("Undoing log2 of data")
        adata.X = np.expm1(adata.X)

    pipeLog("Data has %d samples/observations" % len(adata.obs))
    pipeLog("Data has %d genes/variables" % len(adata.var))

    if conf["doTrimCells"]:
        minGenes = conf["minGenes"]
        pipeLog("Basic filtering: keep only cells with min %d genes" % (minGenes))
        sc.pp.filter_cells(adata, min_genes=minGenes) # adds n_genes to adata

    if conf["doTrimGenes"]:
        minCells = conf["minCells"]
        pipeLog("Basic filtering: keep only gene with min %d cells" % (minCells))
        sc.pp.filter_genes(adata, min_cells=minCells) # adds n_counts to adata.obs?

    pipeLog("After filtering: Data has %d samples/observations and %d genes/variables" % (len(adata.obs), len(adata.var)))

    if len(list(adata.obs_names))==0:
        errAbort("No cells left after filtering. Consider lowering the minGenes/minCells cutoffs in scanpy.conf")
    if len(list(adata.var_names))==0:
        errAbort("No genes left after filtering. Consider lowering the minGenes/minCells cutoffs in scanpy.conf")

    if not "n_counts" in list(adata.obs.columns.values):
        logging.debug("Adding obs.n_counts")
        adata.obs['n_counts'] = np.sum(adata.X, axis=1)

    #### PARAMETERS FOR GATING CELLS (must be changed) #####

    if conf["doFilterMito"]:
        geneIdType = conf["geneIdType"]
        if geneIdType is None or geneIdType=="auto":
            logging.info("'geneIdType' is not specified in config file or set to 'auto'.")
            geneIds = list(adata.var_names)
            geneIdType = guessGeneIdType(geneIds)

        thrsh_mito=conf["mitoMax"]
        pipeLog("Remove cells with more than %f percent of mitochondrial genes" % thrsh_mito)

        pipeLog("Computing percentage of mitochondrial genes")
        mito_genes = [name for name in adata.var_names if name.startswith('MT.') or name.startswith('MT-')]
        if len(mito_genes)>30:
            errAbort("Strange expression matrix - more than 30 mitochondrial genes?")
        if len(mito_genes)==0:
            gencodeMitos = readMitos(geneIdType)
            mito_genes = [name for name in adata.var_names if name.split('.')[0] in gencodeMitos]


        if(len(mito_genes)==0): # no single mitochondrial gene in the expression matrix ?
            pipeLog("WARNING - No single mitochondrial gene was found in the expression matrix.")
            pipeLog("Apoptotic cells cannot be removed - please check your expression matrix")
            doMito = False
            conf["mitosFound"] = False
        else:
            doMito = True

            adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)

            obsFields = list(adata.obs.columns.values)
            if "n_genes" in obsFields: # n_counts is always there
                sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True)
                #fig1 = sc.pl.scatter(adata, x='n_counts', y='n_genes', save="_gene_count")

                #if "n_counts" in obsFields:
                    #fig2 = sc.pl.scatter(adata, x='n_counts', y='percent_mito', save="_percent_mito")

            adata = adata[adata.obs['percent_mito'] < thrsh_mito, :]

    if conf["doFilterGenes"]:
        up_thrsh_genes=conf["filterMaxGenes"]
        low_thrsh_genes=conf["filterMinGenes"]
        pipeLog("Remove cells with less than %d and more than %d genes" % (low_thrsh_genes, up_thrsh_genes))

        #Filtering out cells according to filter parameters
        pipeLog('Keeping only cells with < %d genes' % up_thrsh_genes)
        adata = adata[adata.obs['n_genes'] < up_thrsh_genes, :]
        pipeLog("After filtering: Data has %d samples/observations and %d genes/variables" % (len(adata.obs), len(adata.var)))

        pipeLog('Keeping only cells with > %d genes' % low_thrsh_genes)
        adata = adata[adata.obs['n_genes'] > low_thrsh_genes, :]
        pipeLog("After filtering: Data has %d samples/observations and %d genes/variables" % (len(adata.obs), len(adata.var)))

    if conf["doNormalize"]:
        countsPerCell = conf["countsPerCell"]
        pipeLog('Expression normalization, counts per cell = %d' % countsPerCell)
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=countsPerCell)

    # change June 2019: move up log, use new trim function, make default True
    if conf["doLog"]:
        sc.pp.log1p(adata)
        pipeLog("Did log2'ing of data")

    if conf.get("doTrimVarGenes", True):
        if not conf["doLog"]:
            errAbort("you have set the option doLog=False but doTrimGenes is True. In Scanpy 1.4, this is not allowed anymore. If your dataset is already log2'ed, set doExp=True to reverse this and also doLog=True to re-apply it later, if you want to find variable genes.")

        minMean = conf["varMinMean"]
        maxMean = conf["varMaxMean"]
        minDisp = conf["varMinDisp"]
        #pipeLog('Finding highly variable genes: min_mean=%f, max_mean=%f, min_disp=%f' % (minMean, maxMean, minDisp))
        #filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=minMean, max_mean=maxMean, min_disp=minDisp)
        #sc.pl.filter_genes_dispersion(filter_result)
        #adata = adata[:, filter_result.gene_subset]
        #pipeLog('Number of variable genes identified: %d' % sum(filter_result.gene_subset))
        pipeLog('Finding highly variable genes')
        sc.pp.highly_variable_genes(adata, min_mean=minMean, max_mean=maxMean, min_disp=minDisp)
        sc.pl.highly_variable_genes(adata)
        adata = adata[:, adata.var['highly_variable']]
        pipeLog("After high-var filtering: Data has %d samples/observations and %d genes/variables" % (len(adata.obs), len(adata.var)))

    #Regress out variables nUMI, percent_mito
    if conf["doRegress"]:
        if doMito:
            pipeLog('Regressing out percent_mito and number of UMIs')
            sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
        else:
            pipeLog('Regressing out number of UMIs')
            sc.pp.regress_out(adata, ['n_counts'])

        #Scaling after regression 
        maxValue = conf["regressMax"]
        pipeLog('Scaling data, max_value=%d' % maxValue)
        sc.pp.scale(adata, max_value=maxValue)

    allPcCount = 100
    pipeLog('Performing initial PCA, number of PCs: %d' % allPcCount)
    sc.tl.pca(adata, n_comps=allPcCount)
    #Multiply by -1 to compare with Seurat
    #adata.obsm['X_pca'] *= -1
    #Plot of pca variance ratio to see if formula matches visual determination of pc_nb to use
    sc.pl.pca_variance_ratio(adata, log=True)

    #Computing number of PCs to be used in clustering
    pcCount = conf["pcCount"]
    if pcCount == "auto":
        pipeLog("Estimating number of useful PCs based on Shekar et al, Cell 2016")
        pipeLog("PC weight cutoff used is (sqrt(# of Genes/# of cells) + 1)^2")
        pipeLog("See http://www.cell.com/cell/fulltext/S0092-8674(16)31007-8, STAR methods")
        pc_cutoff= (np.sqrt((adata.n_vars/adata.n_obs))+1)**2
        pc_nb=0
        for i in adata.uns['pca']['variance']:
            if i>pc_cutoff:
                pc_nb+=1
        pipeLog('%d PCs will be used for tSNE and clustering' % pc_nb)
    else:
        pc_nb = pcCount
        pipeLog("Using %d PCs as configured in config" % pcCount)

    if "tsne" in doLayouts:
        pipeLog('Performing tSNE')
        sc.tl.tsne(adata, n_pcs=int(pc_nb), random_state=2, n_jobs=8)

    if not inCluster or "umap" in doLayouts or conf["doLouvain"]:
        neighbors = conf["louvainNeighbors"]
        res = conf["louvainRes"]
        pipeLog('Running knn, using %d PCs and %d neighbors' % (pc_nb, neighbors))
        sc.pp.neighbors(adata, n_pcs=int(pc_nb), n_neighbors=neighbors)

    if not inCluster or conf["doLouvain"]:
        pipeLog('Performing Louvain Clustering, resolution %f' % (res))
        sc.tl.louvain(adata, resolution=res)
        pipeLog("Found %d louvain clusters" % len(adata.obs['louvain'].unique()))

    clusterField = conf.get("inCluster")
    if clusterField is None:
        clusterField = "louvain"

    if "tsne" in doLayouts:
        sc.pl.tsne(adata, color=clusterField)

    #Clustering. Default Resolution: 1
    #res = 1.0
    #pipeLog('Performing Louvain Clustering, resolution = %f' % res)
    #sc.pp.neighbors(adata, n_pcs=int(pc_nb))
    #sc.tl.louvain(adata, resolution=res)
    #sc.pl.tsne(adata, color='louvain')

    if "umap" in doLayouts:
        pipeLog("Performing UMAP")
        sc.tl.umap(adata)

    if "phate" in doLayouts:
        try:
            import phate
            hasPhate = True
        except:
            hasPhate = False

        if not hasPhate:
            pipeLog("Phate is not installed, cannot run PHATE, run 'pip install phate' to install it")
        else:
            pipeLog("Performing PHATE")
            sc.tl.phate(adata)
    obsmKeys = getObsmKeys(adata)

    if "pagaFa" in doLayouts:
        pipeLog("Performing PAGA+ForceAtlas2")
        sc.tl.paga(adata, groups='louvain')
        sc.pl.paga(adata, show=True, layout="fa")
        sc.tl.draw_graph(adata, init_pos='paga', layout="fa")
        if "X_draw_graph_fa" not in obsmKeys:
            logging.warn("After paga, 'X_draw_graph_fa' not found in adata object! Older scanpy version?")
        else:
            adata.obsm["X_pagaFa"] = adata.obsm["X_draw_graph_fa"]

    if "pagaFr" in doLayouts:
        pipeLog("Performing PAGA+Fruchterman Reingold")
        sc.tl.paga(adata, groups='louvain')
        sc.pl.paga(adata, show=True, layout="fr")
        sc.tl.draw_graph(adata, init_pos='paga', layout="fr")
        if "X_draw_graph_fr" not in obsmKeys:
            logging.warn("After paga, 'X_draw_graph_fr' not found in adata object! Older scanpy version?")
        else:
            adata.obsm["X_pagaFr"] = adata.obsm["X_draw_graph_fr"]

    otherLayouts = set(doLayouts).difference(set(["tsne", "phate", "umap", "pagaFa", "pagaFr"]))

    for layoutCode in otherLayouts:
        pipeLog("Performing Layout '%s' = %s" % (layoutCode, coordLabels[layoutCode]))
        sc.tl.draw_graph(adata, layout=layoutCode)

    # Finding Top Markers, sorted by z-score
    if conf["doMarkers"]:
        nGenes = conf["markerCount"]
        pipeLog('Finding top markers for each cluster')

        sc.tl.rank_genes_groups(adata, clusterField)

        if bigDataset:
            pipeLog("Not doing rank_genes plot, big dataset")
        else:
            sc.pl.rank_genes_groups(adata, n_genes=nGenes)

    #Stop Timer
    stop= timeit.default_timer()
    pipeLog("Running time: %f" % (stop-start))

    sys.stdout = sys.stdout.terminal

    import scanpy
    conf["ScanpyVersion"] = scanpy.__version__
    return adata, conf

def mustBePython3():
    if not isPy3:
        print("Unsupported Python version.")
        print("The cellbrowser works on almost any Python version, but Scanpy requires Python3.")
        print("This script has been installed and is running under this Python: %s" % sys.executable)
        print("Most likely it has not been installed with a Python3-pip.")
        print("You will have to install cellbrowser again using a pip command that is")
        print("using Python3 where Scanpy can be installed.")
        print("")
        print("To remove confusion below, it is easier to first remove the current cellbrowser package:")
        print("    pip uninstall cellbrowser")
        print("")
        print("Then use a Python3-pip to re-install the cellbrowser package.")
        print("On mixed Python2/3 systems, the command is sometimes called 'pip3' (see: 'apt-get install python3-pip').")
        print("In conda environments, you can create a Python3 environment and activate it with:")
        print("    conda create -n py36 python=3.6 anaconda && source activate py36")
        print("If using virtualenvs, you can do something like this:")
        print("    virtualenv -p /usr/bin/python3 ~/py3env")
        print("    source ~/py3env/bin/activate")
        print("On OSX and brew, install python3 and adapt the PATH:")
        print("    brew install python3 && export PATH=/usr/local/opt/python/libexec/bin:$PATH")
        print("After this, check that your default python is python3 and pip is using python3 with:")
        print("    python --version")
        print("    pip --version")
        print("")
        print("Then install cellbrowser using the new pip:")
        print("    pip install cellbrowser")
        print("Depending on your OS and setup you may have to remove the old version first and reinstall:")
        print("    pip2 uninstall cellbrowser")
        print("    pip3 install cellbrowser --force-reinstall --upgrade")
        print("and re-run this command")
        print("")
        print("As usual, if you're not root and not on conda/virtualenv, you may need to add --user for pip")
        print("Once this is all done, install scanpy.")
        sys.exit(1)

def generateDataDesc(datasetName, outDir, algParams):
    " write a desc.conf to outDir "
    outFname = join(outDir, "desc.conf")
    c = maybeLoadConfig(outFname)
    #if isfile(outFname):
        #c = loadConfig(outFname)
    #else:
        #c = OrderedDict()

    c["title"] = datasetName
    if not "image" in c:
        c["#image"] = "thumb.png"
    if not "abstract" in c:
        c["abstract"] = "Please edit desc.conf to modify this abstract, then rerun cbBuild"
    if not "methods" in c:
        c["methods"] = "This dataset was created by a generic Scanpy pipeline run through cbScanpy."

    if not "unitDesc" in c:
        c["#unitDesc"] = "count"

    # always overwrite the parameters
    algParams = list(algParams.items())
    c["algParams"] = algParams

    writePyConf(c, outFname)

def copyTsvMatrix(matrixFname, outMatrixFname):
    " copy one file to another, but only if both look like valid input formats for cbBuild "
    if ".mtx" in matrixFname or ".h5" in matrixFname:
        logging.error("Cannot copy %s to %s, as not a text-based file like tsv or csv. You will need to do the conversion yourself manually or with cbTool.")
        return

    if not isfile(outMatrixFname):
        logging.info("File doesn't exist, copying %s to %s" % (matrixFname, outMatrixFname))
        shutil.copyfile(matrixFname, outMatrixFname)
    elif getsize(matrixFname)!=getsize(outMatrixFname):
        logging.info("File size difference, copying %s to %s" % (matrixFname, outMatrixFname))
        shutil.copyfile(matrixFname, outMatrixFname)
    else:
        logging.info("identical and same size, not copying %s to %s" % (matrixFname, outMatrixFname))

def generateQuickGenes(outDir):
    " make a quickGenes.tsv in outDir from markers.tsv "
    outFname = join(outDir, "quickGenes.tsv")
    markerFname = join(outDir, "markers.tsv")
    logging.info("Generating %s from %s" % (outFname, markerFname))

    clusters = parseMarkerTable(markerFname, None)[0]

    genesPerCluster = int(round(18 / len(clusters))) # guess a reasonable number of genes per cluster, ~ 18 genes in total
    quickGenes = defaultdict(list)
    for clusterName, rows in iterItems(clusters):
        for row in rows[:genesPerCluster]:
            sym = row[1]
            quickGenes[sym].append(clusterName)

    ofh = open(outFname, "w")
    for sym, clusterNames in iterItems(quickGenes):
        ofh.write("%s\t%s\n" % (sym, ", ".join(clusterNames)))
    ofh.close()

def cbScanpyCli():
    " command line interface for cbScanpy "
    mustBePython3()

    global options
    args, options = cbScanpy_parseArgs()

    if options.init:
        copyPkgFile("sampleConfig/scanpy.conf")
        sys.exit(1)

    try:
        logging.info("Loading Scanpy libraries")
        import scanpy as sc
    except:
        print("The Python package 'scanpy' is not installed in the current interpreter %s" % sys.executable)
        print("Please install it following the instructions at https://scanpy.readthedocs.io/en/latest/installation.html")
        print("We recommend the miniconda-based installation.")
        print("Then re-run this command.")
        sys.exit(1)

    matrixFname = options.exprMatrix
    metaFname = options.meta
    outDir = options.outDir
    confFname = options.confFname
    inCluster = options.inCluster
    copyMatrix = options.copyMatrix

    if copyMatrix and not matrixFname.endswith(".gz"):
        errAbort("If you use the --copyMatrix option, the input matrix must be gzipped. Please run 'gzip %s' and then re-run cbScanpy" % matrixFname)

    makeDir(outDir)

    figDir = join(outDir, "figs")
    adFname = join(outDir, "anndata.h5ad")
    matrixOutFname = join(outDir, "exprMatrix.tsv.gz")

    logFname = join(outDir, "cbScanpy.log")
    if isfile(logFname):
        os.remove(logFname)

    adata, params = cbScanpy(matrixFname, metaFname, inCluster, confFname, figDir, logFname)


    # anndata in newer versions can't save without the ordering so force an ordering now
    import pandas as pd
    for colName in adata.obs.columns:
        col = adata.obs[colName]
        if col.dtype.kind!="O":
            continue
        logging.debug("Converting column %s to ordered categories" % colName)
        dt = pd.api.types.CategoricalDtype(col.unique(), ordered=True)
        #newCol = col.astype("category", categories=col.unique(), ordered=True)
        newCol = col.astype(dt)
        adata.obs[colName] = newCol

    logging.info("Writing final result as an anndata object to %s" % adFname)
    adata.write(adFname)
    datasetName=options.name

    scanpyToCellbrowser(adata, outDir, datasetName=datasetName,
            clusterField=inCluster, skipMatrix=copyMatrix, useRaw=True)

    if copyMatrix:
        outMatrixFname = join(outDir, "exprMatrix.tsv.gz")
        copyTsvMatrix(matrixFname, outMatrixFname)

    #generateHtmls(datasetName, outDir)
    generateDataDesc(datasetName, outDir, params)

def mtxToTsvGz(mtxFname, geneFname, barcodeFname, outFname, translateIds=False):
    " convert mtx to tab-sep without scanpy. gzip if needed "
    import scipy.io
    import numpy as np
    logging.info("Reading matrix from %s, %s and %s" % (mtxFname, geneFname, barcodeFname))

    genes = []
    sep = sepForFile(geneFname)
    for l in openFile(geneFname):
        row = l.rstrip("\r\n").split(sep) # field 3 is "Gene Expression" in cr3
        if len(row)>1:
            geneId, sym = row[:2]
            genes.append(geneId+"|"+sym)
        else:
            geneId = row[0]
            genes.append(geneId)

    barcodes = [l.strip() for l in openFile(barcodeFname) if l!="\n"]

    if translateIds:
        geneToSym = readGeneSymbols(None, genes)
        genes = [geneToSym[geneId] for geneId in genes]

    logging.info("Read %d genes and %d barcodes" % (len(genes), len(barcodes)))

    logging.info("Reading expression matrix...")
    mat = scipy.io.mmread(mtxFname)

    logging.info("Dimensions of matrix: %d , %d" % mat.shape)
    #geneCount, cellCount = mat.shape

    if mat.shape[0]==len(genes)-1:
        genes = genes[1:]
        logging.info("The genes file seems to have a GEO-like header, removed the first gene, genecount is now %d" % len(genes))

    if mat.shape[1]==len(barcodes)-1:
        barcodes = barcodes[1:]
        logging.info("The barcodes file seems to have a GEO-like header, removed the first barcode, count is now %d" % len(barcodes))

    if mat.shape[0]==len(barcodes):
        logging.info("Matrix looks like it's in transposed format (genes on columns), so transposing matrix now")
        mat = mat.transpose()
        logging.info("New dimensions of matrix: %d , %d" % mat.shape)

    logging.info("Converting matrix to row-based layout...")
    mat = mat.tocsr()

    #print(mat.shape[0])
    #print(len(genes))
    assert(mat.shape[0]==len(genes)) # matrix gene count has to match gene tsv file line count
    assert(mat.shape[1]==len(barcodes)) # matrix cell count has to match barcodes tsv file line count

    logging.info("Writing matrix to text")
    tmpFname = outFname+".tmp"

    # could not find a cross-python way to open ofh for np.savetxt
    # see https://github.com/maximilianh/cellBrowser/issues/73 and numpy ticket referenced therein
    if isPy3:
        ofh = open(tmpFname, "w")
    else:
        ofh = open(tmpFname, "wb")

    ofh.write("gene\t")
    ofh.write("\t".join(barcodes))
    ofh.write("\n")

    for i in range(0, len(genes)):
        if i%1000==0:
            logging.info("%d genes written..." % i)
        ofh.write(genes[i])
        ofh.write("\t")
        arr = mat[i].toarray()
        fmt = "%d"
        if arr.dtype==np.int64:
            np.savetxt(ofh, arr, "%d", "\t", "\n")
        else:
            np.savetxt(ofh, arr, "%g", "\t", "\n")
    ofh.close()

    moveOrGzip(tmpFname, outFname)

    logging.info("Created %s" % outFname)

def getAllFields(ifhs, sep):
    " give a list of file handles, get all non-gene headers and return as a list of names "
    fieldToFname = {}
    allFields = []
    colCounts = []
    for ifh in ifhs:
        fields = ifh.readline().rstrip('\r\n').split(sep)
        fields = [f.strip('"') for f in fields] # R sometimes adds quotes
        colCounts.append(len(fields))
        #assert(fields[0]=="#gene")
        for i, field in enumerate(fields):
            # make sure we have no field name overlaps
            if field in fieldToFname and i!=0:
                raise Exception("field %s seen twice, in %s and %s" % (field, ifh.name, fieldToFname[field]))
            fieldToFname[field] = ifh.name

            # only use first field name from first file
            if i==0:
                if len(allFields)==0:
                    allFields.append(field)
            else:
                allFields.append(field)

    return allFields, colCounts

def generateDownloads(datasetName, outDir):
    htmlFname = join(outDir, "downloads.html")
    if isfile(htmlFname):
        logging.info("%s exists, not overwriting" % htmlFname)
        return

    ofh = open(htmlFname, "w")
    ofh.write("<b>Expression matrix:</b> <a href='%s/exprMatrix.tsv.gz'>exprMatrix.tsv.gz</a><p>\n" % datasetName)

    cFname = join(outDir, "cellbrowser.conf")
    conf = None
    if isfile(cFname):
        conf = loadConfig(cFname)
        if "unit" in conf:
            ofh.write("Unit of expression matrix: %s<p>\n" % conf.get("unit", "unknown"))
    ofh.write("<b>Cell meta annotations:</b> <a href='%s/meta.tsv'>meta.tsv</a><p>" % datasetName)

    markerFname = join(outDir, "markers.csv")
    if not isfile(markerFname):
        markerFname = join(outDir, "markers.tsv")

    if isfile(markerFname):
        baseName = basename(markerFname)
        ofh.write("<b>Cluster Marker Genes:</b> <a href='%s/%s'>%s</a><p>\n" % (datasetName, baseName, baseName))

    if conf:
        coordDescs = conf["coords"]
        for coordDesc in coordDescs:
            coordLabel = coordDesc["shortLabel"]
            cleanName = sanitizeName(coordLabel.replace(" ", "_"))
            coordFname = cleanName+".coords.tsv.gz"
            ofh.write("<b>%s coordinates:</b> <a href='%s/%s'>%s</a><br>\n" % (coordLabel, datasetName, coordFname, coordFname))

    rdsFname = join(datasetName, "seurat.rds")
    if isfile(rdsFname):
        ofh.write("<b>Seurat R data file:</b> <a href='%s'>seurat.rds</a><p>\n" % rdsFname)

    scriptFname = join(datasetName, "runSeurat.R")
    if isfile(scriptFname):
        ofh.write("<b>Seurat R analysis script:</b> <a href='%s'>runSeurat.R</a><p>\n" % scriptFname)

    logFname = join(datasetName, "analysisLog.txt")
    if isfile(logFname):
        ofh.write("<b>Analysis Log File:</b> <a href='%s'>analysisLog.txt</a><p>\n" % logFname)

    h5adFname = join(datasetName, "anndata.h5ad")
    if isfile(h5adFname):
        ofh.write("<b>Scanpy Anndata HDF5 file:</b> <a href='%s'>anndata.h5ad</a><p>\n" % h5adFname)

    ofh.close()
    logging.info("Wrote %s" % ofh.name)

def generateHtmls(datasetName, outDir):
    " generate downloads.html and summary.html in outDir, if they don't exist "
    copyPkgFile("sampleConfig/summary.html", outDir, {"datasetName" :datasetName})
    generateDownloads(datasetName, outDir)

if __name__=="__main__":
    args, options = main_parseArgs()

    cmd = args[0]
    if cmd=="cbServe":
        outDir, port = args[1:3]
        savePid()
        startHttpServer(outDir, int(port))
    else:
        errAbort("Unknown command %s" % cmd)
