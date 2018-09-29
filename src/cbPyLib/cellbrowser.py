#!/usr/bin/env python2

# this library mostly contains functions that convert tab-sep files
# (=single cell expression matrix and meta data) into the binary format that is read by the
# javascript viewer cbWeb/js/cellbrowser.js and cbData.js.
# Helper functions here allow importing data from other tools, e.g. cellranger or scanpy.

# requires at least python2.6, version tested was 2.6.6
# should work with python2.5, not tested
# works on python3, version tested was 3.6.5

import logging, sys, optparse, struct, json, os, string, shutil, gzip, re, unicodedata
import zlib, math, operator, doctest, copy, bisect, array, glob, io, time, subprocess
import hashlib
from distutils import spawn
from collections import namedtuple, OrderedDict
from os.path import join, basename, dirname, isfile, isdir, relpath, abspath, getsize, getmtime

try:
    from collections import defaultdict
    from collections import Counter
except:
    # python2.6 has no defaultdict or Counter yet
    from backport_collections import defaultdict # error? -> pip2 install backport-collections
    from backport_collections import Counter # error? -> pip2 install backport-collections

# We do not require numpy but numpy is around 30-40% faster in serializing arrays
# So use it if it's present
numpyLoaded = True
try:
    import numpy as np
except:
    numpyLoaded = False
    logging.error("Numpy could not be loaded. The script will work, but it will be 30% slower when processing the matrix.")

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
dataDir = join(dirname(__file__), "..", "cbData")

defOutDir = os.environ.get("CBOUT")

# ==== functions =====

def setDebug(options):
    " activate debugging if needed "
    if options.debug:
        logging.basicConfig(level=logging.DEBUG)
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
        logging.getLogger().setLevel(logging.INFO)

def cbBuild_parseArgs(showHelp=False):
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] -i cellbrowser.conf -o outputDir - add a dataset to the single cell viewer directory

    If you have previously built into the same output directory with the same dataset and the
    expression matrix has not changed its filesize, this will be detected and the expression
    matrix will not be copied again. This means that an update of a few meta data attributes
    is quite quick.

    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages")

    parser.add_option("-i", "--inConf", dest="inConf", action="append",
        help="a cellbrowser.conf file that specifies labels and all input files, default %default, can be specified multiple times")

    parser.add_option("-o", "--outDir", dest="outDir", action="store", help="output directory, default can be set through the env. variable CBOUT, current value: %default", default=defOutDir)

    parser.add_option("-p", "--port", dest="port", action="store",
        help="if build is successful, start an http server on this port and serve the result via http://localhost:port", type="int")

    (options, args) = parser.parse_args()

    if showHelp:
        parser.print_help()
        exit(1)

    setDebug(options)

    return args, options

def cbMake_parseArgs():
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("usage: %prog [options] outDir - copy all relevant js/css files into outDir, look for datasets in it and create index.html")

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages")
    parser.add_option("-o", "--outDir", dest="outDir", action="store", help="output directory, default can be set through the env. variable CBOUT, current value: %default", default=defOutDir)

    #parser.add_option("-m", "--meta", dest="meta", action="store",
        #help="meta data tsv file, aka sample sheet. One row per sample, first row has headers, first column has sample name."
        #)
    #parser.add_option("-e", "--matrix", dest="matrix", action="store",
        #help="expression matrix file, one gene per row, one sample per column. First column has gene identifiers (Ensembl or symbol), First row has sample names. ")
    #parser.add_option("-c", "--coords", dest="coords", action="append", help="tab-sep table with cell coordinates, format: metaId, x, y. Can be specified multiple times, if you have multiple coordinate files.")

    (options, args) = parser.parse_args()

    if options.outDir==None:
        parser.print_help()
        exit(1)

    setDebug(options)
    return args, options

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

def lineFileNextRow(inFile, utfHacks=False):
    """
    parses tab-sep file with headers in first line
    yields collection.namedtuples
    strips "#"-prefix from header line
    utfHacks forces all chars to latin1 and removes anything that doesn't fit into latin1
    """

    if isinstance(inFile, str):
        # input file is a string = file name
        fh = openFile(inFile)
        sep = sepForFile(inFile)
    else:
        fh = inFile
        sep = "\t"

    line1 = fh.readline()
    line1 = line1.strip("\n").lstrip("#")
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

    headers = [re.sub("[^a-zA-Z0-9_]","_", h) for h in headers]
    headers = [re.sub("^_","", h) for h in headers] # remove _ prefix
    #headers = [x if x!="" else "noName" for x in headers]
    if headers[0]=="": # R does not name the first column by default
        headers[0]="rowName"

    if "" in headers:
        logging.error("Found empty cells in header line of %s" % inFile)
        logging.error("This often happens with Excel files. Make sure that the conversion from Excel was done correctly. Use cut -f-lastColumn to fix it.")
        assert(False)

    filtHeads = []
    for h in headers:
        if h[0].isdigit():
            filtHeads.append("x"+h)
        else:
            filtHeads.append(h)
    headers = filtHeads


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
        #if isPy3:
            #fields = line.split(sep, maxsplit=len(headers)-1)
        #else:
            #fields = string.split(line, sep, maxsplit=len(headers)-1)
        fields = line.split(sep)

        if sep==",":
            fields = [x.lstrip('"').rstrip('"') for x in fields]

        try:
            rec = Record(*fields)
        except Exception as msg:
            logging.error("Exception occured while parsing line, %s" % msg)
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
    line1 = open(fname).readline().rstrip("\r\n")
    fieldCount = line1.split('\t')
    if "geneId" not in line1:
        d = parseDict(fname)
    # my gencode tables contain a symbol for all genes
    # the old format
    elif line1=="transcriptId\tgeneId\tsymbol":
        for row in lineFileNextRow(fname):
            if row.symbol=="":
                continue
            d[row.geneId.split(".")[0]]=row.symbol
    # my new files are smaller and have headers
    elif line1=="geneId\tsymbol" or fieldCount==2:
        d = {}
        for row in lineFileNextRow(fname):
            if row.symbol=="":
                continue
            d[row.geneId.split(".")[0]]=row.symbol
    else:
        assert(False)
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

def findBins(numVals, breakVals):
    """
    find the right bin index defined by breakVals for every value in numVals.
    Special handling for the last value. The comparison uses "<=". The first
    break is assumed to be the minimum of numVals and is therefore ignored.
    Also returns an array with the count for every bin.
    >>> findBins([1,1,1,2,2,2,3,3,4,4,5,5,6,6], [1, 2,3,5,6])
    ([0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 3, 3], [6, 2, 4, 2])
    """
    breaks = breakVals[1:]
    bArr = []
    binCounts = [0]*len(breaks)
    for x in numVals:
        binIdx = bisect.bisect_left(breaks, x)
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
        for i, x in values:
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

    dArr, binCounts = findBins(numVals, breakVals)
    assert(len(binCounts)==10)
    logging.info("Number of values per decile-bin: %s" % binCounts)

    fieldMeta["binMethod"] = "quantiles"
    fieldMeta["binCounts"] = binCounts
    fieldMeta["breaks"] = breakVals

    return dArr, fieldMeta

def discretizeNumField(numVals, fieldMeta, numType):
    " given a list of numbers, add attributes to fieldMeta that describe the binning scheme "
    #digArr, fieldMeta = discretizeArr_uniform(numVals, fieldMeta)
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

def guessFieldMeta(valList, fieldMeta, colors, forceEnum):
    """ given a list of strings, determine if they're all int, float or
    strings. Return fieldMeta, as dict, and a new valList, with the correct python type
    - 'type' can be: 'int', 'float', 'enum' or 'uniqueString'
    - if int or float: 'deciles' is a list of the deciles
    - if uniqueString: 'maxLen' is the length of the longest string
    - if enum: 'values' is a list of all possible values
    - if colors is not None: 'colors' is a list of the default colors
    """
    intCount = 0
    floatCount = 0
    valCounts = defaultdict(int)
    #maxVal = 0
    for val in valList:
        fieldType = "string"
        try:
            newVal = int(val)
            intCount += 1
            floatCount += 1
            #maxVal = max(newVal, val)
        except:
            try:
                newVal = float(val)
                floatCount += 1
                #maxVal = max(newVal, val)
            except:
                pass

        valCounts[val] += 1

    valToInt = None

    if floatCount==len(valList) and intCount!=len(valList) and len(valCounts) > 10 and not forceEnum:
        # field is a floating point number: convert to decile index
        numVals = [float(x) for x in valList]

        newVals, fieldMeta = discretizeNumField(numVals, fieldMeta, "float")

        fieldMeta["type"] = "float"
        #fieldMeta["maxVal"] = maxVal

    elif intCount==len(valList) and not forceEnum:
        # field is an integer: convert to decile index
        numVals = [int(x) for x in valList]
        newVals, fieldMeta = discretizeNumField(numVals, fieldMeta, "int")
        fieldMeta["type"] = "int"
        #fieldMeta["maxVal"] = maxVal

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
                    colArr.append("DDDDDD") # wonder if I should not stop here
            if foundColors > 0:
                fieldMeta["colors"] = colArr
                if len(notFound)!=0:
                    logging.warn("No default color found for field values %s" % notFound)

        fieldMeta["valCounts"] = valCounts
        fieldMeta["arrType"], fieldMeta["_fmt"] = bytesAndFmt(len(valArr))
        valToInt = dict([(y[0],x) for (x,y) in enumerate(valCounts)]) # dict with value -> index in valCounts
        newVals = [valToInt[x] for x in valList] #

    #fieldMeta["valCount"] = len(valList)
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

def runGzip(fname, finalFname=None):
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
    return finalFname

def metaToBin(inConf, outConf, fname, colorFname, outDir, enumFields):
    """ convert meta table to binary files. outputs fields.json and one binary file per field.
    adds names of metadata fields to outConf and returns outConf
    """
    logging.info("Converting to numbers and compressing meta data fields")
    makeDir(outDir)

    colData = parseIntoColumns(fname)

    colors = parseColors(colorFname)

    fieldInfo = []
    for colIdx, (fieldName, col) in enumerate(colData):
        logging.info("Meta data field index %d: '%s'" % (colIdx, fieldName))

        forceEnum = False
        if enumFields!=None:
            forceEnum = (fieldName in enumFields)
        cleanFieldName = cleanString(fieldName)
        binName = join(outDir, cleanFieldName+".bin")

        fieldMeta = OrderedDict()
        fieldMeta["name"] = cleanFieldName
        fieldMeta["label"] = fieldName

        if fieldName=="cluster" or fieldName=="Cluster":
            forceEnum=True
        fieldMeta, binVals = guessFieldMeta(col, fieldMeta, colors, forceEnum)
        fieldType = fieldMeta["type"]

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

        del fieldMeta["_fmt"]
        fieldInfo.append(fieldMeta)
        if "type" in fieldMeta:
            logging.info(("Type: %(type)s, %(diffValCount)d different values" % fieldMeta))
        else:
            logging.info(("Type: %(type)s, %(diffValCount)d different values, max size %(maxSize)d " % fieldMeta))

    return fieldInfo

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
    " open a .tsv file and yield rows via iterRows. gz and csv OK."

    def __init__(self, geneToSym=None):
        self.geneToSym = geneToSym

    def open(self, fname, matType=None):
        " return something for iterMatrixTsv "
        logging.debug("Opening %s" % fname)
        self.fname = fname
        if fname.endswith(".gz"):
            #ifh = gzip.open(fname)
            self.ifh = subprocess.Popen(
                ['gunzip', '-c', fname],
                stdout=subprocess.PIPE,
                encoding='utf-8',
            ).stdout # faster, especially with two CPUs
        else:
            self.ifh = open(fname, "rU")

        self.sep = "\t"
        if ".csv" in fname.lower():
            self.sep = ","
            logging.debug("Field separator is %s" % repr(self.sep))

        headLine = self.ifh.readline().rstrip("\r\n")
        self.sampleNames = headLine.split(self.sep)[1:]
        self.sampleNames = [x.strip('"') for x in self.sampleNames]
        assert(len(self.sampleNames)!=0)
        logging.debug("Read %d sampleNames, e.g. %s" % (len(self.sampleNames), self.sampleNames[0]))

        if matType is None:
            self.matType = self.autoDetectMatType(10)
            logging.info("Numbers in matrix are of type '%s'", self.matType)
        else:
            self.matType = matType

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

        skipIds = 0
        doneGenes = set()
        lineNo = 0
        sep = self.sep
        sampleCount = len(self.sampleNames)
        geneToSym = self.geneToSym
        for line in self.ifh:
            self.lineLen = len(line)

            if isPy3:
                gene, rest = line.rstrip("\r\n").split(sep, maxsplit=1)
            else:
                gene, rest = string.split(line.rstrip("\r\n"), sep, maxsplit=1)

            if numpyLoaded:
                arr = np.fromstring(rest, dtype=npType, sep=sep, count=sampleCount)
            else:
                if self.matType=="int":
                    #a = [int(x) for x in rest.split(sep)]
                    arr = map(int, rest.split(sep))
                else:
                    #a = [float(x) for x in rest.split(sep)]
                    arr = map(float, rest.split(sep))

            if "|" in gene:
                gene, symbol = gene.split("|")
            else:
                if geneToSym is None:
                    symbol = gene
                else:
                    symbol = geneToSym.get(gene)
                    if symbol is None:
                        skipIds += 1
                        logging.warn("line %d: %s is not a valid Ensembl gene ID, check geneIdType setting in cellbrowser.conf" % (lineNo, gene))
                        continue
                    if symbol.isdigit():
                        logging.warn("line %d in gene matrix: gene identifier %s is a number. If this is indeed a gene identifier, you can ignore this warning." % (lineNo, symbol))

            if symbol in doneGenes:
                logging.warn("line %d: Gene %s/%s is duplicated in matrix, using only first occurence" % (lineNo, gene, symbol))
                skipIds += 1
                continue

            doneGenes.add(gene)

            lineNo += 1

            yield gene, symbol, arr

        if skipIds!=0:
            logging.warn("Skipped %d expression matrix lines because of duplication/unknown ID" % skipIds)

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

def matrixToBin(fname, geneToSym, binFname, jsonFname, discretBinFname, discretJsonFname):
    """ convert gene expression vectors to vectors of deciles
        and make json gene symbol -> (file offset, line length)
    """
    logging.info("converting %s to %s and writing index to %s" % (fname, binFname, jsonFname))
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
    matReader.open(fname)
    matType = matReader.getMatType()
    sampleNames = matReader.getSampleNames()

    geneCount = 0
    for geneId, sym, exprArr in matReader.iterRows():
        geneCount += 1

        #if maxVal(exprArr) > 200:
            #highCount += 1

        logging.debug("Processing %s, symbol %s" % (geneId, sym))
        exprStr = exprEncode(geneId, exprArr, matType)
        exprIndex[sym] = (ofh.tell(), len(exprStr))
        ofh.write(exprStr)

        if geneCount % 1000 == 0:
            logging.info("Wrote expression values for %d genes" % geneCount)

    discretOfh.close()
    ofh.close()

    #if highCount==0:
        #logging.warn("No single value in the matrix is > 200. It looks like this matrix has been log'ed before. Our recommendation for visual inspection is to not transform matrices, but that is of course up to you.")
        #logging.error("Rerun with --skipLog.")
        #sys.exit(1)

    if len(exprIndex)==0:
        errAbort("No genes from the expression matrix could be mapped to symbols."
            "Are you sure these are Ensembl IDs? Adapt geneIdType in cellbrowser.conf. Example ID: %s" % geneId)

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
    if ".csv" in fname:
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

        if isPy3:
            row = line.split(sep, maxsplit=1)
        else:
            row = string.split(line, sep, maxsplit=1)
        field1 = row[0]

        lineLen = end - start
        assert(lineLen!=0)
        assert(lineLen<65535) # meta data line cannot be longer than 2 bytes
        ofh.write(struct.pack("<L", start))
        ofh.write(struct.pack("<H", lineLen))
    ofh.close()

def testMetaIndex(outDir):
    # test meta index
    fh = open(join(outDir, "meta.index"))
    #fh.seek(10*6)
    o = fh.read(4)
    s = fh.read(2)
    offset = struct.unpack("<L", o) # little endian
    l = struct.unpack("<H", s)
    #print "offset, linelen:", offset, l

    #fh = open(join(outDir, "meta/meta.tsv"))
    #fh.seek(offset[0])
    #print fh.read(l[0])

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
    Flip the y coordinates to make it more look like plots in R, for people transitioning from R.
    """
    logging.info("Parsing coordinates from %s. FlipY=%s, useTwoBytes=%s" % (fname, flipY, useTwoBytes))
    coords = []
    maxY = 0
    minX = 2^32
    minY = 2^32
    maxX = -2^32
    maxY = -2^32
    skipCount = 0

    # parse and find the max values
    for row in lineFileNextRow(fname):
        assert(len(row)==3) # coord file has to have three rows (cellId, x, y), we just ignore the headers
        cellId = row[0]
        x = float(row[1])
        y = float(row[2])
        minX = min(x, minX)
        minY = min(y, minY)
        maxX = max(x, maxX)
        maxY = max(y, maxY)
        coords.append( (cellId, x, y) )

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

def metaReorder(matrixFname, metaFname, fixedMetaFname):
    """ check and reorder the meta data, has to be in the same order as the
    expression matrix, write to fixedMetaFname """

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

    matNotMeta = meta - mat
    metaNotMat = mat - meta
    stop = False
    mustFilterMatrix = False
    if len(metaNotMat)!=0:
        logging.warn("%d samples names are in the meta data, but not in the expression matrix. Examples: %s" % (len(metaNotMat), list(metaNotMat)[:10]))
        logging.warn("These samples will be removed from the meta data")
        matrixSampleNames = [x for x in matrixSampleNames if x in meta]
        mustFilterMatrix = True

    if len(matNotMeta)!=0:
        logging.warn("%d samples names are in the expression matrix, but not in the meta data. Examples: %s" % (len(matNotMeta), list(matNotMeta)[:10]))
        logging.warn("These samples will be removed from the expression matrix")

    # filter the meta data file
    logging.info("Data contains %d samples/cells" % len(matrixSampleNames))

    # slurp in the whole meta data
    tmpFname = fixedMetaFname+".tmp"
    ofh = open(tmpFname, "w")
    metaToRow = {}
    sep = sepForFile(metaFname)
    for lNo, line in enumerate(open(metaFname)):
        row = line.rstrip("\r\n").split(sep)
        if lNo==0:
            # copy header over
            ofh.write("\t".join(row))
            ofh.write("\n")
            continue
        row = line.rstrip("\r\n").split(sep)
        metaToRow[row[0]] = row

    # and write it in the right order
    for matrixName in matrixSampleNames:
        ofh.write("\t".join(metaToRow[matrixName]))
        ofh.write("\n")
    ofh.close()
    os.rename(tmpFname, fixedMetaFname)

    return matrixSampleNames, mustFilterMatrix

def writeCoords(coordName, coords, sampleNames, coordBinFname, coordJson, useTwoBytes, coordInfo, textOutName):
    """ write coordinates given as a dictionary to coordBin and coordJson, in the order of sampleNames
    Also return as a list.
    """
    tmpFname = coordBinFname+".tmp"
    logging.info("Writing coordinates to %s and %s" % (coordBinFname, coordJson))
    binFh = open(tmpFname, "wb")

    minX = 2^32
    minY = 2^32
    maxX = -2^32
    maxY = -2^32
    xVals = []
    yVals = []

    textOutTmp = textOutName+".tmp"
    textOfh = open(textOutTmp, "w")

    for sampleName in sampleNames:
        coordTuple = coords.get(sampleName)
        if coordTuple is None:
            logging.warn("sample name %s is in meta file but not in coordinate file %s, setting to (0,0)" % (sampleName, coordName))
            x = 0
            y = 0
        else:
            x, y = coordTuple
            textOfh.write("%s\t%f\t%f\n" % (sampleName, x, y))
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

    binFh.close()
    os.rename(tmpFname, coordBinFname)

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

    logging.info("Wrote %d coordinates to %s and %s" % (len(sampleNames), coordBinFname, textOutName))
    return coordInfo, xVals, yVals

def runCommand(cmd):
    " run command "
    logging.debug("Running %s" % cmd)
    err = os.system(cmd)
    if err!=0:
        errAbort("Could not run: %s" % cmd)
    return 0

def copyMatrixTrim(inFname, outFname, filtSampleNames, doFilter):
    " copy matrix and compress it. If doFilter is true: keep only the samples in filtSampleNames"
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
        return

    sep = "\t"

    logging.info("Copying+reordering+trimming %s to %s, keeping only the %d columns with a sample name in the meta data" % (inFname, outFname, len(filtSampleNames)))

    matIter = MatrixTsvReader()
    matIter.open(inFname)

    sampleNames = matIter.getSampleNames()

    keepFields = set(filtSampleNames)
    keepIdx = []
    for i, name in enumerate(sampleNames):
        if name in keepFields:
            keepIdx.append(i)

    tmpFname = outFname+".tmp"

    ofh = openFile(tmpFname, "w")
    ofh.write("gene\t")
    ofh.write("\t".join(filtSampleNames))
    ofh.write("\n")

    count = 0
    for geneId, sym, exprArr in matIter.iterRows():
        newRow = [geneId]
        for idx in keepIdx:
            newRow.append(str(exprArr[idx]))
        ofh.write("\t".join(newRow))
        ofh.write("\n")
        count += 1
        if count%1000==0:
            logging.info("Wrote %d rows" % count)
    ofh.close()

    #tmpFnameGz = outFname+".tmp.gz"
    #runCommand("gzip -c %s > %s " % (tmpFname, tmpFnameGz))
    #os.remove(tmpFname)
    #os.rename(tmpFnameGz, outFname)
    runGzip(tmpFname, outFname)

def convIdToSym(geneToSym, geneId):
    if geneToSym is None:
        return geneId
    else:
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
    logging.debug("Sanitizing %s -> %s" % (repr(name), newName))
    assert(len(newName)!=0)
    return newName

def splitMarkerTable(filename, geneToSym, outDir):
    """ split .tsv on first field and create many files in outDir with columns 2-end.
    """
    if filename is None:
        return
    logging.info("Splitting cluster markers from %s into directory %s" % (filename, outDir))
    #logging.debug("Splitting %s on first field" % filename)
    ifh = openFile(filename)

    headers = ifh.readline().rstrip("\r\n").split('\t')
    otherHeaders = headers[2:]

    data = defaultdict(list)
    columns = defaultdict(list)
    for line in ifh:
        row = line.rstrip("\r\n").split('\t')
        for colIdx, val in enumerate(row[1:]):
            columns[colIdx].append(val)

        clusterName = row[0]
        geneId = row[1]
        scoreVal = float(row[2])
        otherFields = row[3:]

        #geneSym = convIdToSym(geneToSym, geneId)
        geneSym = geneId # let's assume for now that the marker table already has symbols

        newRow = []
        newRow.append(geneId)
        newRow.append(geneSym)
        newRow.append(scoreVal)
        newRow.extend(otherFields)

        data[clusterName].append(newRow)

    colTypes = {}
    for colIdx, vals in iterItems(columns):
        colTypes[colIdx] = typeForStrings(vals)

    headersWithType = []
    for colIdx, header in enumerate(otherHeaders):
        if colTypes[colIdx+1]!="string":
            header = header+"|"+colTypes[colIdx+1]
        headersWithType.append(header)

    newHeaders = ["id", "symbol"]
    newHeaders.extend(headersWithType)

    fileCount = 0
    sanNames = set()
    for clusterName, rows in iterItems(data):
        #rows.sort(key=operator.itemgetter(2), reverse=True) # rev-sort by score (fold change)
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
            scoreVal = row[2]
            row[2] = "%0.5E" % row[2] # limit to 5 digits
            ofh.write("\t".join(row))
            ofh.write("\n")
        ofh.close()

        runGzip(outFname)

        fileCount += 1
    logging.info("Wrote %d .tsv.gz files into directory %s" % (fileCount, outDir))

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

def loadConfig(fname):
    """ parse python in fname and return variables as dictionary.
    add the directory of fname to the dict as 'inDir'.
    """
    logging.debug("Loading config from %s" % fname)
    g = {}
    l = OrderedDict()
    execfile(fname, g, l)

    conf = l

    if not "coords" in conf:
        errAbort("The input configuration has to define the 'coords' statement")
    if not "meta" in conf:
        errAbort("The input configuration has to define the 'meta' statement")
    if not "exprMatrix" in conf:
        errAbort("The input configuration has to define the 'exprMatrix' statement")
    if "tags" in conf and type(conf["tags"])!=type([]):
        errAbort("'tags' in input config file must be a list")

    conf["inDir"] = dirname(fname)

    return conf

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
    " copy dataset description html files to output directory "
    filesToCopy = []

    outConf["desc"] = {}

    for fileBase in ["summary.html", "methods.html", "downloads.html", "thumb.png"]:
        inFname = makeAbs(inDir, fileBase)
        if not isfile(inFname):
            logging.info("%s does not exist" % inFname)
        else:
            #copyFiles.append( (fname, "summary.html") )
            outPath = join(datasetDir, fileBase)
            logging.debug("Copying %s -> %s" % (inFname, outPath))
            shutil.copy(inFname, outPath)

            fileDesc = fileBase.split(".")[0]
            outConf["desc"][fileDesc] = fileBase

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
    logging.info("Calculating cluster midpoints")
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
        clusterXVals[clusterIdx].append(xVals[i])
        clusterYVals[clusterIdx].append(yVals[i])

    midInfo = []
    for clustIdx, xList in enumerate(clusterXVals):
        yList = clusterYVals[clustIdx]
        # get the midpoint of this cluster
        midX = sum(xList) / float(len(xList))
        midY = sum(yList) / float(len(yList))

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

        clusterName = labelVals[clustIdx]
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
    ifh = openFile(fname, "rt")
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

    geneInfo = []
    hasDesc = None
    hasPmid = None
    for row in lineFileNextRow(fname):
        if hasDesc == None:
            if "desc" in row._fields:
                hasDesc = True
        if hasPmid == None:
            if "pmid" in row._fields:
                hasPmid = True
        sym = row.symbol
        if validSyms is not None and sym not in validSyms:
            logging.error("'%s' is not a valid gene gene symbol, skipping it" % sym)
            continue

        info = [sym]
        if hasDesc:
            info.append(row.desc)
        if hasPmid:
            info.append(row.pmid)
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
            logging.error("invalid sample name - line %d in %s: sample name (first field) is empty" % (i, fname))
            sys.exit(1)
        if metaName in doneNames:
            logging.error("sample name duplicated - line %d in %s: sample name %s (first field) has been seen before" % (i, fname, metaName))
            sys.exit(1)

        doneNames.add(metaName)
        sampleNames.append(row[0])
        i+=1
    logging.debug("Found %d sample names, e.g. %s" % (len(sampleNames), sampleNames[:3]))
    return sampleNames

def convertExprMatrix(inConf, outMatrixFname, outConf, metaSampleNames, geneToSym, outDir, needFilterMatrix):
    """ trim a copy of the expression matrix for downloads, also create an indexed
    and compressed version
    """

    # step1: copy expression matrix, so people can download it, potentially
    # removing those sample names that are not in the meta data
    matrixFname = getAbsPath(inConf, "exprMatrix")
    outConf["fileVersions"]["inMatrix"] = getFileVersion(matrixFname)
    copyMatrixTrim(matrixFname, outMatrixFname, metaSampleNames, needFilterMatrix)

    # step2: discretize expression matrix for the viewer, compress and index to file
    #logging.info("quick-mode: Not compressing matrix, because %s already exists" % binMat)
    binMat = join(outDir, "exprMatrix.bin")
    binMatIndex = join(outDir, "exprMatrix.json")
    discretBinMat = join(outDir, "discretMat.bin")
    discretMatrixIndex = join(outDir, "discretMat.json")

    matType = matrixToBin(outMatrixFname, geneToSym, binMat, binMatIndex, discretBinMat, discretMatrixIndex)

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
    useTwoBytes = inConf.get("useTwoBytes", False)

    hasLabels = False
    if "labelField" in inConf and inConf["labelField"] is not None:
        hasLabels = True
        clusterLabelField = inConf["labelField"]
        labelVec, labelVals = parseTsvColumn(outMeta, clusterLabelField)
        outConf["labelField"] = clusterLabelField

    newCoords = []
    for coordIdx, coordInfo in enumerate(coordFnames):
        coordFname = coordInfo["file"]
        coordLabel = coordInfo["shortLabel"]
        coords = parseScaleCoordsAsDict(coordFname, useTwoBytes, flipY)
        coordName = "coords_%d" % coordIdx
        coordDir = join(outDir, "coords", coordName)
        makeDir(coordDir)
        coordBin = join(coordDir, "coords.bin")
        coordJson = join(coordDir, "coords.json")
        coordInfo = OrderedDict()
        coordInfo["name"] = coordName
        coordInfo["shortLabel"] = coordLabel
        cleanName = sanitizeName(coordLabel.replace(" ", "_"))
        textOutName = join(outDir, cleanName+".coords.tsv.gz")
        coordInfo, xVals, yVals = writeCoords(coordLabel, coords, sampleNames, coordBin, coordJson, useTwoBytes, coordInfo, textOutName)
        newCoords.append( coordInfo )

        if hasLabels:
            clusterMids = makeMids(xVals, yVals, labelVec, labelVals, coordInfo)

            midFname = join(coordDir, "clusterLabels.json")
            midFh = open(midFname, "w")
            json.dump(clusterMids, midFh, indent=2)
            logging.info("Wrote cluster labels and midpoints to %s" % midFname)

    outConf["coords"] = newCoords
    copyConf(inConf, outConf, "labelField")
    copyConf(inConf, outConf, "useTwoBytes")

def readAcronyms(inConf, outConf):
    " read the acronyms and save them into the config "
    inDir = inConf["inDir"]
    fname = inConf.get("acroFname")
    if fname is not None:
        fname = makeAbs(inDir, fname)
        if not isfile(fname):
            logging.warn("%s specified in config file, but does not exist, skipping" % fname)
        else:
            acronyms = parseDict(fname)
            logging.info("Read %d acronyms from %s" % (len(acronyms), fname))
            outConf["acronyms"] = acronyms

def convertMarkers(inConf, outConf, geneToSym, outDir):
    " split the marker tables into one file per cluster "
    markerFnames = []
    if "markers" in inConf:
        markerFnames = makeAbsDict(inConf, "markers")

    newMarkers = []
    for markerIdx, markerInfo in enumerate(markerFnames):
        markerFname = markerInfo["file"]
        markerLabel = markerInfo["shortLabel"]

        clusterName = "markers_%d" % markerIdx # use sha1 of input file ?
        markerDir = join(outDir, "markers", clusterName)
        makeDir(markerDir)

        splitMarkerTable(markerFname, geneToSym, markerDir)

        newMarkers.append( {"name" : sanitizeName(clusterName), "shortLabel" : markerLabel})
    outConf["markers"] = newMarkers

def readQuickGenes(inConf, geneToSym, outConf):
    quickGeneFname = inConf.get("quickGenesFile")
    if quickGeneFname:
        fname = getAbsPath(inConf, "quickGenesFile")
        quickGenes = parseGeneInfo(geneToSym, fname)
        outConf["quickGenes"] = quickGenes
        logging.info("Read %d quick genes from %s" % (len(quickGenes), fname))

def getFileVersion(fname):
    metaVersion = {}
    metaVersion["fname"] = fname
    hexHash = md5ForFile(fname).decode("ascii")
    metaVersion["md5"] = hexHash
    metaVersion["size"] = getsize(fname)
    metaVersion["mtime"] = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(getmtime(fname)))
    return metaVersion

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
    fieldConf = metaToBin(inConf, outConf, finalMetaFname, colorFname, metaDir, enumFields)
    outConf["metaFields"] = fieldConf

    indexMeta(finalMetaFname, metaIdxFname)

    logging.info("Kept %d cells present in both meta data file and expression matrix" % len(sampleNames))

    outConf["fileVersions"]["outMeta"] = getFileVersion(finalMetaFname)

    return sampleNames, needFilterMatrix, finalMetaFname

def readGeneSymbols(inConf):
    " return geneToSym, based on gene tables "
    geneIdType = inConf.get("geneIdType")
    if geneIdType==None:
        logging.warn("'geneIdType' is not set in input config. Gene IDs will not be converted to symbols. Assuming that the matrix already has symbols. ")
        geneIdType = "symbols"

    if geneIdType.startswith('symbol'):
        return None

    searchMask = join(dataDir, "genes", geneIdType+".symbols.tsv")
    fnames = glob.glob(searchMask)
    assert(len(fnames)<=1)
    if(len(fnames)==0):
        errAbort("Could not find any files matching %s. Possible files that were found: %s" % (searchMask, searchMask))
    geneIdTable = fnames[0]
    geneToSym = readGeneToSym(geneIdTable)
    return geneToSym

def readMitos(org):
    ' return the gene IDs of all mitochondrial genes. 37 for human for all gencode versions '
    if org=="human":
        geneToSym = readGeneSymbols({'geneIdType':"gencode22"})
    else:
        assert(False) # not doing mouse just yet

    mitos = []
    for geneId, sym in iterItems(geneToSym):
        if sym.startswith("MT-"):
            mitos.append(geneId)
    logging.debug("Found %d mitochondrial genes for %s, e.g. %s" % (len(mitos), org, mitos[0]))
    return mitos

def getAbsPath(conf, key):
    " get assume that value of key in conf is a filename and use the inDir value to make it absolute "
    return abspath(join(conf["inDir"], conf[key]))

def getMd5Using(md5Cmd, fname):
    " posix command line tool is much faster than python "
    logging.debug("Getting md5 of %s using %s command line tool" % (fname, md5Cmd))
    cmd = [md5Cmd, fname]
    logging.debug("Cmd: %s" % cmd)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    md5 = proc.stdout.readline()
    proc.stdout.close()
    stat = os.waitpid(proc.pid, 0)
    err = stat[1]
    assert(err==0)
    return md5

def md5ForFile(fname):
    " return the md5sum of a file. Use a command line tool, if possible. "
    logging.info("Getting md5 of %s" % fname)
    if spawn.find_executable("md5sum")!=None:
        md5 = getMd5Using("md5sum", fname).split()[0]
    elif spawn.find_executable("md5")!=None:
        md5 = getMd5Using("md5", fname).split()[-1]
    else:
        hash_md5 = hashlib.md5()
        with open(fname, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                hash_md5.update(chunk)
        md5 = hash_md5.hexdigest()
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

    # this obscure command gets the cell identifiers in the dataset directory
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

def convertDataset(inConf, outConf, datasetDir):
    """ convert everything needed for a dataset to datasetDir, write config to outConf.
    If the expression matrix has not changed since the last run, and the sampleNames are the same,
    it won't be converted again.
    """
    copyDatasetHtmls(inConf["inDir"], outConf, datasetDir)

    # some config settings are passed through unmodified to the javascript
    for tag in ["name", "shortLabel", "radius", "alpha", "priority", "tags",
        "clusterField", "hubUrl", "showLabels", "ucscDb", "unit"]:
        copyConf(inConf, outConf, tag)

    if " " in inConf["name"]:
        errAbort("Sorry, please no whitespace in the dataset name in the .conf file")

    # convertMeta also compares the sample IDs between meta and matrix
    # outMeta is a reordered/trimmed tsv version of the meta table
    sampleNames, needFilterMatrix, outMeta = convertMeta(inConf, outConf, datasetDir)

    geneToSym = None

    outMatrixFname = join(datasetDir, "exprMatrix.tsv.gz")
    geneToSym = -1 # None would mean "there are no gene symbols to map to"
    inMatrixFname = getAbsPath(inConf, "exprMatrix")
    doMatrix = matrixOrSamplesHaveChanged(datasetDir, inMatrixFname, outMatrixFname, outConf)

    if doMatrix:
        geneToSym = readGeneSymbols(inConf)
        convertExprMatrix(inConf, outMatrixFname, outConf, sampleNames, geneToSym, datasetDir, needFilterMatrix)
        # in case script crashes after this, keep the current state of the config
        writeConfig(inConf, outConf, datasetDir)
    else:
        logging.info("Matrix and meta sample names have not changed, not indexing matrix again")

    convertCoords(inConf, outConf, sampleNames, outMeta, datasetDir)

    if geneToSym==-1:
        geneToSym = readGeneSymbols(inConf)
    convertMarkers(inConf, outConf, geneToSym, datasetDir)

    readAcronyms(inConf, outConf)

    readQuickGenes(inConf, geneToSym, outConf)

def writeAnndataCoords(anndata, fieldName, outDir, filePrefix, fullName, desc):
    #if 'X_draw_graph_fa' in anndata.obsm.dtype.names:
    import pandas as pd
    fileBase = filePrefix+"_coords.tsv"
    fname = join(outDir, fileBase)
    if fieldName in anndata.obsm.dtype.names:
        logging.info("Writing %s coords to %s" % (fullName, fname))
        fa2_coord=pd.DataFrame(anndata.obsm[fieldName],index=anndata.obs.index)
        fa2_coord.columns=['x','y']
        fa2_coord.to_csv(fname,sep='\t')
        desc.append( {'file':fileBase, 'shortLabel': fullName} )
    else:
        logging.warn('Couldnt find %s coordinates' % fullName)

def writeCellbrowserConf(name, coordsList, fname, args={}):
    for c in name:
        assert(c.isalnum() or c in ["-", "_"]) # only digits and letters are allowed in dataset names

    metaFname = args.get("meta", "cell_to_cluster.tsv")
    clusterField = args.get("clusterField", "Louvain CLuster")

    conf = """
name='%(name)s'
shortLabel='%(name)s'
exprMatrix='exprMatrix.tsv.gz'
#tags = ["10x", 'smartseq2']
meta='%(metaFname)s'
geneIdType='symbols'
clusterField='%(clusterField)s'
labelField='%(clusterField)s'
enumFields=['%(clusterField)s']
markers = [{"file": "markers.tsv", "shortLabel":"Cluster Markers"}]
coords=%(coordsList)s
radius=5
alpha=0.6
""" % locals()

    if "geneToSym" in args:
        conf += "geneToSym='%s'\n" % args["geneToSym"]

    #fname = join(outDir, 'cellbrowser.conf')
    if isfile(fname):
        logging.info("Not overwriting %s, file already exists." % fname)
        return

    ofh = open(fname, "w")
    ofh.write(conf)
    ofh.close()
    logging.info("Wrote %s" % ofh.name)

def anndataToTsv(anndata, matFname):
    " write anndata to .tsv file and gzip it "
    logging.info("Writing matrix to %s" % matFname)
    tmpFname = matFname+".tmp"
    import pandas as pd
    adT = anndata.T
    data_matrix=pd.DataFrame(adT.X, index=adT.obs.index.tolist(), columns=adT.var.index.tolist())
    data_matrix.to_csv(tmpFname, sep='\t', index=True)
    os.rename(tmpFname, matFname)

def scanpyToTsv(anndata, path, datasetName, meta_option=None, nb_marker=50):
    """
    Written by Lucas Seninge, lucas.seninge@etu.unistra.fr

    Given a scanpy object, write dataset to output directory under path.

    This function export files needed for the ucsc cells viewer from the Scanpy Anndata object
    :param anndata: Scanpy AnnData object where information are stored
    :param path : Path to folder where to save data (tsv tables)
    :param meta_option: list of metadata names (string) present
    in the AnnData objects(other than 'louvain' to also save (eg: batches, ...))
    :param nb_marker: number of cluster markers to store. Default: 100

    """
    confName = join(path, "cellbrowser.conf")
    if isfile(confName):
        errAbort("File %s already exists. Cowardly refusing to overwrite it. Please move the file and re-run this command" % confName)

    import numpy as np
    import pandas as pd
    import scanpy.api as sc

    ##Save data matrix to tsv
    #if "raw" in dir(anndata):
    #    adT = anndata.raw.T
    #else:
    matFname = join(path, 'exprMatrix.tsv')
    anndataToTsv(anndata, matFname)
    matFname = runGzip(matFname)

    coordDescs = []
    writeAnndataCoords(anndata, "X_tsne", path, "tsne", "T-SNE", coordDescs)
    writeAnndataCoords(anndata, "X_umap", path, "umap", "UMAP", coordDescs)
    writeAnndataCoords(anndata, "X_draw_graph_fa", path, "fa2", "ForceAtlas2", coordDescs)
    writeAnndataCoords(anndata, "X_pagaFa2", path, "pagaFa2", "PAGA+ForceAtlas2", coordDescs)
    writeAnndataCoords(anndata, "X_pagaUmap", path, "pagaUmap", "PAGA+UMAP", coordDescs)
    writeAnndataCoords(anndata, "X_phate", path, "phate", "PHATE", coordDescs)

    ##Check for louvain clustering
    if 'louvain' in anndata.obs:
        #Export cell <-> cluster identity
        fname = join(path, 'meta.tsv')
        # add prefix to make sure that it's not treated as a number
        #anndata.obs[['louvain']]['louvain'] = "cluster "+anndata.obs[['louvain']]['louvain'].astype(str)
        anndata.obs[['louvain']].to_csv(fname,sep='\t', header=["Louvain Cluster"])
    else:
        errAbort('Couldnt find clustering information')

    ##Check for cluster markers
    if 'rank_genes_groups' in anndata.uns:
        top_score=pd.DataFrame(anndata.uns['rank_genes_groups']['scores']).loc[:nb_marker]
        top_gene=pd.DataFrame(anndata.uns['rank_genes_groups']['names']).loc[:nb_marker]
        marker_df= pd.DataFrame()
        for i in range(len(top_score.columns)):
            concat=pd.concat([top_score[[str(i)]],top_gene[[str(i)]]],axis=1,ignore_index=True)
            concat['cluster_number']=i
            col=list(concat.columns)
            col[0],col[-2]='z_score','gene'
            concat.columns=col
            marker_df=marker_df.append(concat)
    else:
        errAbort ('Couldnt find cluster markers list')

    #Rearranging columns -> Cluster, gene, score
    cols=marker_df.columns.tolist()
    cols=cols[::-1]
    marker_df=marker_df[cols]
    #Export
    fname = join(path, "markers.tsv")
    pd.DataFrame.to_csv(marker_df,fname,sep='\t',index=False)

    ##Save more metadata
    if meta_option != None:
        meta_df=pd.DataFrame()
        for element in meta_option:
            if element not in anndata.obs:
                print(str(element) + ' field is not present in the AnnData object')
            else:
                temp=anndata.obs[[element]]
                meta_df=pd.concat([meta_df,temp],axis=1)
        fname = join(path, "meta.tsv")
        meta_df.to_csv(fname,sep='\t')

    writeCellbrowserConf(datasetName, coordDescs, confName)
    #ofh = open(confName, "w")
    #ofh.write("coords = %s\n" % repr(coordDescs))
    #ofh.write("meta = 'meta.tsv'\n")
    #ofh.write("name = %s" % repr(datasetName))
    #ofh.write("exprMatrix = 'exprMatrix.tsv'")
    #ofh.close()

def writeJson(data, outFname):
    """ https://stackoverflow.com/a/37795053/233871 """
    # Make it work for Python 2+3 and with Unicode
    try:
        to_unicode = unicode
    except NameError:
        to_unicode = str

    # Write JSON file
    tmpName = outFname+".tmp"
    with io.open(tmpName, 'w', encoding='utf8') as outfile:
        #str_ = json.dumps(data, indent=2, sort_keys=True,separators=(',', ': '), ensure_ascii=False)
        str_ = json.dumps(data, indent=2, separators=(',', ': '), ensure_ascii=False)
        outfile.write(to_unicode(str_))
    os.rename(tmpName, outFname)
    logging.info("Wrote %s" % outFname)

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
    import RangeHTTPServer
    try:
        # py3
        import http.server as SimpleHTTPServer
        from http.server import HTTPServer
    except:
        # py2
        import SimpleHTTPServer
        from BaseHTTPServer import HTTPServer

    #server_address = ('localhost', port)
    server_address = ('', port)
    HandlerClass = RangeHTTPServer.RangeRequestHandler
    HandlerClass.protocol_version = "HTTP/1.0"
    httpd = HTTPServer(server_address, HandlerClass)

    sa = httpd.socket.getsockname()
    os.chdir(outDir)
    print("Serving "+outDir+". Press Ctrl-C to exit.")
    print("Point your internet browser to http://"+sa[0]+":"+str(sa[1])+" (or the address of this server)")
    sys.stderr = open("/dev/null", "w") # don't show http status message on console
    httpd.serve_forever()

def convertAndCopy(confFnames, outDir, port):
    " build browser from config files confFnames into directory outDir and serve on port "
    for inConfFname in confFnames:
        inConf = loadConfig(inConfFname)
        datasetDir = join(outDir, inConf["name"])
        makeDir(datasetDir)

        outConfFname = join(outDir, "dataset.conf")
        #if onlyMeta:
            #outConf = json.parse(open(outConfFname)) # reuse the old config
        #else:
        outConf = OrderedDict()

        convertDataset(inConf, outConf, datasetDir)

        writeConfig(inConf, outConf, datasetDir)

    cbMake(outDir)

    if port:
        startHttpServer(outDir, port)

def convertAndCopyCli():
    " command line interface for dataset converter, also copies the html/js/etc files "
    args, options = cbBuild_parseArgs()

    confFnames = options.inConf
    if confFnames==None:
        confFnames = ["cellbrowser.conf"]

    for fname in confFnames:
        if not isfile(fname):
            logging.error("File %s does not exist." % fname)
            cbBuild_parseArgs(showHelp=True)
    if options.outDir is None:
        logging.error("You have to specify at least the output directory or set the env. variable CBOUT.")
        cbBuild_parseArgs(showHelp=True)

    outDir = options.outDir
    #onlyMeta = options.onlyMeta
    port = options.port

    convertAndCopy(confFnames, outDir, port)

def cbCellrangerCli_parseArgs(showHelp=False):
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] -i cellRangerDir -o outputDir - convert the cellranger output to cellbrowser format and create a cellranger.conf file

    """)

    parser.add_option("-d", "--debug", dest="debug", action="store_true",
        help="show debug messages")

    parser.add_option("-i", "--inDir", dest="inDir", action="store", help="input folder with the cellranger analysis output. This is the directory with the .h5 file.")
    parser.add_option("-o", "--outDir", dest="outDir", action="store", help="output directory")
    #parser.add_option("-g", "--geneSet", dest="geneSet", action="store", help="geneset, e.g. gencode28 or gencode-m13 or similar. Default: %default", default="gencode24")
    parser.add_option("-n", "--name", dest="datasetName", action="store", help="name of the dataset, default is %default", default="cellrangerImport")

    (options, args) = parser.parse_args()

    if showHelp:
        parser.print_help()
        exit(1)

    setDebug(options)

    return args, options

def cbCellrangerCli():
    args, options = cbCellrangerCli_parseArgs()

    if options.outDir is None or options.inDir is None:
        logging.error("You have to specify at least an input and an output directory.")
        cbCellrangerCli_parseArgs(showHelp=True)

    crangerToCellbrowser(options.datasetName, options.inDir, options.outDir)

def readMatrixAnndata(matrixFname, samplesOnRows=False):
    " read an expression matrix and return an adata object. Supports .mtx, .h5 and .tsv (not .tsv.gz) "
    import scanpy.api as sc
    #adata = sc.read(matFname)
    if matrixFname.endswith(".mtx"):
        import pandas as pd
        logging.info("Loading expression matrix: mtx format")
        adata = sc.read(matrixFname, cache=False).T

        mtxDir = dirname(matrixFname)
        adata.var_names = pd.read_csv(join(mtxDir, 'genes.tsv'), header=None, sep='\t')[1]
        adata.obs_names = pd.read_csv(join(mtxDir, 'barcodes.tsv'), header=None)[0]

    else:
        logging.info("Loading expression matrix: tab-sep format")
        adata = sc.read(matrixFname, cache=False , first_column_names=True)
        if not samplesOnRows:
            info("Transposing the expression matrix")
            adata = adata.T

    return adata

def mtxToTsvGz(mtxFname, geneFname, barcodeFname, outFname):
    " convert mtx to tab-sep without scanpy. gzip if needed "
    from scipy import io
    import numpy as np
    logging.info("Reading matrix from %s, %s and %s" % (mtxFname, geneFname, barcodeFname))
    mat = io.mmread(mtxFname)
    genes = [l.strip() for l in open(geneFname)]
    genes = [g.replace("\t", "|") for g in genes]
    barcodes = [l.strip() for l in open(barcodeFname)]
    mat = mat.tocsr()

    logging.info("Writing matrix to text")
    tmpFname = outFname+".tmp"
    ofh = open(tmpFname, "w")

    ofh.write("gene\t")
    ofh.write("\t".join(barcodes))
    ofh.write("\n")

    geneCount, cellCount = mat.shape

    assert(geneCount==len(genes))
    assert(cellCount==len(barcodes))

    for i in range(0, geneCount):
        ofh.write(genes[i])
        ofh.write("\t")
        arr = mat[i].toarray()
        fmt = "%d"
        assert(arr.dtype==np.int64) # float not supported yet. Email me or open a ticket.
        # todo for float: when using float, we need to find a way to store 0 as 0 and not as 0+0000
        np.savetxt(ofh, arr, "%d", "\t", "\n")
    ofh.close()

    logging.info("Compressing expression matrix...")
    runGzip(tmpFname, outFname)
    logging.info("Wrote %s" % outFname)

def crangerToCellbrowser(datasetName, inDir, outDir):
    " convert cellranger output to a cellbrowser directory "
    # copy over the clusters
    clustFname = join(inDir, "analysis/clustering/graphclust/clusters.csv")
    metaFname = join(outDir, "meta.csv")
    shutil.copy(clustFname, metaFname)

    # copy over the t-SNE coords
    tsneFname = join(inDir, "analysis/tsne/2_components/projection.csv")
    coordFname = join(outDir, "tsne.coords.csv")
    shutil.copy(tsneFname, coordFname)

    # copy over the markers
    dgeFname = join(inDir, "analysis/diffexp/graphclust/differential_expression.csv")
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
    crangerWriteDownloads(datasetName, outDir)

def crangerWriteMethods(inDir, outDir, matFname):
    htmlFname = join(outDir, "methods.html")
    if isfile(htmlFname):
        logging.info("%s exists, not overwriting" % htmlFname)
        return

    import csv
    csvMask = join(inDir, "*_metrics_summary.csv")
    csvFnames = glob.glob(csvMask)
    assert(len(csvFnames)==1)
    qcVals = list(csv.DictReader(open(csvFnames[0])))[0]

    ofh = open(htmlFname, "w")
    ofh.write("This dataset was imported from a CellRanger analysis directory with cbCellranger.<p><p>")
    ofh.write("<p><b>QC metrics reported by CellRanger:</b></p>\n")

    for key, value in iterItems(qcVals):
        ofh.write("%s: %s<br>\n" % (key, value))
    ofh.close()
    logging.info("Wrote %s" % ofh.name)

def crangerWriteDownloads(datasetName, outDir):
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

def crangerSignMarkers(dgeFname, markerFname, geneFname, maxPval, maxGenes):
    " convert cellranger diff exp file to markers.tsv file "
    ofh = open(markerFname, "w")
    ofh.write("cluster\tgene\tAdj. P-Value\tLog2 fold change\tMean UMI Counts\n")

    clusterToGenes = defaultdict(list)

    # read the significant markers and their p-Values
    for line in open(dgeFname):
        if line.startswith("Gene ID"):
            continue
        row = line.rstrip("\n\r").split(",")
        geneId = row[0]
        sym = row[1]
        clusterCount = int(len(row) / 3)
        for clusterIdx in range(0, clusterCount):
            startField = (clusterIdx*3)+2
            mean = float(row[startField])
            fc = float(row[startField+1])
            pVal = float(row[startField+2])
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

        datasetDesc = json.load(open(fname))
        assert("name" in datasetDesc) # every dataset has to have a name

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
        if isdir(filename):
            continue
        logging.debug("Copying %s to %s" % (filename, outDir))
        shutil.copy(filename, outDir)

def copyStatic(baseDir, outDir):
    " copy all js, css and img files to outDir "
    logging.info("Copying js, css and img files to %s" % outDir)
    imgDir = join(outDir, "img")

    copyAllFiles(baseDir, "ext/images", outDir)
    copyAllFiles(baseDir, "img", outDir)
    copyAllFiles(baseDir, "ext", outDir)
    copyAllFiles(baseDir, "js", outDir)
    copyAllFiles(baseDir, "css", outDir)

def makeIndexHtml(baseDir, datasets, outDir):
    dsList = []
    for ds in datasets:
        summDs = {
                "shortLabel" : ds["shortLabel"],
                "sampleCount" : ds["sampleCount"],
                "name" : ds["name"]
                }

        if "tags" in ds:
            assert(type(ds["tags"])==type([])) # "tags" have to be a list, not a string or dict
            summDs["tags"] = ds["tags"]

        dsList.append(summDs)

    indexFname = join(baseDir, "html", "index.html")
    indexStr = open(indexFname).read()
    old = "datasetList = null"
    new = "datasetList = "+json.dumps(dsList, sort_keys=True, indent=4, separators=(',', ': '))
    newIndexStr = indexStr.replace(old, new)
    assert(newIndexStr!=indexStr)

    newFname = join(outDir, "index.html")
    ofh = open(newFname, "w")
    ofh.write(newIndexStr)
    ofh.close()

    datasetLabels = [x["name"] for x in dsList]
    logging.info("Wrote %s, added datasets: %s" % (newFname, " - ".join(datasetLabels)))

def cbMake(outDir):
    " create index.html in outDir and copy over all other static files "
    baseDir = dirname(__file__) # = directory of this script
    webDir = join(baseDir, "..", "cbWeb")
    copyStatic(webDir, outDir)
    datasets = findDatasets(outDir)
    makeIndexHtml(webDir, datasets, outDir)

def cbMake_cli():
    " command line interface for copying over the html and js files "
    args, options = cbMake_parseArgs()
    outDir = options.outDir

    if outDir is None:
        errAbort("You have to specify at least the output directory or set the environment variable CBOUT.")

    cbMake(outDir)

def parseGeneLocs(geneType):
    """
    return dict with geneId -> list of bedRows
    bedRows have (chrom, start, end, geneId, score, strand)
    """
    fname = join(dataDir, "genes", geneType+".genes.bed")
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
    fname = join(dataDir, "genomes", genome+".sizes")
    assert(isfile(fname))
    return fname

def makeBarGraphBigBed(genome, inMatrixFname, outMatrixFname, geneType, clusterToCells, \
        clusterOrder, clusterFname, bbFname):
    """ create a barGraph bigBed file for an expression matrix
    clusterToCells is a dict clusterName -> list of cellIDs
    clusterOrder is a list of the clusterNames in the right order
    """
    logging.info("*** Creating barChartGraph bigbed file")
    if geneType.startswith("symbol"):
        # create a mapping from symbol -> gene locations
        if "/" in geneType:
            defGenes = geneType.split("/")[1]
        elif genome=="hg38":
            defGenes = "gencode24"
        elif genome=="hg19":
            defGenes = "gencode19"
        elif genome=="mm10":
            defGenes = "gencode-m13"
        else:
            errAbort("Unclear how to map symbols to genome for db %s. Please adapt cellbrowser.py" % genome)

        logging.info("Using %s to map symbols to genome" % defGenes)

        geneToSym = readGeneSymbols({'geneIdType':defGenes})
        geneLocsId = parseGeneLocs(defGenes)
        geneLocs = {}
        for geneId, locs in iterItems(geneLocsId):
            sym = geneToSym[geneId]
            geneLocs[sym] = locs
    else:
        geneToSym = readGeneSymbols({'geneIdType':geneType})
        geneLocs = parseGeneLocs(geneType)

    matOfh = open(outMatrixFname, "w")
    clustOfh = open(clusterFname, "w")

    mr = MatrixTsvReader()
    mr.open(inMatrixFname)
    matType, cellNames = mr.matType, mr.sampleNames

    cellIds = range(0, len(cellNames))
    cellNameToId = dict(zip(cellNames, cellIds))

    # make a list of lists of cellIds, one per cluster, in the right order
    clusterCellIds = [] # list of tuples with cell-indexes, one per cluster
    allCellNames = [] # list for cellIds, with a matrix, meta and with bam file
    allCellIndices = [] # position of all cellIds in allCellNames
    for clusterName in clusterOrder:
        cellIdxList = []
        for cellName in clusterToCells[clusterName]:
            if cellName not in cellNameToId:
                logging.warn("%s is in meta but not in expression matrix." % cellName)
                continue
            idx = cellNameToId[cellName]
            cellIdxList.append(idx)
            allCellNames.append(cellName)
            allCellIndices.append(idx)
            sanClusterName = clusterName.replace(" ", "_")
            clustOfh.write("%s\t%s\n" % (cellName, sanClusterName))

        if len(cellIdxList)==0:
            logging.warn("No cells assigned to cluster %s" % clusterName)

        clusterCellIds.append(tuple(cellIdxList))
    clustOfh.close()

    # write header line
    matOfh.write("#gene\t")
    matOfh.write("\t".join(allCellNames))
    matOfh.write("\n")

    # make the barchart bed file. format:
    # chr14 95086227 95158010 ENSG00000100697.10 999 - DICER1 5 10.94,11.60,8.00,6.69,4.89 93153 26789
    #bedFname = join(outDir, "barchart.bed")
    bedFname = bbFname.replace(".bb", ".bed")
    assert(bedFname!=bbFname)

    bedFh = open(bedFname, "w")

    skipCount = 0
    for geneId, sym, exprArr in mr.iterRows():
        logging.debug("Writing BED and matrix line for %s" % geneId)

        # write the new matrix row
        offset = matOfh.tell()
        rowHeader = "%s\t" % (geneId)
        matOfh.write(rowHeader)

        newRow = []
        for idx in allCellIndices:
            newRow.append(str(exprArr[idx]))
        newLine = "\t".join(newRow)
        matOfh.write(newLine)
        matOfh.write("\n")
        lineLen = len(geneId)+len(newLine)+2 # include tab and newline

        medianList = []

        for cellIds in clusterCellIds:
            exprList = []
            for cellId in cellIds:
                exprList.append(exprArr[cellId])
            n = len(cellIds)
            if len(exprList)==0:
                median = 0
            else:
                median = sorted(exprList)[n//2] # approx OK, no special case for even n's
            medianList.append(str(median))
            bedScore = len([x for x in exprList if x!=0]) # score = non-zero medians
            bedScore = min(1000, bedScore)

        if geneId not in geneLocs:
            geneId2 = geneId.replace(".", "-", 1) # does this make sense? (for R)
            if geneId2 not in geneLocs:
                logging.warn("Cannot place gene '%s' onto genome, dropping it" % geneId)
                skipCount += 1
                continue
            else:
                geneId = geneId2

        bedRows = geneLocs[geneId]

        # one geneId may have multiple placements, e.g. Ensembl's rule for duplicate genes
        for bedRow in bedRows:
            sym = geneToSym.get(geneId, geneId)
            bedRow[4] = str(bedScore) # 4 = score field
            bedRow.append(sym)
            bedRow.append(str(len(medianList)))
            bedRow.append(",".join(medianList))
            bedRow.append(str(offset))
            bedRow.append(str(lineLen))

            bedFh.write("\t".join(bedRow))
            bedFh.write("\n")

    bedFh.close()

    if skipCount != 0:
        logging.info("Could not place %d genes, these were skipped" % skipCount)

    bedFname2 = bedFname.replace(".bed", ".sorted.bed")
    cmd = "LC_COLLATE=C sort -k1,1 -k2,2n %s > %s" % (bedFname, bedFname2)
    runCommand(cmd)

    # convert to .bb using .as file
    # from https://genome.ucsc.edu/goldenpath/help/examples/barChart/barChartBed.as
    asFname = join(dataDir, "genomes", "barChartBed.as")
    sizesFname = getSizesFname(genome)

    cmd = "bedToBigBed -as=%s -type=bed6+5 -tab %s %s %s" % (asFname, bedFname2, sizesFname, bbFname)
    runCommand(cmd)

#if __name__=="__main__":
    #main()
    #import scanpy.api as sc
    #ad = sc.read("sampleData/quakeBrainGeo1.old/geneMatrix.tsv")
    #ad = ad.T
    #convScanpy(ad, "temp", "./")
