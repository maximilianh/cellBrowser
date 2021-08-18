
from collections import defaultdict 
import json, gzip, hashlib
from os.path import *

fileInfo = {}

fileInfoFname = "files.json"
if isfile(fileInfoFname):
    fileInfo = json.load(open(fileInfoFname))

def geneTsvToJson(bedFname, code, label, jsonFname):
    global fileInfo
    " convert BED file to more compact json file: chrom -> list of (start, end, strand, gene) "
    print("Processing %s" % fname)
    bySym = defaultdict(dict)
    for line in open(bedFname):
        row = line.rstrip("\n").split("\t")
        #ENST00000456328.2_1	chr1	+	11868	14409	DDX11L1
        #transId, chrom, strand, start, end, sym = row
        chrom, start, end, score, strand = row[6]
        geneSym = row[-1]
        start = int(start)
        end = int(end)
        transLen = end-start
        bySym[sym].setdefault(chrom, []).append( (transLen, start, end, strand, transId) )

    symLocs = defaultdict(list)
    for sym, chromDict in bySym.items():
        #if len(chromDict) > 1:
            #print("sym", sym)
        for chrom, transList in chromDict.items():
            transList.sort(reverse=True) # take longest transcript per chrom
            _, start, end, strand, transId = transList[0]
            symLocs[chrom].append( (start, end, strand, sym) )
    
    sortedLocs = {}
    for chrom, geneList in symLocs.items():
        geneList.sort()
        sortedLocs[chrom] = geneList

    ofh = gzip.open(jsonFname, "wt")
    outs = json.dumps(sortedLocs)
    md5 = hashlib.md5(outs.encode("utf8")).hexdigest()[:10]
    ofh.write(outs)

    ofh.close()

    fileInfo[code] = {"label":label, "file" : jsonFname, "md5" :md5}

# hgsql hg19 -NB -e 'select * from wgEncodeGencodeBasicV34lift37' | cut -f2- | cut -f1,2,3,4,5,12 > hg19.gc34.tsv
# hgsql hg38 -NB -e 'select * from wgEncodeGencodeBasicV34' | cut -f2- | cut -f1,2,3,4,5,12 > hg38.gc34.tsv
# hgsql mm10 -NB -e 'select * from wgEncodeGencodeBasicVM25' | cut -f2- | cut -f1,2,3,4,5,12 > mm10.vm25.tsv

geneTsvToJson("hg19.gc34.tsv", "hg19gc34", "hg19 Gencode34backlift", "hg19.json.gz")
geneTsvToJson("hg38.gc34.tsv", "hg38gc38", "hg38 Gencode34", "hg38.json.gz")
geneTsvToJson("mm10.vm25.tsv", "mm10vm25", "mm10 GencodeVM25", "mm10.json.gz")

json.dump(fileInfo, open(fileInfoFname, "wt"), indent=4)
