# create a track hub from an expression matrix

import subprocess, colorsys
import logging, sys, optparse, re, unicodedata, string, glob, distutils.spawn, gzip
from collections import defaultdict, namedtuple
from os.path import join, basename, dirname, isfile, isdir, splitext
import os

import sys
import cellbrowser

CBEMAIL = os.getenv("CBEMAIL", "unknown")

# ==== functions =====
    
def cbHub_parseArgs():
    " setup logging, parse command line arguments and options. -h shows auto-generated help page "
    parser = optparse.OptionParser("""usage: %prog [options] - create a UCSC track hub for a single cell dataset

    Run the program with "--init" to write a sample cellbrowser.conf into the current directory.
    Adapt this file to your needs. For %prog, only the parts under "hub" are relevant.

    Call the program without any arguments if there is cellbrowser.conf in the current dir:
    %prog
    It will read cellbrowser.conf, create hub.txt and write the job scripts.
    You can then run the job scripts as described below.

    This will create one shell script per cell cluster, with names like 'job-xxx.sh',
    in the current directory. To run these scripts in parallel on a big server,
    e.g. with 10 processes use the Unix command "parallel":
      ls *.sh | parallel -j 10 bash

    Each job will write a log of stdout/stderr to files named e.g. job-1.log

    Requirements that have to be installed on this machine with binaries in PATH:
    - samtools https://github.com/samtools/samtools
    - wiggletools https://github.com/Ensembl/WiggleTools
    - intronProspector https://github.com/diekhans/intronProspector
    - UCSC tools wigToBigWig, bedToBigBed from http://hgdownload.soe.ucsc.edu/admin/exe/

    Run the UCSC tool hubCheck on the resulting hub.txt to make sure that all jobs
    have successfully completed.

    You can override some settings in cellbrowser.conf with the command line options:
    %prog -m metaFname -c clusterField -m exprMatrix -o hubDir
    """)

    parser.add_option("", "--init", dest="init", action="store_true", \
            help="write a sample cellbrowser.conf to the current directory")

    parser.add_option("-i", "--inConf", dest="inConf", action="store", \
            help="a cellbrowser.conf input file to read all options from, default %default", \
            default = "cellbrowser.conf")

    #parser.add_option("-g", "--genome", dest="genome", action="store", \
            #help="a ucsc assembly identifier, like hg19, hg38 or mm10")

    parser.add_option("-m", "--metaFname", dest="meta", action="store", \
            help="a csv or tsv matrix, one row per cell")

    parser.add_option("-e", "--exprMatrix", dest="exprMatrix", action="store", \
            help="exprMatrix is a tsv or csv expression matrix, one line per cell")

    parser.add_option("-c", "--clusterField", dest="clusterField", action="store", \
            help="field in expr matrix that contains the cluster name")

    parser.add_option("-o", "--hubDir", dest="hubDir", action="store", \
            help="the output directory for the hub, default is %default")

    parser.add_option("-d", "--debug", dest="debug", action="store_true", help="show debug messages")
    #parser.add_option("", "--fixDot", dest="fixDot", action="store_true", help="replace dots in cell meta IDs with dashes (for R)")
    #parser.add_option("-t", "--geneType", dest="geneType", help="type of gene IDs in expression matrix. values like 'symbols', or 'gencode22', 'gencode28' or 'gencode-m13'.")
    #parser.add_option("", "--bamDir", dest="bamDir", help="directory with BAM files, one per cell. Merges small BAM files into one per cell cluster.")
    parser.add_option("", "--clusterOrder", dest="clusterOrder", help="file with cluster names in the order that they should appear in the track. default is alphabetical order.")
    parser.add_option("-s", "--skipBarchart", dest="skipBarchart", help="do not create the bar chart graph", action="store_true")
    #parser.add_option("", "--name", dest="name", help="name of track hub.")
    #parser.add_option("", "--email", dest="email", help="contact email for track hub. Default is %default, taken from the env. variable CBEMAIL", default=CBEMAIL)
    #parser.add_option("-f", "--file", dest="file", action="store", help="run on file") 
    #parser.add_option("", "--test", dest="test", action="store_true", help="do something") 
    (options, args) = parser.parse_args()

    if not options.exprMatrix and not isfile(options.inConf) and not options.init:
        parser.print_help()
        exit(1)

    cellbrowser.setDebug(options)
    return args, options

def parseClustersFromMeta(metaFname, clusterFieldName, fixDot):
    " parse cluster -> cellId assignment from meta file, return as dict clusterName -> list of cellIds "
    logging.info("Parsing and using first field as the cell ID in file %s" % metaFname)
    clusterToCells = defaultdict(list)
    metaCellIds = set()
    skipCount = 0
    for row in cellbrowser.lineFileNextRow(metaFname):
        clusterName = row._asdict()[clusterFieldName]
        cellId = row[0]
        if fixDot:
            cellId = cellId.replace(".", "-")

        # skip all cells assigned to the "" cluster
        if clusterName=="":
            skipCount +=1
            continue

        if cellId in metaCellIds:
            errAbort("duplicated cell ID %s in meta data" % cellId)

        metaCellIds.add(cellId)
        clusterToCells[clusterName].append(cellId)

    logging.info("Got %d clusters and %d cell IDs assigned to them" % (len(clusterToCells), len(metaCellIds)))
    if skipCount!=0:
        logging.info("Skipped %d meta rows with empty cluster names" % skipCount)
    return clusterToCells

# ----------- main --------------

def writeHubStanza(tfh, hubName, db, email):
    "  write a single-file hub.txt "
    tfh.write("""hub scHub
shortLabel %s
longLabel %s
useOneFile on
email %s
descriptionUrl hub.txt

genome %s

""" % (hubName, hubName, email, db))

def writeDescPages(hubDir, hubName, hubUrl, refHtmlFname):
    saneHubName = sanitizeName(hubName)
    readDesc = join(hubDir, saneHubName+".html")
    readFh = open(readDesc, "w")
    hubBase = hubUrl.rsplit("/")[0]
    refHtml = ""
    if refHtmlFname is not None:
        refHtml = open(refHtmlFname).read()

    readFh.write("""
<h2>Description</h2>
<p>This track shows the alignment of sequencing reads to the genome, their coverage and the splice junctions in them. This can be helpful for quality checking of cluster markers, when designing in-situ probes or when trying to determine how specific a transcript is for a given cell type.

<h2>Display Conventions and Configuration</h2>

<p>This track has multiple "Views", which can be configured independently. See
the <a href="/goldenPath/help/multiView.html">documentation on Multi-View
tracks</a>. There are three different types of view sub tracks in this track:</a>

<p><b>BAM Reads:</b> Alignable regions are shown in black, unalignable regions
between two alignable ones shown with thin lines. For configuration options,
see <a href="https://genome.ucsc.edu/goldenpath/help/hgBamTrackHelp.html">the
BAM tracks help page</a>. By default, these tracks are hideen, set any of them to "squish" or "pack" (=clickable) on the track configuration page to see the reads. Click reads to show the alignment and read group.</p>

<p><b>Coverage:</b> bar graphs indicate the number of reads at this base pair.
You may want to switch on auto-scaling of the y axis. For configuration
options, see <a
href="https://genome.ucsc.edu/goldenpath/help/hgWiggleTrackHelp.html">the graph
tracks configuration help page</a>. These tracks are shown in "dense" by default, set any of the tracks to "full" to see the detailed coverage plot.</p>

<p><b>Splice Junctions:</b> thick rectangles show exons around a splice site,
connected by a line that indicates the intron. These gaps are shown and are
annotated with the number of reads, in the 'score' field. You can use the
'score' filter on the track configuration page to show only introns with a
certain number of supporting reads. The maximum number of reads that are shown
is 1000, even if more reads support an intron. These tracks are shown in dense by default, set this track to "pack" to see. Then click the splice junctions to see their score.</p>

<h2>Methods</h2>
<p>BAM files were provided by the data submitters, one (single end) or two files (paired end) per cell. The BAM alignments were used as submitted. They were merged with "samtools merge" into a single BAM file per cluster. The readgroup (RG) BAM tag indicates the original cell.</p>

<p>From the resulting merged BAM file, coverage was obtained using "wiggletools coverage" a tool written by Daniel Zerbino and the result was converted with the UCSC tool "wigToBigWig".</p>

<p>Also on the merged BAM file, the software IntronProspector was run with default settings. It retains reads with a gap longer than 70 bp and shorter than 500 kbp and merges them into annotated splice junctions.</p>

<h2>Data Access</h2>
<p>The merged BAM files, coverage bigWig files and splice junctions in bigBed format can be downloaded from the <a href='{hubBase}'>track hub directory</a>.</p>

<p>Since the splice junction .bigBed files have their scores capped at 1000, the original IntronProspector .bed files can also be downloaded from the <a href='{hubBase}'>track hub directory</a>. You can also find there *.calls.tsv files with more details about each junction, e.g. the number of uniquely mapping reads.</p>

<h2>Credits</h2>
<p>WiggleTools was written by Daniel Zerbino, IntronProspector was written by Mark Diekhans, track hubs were written to a large extent by Brian Raney.</p>

<h2>References</h2>
<p>
Zerbino DR, Johnson N, Juettemann T, Wilder SP, Flicek P.
<a href="https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btt737"
target="_blank">
WiggleTools: parallel processing of large collections of genome-wide datasets for visualization and
statistical analysis</a>.
<em>Bioinformatics</em>. 2014 Apr 1;30(7):1008-9.
PMID: <a href="https://www.ncbi.nlm.nih.gov/pubmed/24363377" target="_blank">24363377</a>; PMC: <a
href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3967112/" target="_blank">PMC3967112</a>
</p>

<p>
Mark Diekhans,
<a href="https://github.com/diekhans/intronProspector"
target="_blank">
IntronProspector GitHub Repository </a>.
<em>Github</em> 2018
</p>

<p>
{refHtml}
</p>""".format(**locals()))
    readFh.close()

def writeParentStanzas(tfh, saneHubName, hubName, cellCount):
    " write tdb entries for the composite parent tracks "
    tfh.write("""track %s\n""" % (saneHubName))
    tfh.write("""shortLabel %s\n""" % (hubName))
    tfh.write("""longLabel %s - %d cells\n""" % (hubName, cellCount))
    tfh.write("compositeTrack on\n")
    tfh.write("subGroup1 view Views reads=Reads cov=Coverage junc=Junctions\n")
    tfh.write("type bed 3\n")
    tfh.write("visibility dense\n")
    tfh.write("\n")

    tfh.write("""track %sViewReads\n""" % (saneHubName))
    tfh.write("""shortLabel Reads\n""")
    tfh.write("""longLabel %s Reads - %d cells\n""" % (hubName, cellCount))
    tfh.write("parent %s\n" % saneHubName)
    tfh.write("type bam\n")
    tfh.write("visibility hide\n")
    tfh.write("view reads\n")
    tfh.write("\n")

    tfh.write("""track %sViewCov\n""" % (saneHubName))
    tfh.write("""shortLabel Coverage\n""")
    tfh.write("""longLabel %s Coverage - %d cells\n""" % (hubName, cellCount))
    tfh.write("parent %s\n" % saneHubName)
    tfh.write("type bigWig\n")
    tfh.write("visibility dense\n")
    tfh.write("view cov\n")
    tfh.write("\n")

    tfh.write("""track %sViewJunc\n""" % (saneHubName))
    tfh.write("""shortLabel Splicing\n""")
    tfh.write("""longLabel %s Splice Junctions - %d cells\n""" % (hubName, cellCount))
    tfh.write("parent %s\n" % saneHubName)
    tfh.write("spectrum on\n")
    tfh.write("type bigBed 12\n")
    tfh.write("visibility dense\n")
    tfh.write("view junc\n")
    tfh.write("\n")

    #tfh.write("""track %sCov\n""" % (saneHubName))
    #tfh.write("""shortLabel %s Coverage\n""" % (hubName))
    #tfh.write("""longLabel %s Read Coverage - %d cells\n""" % (hubName, cellCount))
    #tfh.write("compositeTrack on\n")
    #tfh.write("type bigWig\n")
    #tfh.write("visibility dense\n")
    #tfh.write("\n")

    #tfh.write("""track %sJunc\n""" % (saneHubName))
    #tfh.write("""shortLabel %s Introns\n""" % (hubName))
    #tfh.write("""longLabel %s Intron Support - %d cells\n""" % (hubName, cellCount))
    #tfh.write("compositeTrack on\n")
    #tfh.write("type bigBed 12 +\n")
    #tfh.write("visibility dense\n")
    #tfh.write("\n")

def findProgram(name):
    " search PATH for binary and return full path "
    progPath = distutils.spawn.find_executable(name)
    if progPath==None:
        errAbort("Cannot find program '%s' in PATH." % name)
    return progPath

def writeJobScript(jobFn, cmds):
    " write commands to a shell script with file name "
    jobFh = open(jobFn, "w")
    logFname = splitext(jobFn)[0]+'.log'

    jobFh.write("#!/bin/bash\n")
    jobFh.write("set -e # = abort on error\n")
    jobFh.write("set -o pipefail # = even in pipes\n")
    jobFh.write("echo Job %s start, stdout/stderr go to %s\n" % (jobFn, logFname))
    jobFh.write("exec 1> %s\n" % logFname)
    jobFh.write("exec 2>&1\n");
    jobFh.write("echo PID=$$\n");
    for cmd in cmds:
        jobFh.write(cmd+"\n")
    jobFh.close()

    logging.info("Wrote job script %s" % jobFn)

def writeTdbStanzas(tfh, clusterName, saneHubName, cellCount):
    " write the track db statements for a one cluster "
    #saneClusterName = filter(str.isalnum, clusterName) # remove non alpha chars
    saneClusterName = sanitizeName(clusterName)

    tfh.write("""track %sReads\n""" % saneClusterName)
    tfh.write("""shortLabel %s BAM\n""" % clusterName)
    tfh.write("""longLabel %s sequencing reads from BAM (%d cells)\n""" % (clusterName, cellCount))
    tfh.write("type bam\n")
    tfh.write("view reads\n")
    tfh.write("parent %sViewReads\n" % saneHubName)
    tfh.write("subGroups view=reads\n")
    tfh.write("bigDataUrl %s\n" % (saneClusterName+".bam"))
    tfh.write("visibility dense\n")
    tfh.write("\n")

    tfh.write("""track %sCov\n""" % saneClusterName)
    tfh.write("""shortLabel %s Cov\n""" % (clusterName))
    tfh.write("""longLabel %s Coverage (%d cells)\n""" % (clusterName, cellCount))
    tfh.write("type bigWig\n")
    tfh.write("view cov\n")
    tfh.write("parent %sViewCov\n" % saneHubName)
    tfh.write("subGroups view=cov\n")
    tfh.write("autoScale on\n")
    tfh.write("bigDataUrl %s\n" % (saneClusterName+".bw"))
    tfh.write("visibility dense\n")
    tfh.write("\n")

    tfh.write("""track %sJunc\n""" % saneClusterName)
    tfh.write("""shortLabel %s Junc\n""" % (clusterName))
    tfh.write("""longLabel %s Splice Junctions (%d cells)\n""" % (clusterName, cellCount))
    tfh.write("type bigBed 12 +\n")
    tfh.write("spectrum on\n")
    tfh.write("scoreLabel Number of supporting reads (max: 1000)\n")
    tfh.write("parent %sViewJunc\n" % saneHubName)
    tfh.write("subGroups view=junc\n")
    tfh.write("bigDataUrl %s\n" % (saneClusterName+".junctions.bb"))
    tfh.write("visibility dense\n")
    tfh.write("\n")

def cellIdsForBam(bamPath):
    # we try different ways to map the bam file name to a cellID:
    # 1) remove the file ext, e.g. 123.3.bam -> cellId 123.3
    # 2) remove everything after the first dot, 123.3.bam -> cellId 123
    # 3 and 4) like 1) and 2) but with dashes replaced with dots (->R)
    basePath = basename(bamPath)
    cellId = splitext(basePath)[0]
    yield cellId.replace("_R1", "").replace("_R2", "")

    cellId2 = basePath.split(".")[0]
    yield cellId2.replace("_R1", "").replace("_R2", "")

    cellId3 = cellId.replace("-", ".")
    yield cellId3.replace("_R1", "").replace("_R2", "")

    cellId4 = cellId2.replace("-", ".")
    yield cellId4.replace("_R1", "").replace("_R2", "")

def findBams(bamDir, clusterToCells):
    """
    create dict with cluster -> (cellIds, bamFnames)
    """
    metaCellIds = set()
    for cluster, cellIds in clusterToCells.iteritems():
        for cellId in cellIds:
            metaCellIds.add(cellId)

    bamMask =join(bamDir, "*.bam")
    logging.info("Getting list of BAM files matching %s" % bamMask)
    bamPaths = glob.glob(bamMask)
    bamPaths.sort()
    logging.info("Found %d BAM files for %d meta cell IDs" % (len(bamPaths), len(metaCellIds)))

    # make map cellId -> bam files
    cellIdToBamFnames = defaultdict(list)
    fileCount = 0
    for bamPath in bamPaths:
        logging.debug("BAM file name: %s" % bamPath)

        cellId = None
        for tryCellId in cellIdsForBam(bamPath):
            logging.debug("Trying cell ID %s" % tryCellId)

            if tryCellId in metaCellIds:
                cellId = tryCellId
                break

            logging.debug("%s is not in meta data" % tryCellId)

        if cellId==None:
            logging.debug("Could not resolve %s to any cellId in meta" % (bamPath))
        else:
            logging.debug("Found file %s -> cellId=%s" % (bamPath, cellId))
            fileCount += 1
            cellIdToBamFnames[cellId].append(bamPath)

    clusterBams = dict()

    #logging.info("Got %s BAM files and %s cell ids in the meta data" % (len(metaCellIds), len(cellIdToBamFnames)))

    emptyClusterCount = 0
    for clusterName, cellIds in clusterToCells.iteritems():

        bamFnames = []
        for cellId in cellIds:
            cellBams = cellIdToBamFnames[cellId]
            bamFnames.extend(cellBams)

        if len(bamFnames)==0:
            logging.warn("No single BAM file found for cluster %s!" % clusterName)
            emptyClusterCount += 1
            continue
        elif len(bamFnames)==1:
            logging.warn("Only 1 BAM file for cluster %s?" % clusterName)

        logging.info("Cluster: %s, %d cell IDs in meta, found %d BAM files, example %s" % (clusterName, len(cellIds), len(bamFnames), bamFnames[0]))

        clusterBams[clusterName] = (cellIds, bamFnames)

    # now do some sanity checks
    missMeta = set(cellIdToBamFnames) - metaCellIds
    missBam = set(metaCellIds) - set(cellIdToBamFnames)
    allCellIds = set(cellIdToBamFnames).intersection(metaCellIds)
    logging.info("%d BAM cell ids have no meta data. Examples: %s" % (len(missMeta), " ".join(list(missMeta)[:10])))
    logging.info("%d meta cell ids have no BAM file. Examples: %s" % (len(missBam), " ".join(list(missBam)[:10])))

    return clusterBams, metaCellIds, cellIdToBamFnames

def writeDebugReport(allMetaCellIds, cellIdToBams, clusterBams, idReportFname):
    """ write a three-column table for debugging BAM file name <-> meta differences.
    Return number of cells with both at least one BAM file and meta data as a number.
    """
    logging.info("Writing %s to help track down BAM file problems" % idReportFname)
    cellCount = 0
    ofh = open(idReportFname, "w")
    ofh.write("#cellId\thasMeta\tcluster\tbamFname1\tbamFname2OrMore\n")

    allMetaCellIds = set(allMetaCellIds)

    allCellIds = set(allMetaCellIds).union(cellIdToBams)

    cellToCluster = {}
    for clusterName, (cellIds, bamFnames) in clusterBams.iteritems():
        for cellId in cellIds:
            cellToCluster[cellId] = clusterName

    for cellId in allCellIds:
        ofh.write(cellId+"\t")

        if cellId in allMetaCellIds:
            ofh.write("Yes\t")
        else:
            ofh.write("No\t")

        clusterName = cellToCluster.get(cellId, "!noClusterFound")
        ofh.write("%s\t" % clusterName)

        bamFnames = cellIdToBams[cellId]
        if len(bamFnames)==0:
            ofh.write("!noBamFound\t!noBamFound\n")
            continue

        if cellId in allMetaCellIds:
            cellCount += 1

        ofh.write(bamFnames[0])
        ofh.write("\t")

        if len(bamFnames)>1:
            ofh.write(",".join(bamFnames[1:]))
        ofh.write("\n")

    ofh.close()
    return cellCount


def mergeBams(hubName, db, tfh, bamDir, clusterToCells, outDir):
    logging.info("*** Merging and analyzing BAM files")

    clusterBams, allMetaCellIds, cellIdToBams = findBams(bamDir, clusterToCells)

    idReportFname = join(outDir, "metaBamMatch.txt")
    cellCount = writeDebugReport(allMetaCellIds, cellIdToBams, clusterBams, idReportFname)

    jlFh = open("jobList", "w")

    #cellCount = 0
    for clusterName, (cellIds, bamFnames) in clusterBams.iteritems():
        uniqueCellIds = set(cellIds)
        #cellCount += len(uniqueCellIds)

    logging.info("Merging BAM files and writing hub")
    saneHubName = sanitizeName(hubName)
    writeParentStanzas(tfh, saneHubName, hubName, cellCount)

    chromSizes = cellbrowser.getSizesFname(db)

    jobNo = 0
    emptyClusterCount = 0
    for clusterName, (cellIds, bamFnames) in clusterBams.iteritems():
        saneClusterName = sanitizeName(clusterName)
        logging.info("Processing cluster %s, %d cellIds/BAM files, examples: %s" % (clusterName, len(cellIds), cellIds[0]))
        cmds = []

        outBam = join(outDir, saneClusterName+".bam")
        outStat = join(outDir, saneClusterName+".stats.txt")
        outCalls = join(outDir, saneClusterName+".calls.tsv")
        junctionBed = join(outDir, saneClusterName+".junctions.bed")
        intronBed = join(outDir, saneClusterName+".introns.bed")

        junctionBedSorted = junctionBed.replace(".bed", ".sorted.bed")
        intronBedSorted = intronBed.replace(".bed", ".sorted.bed")

        junctionBb = junctionBed.replace(".bed", ".bb")
        intronBb = intronBed.replace(".bed", ".bb")

        cmd = "echo merging %d bam files for %d cell IDs" % (len(bamFnames), len(set(cellIds)))
        cmds.append(cmd)

        intronProspector = findProgram("intronProspector")
        samtools = findProgram("samtools")
        if len(bamFnames)==1:
            # samtools merge aborts if only one input BAM file
            cmd = "cat %s " % (bamFnames[0])
        else:
            # -r construct read groups
            # -f overwrite output files if needed
            cmd = samtools+" merge -f -r - %s " % (" ".join(bamFnames))

        cmd += "| tee %s | %s -p /dev/stdout -U -c %s -j %s -n %s " % \
            (outBam, intronProspector, outCalls, junctionBed, intronBed)
        cmd += "| samtools flagstat - > %s" % outStat
        cmds.append(cmd)
        
        #cmd = "echo running intronprospector on merged BAM file"
        #cmds.append(cmd)

        # XX could be part of the samtools merge command, in a pipe
        #cmds.append(cmd)

        cmd = "samtools index %s" % outBam
        cmds.append(cmd)

        cmd = """export LC_COLLATE=C; cat %s | awk '($5>1000) {$5=1000;} { OFS="\\t"; print}'  | sort -k1,1 -k2,2n > %s""" % (junctionBed, junctionBedSorted)
        cmds.append(cmd)

        cmd = """export LC_COLLATE=C; cat %s | awk '($5>1000) {$5=1000;} { OFS="\\t"; print}' | sort -k1,1 -k2,2n > %s""" % (intronBed, intronBedSorted)
        cmds.append(cmd)

        bigBed = findProgram("bedToBigBed")
        cmd = bigBed+" %s %s %s -type=bed6 -tab" % (intronBedSorted, chromSizes, intronBb)
        cmds.append(cmd)

        cmd = bigBed+" %s %s %s -type=bed12 -tab" % (junctionBedSorted, chromSizes, junctionBb)
        cmds.append(cmd)

        # alternative, 4x slower, silently drops ERCC chroms(?)
        # ~/bin/x86_64/bamToPsl in.bam stdout -chromAlias=mm10.chromAlias.tab -nohead | pslToBed stdin stdout | bedItemOverlapCount mm10 stdin -bed12 > hub/largeintestinegobletcell.bedgraph
        cmd = "echo creating coverage for merged BAM file"
        cmds.append(cmd)

        wiggletools = findProgram("wiggletools")
        outWig = join(outDir, saneClusterName+".wig")
        outBw = join(outDir, saneClusterName+".bw")
        # fix up the chrom field and the chrom=xxx stepSize fields
        # For wiggle files with Ensembl chrom names, we need UCSC chrom names
        cmd = wiggletools+ """ coverage %s | awk '($1 ~ /^[0-9XYM]/ && (NF==4)) { if ($1=="MT") {$1="M"}; $1="chr"$1 } ($2 ~ /^chrom=[0-9XYM]/) {split($2,a, "="); if (a[2]=="MT") {a[2]="M"}; $2="chrom=chr"a[2]} // {OFS=" "; print}' > %s""" % (outBam, outWig)
        cmds.append(cmd)

        wigToBigWig = findProgram("wigToBigWig")
        # -clip also means to not check the sequence names against chrom.sizes anymore
        # clip is needed so GL and KI seqs in hg38 don't mess up the conversion
        cmd = wigToBigWig+ " %s %s %s -clip" % (outWig, chromSizes, outBw)
        cmds.append(cmd)

        cmd = "rm -f %s %s %s" % (outWig, intronBedSorted, junctionBedSorted)
        cmds.append(cmd)

        #if False:
        if isfile(outBw):
            logging.info("Not running anything for %s, file %s already exists" % (saneClusterName, outBw))
        else:
            jobFn = "job-%d.sh" % jobNo
            writeJobScript(jobFn, cmds)
            jlFh.write("/bin/bash %s {check out exists %s}\n" % (jobFn, outBw))
        jobNo+=1

        writeTdbStanzas(tfh, clusterName, saneHubName, len(cellIds))

    if emptyClusterCount==len(clusterToCells):
        logging.error("No single BAM file could be linked to the meta data. Please check the identifiers in the meta data table against the file names in %s." % (bamDir))
        sys.exit(1)

    jlFh.close()
    logging.info("Wrote jobList file %s" % jlFh.name)

    logging.info("Wrote %s" % tfh.name)

def clamp(x):
    return max(0, min(int(x), 255))

def toHex(rgb):
    r = rgb[0]
    g = rgb[1]
    b = rgb[2]
    return "#{0:02x}{1:02x}{2:02x}".format(clamp(r*255), clamp(g*255), clamp(b*255))

def writeBarChartTdb(tfh, bbFname, clusterNames, unitName):
    " write the barChart tdb stanza "
    stepSize = 1.0 / len(clusterNames)
    colorCodes = []
    x = 0
    while x < 1.0:
        colorCodes.append(colorsys.hsv_to_rgb(x, 1.0, 1.0))
        x+=stepSize

    hexCodes = [toHex(x) for x in colorCodes]

    tfh.write('track barChart\n')
    tfh.write('type bigBarChart\n')
    tfh.write('visibility full\n')
    tfh.write('shortLabel Cluster expression\n')
    tfh.write('longLabel Median Cluster expression\n')

    # only remove spaces, this is more readable than sanitizeName
    saneClusterNames = []
    for cn in clusterNames:
        cn = cn.replace(" ", "_")
        #cn = ''.join(ch for ch in newName if (ch.isalnum() or ch=="_"))
        saneClusterNames.append(cn)

    tfh.write('barChartBars %s\n' % " ".join(saneClusterNames))
    tfh.write('barChartColors %s\n' % " ".join(hexCodes))
    tfh.write('barChartMetric median\n')
    tfh.write('barChartUnit %s\n' % unitName)
    tfh.write('barChartMatrixUrl exprMatrix.tsv\n')
    tfh.write('barChartSampleUrl clusters.tsv\n')
    tfh.write('defaultLabelFields name2\n')
    tfh.write('labelFields name2,name\n')
    tfh.write('bigDataUrl barChart.bb\n\n')

def to_camel_case(snake_str):
    components = snake_str.split('_')
    # We capitalize the first letter of each component except the first one
    # with the 'title' method and join them together.
    return components[0] + ''.join(x.title() for x in components[1:])

def sanitizeName(name):
    " remove all nonalpha chars and camel case name "
    newName = to_camel_case(name.replace(" ", "_"))
    newName = ''.join(ch for ch in newName if (ch.isalnum() or ch=="_"))
    return newName

def writeCatFile(clusterToCells, catFname):
    " write file with cellId<tab>clusterName "
    logging.info("Writing %s" % catFname)
    ofh = open(catFname, "w")
    for clusterName, cellIds in clusterToCells.iteritems():
        for cellId in cellIds:
            ofh.write("%s\t%s\n" % (cellId, clusterName))
    ofh.close()

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

        geneToSym = cellbrowser.readGeneSymbols({'geneIdType':defGenes})
        geneLocsId = cellbrowser.parseGeneLocs(defGenes)
        geneLocs = {}
        for geneId, locs in iterItems(geneLocsId):
            sym = geneToSym[geneId]
            geneLocs[sym] = locs
    else:
        geneToSym = cellbrowser.readGeneSymbols({'geneIdType':geneType})
        geneLocs = cellbrowser.parseGeneLocs(geneType)

    matOfh = open(outMatrixFname, "w")
    clustOfh = open(clusterFname, "w")

    mr = cellbrowser.MatrixTsvReader()
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
    cellbrowser.runCommand(cmd)

    # convert to .bb using .as file
    # from https://genome.ucsc.edu/goldenpath/help/examples/barChart/barChartBed.as
    #asFname = join(dataDir, )
    asFname = cellbrowser.getStaticFile(["genomes", "barChartBed.as"])
    sizesFname = cellbrowser.getSizesFname(genome)

    cmd = "bedToBigBed -as=%s -type=bed6+5 -tab %s %s %s" % (asFname, bedFname2, sizesFname, bbFname)
    cellbrowser.runCommand(cmd)

def buildTrackHub(db, inMatrixFname, metaFname, clusterFieldName, clusterOrderFile, hubName, bamDir, fixDot, geneType, unitName, email, refHtmlFname, hubUrl, skipBarchart, outDir):

    clusterToCells = parseClustersFromMeta(metaFname, clusterFieldName, fixDot)

    if clusterOrderFile is None:
        clusterOrder = list(sorted(clusterToCells.keys()))
        logging.info("No cluster order specified, using clusters in alphabetical order")
    else:
        logging.info("Reading cluster order from %s" % clusterOrderFile)
        clusterOrder = open(clusterOrderFile).read().splitlines()

    assert(set(clusterOrder)==set(clusterToCells)) # cluster order has to match actual cluster names

    tdbFname = join(outDir, "hub.txt")
    logging.info("Creating %s" % tdbFname)
    tfh = open(tdbFname, "w")

    writeHubStanza(tfh, hubName, db, email)

    matrixFname = join(outDir, "exprMatrix.tsv")

    bbFname = join(outDir, 'barChart.bb')
    catFname = join(outDir, 'clusters.tsv')
    if skipBarchart:
        logging.info("Not creating barChart file, got command line option")
    else:
        writeBarChartTdb(tfh, bbFname, clusterOrder, unitName)
        if isfile(bbFname):
            logging.info("Not creating barChart file, %s already exists" % bbFname)
        else:
            makeBarGraphBigBed(db, inMatrixFname, matrixFname, geneType, clusterToCells, clusterOrder, catFname, bbFname)

    if bamDir:
        mergeBams(hubName, db, tfh, bamDir, clusterToCells, outDir)

    tfh.close()

    writeDescPages(outDir, hubName, hubUrl, refHtmlFname)

def cbTrackHub(options):
    " make track hub given meta file and directory with bam files "
    if options.init:
        cellbrowser.copyPkgFile("sampleConfig/cellbrowser.conf")
        sys.exit(0)

    if isfile(options.inConf):
        conf = cellbrowser.loadConfig(options.inConf)

        db = conf["ucscDb"]
        inMatrixFname = conf["exprMatrix"]
        metaFname = conf["meta"]
        clusterFieldName = conf["clusterField"]
        clusterOrderFile = conf.get("clusterOrder")
        bamDir = conf.get("bamDir", "bam")
        fixDot = conf.get("fixDot", False)
        email = conf.get("hubEmail", CBEMAIL)
        geneType = conf["geneIdType"]
        outDir = conf["hubDir"]
        unitName = conf.get("unit", "TPM")
        hubUrl = conf.get("hubUrl", "")
        refHtmlFname = conf.get("refHtml", None)

        # use name, shortLabel or hubName from conf
        hubName = conf.get("hubName", conf.get("shortLabel", conf["name"]))

    if options.hubDir:
        outDir = options.hubDir
    if options.exprMatrix:
        inMatrixFname = options.exprMatrix
    if options.meta:
        metaFname = options.meta
    if options.clusterField:
        clusterFieldName = options.clusterField
    if options.clusterOrder:
        clusterOrderFile = options.clusterOrder
    if options.hubDir:
        outDir = options.hubDir
    skipBarchart = options.skipBarchart

    if not isdir(outDir):
        logging.info("Making %s" % outDir)
        os.makedirs(outDir)

    buildTrackHub(db, inMatrixFname, metaFname, clusterFieldName, clusterOrderFile, hubName, bamDir, fixDot, geneType, unitName, email, refHtmlFname, hubUrl, skipBarchart, outDir)

def cbHubCli():
    args, options = cbHub_parseArgs()
    cbTrackHub(options, *args)
