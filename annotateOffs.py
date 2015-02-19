# annotate guideseq offtargets with all possible scores we have
import glob, copy, sys, math
from collections import defaultdict
from os.path import basename

def parseBeds(dirName):
    " parse all beds in dir and return as baseFname -> set of (seq,score, count,guideName)"
    ret = {}
    for fname in glob.glob(dirName+"/*.bed"):
        counts = defaultdict(int)
        for line in open(fname):
            fs = line.rstrip("\n").split("\t")
            seq = fs[3]
            score = float(fs[4])
            counts[(seq,score)]+=1

        guideName = basename(fname).split('.')[0]
        seqs = set()
        for data, count in counts.iteritems():
            seq, score = data
            seqs.add( (seq, score, count, guideName) )

        ret[guideName]= seqs
    return ret

def parseFlanks(faDir):
    seqToFlank = {}
    for fname in glob.glob(join(faDir, "*.fa")):
        for shortSeq, seq in parseFasta(open(fname)).iteritems():
            assert(shortSeq not in seqToFlank)
            seqToFlank[shortSeq] = seq
    return seqToFlank

def main():
    gSeq = parseBeds("guideSeq")
    flanks = parseFlanks("guideSeq/flankSeq/")

main()
