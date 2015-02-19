# compare offtargets between different algorithms and against guideSeq
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

def gcContent(seq):
    c = 0
    for x in seq:
        if x in ["G","C"]:
            c+= 1
    return 100.0*(float(c)/len(seq))
            
# DOENCH SCORING 
params = [
# pasted/typed table from PDF and converted to zero-based positions
(1,'G',-0.2753771),(2,'A',-0.3238875),(2,'C',0.17212887),(3,'C',-0.1006662),
(4,'C',-0.2018029),(4,'G',0.24595663),(5,'A',0.03644004),(5,'C',0.09837684),
(6,'C',-0.7411813),(6,'G',-0.3932644),(11,'A',-0.466099),(14,'A',0.08537695),
(14,'C',-0.013814),(15,'A',0.27262051),(15,'C',-0.1190226),(15,'T',-0.2859442),
(16,'A',0.09745459),(16,'G',-0.1755462),(17,'C',-0.3457955),(17,'G',-0.6780964),
(18,'A',0.22508903),(18,'C',-0.5077941),(19,'G',-0.4173736),(19,'T',-0.054307),
(20,'G',0.37989937),(20,'T',-0.0907126),(21,'C',0.05782332),(21,'T',-0.5305673),
(22,'T',-0.8770074),(23,'C',-0.8762358),(23,'G',0.27891626),(23,'T',-0.4031022),
(24,'A',-0.0773007),(24,'C',0.28793562),(24,'T',-0.2216372),(27,'G',-0.6890167),
(27,'T',0.11787758),(28,'C',-0.1604453),(29,'G',0.38634258),(1,'GT',-0.6257787),
(4,'GC',0.30004332),(5,'AA',-0.8348362),(5,'TA',0.76062777),(6,'GG',-0.4908167),
(11,'GG',-1.5169074),(11,'TA',0.7092612),(11,'TC',0.49629861),(11,'TT',-0.5868739),
(12,'GG',-0.3345637),(13,'GA',0.76384993),(13,'GC',-0.5370252),(16,'TG',-0.7981461),
(18,'GG',-0.6668087),(18,'TC',0.35318325),(19,'CC',0.74807209),(19,'TG',-0.3672668),
(20,'AC',0.56820913),(20,'CG',0.32907207),(20,'GA',-0.8364568),(20,'GG',-0.7822076),
(21,'TC',-1.029693),(22,'CG',0.85619782),(22,'CT',-0.4632077),(23,'AA',-0.5794924),
(23,'AG',0.64907554),(24,'AG',-0.0773007),(24,'CG',0.28793562),(24,'TG',-0.2216372),
(26,'GT',0.11787758),(28,'GG',-0.69774)]

intercept =  0.59763615
gcHigh    = -0.1665878
gcLow     = -0.2026259

def calcDoenchScore(seq):
    assert(len(seq)==30)
    score = intercept

    guideSeq = seq[4:24]
    gcCount = guideSeq.count("G") + guideSeq.count("C")
    if gcCount <= 10:
        gcWeight = gcLow
    if gcCount > 10:
        gcWeight = gcHigh
    score += abs(10-gcCount)*gcWeight

    for pos, modelSeq, weight in params:
        subSeq = seq[pos:pos+len(modelSeq)]
        if subSeq==modelSeq:
            score += weight
    return 1.0/(1.0+math.exp(-score))

def rangeIntersection(start1, end1, start2, end2):
    """ return amount that two ranges intersect, <0 if no intersection 
    >>> rangeIntersection(1,10, 9, 20)
    1
    """
    s = max(start1,start2);
    e = min(end1,end2);
    return e-s;

def overlapSeqs(beds1, beds2):
    seqs2 = set()
    for b in beds2:
        seqs2.add(b[0])
    
    overlapBeds = []
    nonOverlapBeds = []
    for b in beds1:
        if b[0] in seqs2:
            overlapBeds.append(b)
        else:
            nonOverlapBeds.append(b)
    return overlapBeds, nonOverlapBeds
        

def overlapBedsCoords(beds1, beds2):
    """ given two lists of bed features, return only those from beds1 that overlap a match in beds2
    
    modified: USING SEQUENCE, not coordinates
    """
    #starts2 = set()
    #for b in beds2:
        #starts2.add(b[1])
    
    bedRet = []
    notOverlapBeds = []
    for b1 in beds1:
        #if b1[1] in starts2:
            #bedRet.append(b)
        found = False
        for b2 in beds2:
            if rangeIntersection(b1[1], b1[2], b2[1], b2[2])>0:
                bedRet.append(b1)
                found = True
                break
        if not found:
            notOverlapBeds.append(b1)
                
            #bedRet.append(b)
    return bedRet, notOverlapBeds
            
def filterBedScore(fts, minScore):
    " remove all bed fts if they dont' have at least minScore in their score column "
    newFts = []
    for b in fts:
        if b[4] > minScore:
            newFts.append(b)
    return newFts

def countMm(str1, str2):
    mmCount = 0
    for x,y in zip(str1, str2):
        if x!=y:
            mmCount += 1
    return mmCount

def parseGuides():
    " return guides as dict name -> seq "
    guides = {}
    guidesExt = {}
    for l in open("guides.txt"):
        name, seq, extSeq, pos = l.strip().split()
        guides[name] = seq
        guidesExt[name] = extSeq
    return guides, guidesExt
        
hitScoreM = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]

def calcHitScore(string1,string2):
    " see 'Scores of single hits' on http://crispr.mit.edu/about "
    # The Patrick Hsu weighting scheme
    #print string1, string2
    assert(len(string1)==len(string2)==20)

    dists = [] # distances between mismatches, for part 2
    mmCount = 0 # number of mismatches, for part 3
    lastMmPos = None # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(0, len(string1)):
        if string1[pos]!=string2[pos]:
            mmCount+=1
            if lastMmPos!=None:
                dists.append(pos-lastMmPos)
            score1 *= 1-hitScoreM[pos]
            lastMmPos = pos
    # 2nd part of the score
    if mmCount<2: # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists)/len(dists)
        score2 = 1.0 / (((19-avgDist)/19.0) * 4 + 1)
    # 3rd part of the score
    if mmCount==0: # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount**2)

    score = score1 * score2 * score3 * 100
    return score

def compSeq(str1, str2):
    " return a string that marks mismatches between str1 and str2 with * "
    s = []
    for x, y in zip(str1, str2):
        if x==y:
            s.append(".")
        else:
            s.append("*")
    return "".join(s)
            
def annotateOtData(beds, guideSeqs, guideCounts=None, guideExtSeqs=None):
    " add guide name, guide seq, mismatches, hitscore to tuples. Return a copy of beds. "
    newBeds = []
    for b in beds:
        b = list(b)
        otSeq, otScore, otCopyCount, guideName = b
        guideSeq = guideSeqs[guideName]
        b.append(guideSeq)
        b.append(gcContent(guideSeq))
        b.append(compSeq(guideSeq, otSeq))
        b.append(countMm(guideSeq, otSeq))
        b.append(calcHitScore(guideSeq[:20], otSeq[:20]))
        if guideCounts:
            b.append( float(otScore) / guideCounts[guideName] )
        if guideExtSeqs:
            b.append( calcDoenchScore(guideExtSeqs[guideName]) )
        newBeds.append(b)
    return newBeds

def writeTab(rows, fname):
    ofh = open(fname, "w")
    for b in rows:
        ofh.write("\t".join([str(x) for x in b]))
        ofh.write("\n")
    ofh.close()

def sumReads(gSeq):
    " return dict with name -> total number of reads obtained "
    ret = {}
    for name, beds in gSeq.iteritems():
        total = 0
        for bed in beds:
            score = bed[1]
            total+= score
        ret[name] = total
    return ret
    
def main():
    crisporDir = sys.argv[1]

    gSeq = parseBeds("guideSeq")
    mit = parseBeds("mit")
    crispor = parseBeds(crisporDir)
    #crispor2 = parseBeds("crispor2")
    guideSeqs, guideExtSeqs = parseGuides()

    headers = ["guideName", "guideSeq_otCount", "mit_otCount", "crispor_count", "guideSeq_mit_overlap", "guideSeq_crispor_overlap", "mit_share", "crispor_share"]
    print "\t".join(headers)
    totalCrispor = 0
    totalMit = 0
    totalGseq = 0
    #totalCrispor2 = 0
    allCrisporMissed = []
    allGuideSeqMissed = []
    for name, gSeqOts in gSeq.iteritems():
        #gSeqOts = filterBedScore(gSeqOts, 500)
        mitOts = mit[name]
        crisporOts = crispor.get(name, [])

        mitgSeq, mitMissed = overlapSeqs(gSeqOts, mitOts)
        crisporgSeq, crisporMissed = overlapSeqs(gSeqOts, crisporOts)
        _, guideSeqMissed = overlapSeqs(crisporOts, gSeqOts)

        mitShare = 100 * float(len(mitgSeq)) / len(gSeqOts)
        crisporShare = 100 * float(len(crisporgSeq)) / len(gSeqOts)

        allCrisporMissed.extend(crisporMissed)
        allGuideSeqMissed.extend(guideSeqMissed)

        row = [name, len(gSeqOts), len(mitOts), len(crisporOts), len(mitgSeq), len(crisporgSeq), "%0.2f" % mitShare, "%0.2f" % crisporShare]
        row = [str(x) for x in row]
        print "\t".join(row)

        totalCrispor += len(crisporgSeq)
        totalMit += len(mitgSeq)
        totalGseq += len(gSeqOts)

    print "total guide-seq %d, total found MIT: %0.2f%%, total found CRISPOR: %0.2f%%" % (totalGseq, 100.0 * totalMit / totalGseq, 100.0*totalCrispor/totalGseq)

    guideCounts = sumReads(gSeq)
    writeTab(annotateOtData(allCrisporMissed, guideSeqs, guideCounts, guideExtSeqs), "crisporMiss.tab")
    writeTab(annotateOtData(allGuideSeqMissed, guideSeqs, guideCounts, guideExtSeqs), "gSeqMissed.tab")
    #writeTab(gSeq.values(), "gSeqOfftargets.tab", guideSeqs, guideCounts, guideExtSeqs)

    # also annotate all guideSeq offtarget sequences
    allGuideSeqs = []
    for name, gSeqOts in gSeq.iteritems():
        annotBeds = annotateOtData(gSeqOts, guideSeqs, guideCounts, guideExtSeqs)
        allGuideSeqs.extend(annotBeds)
    writeTab(allGuideSeqs, "allGuideSeqOts.tab")

main()
