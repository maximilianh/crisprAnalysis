guides.txt: the guide sequences, their 30mer extension and the position in the genome


guidesGc.tab: the list of guides with their offtarget counts and GC content
allGuideSeqOts.tab: list of off-targets with sequence len=20 from the guideSeq paper, annotated with
various measures:

fields:
1 - offtarget sequence
2 - readcount (for guide-seq data) or off-target score
3 - number of copies of this sequence in input file
4 - name of guide sequence
5 - guide sequence
6 - GC content of guide sequence
7 - mismatch locations relative to guide seq
8 - mismatch count relative to guide seq
9 - MIT offtarget score of offtarget seq rel. to guide seq.
10- readcount share relative to all reads of this guide
11- doench score of extended guide

