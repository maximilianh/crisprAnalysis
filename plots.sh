# scatter plot offtarget score versus read count 
echo 1
less allGuideSeqOts.bed | grep -v _site2 | grep -v sgRNA4 | gawk '($11<=4)' | cut -f9,12 > plots/lowGcMm4.tab
quick_plot plots/lowGcMm4.tab --out /var/www/crisporMax/downloads/guideSeq-score-readCount-Mm4lowGc --mode scatter --alpha 0.2 --title "guideSeq: readCount vs MIT offtarget score (MM<=4,low-GC guides)" --dot_size 3 --height 12 --logy --logx

echo 2
less allGuideSeqOts.bed | cut -f9,12 > plots/all.tab
quick_plot plots/all.tab --out /var/www/crisporMax/downloads/guideSeq-score-readCount --mode scatter --alpha 0.2 --title "guideSeq: readCount vs MIT offtarget score (all offtargets)" --dot_size 3 --height 12 --logy --logx

echo 3
less allGuideSeqOts.bed | grep "GG	" | grep -v _site2 | grep -v sgRNA4 | gawk '($11<=4)' | cut -f9,12 > plots/lowGcMm4NGG.tab
quick_plot plots/lowGcMm4NGG.tab --out /var/www/crisporMax/downloads/guideSeq-score-readCount-Mm4lowGcMainPam --mode scatter --alpha 0.2 --title "guideSeq: readCount vs MIT offtarget score (MM<=4, low-GC guides, only NGG)" --dot_size 3 --height 12 --logy --logx

echo 4
less allGuideSeqOts.bed | gawk '(substr($8,22,2)!="GG")' | grep -v _site2 | grep -v sgRNA4 | gawk '($11<=4)' | cut -f9,12 > plots/lowGcMm4NonNGG.tab
quick_plot plots/lowGcMm4NonNGG.tab --out /var/www/crisporMax/downloads/guideSeq-score-readCount-Mm4lowGcAltPam --mode scatter --alpha 0.2 --title "guideSeq: readCount vs MIT offtarget score (MM<=4, low-GC guides, non-NGG)" --dot_size 3 --height 12 --logy --logx

#cat gSeqMissed.tab | gawk '(substr($8,22,2)=="GG")' | gawk '($3<75 && $11<=4)' | cut -f9,12 | gawk  '{for (i=0;i<$1;i++) {print $2}}' > plots/missedScores.tab
#cat gSeqMissed.tab | gawk '($3<75 && $11<=4)' | cut -f12 > plots/crisporNotGuideSeq.tab
#cat allGuideSeqOts.tab | gawk '($3<75 && $11<=4)' | cut -f9,12 | gawk  '{for (i=0;i<$1;i++) {print $2}}' > plots/gSeqScores.tab
#cat allGuideSeqOts.tab | gawk '($3<75 && $11<=4)' | cut -f12 > plots/guideSeq.tab
cat allGuideSeqOts.tab | gawk '($6<75 && $8<=4) {print $9}' > plots/guideSeq.tab 
cat gSeqMissed.tab | gawk '($6<75 && $8<=4) {print $9}' > plots/crisporNotGuideSeq.tab
# weight scores by number of reads

cat allGuideSeqOts.tab | gawk '($6<75 && $8<=4) {for (i=0;i<=$2;i++) {print $9}}' > plots/guideSeqWeighted.tab
plotHist plots/guideSeq.tab plots/crisporNotGuideSeq.tab --out /var/www/crisporMax/downloads/histGseqVsMissZoom.pdf --max 10 --binCount 2000 --maxY 0.30
