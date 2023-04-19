#!/bin/bash

set -e
set -x

getorf -sequence novel.fa -outseq novel_ORF.fa -find 3 -noreverse

./getorf2genomicGFF.awk novel_ORF.fa > novel_ORF.gff

./GTFoverlap2.awk CDS.gff novel_ORF.gff > novel_ORF_frames.txt

./framecounts.awk novel_ORF_frames.txt > framecounts.txt

bedtools intersect -a CDS.gff -b novel.gff > overlappingRNAregion.gff

mawk 'BEGIN{FS = OFS ="\t"} $10==100{$3 = "frame-"$9; print $0}' novel_ORF_frames.txt | bedtools getfasta -fi GCF_000146045.2_R64_genomic.fna -bed /dev/stdin -s -name > antisense_ORFs.fa

