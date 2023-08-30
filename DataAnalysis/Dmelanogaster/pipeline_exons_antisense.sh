#!/bin/bash

set -e
set -x

i=$1

echo -e "\n Obtaining antisense RNA coordinates"
./antisenseExons.awk "$i"_GA.gtf "$i"_transcriptome.gtf > "$i"_antisenseExons.gtf
####
echo -e "\n Extracting antisense RNA sequences"
bedtools getfasta -fi "$i"_Genome.fa -bed "$i"_antisenseExons.gtf -name -s > "$i"_antisenseExons.fa
####
echo -e "\n Extracting ORFs from antisense exons"
getorf -sequence "$i"_antisenseExons.fa -outseq "$i"_antisense_exonic_ORF.fa -find 3 -noreverse -sformat pearson
####
echo -e "\n Extracting antisense ORF coordinates and frame"
getorf2bed.awk "$i"_antisense_exonic_ORF.fa > "$i"_antisense_exonic_ORF.gtf
./GTFoverlap3.awk "$i"_GA.gtf "$i"_antisense_exonic_ORF.gtf > "$i"_antisense_exonic_ORF2.gtf
####
echo -e "\n Extracting true antisense ORF sequences"
bedtools getfasta -fi "$i"_Genome.fa -bed "$i"_antisense_exonic_ORF2.gtf -s -name > "$i"_realasORFs.fa "$i"_antisense_exonic_ORF.fa
###
echo -e "\n Calculating frame statistics"
###
./framecounts.awk "$i"_realasORFs.gtf > "$i"_antisense_exonic_ORF2.gtf
###
mawk '/>/{f=substr($0,8,1);next} {for(i=1;i<=length($0)-30;i+=3){if(substr($0,i,3)=="ATG") x[f]++}} END{for(i in x) print substr(FILENAME,1,2), i,x[i]}' "$i"_antisense_exonic_ORF.fa




