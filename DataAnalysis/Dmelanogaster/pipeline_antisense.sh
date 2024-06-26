#!/bin/bash

set -e
set -x

i=$1
#echo -e "Generating Genome BLAST database" 
#makeblastdb -dbtype nucl -title $i -out $i"_blast" -in $i"_Genome.fa"; done
###
#echo -e "\n Obtaining antisense RNA coordinates"
#./antisenseRNA.awk "$i"_GA.gtf "$i"_transcriptome.gtf | sort -k1,1 -k4,4g > "$i"_antisenseRNA.gtf
####
echo -e "\n Extracting antisense RNA sequences"
gffread -g "$i"_Genome.fa -w "$i"_antisenseRNA.fa "$i"_antisenseRNA.gtf
####
echo -e "\n Extracting ORFs from antisense RNAs"
getorf -sequence "$i"_antisenseRNA.fa -outseq "$i"_antisense_ORFs.fa -find 3 -noreverse
./removespuriousORFs.awk "$i"_antisenseRNA.fa "$i"_antisense_ORFs.fa > temporfs"$i"
mv temporfs"$i" "$i"_antisense_ORFs.fa
####
echo -e "\n Mapping the coordinates of ORFs from antisense RNAs"
blastn -query "$i"_antisense_ORFs.fa -db "$i"_blast -outfmt "6 qseqid sseqid sstart send qcovs qcovhsp pident" -qcov_hsp_perc 100 -perc_identity 100 | ./blast2gtf.awk  > "$i"_antisense_ORFs_blastn.gtf
####
echo -e "\n Extracting true antisense ORF coordinates"
./GTFoverlap2.awk "$i"_GA.gtf "$i"_antisense_ORFs_blastn.gtf > "$i"_realasORFs.gtf
####
echo -e "\n Calculating frame statistics"
./framecounts.awk "$i"_realasORFs.gtf > "$i"_framecounts.txt
###
echo -e "\n Extracting true antisense ORF sequences"
bedtools getfasta -fi "$i"_Genome.fa -bed "$i"_realasORFs.gtf -s -name > "$i"_realasORFs.fa
###

echo -e "\n Obtaining antisense exon coordinates"
./antisenseExons.awk "$i"_GA.gtf "$i"_transcriptome.gtf > "$i"_antisenseExons.gtf
####
echo -e "\n Extracting antisense exon sequences (overlapping regions)"
bedtools getfasta -fi "$i"_Genome.fa -bed "$i"_antisenseExons.gtf -name -s -tab > "$i"_antisenseExons.tsq
####

mawkt '$6==100' "$i"_realasORFs.gtf | bedtools intersect -a "$i"_GA.gtf -b /dev/stdin > "$i"_OLR_100.gtf

./Exp_vs_Obs_dmel-asORF.awk "$i"_antisenseExons.tsq "$i"_realasORFs_100.tsq > asORFstats.txt