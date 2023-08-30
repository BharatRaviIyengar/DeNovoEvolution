#!/bin/bash

set -x
set -e

echo -e "Annotating ncRNAs in the genome assemblies"

echo -e "0.1 Combining assemblies"
for prefix in AK SW GI YE DK UM ZB; do sed "/>/s/$/_$prefix/" "$prefix"_Genome.fa; done > allgenomes.fa

echo -e "0.2 BLASTing dmel-ncRNAs to genome assemblies"
makeblastdb -dbtype nucl -title allgenomes -in allgenomes.fa -out allgenomes
blastn -db allgenomes -query dmel-6.53-ncRNAs.fa -outfmt "6 qseqid sseqid score bitscore evalue sstart send pident positive length qlen qcovs qcovhsp" -qcov_hsp_perc 99 -perc_identity 90 -num_threads 5 -out dmel-6.53-ncRNAs_annaGenomes_BLASTN.tsv

echo -e "0.3 Converting BLAST results to GTF"
./consolidateBLAST.awk dmel-6.53-ncRNAs_annaGenomes_BLASTN.tsv

echo -e "0.4 Generating genome bed file"
for prefix in AK SW GI YE DK UM ZB
	do
	mawk 'BEGIN{FS=OFS="\t"} />/{h=substr($0,2);next} {l[h]+=length($0)} END{for(j in l) print j, l[j]}' "$prefix"_Genome.fa | sort -k1,1 > "$prefix"_genome.bed
done

for prefix in AK SW GI YE DK UM ZB
	do
	echo -e "PROCESSING "$prefix
	echo -e "1. Creating sorted GTF"
	#if [ ! -f $prefix"_sorted.gtf" ]
	#then
		cat $prefix"_GA.gtf" $prefix"_dmel-ncRNAs.gtf" | grep -v "#" | sort -k1,1 -k4,4n -k5,5n > $prefix"_sorted.gtf"
	#fi
	echo -e "2. Extracting intergenic regions"
	#if [ ! -f $prefix"_intergenic.bed" ]
	#then
		bedtools complement -i $prefix"_sorted.gtf" -g $prefix"_genome.bed" > $prefix"_intergenic.bed"
		bedtools getfasta -fi $prefix"_Genome.fa" -name -bed $prefix"_intergenic.bed" | sed 's/:://;s/:/_/' > $prefix"_intergenic.fa"
	#fi
	echo -e "3. getorf"
	getorf -sequence $prefix"_intergenic.fa" -outseq $prefix"_IGORF1.fa" -minsize 30 -find 3
	#echo -e "4. getorf to bed"
	#../../getorf2bed.awk $prefix"_IGORF1.fa" > $prefix"_IGORF2.bed"
	#echo -e "5. getfasta from bed"
	#bedtools getfasta -fi $prefix"_Genome.fa" -fo $prefix"_IGORF2.fa" -name -s -bed $prefix"_IGORF2.bed"
	#echo -e $prefix" = DONE \n ======= \n"
done


