#!/bin/sh

DIR="$HOME/1k_genomes"

#---
# SNP analysis
#---

#time java -jar snpEff.jar -v -stats $DIR/1000_Genomes_snpEff_summary.html -vcf4 hg37.60 $DIR/all.vcf > /dev/null

#time java -jar snpEff.jar \
#	-v \
#	-stats $DIR/1000_Genomes_snpEff_summary.html \
#	-vcf4 \
#	hg37.60 $DIR/all.vcf \
#	| gzip -c > $DIR/1000_Genomes_snpEff.txt.gz

#---
# Inde analysis
#---

#time java -jar snpEff.jar \
#	-v \
#	-vcf4 \
#	-stats $DIR/1000_Genomes_snpEff_summary_indels.html \
#	-upDownStreamLen 0 \
#	-no-downstream -no-intergenic -no-intron -no-upstream -no-utr \
#	hg37.60 $DIR/all_indels.vcf \
#	| gzip -c > $DIR/1000_Genomes_snpEff_indels.txt.gz




time java -jar snpEff.jar -v -vcf4 -stats $DIR/1000_Genomes_snpEff_summary_indels.html -upDownStreamLen 0 -no-downstream -no-intergenic -no-intron -no-upstream -no-utr hg37.60 $DIR/all_indels.vcf > /dev/null
