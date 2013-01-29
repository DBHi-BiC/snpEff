#!/bin/sh

BASE=reads_gC2T
REF=fake_C2T.fasta
BASE=reads_gG2A
REF=fake_G2A.fasta

BWA="$HOME/tools/bwa/bwa"
SAMTOOLS="$HOME/tools/samtools/samtools"
BCFTOOLS="$HOME/tools/samtools/bcftools/bcftools"

# Delelte old files
rm -vf $BASE.bam $BASE"_sort.bam" $BASE.txt.mpileup  $BASE.mpileup $BASE.bcf $BASE.vcf

echo Create BAM
$SAMTOOLS view -S -b $BASE.sam > $BASE.bam

echo Sort BAM
SORT=$BASE"_sort"
$SAMTOOLS sort $BASE.bam $SORT

echo Create mpileup and vcf
$SAMTOOLS mpileup -Bf $REF $SORT.bam > $BASE.txt.mpileup
$SAMTOOLS mpileup -ugBf $REF $SORT.bam > $BASE.mpileup
$BCFTOOLS view -cvg $BASE.mpileup > $BASE.vcf


