#!/bin/sh

echo Create fake genome
./scripts/createGenome.pl > fake.fasta 
cat fake.fasta | grep -v "^>" | tr -d "\n" > fake.txt

echo Convert genome C2T
cat fake.fasta | tr "cC" "TT" > fake_C2T.fasta 

echo Convert genome G2A
cat fake.fasta | tr "gG" "AA" > fake_G2A.fasta 

echo Create fake reads
cat fake.txt | ./scripts/createReadsBs.pl > reads.fastq

echo Create CT reads
cat reads.fastq | tr "cC" "TT" > reads_C2T.fastq

echo Create GA reads
cat reads.fastq | tr "gG" "AA" > reads_G2A.fastq

BWA="/home/pcingola/tools/bwa/bwa"
SAMTOOLS="/home/pcingola/tools/samtools/samtools"
BCFTOOLS="/home/pcingola/tools/samtools/bcftools/bcftools"

echo Create index
$BWA index fake_C2T.fasta
$BWA index fake_G2A.fasta

echo Align to SAI
$BWA aln fake_C2T.fasta reads_C2T.fastq > reads_C2T_gC2T.sai
$BWA aln fake_C2T.fasta reads_G2A.fastq > reads_G2A_gC2T.sai
$BWA aln fake_G2A.fasta reads_C2T.fastq > reads_C2T_gG2A.sai
$BWA aln fake_G2A.fasta reads_G2A.fastq > reads_G2A_gG2A.sai

echo Create SAM
$BWA samse fake_C2T.fasta reads_C2T_gC2T.sai reads_C2T.fastq > reads_C2T_gC2T.sam
$BWA samse fake_C2T.fasta reads_G2A_gC2T.sai reads_G2A.fastq > reads_G2A_gC2T.sam
$BWA samse fake_G2A.fasta reads_C2T_gG2A.sai reads_C2T.fastq > reads_C2T_gG2A.sam
$BWA samse fake_G2A.fasta reads_G2A_gG2A.sai reads_G2A.fastq > reads_G2A_gG2A.sam

#echo Create BAM
#$SAMTOOLS view -S -b fake_C2T.sam > fake_C2T.bam
#
#echo Sort BAM
#$SAMTOOLS sort fake_C2T.bam fake_C2T_sort
#
#echo Create mpileup and vcf
#$SAMTOOLS mpileup -f fake_C2T.fasta fake_C2T_sort.bam > fake_C2T.mpileup
#$SAMTOOLS mpileup -ugf fake_C2T.fasta fake_C2T_sort.bam | $BCFTOOLS view - > fake_C2T.raw.vcf 

