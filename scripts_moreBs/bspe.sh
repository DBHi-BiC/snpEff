#!/bin/sh

echo Create fake genome
./scripts/createGenome.pl > fake.fasta 
cat fake.fasta | grep -v "^>" | tr -d "\n" > fake.txt

echo Convert genome C2T
cat fake.fasta | tr "cC" "TT" > fake_C2T.fasta 

echo Convert genome G2A
cat fake.fasta | tr "gG" "AA" > fake_G2A.fasta 

echo Create fake reads
cat fake.txt | ./scripts/createReadsBsPe.pl reads1.fastq reads2.fastq

echo Create CT reads
cat reads1.fastq | tr "cC" "TT" > reads1_C2T.fastq
cat reads2.fastq | tr "cC" "TT" > reads2_C2T.fastq

echo Create GA reads
cat reads1.fastq | tr "gG" "AA" > reads1_G2A.fastq
cat reads2.fastq | tr "gG" "AA" > reads2_G2A.fastq

BWA="/home/pcingola/tools/bwa/bwa"
SAMTOOLS="/home/pcingola/tools/samtools/samtools"
BCFTOOLS="/home/pcingola/tools/samtools/bcftools/bcftools"

echo Create index
$BWA index fake_C2T.fasta
$BWA index fake_G2A.fasta

echo Align to SAI
$BWA aln fake_C2T.fasta reads1_C2T.fastq > reads1_C2T_gC2T.sai
$BWA aln fake_C2T.fasta reads1_G2A.fastq > reads1_G2A_gC2T.sai
$BWA aln fake_G2A.fasta reads1_C2T.fastq > reads1_C2T_gG2A.sai
$BWA aln fake_G2A.fasta reads1_G2A.fastq > reads1_G2A_gG2A.sai

$BWA aln fake_C2T.fasta reads2_C2T.fastq > reads2_C2T_gC2T.sai
$BWA aln fake_C2T.fasta reads2_G2A.fastq > reads2_G2A_gC2T.sai
$BWA aln fake_G2A.fasta reads2_C2T.fastq > reads2_C2T_gG2A.sai
$BWA aln fake_G2A.fasta reads2_G2A.fastq > reads2_G2A_gG2A.sai

#echo Create SAM
$BWA samse fake_C2T.fasta reads1_C2T_gC2T.sai reads1_C2T.fastq > reads1_C2T_gC2T.sam
$BWA samse fake_C2T.fasta reads1_G2A_gC2T.sai reads1_G2A.fastq > reads1_G2A_gC2T.sam
$BWA samse fake_G2A.fasta reads1_C2T_gG2A.sai reads1_C2T.fastq > reads1_C2T_gG2A.sam
$BWA samse fake_G2A.fasta reads1_G2A_gG2A.sai reads1_G2A.fastq > reads1_G2A_gG2A.sam
                                                                       
$BWA samse fake_C2T.fasta reads2_C2T_gC2T.sai reads2_C2T.fastq > reads2_C2T_gC2T.sam
$BWA samse fake_C2T.fasta reads2_G2A_gC2T.sai reads2_G2A.fastq > reads2_G2A_gC2T.sam
$BWA samse fake_G2A.fasta reads2_C2T_gG2A.sai reads2_C2T.fastq > reads2_C2T_gG2A.sam
$BWA samse fake_G2A.fasta reads2_G2A_gG2A.sai reads2_G2A.fastq > reads2_G2A_gG2A.sam

