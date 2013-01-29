#!/bin/sh

echo Create fake genome
./scripts/createGenome.pl > fake.fasta 
cat fake.fasta | grep -v "^>" | tr -d "\n" > fake.txt

echo Create fake reads
cat fake.txt | ./scripts/createReadsPe.pl reads1.fastq reads2.fastq

echo Create index
bwa index fake.fasta

echo Align to SAI
bwa aln fake.fasta reads1.fastq > reads1.sai
bwa aln fake.fasta reads2.fastq > reads2.sai

echo Create SAM
bwa samse fake.fasta reads1.sai reads1.fastq > reads1.sam
bwa samse fake.fasta reads2.sai reads2.fastq > reads2.sam
                                                                       

