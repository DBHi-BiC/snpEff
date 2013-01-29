#!/bin/sh

echo Create fake genome
./scripts/createGenome.pl > fake.fasta 
cat fake.fasta | grep -v "^>" | tr -d "\n" > fake.txt

echo Create fake reads
cat fake.txt | ./scripts/createReads.pl > reads.fastq

echo Create index
bwa index fake.fasta

echo Align to SAI
bwa aln fake.fasta reads.fastq > reads.sai

echo Create SAM
bwa samse fake.fasta reads.sai reads.fastq > reads.sam

