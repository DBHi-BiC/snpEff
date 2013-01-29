#!/bin/sh

WGSIM="$HOME/tools/samtools/misc/wgsim"

echo Create fake genome
./scripts/createGenome.pl > fake.fasta 
cat fake.fasta | grep -v "^>" | tr -d "\n" > fake.txt
cat fake.fasta | tr "cC" "TC" > fake_CT.fasta
cat fake.fasta | tr "gG" "AG" > fake_GA.fasta


echo Create fake reads
#cat fake.txt | ./scripts/createReadsBsPe.pl reads1.fastq reads2.fastq
#cat fake.txt | ./scripts/createReadsPe.pl reads1_nonBs.fastq reads2_nonBs.fastq

$WGSIM -d 300 -N 5000 fake_CT.fasta rct1 rct2
$WGSIM -d 300 -N 5000 fake_GA.fasta rga1 rga2
cat rct1 rga1 > reads1.fastq
cat rct2 rga2 > reads2.fastq
rm -vf rct1 rct2 rga1 rga2
