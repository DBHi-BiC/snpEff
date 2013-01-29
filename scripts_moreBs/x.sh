#!/bin/sh

# Index
./scripts/moreBs.sh index -g fake.fasta -pathBwa $HOME/tools/bwa

# Map reads
./scripts/moreBs.sh map \
	-g fake.fasta \
	-pathBwa $HOME/tools/bwa \
	-pathSam $HOME/tools/samtools \
	-se reads1.fastq \
	-t 4

