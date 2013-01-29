# BLOSUM Makefile for UNIX BSD with Standard C
#   Type  "make blosums" to compile all blosum programs
#   10/17/92
#
CC = cc
CFLAGS = -c
CL = cc
LFLAGS = -lm -o
#
blosums:	blosum matrix matalin matblas matfas matspace
#
#				BLOSUM programs
#	blosum: Make clustered blosum matrices from a blocks database
blosum:		blosum.o motmisc.o
		$(CL) blosum.o motmisc.o $(LFLAGS) blosum
#	matrix: Compare two matrices
matrix:		matrix.o motmisc.o
		$(CL) matrix.o motmisc.o $(LFLAGS) matrix
#	matalin: Reformat matrix for multalin
matalin:	matalin.o motmisc.o
		$(CL) matalin.o motmisc.o $(LFLAGS) matalin
#	matblas: Reformat matrix for blast
matblas:	matblas.o motmisc.o
		$(CL) matblas.o motmisc.o $(LFLAGS) matblas
#	matfas:  Reformat matrix for fasta/ssearch
matfas:		matfas.o motmisc.o
		$(CL) matfas.o motmisc.o $(LFLAGS) matfas
#	matspace: Reformat matrix for spacer
matspace:	matspace.o motmisc.o
		$(CL) matspace.o motmisc.o $(LFLAGS) matspace
#
#
blosum.o:	blosum.c motifj.h
		$(CC) $(CFLAGS) blosum.c
matrix.o:	matrix.c motifj.h
		$(CC) $(CFLAGS) matrix.c
matalin.o:	matalin.c motifj.h
		$(CC) $(CFLAGS) matalin.c
matblas.o:	matblas.c motifj.h
		$(CC) $(CFLAGS) matblas.c
matfas.o:	matfas.c motifj.h
		$(CC) $(CFLAGS) matfas.c
matspace.o:	matspace.c motifj.h
		$(CC) $(CFLAGS) matspace.c
#
motmisc.o:	motmisc.c motifj.h
		$(CC) $(CFLAGS) motmisc.c
#
rmos:
	rm blosum.o
	rm matrix.o
	rm matalin.o
	rm matblas.o
	rm matfas.o
	rm matspace.o
