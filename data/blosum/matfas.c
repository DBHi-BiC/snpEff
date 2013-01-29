/*=======================================================================
 MATFAS. (C) Copyright 1992 Fred Hutchinson Cancer Research Center.
      matfas.c reads a file containing a scoring matrix and
      outputs a fasta matrix.  Assumes alphabet for file is
      listed on first non-blank, non-comment line.
      Lines beginning with # or > are assumed to be comment lines.
      Always sets gap penalties to -12, -4.

          matfas blosum62.iij >blosum62.fas

   4/5/92  J. Henikoff
  10/19/92 Added comments.
=========================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define NO 0
#define YES 1
#define AAS 23

void read_file();
void histogram();

/*---- Global scoring matrix , order is :
  A R N D C Q E G H I L K M F P S T W Y V B Z X   -----------*/

char Alphabet[23]={'A','R','N','D','C','Q','E','G','H','I','L','K',
		   'M','F','P','S','T','W','Y','V','B','Z','X'};
/*=======================================================================*/
void main(argc, argv)
int argc;
char *argv[];
{
   char filename[50], ctemp[3];
   FILE *fin;
   int score1[AAS][AAS], row, col;

   if (argc > 1) strcpy(filename, argv[1]);
   else
   {
      printf("\nEnter name of file containing scoring matrix: " );
      gets(filename);
   }
   if ((fin = fopen(filename, "rt")) == NULL)
   {
      printf("\nCould not open %s\n", filename);
      exit(-1);
   }
   read_file(fin, score1);
   fclose(fin);
/*--------------------Print fasta matrix ---------------------------------*/
   printf(";P %s\n", filename);
   printf(" 4 27 200 5 2 50 2\n");
   printf(" -12 -4\n");
   printf("@ *\n");
   printf("A R N D C Q E G H I L  K  M  F  P  S  T  W  Y  V  B Z X J\n");
   printf("0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 3 6 0 0\n");
   for (row=0; row<AAS; row++)
   {
     strncpy(ctemp, Alphabet+row, 1); ctemp[1] = '\0';
     for (col=0; col<=row; col++)
     {
	if ((row==22 || col==22) && score1[row][col]==-999)
	   score1[row][col] = -1;		/* X -> -1 */
	printf("%2d ", score1[row][col]);
     }
     printf("\n");
   }
   for (col=0; col < AAS; col++)
      printf("%2d ", score1[22][col]);
   printf("%2d\n", score1[22][22]);

   exit(0);
}
/*======================================================================*/
void read_file(fin, scores)
FILE *fin;
int scores[AAS][AAS];
{
   char line[132], *ptr, ctemp[3];
   int alpha[AAS], nrows, ncols, row, col, i;

/*----------Read file until first non-blank line --------------*/
   line[0] = '\0';
   while (!feof(fin) && fgets(line, sizeof(line), fin) != NULL &&
          (line[0] == '#' || line[0] == '>') )
	    ;
/*------See if the first line has characters on it ------------*/
   for (col=0; col < AAS; col++) alpha[col] = -1;
   if (strstr(line, "A") != NULL)	/* This line has characters */
   {
      row = 0;	/* # of alphabetic characters on the line */
      for (i=0; i<strlen(line); i++)
      {
	 col = -1;
	 if (line[i] == 'A') col = 0;
	 if (line[i] == 'R') col = 1;
	 if (line[i] == 'N') col = 2;
	 if (line[i] == 'D') col = 3;
	 if (line[i] == 'C') col = 4;
	 if (line[i] == 'Q') col = 5;
	 if (line[i] == 'E') col = 6;
	 if (line[i] == 'G') col = 7;
	 if (line[i] == 'H') col = 8;
	 if (line[i] == 'I') col = 9;
	 if (line[i] == 'L') col = 10;
	 if (line[i] == 'K') col = 11;
	 if (line[i] == 'M') col = 12;
	 if (line[i] == 'F') col = 13;
	 if (line[i] == 'P') col = 14;
	 if (line[i] == 'S') col = 15;
	 if (line[i] == 'T') col = 16;
	 if (line[i] == 'W') col = 17;
	 if (line[i] == 'Y') col = 18;
	 if (line[i] == 'V') col = 19;
	 if (line[i] == 'B') col = 20;
	 if (line[i] == 'Z') col = 21;
	 if (line[i] == 'J' || line[i] == 'O' || line[i] == 'X') col = 22;
	 if (col >= 0)
	 {
	    alpha[row] = col;
	    row++;
	 }
	 else if (isalpha(line[i])) row++;
      }
   }
/*-------Get the data values now ------------*/
   for (row=0; row<AAS; row++)
     for (col=0; col<AAS; col++)
	scores[row][col] = -999;		/* Null value */
   nrows = 0;
   line[0] = '\0';
   while (fgets(line, sizeof(line), fin) != NULL)
   {
      if (strlen(line) > 1)
      {
	 if (alpha[nrows] >= 0 && alpha[nrows] < AAS)
	 {
	    row = alpha[nrows]; ncols = 0;
	    ptr = strtok(line, " ,\n");
	    while (ptr != NULL)
	    {
	       if (strspn(ptr, "+-0123456789") == strlen(ptr))
	       {
		  col = alpha[ncols];
		  if (col >= 0 && col < AAS)
		     scores[row][col] = atoi(ptr);
		  ncols++;
	       }
	       ptr = strtok(NULL, " ,\n");
	    }
	 }
	 nrows++;
      }
   }

/*-------If some entries are still missing, assume symmetry ---------*/
/*   printf("   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X\n");*/
   for (row=0; row<AAS; row++)
   {
     strncpy(ctemp, Alphabet+row, 1); ctemp[1] = '\0';
/*     printf("%1s ", ctemp); */
     for (col=0; col<AAS; col++)
     {
	if (scores[row][col] == -999) scores[row][col] = scores[col][row];
	if (row==20 && scores[row][col]==-999)	/*  B -> D */
	   scores[row][col] = scores[3][col];
	if (row==21 && scores[row][col]==-999)	/* Z -> E */
	   scores[row][col] = scores[6][col];
	if (col==20 && scores[row][col]==-999)	/*  B -> D */
	   scores[row][col] = scores[row][3];
	if (col==21 && scores[row][col]==-999)	/* Z -> E */
	   scores[row][col] = scores[row][6];
/*	printf("%2d ", scores[row][col]);  blast */
     }
/*     printf("\n"); */
   }
/*   printf("   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X\n");*/
}  /* end of read_file */
/*======================================================================*/
void histogram(scores)
int scores[AAS][AAS];
{
   int row, col, n, tot, nbucket, bucket[30], median, ntot;
   int max=17, min=-8;

   nbucket = max-min+1;
   for (n=0; n<nbucket; n++)  bucket[n] = 0;
   tot = 0;

   for (row=0; row<20; row++)
     for (col=row; col<20; col++)		/* Half matrix */
     {
	n = scores[row][col] - min;
	if (n >= 0 && n <= max-min)
	{
	   bucket[n]++;  tot++;
	}
     }

   printf("\nHalf matrix histogram, number of scores=%d", tot);
   median=-1;  ntot=0;
   for (n=0; n<nbucket; n++)
   {
      printf("\n% 3d=%3d ", n+min, bucket[n]);
      for (col=0; col<bucket[n]; col++)
	 printf("X");
      ntot += bucket[n];
      if (median < 0 && ntot >= (int) tot/2) median=n;
   }
   printf("\nMedian score=%d", median+min);
}  /* end of histogram */
