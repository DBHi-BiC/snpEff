/*=======================================================================
 MATBLAS. (C) Copyright 1992 Fred Hutchinson Cancer Research Center.
      matblas.c reads a file containing a scoring matrix and
      outputs a blast matrix.  Assumes alphabet for file is
      listed on first non-blank line.

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
#define AAS 24
#define MISSING -999			/* Null matrix value */

void read_file();

/*---- Global scoring matrix , order is :
  A R N D C Q E G H I L K M F P S T W Y V B Z X *  -----------*/

char Alphabet[AAS]={'A','R','N','D','C','Q','E','G','H','I','L','K',
		   'M','F','P','S','T','W','Y','V','B','Z','X','*'};
/*=======================================================================*/
void main(argc, argv)
int argc;
char *argv[];
{
   char filename[50];
   FILE *fin;
   int score1[AAS][AAS];

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
   else
   {
      printf("#  Matrix made by matblas from %s\n", filename);
      printf("#  * column uses minimum score\n");
   }

   read_file(fin, score1);
   fclose(fin);
   exit(0);
}
/*======================================================================*/
void read_file(fin, scores)
FILE *fin;
int scores[AAS][AAS];
{
   char line[132], *ptr;
   int alpha[AAS], nrows, ncols, row, col, i;
   int minscore, maxscore;
   double x, xx, sumfreq, dtemp, dtemp1;

/*----------Read file until first non-blank, non-comment line -----------*/
   line[0] = '\0';
   while(!feof(fin) && fgets(line, sizeof(line), fin) != NULL &&
          (line[0]=='#' || line[0]=='>') )
             printf("%s", line);
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
	 if (line[i] == '*') col = 23;
	 if (col >= 0)
	 {
	    alpha[row] = col;
	    row++;
	 }
	 else if (isalpha(line[i])) row++;
      }
   }
   else
   {
      printf("\nFirst non-comment line must contain alphabet.\n");
      exit(-1);
   }
/*-------Get the data values now ------------*/
   for (row=0; row<AAS; row++)
     for (col=0; col<AAS; col++)
	scores[row][col] = MISSING;		/* Null value */
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
		  if (col >= 0 && col < AAS) scores[row][col] = atoi(ptr);
		  ncols++;
	       }
	       ptr = strtok(NULL, " ,\n");
	    }
	 }
	 nrows++;
      }
   }
/*-----------------Fill in symmetric entries if missing ------------*/
   minscore = 0 - MISSING;
   maxscore = MISSING;
   for (row=0; row<AAS; row++)
     for (col=0; col<AAS; col++)
     {
	if (scores[row][col] == MISSING)
		 scores[row][col] = scores[col][row];
	if (scores[row][col] != MISSING && scores[row][col] < minscore)
		 minscore=scores[row][col];
	if (scores[row][col] != MISSING && scores[row][col] > maxscore)
		 maxscore=scores[row][col];
     }

/*-------If some entries are still missing, assume symmetry ---------*/
   printf(" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n");
   for (row=0; row<AAS; row++)
   {
     for (col=0; col<AAS; col++)
     {
	if (scores[row][col] == MISSING) scores[row][col] = scores[col][row];
	if ((row==23 || col==23) && row != col && scores[row][col]==MISSING)
	   scores[row][col] = minscore;		/* * -> min score */
	if (scores[23][23]==MISSING) scores[23][23] = 1;
	printf("%2d ", scores[row][col]);
     }
    printf("\n");
   }
}  /* end of read_file */
