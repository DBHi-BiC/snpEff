/*=======================================================================
     MATSPACE. Copyright (C) 1992, Fred Hutchinson Cancer Research Center

      matspace.c reads a file containing a scoring matrix with the
      alphabet listed on the first non-blank line
      and outputs a matrix where each element is subtracted from the
      maximum element for input to Des Higgins' "spacer" program.
      Note: spacer requires a separate "symbs" file containing the
      alphabet for the output from this program which should contain
      the single line:
	   A R N D C Q E G H I L K M F P S T W Y V

   5/26/92  J. Henikoff
  10/19/92  Added comments.
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

/*---- Global scoring matrix , order is :
  A R N D C Q E G H I L K M F P S T W Y V B Z X   -----------*/

char Alphabet[23]={'A','R','N','D','C','Q','E','G','H','I','L','K',
		   'M','F','P','S','T','W','Y','V','B','Z','X'};
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
/*       printf("Reading %s ...\n", filename);*/

   read_file(fin, score1);
   fclose(fin);
   exit(0);
}
/*======================================================================*/
void read_file(fin, scores)
FILE *fin;
int scores[AAS][AAS];
{
   char line[132], *ptr, ctemp[3];
   int alpha[AAS], nrows, ncols, row, col, i, maxval;

   maxval = -999;
/*----------Read file until first non-blank, non-comment line --------*/
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
                  if (scores[row][col] > maxval) maxval = scores[row][col];
		  ncols++;
	       }
	       ptr = strtok(NULL, " ,\n");
	    }
	 }
	 nrows++;
      }
   }

/*-------If some entries are still missing, assume symmetry ---------*/
   printf("20\n");
   for (row=0; row<20; row++)
   {
     strncpy(ctemp, Alphabet+row, 1); ctemp[1] = '\0';
     for (col=0; col<=row; col++)
     {
	if (scores[row][col] == -999) scores[row][col] = scores[col][row];
	printf("%2d ", maxval-scores[row][col]);
     }
     printf("\n");
   }
}  /* end of read_file */
