/*=======================================================================
    MATRIX: (C) Copyright 1992, Fred Hutchinson Cancer Research Center
      matrix.c reads two files, each containing a scoring matrix and
      compares them.  Assumes alphabet for a file is
      listed on first non-blank line after comments, and that
      the matrix, either half or full, follows on subsequent lines.
      Comment lines must appear at the beginning of the file & must
      have the character # or > or ; in the first column.

           matrix file1 file2

    Tries to infer the entropy of the matrix, but always assumes the
    matrix is scaled in half bits, so the result may be incorrect.
--------------------------------------------------------------------------
  1/12/92 J. Henikoff
  6/29/92 Added Jones and BLUS5065 frequencies.
  9/2/92  Changed to ignore lines starting with #, > or ;.
          Ignore * columns.
 10/17/92 Added comments.
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
void display();
void display2();
void entropy();

/*---- Global scoring matrix , order is :
  A R N D C Q E G H I L K M F P S T W Y V B Z J   -----------*/

/*------ Dayhoff AA frequencies ----------------------*/
double Dayhoff[20] = {0.087, 0.041, 0.040, 0.047, 0.033, 0.038, 0.050,
	0.089, 0.034, 0.037, 0.085, 0.081, 0.015, 0.040, 0.051, 0.070,
	0.058, 0.010, 0.030, 0.065};
/*---------Jones AA frequencies -------------------------*/
double Jones[20] =   {0.077, 0.051, 0.043, 0.052, 0.020, 0.041, 0.062,
        0.074, 0.023, 0.053, 0.091, 0.059, 0.024, 0.040, 0.051, 0.069,
        0.059, 0.014, 0.032, 0.066};
/*--------Blosum50 62% clustering marginal AA frequencies--------------*/
double Blosum[20] = {0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054,
	  0.074, 0.026, 0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057,
	  0.051, 0.013, 0.032, 0.073};
char Alphabet[23]={'A','R','N','D','C','Q','E','G','H','I','L','K',
		   'M','F','P','S','T','W','Y','V','B','Z','J'};
/*-----------  Dayhoff row order ----------------------*/
int DayRow[20] = {4,15,16,14,0,7,2,3,6,5,8,1,11,12,9,10,19,13,18,17};
/*=======================================================================*/
void main(argc, argv)
int argc;
char *argv[];
{
   char filename1[80], filename2[80];
   FILE *fin;
   int score1[AAS][AAS], score2[AAS][AAS], row, col;
   int score3[AAS][AAS];

   printf("\nMATRIX: (C) Copyright 1992, Fred Hutchinson Cancer");
   printf(" Research Center\n");
   if (argc > 1) strcpy(filename1, argv[1]);
   else
   {
      printf("\nEnter name of first file containing scoring matrix: " );
      gets(filename1);
   }
   if ((fin = fopen(filename1, "rt")) == NULL)
   {
      printf("\nCould not open %s\n", filename1);
      exit(-1);
   }
   else
      printf("Reading %s ...\n", filename1);
   read_file(fin, score1);
   fclose(fin);
   display(score1);
   printf("Dayhoff frequencies.\n"); entropy(score1, Dayhoff);
   printf("Jones frequencies.\n"); entropy(score1, Jones);
   printf("Blosum 62 frequencies.\n"); entropy(score1, Blosum);
   histogram(score1);

   if (argc > 2) strcpy(filename2, argv[2]);
   else
   {
      printf("\nEnter name of second file containing scoring matrix: " );
      gets(filename2);
   }
   if ((fin = fopen(filename2, "rt")) == NULL)
   {
      printf("\nCould not open %s\n", filename2);
      exit(-1);
   }
   else
      printf("\nReading %s ...\n", filename2);
   read_file(fin, score2);
   fclose(fin);
   display(score2);
   printf("Dayhoff frequencies.\n"); entropy(score2, Dayhoff);
   printf("Jones frequencies.\n"); entropy(score2, Jones);
   printf("Blosum 62 frequencies.\n"); entropy(score2, Blosum);
   histogram(score2);

   printf("\n%s below, %s above...", filename1, filename2);
   display2(score1, score2);	/* score1 below, score2 above */

/*--------------------Print half matrix  difference-----------------------*/
   printf("\nDifference %s-%s...\n", filename1, filename2);
   for (row=0; row<AAS; row++)
      for (col=0; col<AAS; col++)
	 score3[row][col] = score1[row][col] - score2[row][col];
   display(score3);
   histogram(score3);

   printf("\n%s below, %s-%s above...", filename1, filename1, filename2);
   display2(score1, score3);

   exit(0);
}
/*======================================================================*/
void read_file(fin, scores)
FILE *fin;
int scores[AAS][AAS];
{
   char line[132], *ptr;
   int alpha[AAS], nrows, ncols, row, col, i;

/*----------Read file until first non-blank line --------------*/
/* Skip comments at beginning of file - 1st char = #, > or ;   */
   line[0] = '\0';
   while ((strlen(line) < 1 || line[0]=='#' || line[0]=='>' || line[0]==';')
          && fgets(line, sizeof(line), fin) != NULL)
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
	 if (line[i]=='J' || line[i]=='O' || line[i]=='X') 
             col = 22;
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
      if (strlen(line) > 1 && nrows < AAS)
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
   for (row=0; row<AAS; row++)
   {
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
     }
   }

}  /* end of read_file */
/*========================================================================*/
void display(scores)
int scores[AAS][AAS];
{
   char ctemp[3];
   int row, col, i, j;
/*-------Print full matrix --------------------------------------------*/
   printf("   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  J\n");
   for (row=0; row<AAS; row++)
   {
     strncpy(ctemp, Alphabet+row, 1); ctemp[1] = '\0';
     printf("%1s ", ctemp);
     for (col=0; col<AAS; col++)
	printf("%2d ", scores[row][col]);
     printf("\n");
   }
   printf("   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  J\n");

/*--------------------Print half matrix ---------------------------------*/
   printf("\n   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  J\n");
   for (row=0; row<AAS; row++)
   {
     strncpy(ctemp, Alphabet+row, 1); ctemp[1] = '\0';
     printf("%1s ", ctemp);
     for (col=0; col<=row; col++)
	printf("%2d ", scores[row][col]);
     printf("\n");
   }
   printf("   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  J\n");
/*-------------Print half matrix in Dayhoff order -----------------------*/
   printf("\n   C  S  T  P  A  G  N  D  E  Q  H  R  K  M  I  L  V  F  Y  W\n");
   for (i=0; i<20; i++)
   {
      row = DayRow[i];
      strncpy(ctemp, Alphabet+row, 1); ctemp[1] = '\0';
      printf("%1s ", ctemp);
      for (j=0; j<=i; j++)
      {
	 col = DayRow[j];
	 printf("%2d ", scores[row][col]);
      }
      printf("\n");
   }
   printf("   C  S  T  P  A  G  N  D  E  Q  H  R  K  M  I  L  V  F  Y  W\n");
}  /* end of display */
/*=====================================================================*/
/*-------Altschul's entropy, inferred from Dayhoff frequencies -------*/
/*   Assumes matrix is in half bits, lambda = ln2, entropy is in bits
      qij = pi*pj*exp(lambda*sij), entropy=sum(qij*sij),
	      expected=sum(pi*pj*sij) 
      Since half bits: Use sij/2 for sij */
/*=======================================================================*/
void entropy(scores, freqs)
int scores[AAS][AAS];
double freqs[20];
{
   int row, col;
   double entropy, lambda, eij, dtemp2, expected, qij, sumqij, sij;

   printf("NOTE: Entropy calculations assume matrix is scaled in half bits\n");
   lambda = log(2.0);

   entropy = expected = sumqij = 0.0;
   for (row=0; row<20; row++)
   {
     for (col=0; col<=row; col++)
     {
	/* eij = pi*pj*sij/2          */
	/* dtemp2 = exp(lambda*sij/2)     */

	sij = (double) scores[row][col] / 2.0;
	eij = (double) freqs[row]*freqs[col]*sij;
	dtemp2 = (double) lambda * sij;
	dtemp2 = exp (dtemp2);
	qij = freqs[row]*freqs[col] * dtemp2;
	if (row==col) sumqij += qij;
	else sumqij += 2.0*qij;
	if (row==col)
	{
           entropy += qij * sij;
	   expected += eij;
	}
	else					/* matrix is symmetric */
	{
           entropy += 2.0 * qij * sij;
	   expected += 2.0*eij;
	}
     }
   }
   printf("sumqij=%.4f\n", sumqij);
   printf("Entropy, expected: = %.4f, %.4f\n", entropy, expected);
}
/*  end of entropy */
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
/*========================================================================*/
void display2(score1, score2)
int score1[AAS][AAS], score2[AAS][AAS];
{
   char ctemp[3];
   int row1, col1, row2, col2, i, j;

/*-------------Print both half matrices in Dayhoff order ----------------*/
   printf("\n   C  S  T  P  A  G  N  D  E  Q  H  R  K  M  I  L  V  F  Y  W\n");
   for (i=0; i<22; i++)
   {
      row1=row2= -1;
      if (i >=2 && i<20)
      {  row1=DayRow[i-2]; row2=DayRow[i]; }
      else if (i < 2) row2=DayRow[i];
      else if (i >= 20) row1=DayRow[i-2];
      if (row1 >= 0)
      {
	  strncpy(ctemp, Alphabet+row1, 1); ctemp[1] = '\0';
	  printf("%1s ", ctemp);
      }
      for (j=0; j<20; j++)
      {
	 if (row1 >= 0 && j <= i-2)		/* lower half matrix */
	 {
	    col1 = DayRow[j];
	    printf("%2d ", score1[row1][col1]);
	 }
	 if (j==i-2) printf("   ");			/* separation */
	 else if (i==0 && j==0) printf("  ");
	 else if (i==1 && j==0) printf("     ");
	 if (row2 >=0 && row2 < 20 && j >= i)		/* upper half */
	 {
	    col2 = DayRow[j];
	    printf("%2d ", score2[col2][row2]);
	 }
      }
      if (row2 >=0 && row2 < 20)
      {
	  strncpy(ctemp, Alphabet+row2, 1); ctemp[1] = '\0';
	  printf("%1s ", ctemp);
      }
      printf("\n");
   }
   printf("   C  S  T  P  A  G  N  D  E  Q  H  R  K  M  I  L  V  F  Y  W\n");
}  /* end of display2 */
