/*=======================================================================
(C) Copyright 1991, Fred Hutchinson Cancer Research Center
     motmisc.c    Miscellaneous PROTOMAT routines.
-------------------------------------------------------------------------
   6/20/91
   6/27/91  Removed printfs from getscore().
   8/26/91  Made Score & HighPass part of score structure & removed
	    external variables.
   9/5/91   Moved init_dbs() from uextract
   9/12/91  Added kr_itoa()
   9/23/91  Added split_names()
   11/20/91 Added B, Z and X to aachar_to_num().
   >>>>>>>>>>>>>>>>>>>>>>>>   4.0   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   1/27/92  Added dir_dos() & dir_unix().
	    Added makedbid(), get_ids().
   2/10/92  Added extract_seqs() from uextract/excluded.
	    Modified makedbid() & get_ids() for modified struct db_id in
	    motifj.h.
   5/12/92  Corrected problem not recognizing ":" in split_names() for DOS.
	    Corrected problem in dir_dos() & dir_unix() if no directory
	    line in .lis file; changed type of argument for these.
   >>>>>>>>>>>>>>>>>>>>>>>>>>>   5.0   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   8/17/92  Removed $ from title and * from end of sequence in extract_seqs().
  10/12/92  Added pvalue field to struct db_id.
=========================================================================*/
/*DOS #include <ctype.h>
#include <io.h>
#include <dir.h>
#include <process.h>
DOS*/
#include <sys/types.h>
#include <sys/dir.h>
#include "motifj.h"

/*---- Global scoring matrix , order is :
		  A R N D C Q E G H I L K M F P S T W Y V J   -----------*/

/* PAM250 matrix with each cell offset by 8 to make all scores non-negative */
char pam250_matrix[20][20]={
10, 6, 8, 8, 6, 8, 8, 9, 7, 7, 6, 7, 7, 4, 9, 9, 9, 2, 5, 8, /* A */
 6,14, 8, 7, 4, 9, 7, 5,10, 6, 5,11, 8, 4, 8, 8, 7,10, 4, 6, /* R */
 8, 8,10,10, 4, 9, 9, 8,10, 6, 5, 9, 6, 4, 7, 9, 8, 4, 6, 6, /* N */
 8, 7,10,12, 3,10,11, 9, 9, 6, 4, 8, 5, 2, 7, 8, 8, 1, 4, 6, /* D */
 6, 4, 4, 3,20, 3, 3, 5, 5, 6, 2, 3, 3, 4, 5, 8, 6, 0, 8, 6, /* C */
 8, 9, 9,10, 3,12,10, 7,11, 6, 6, 9, 7, 3, 8, 7, 7, 3, 4, 6, /* Q */
 8, 7, 9,11, 3,10,12, 8, 9, 6, 5, 8, 6, 3, 7, 8, 8, 1, 4, 6, /* E */
 9, 5, 8, 9, 5, 7, 8,13, 6, 5, 4, 6, 5, 3, 7, 9, 8, 1, 3, 7, /* G */
 7,10,10, 9, 5,11, 9, 6,14, 6, 6, 8, 6, 6, 8, 7, 7, 5, 8, 6, /* H */
 7, 6, 6, 6, 6, 6, 6, 5, 6,13,10, 6,10, 9, 6, 7, 8, 3, 7,12, /* I */
 6, 5, 5, 4, 2, 6, 5, 4, 6,10,14, 5,12,10, 5, 5, 6, 6, 7,10, /* L */
 7,11, 9, 8, 3, 9, 8, 6, 8, 6, 5,13, 8, 3, 7, 8, 8, 5, 4, 6, /* K */
 7, 8, 6, 5, 3, 7, 6, 5, 6,10,12, 8,14, 8, 6 ,6 ,7, 4, 6,10, /* M */
 4, 4, 4, 2, 4, 3, 3, 3, 6, 9,10, 3, 8,17, 3, 5 ,5 ,8,15, 7, /* F */
 9, 8, 7, 7, 5, 8, 7, 7, 8, 6, 5, 7, 6, 3,14, 9, 8, 2, 3, 7, /* P */
 9, 8, 9, 8, 8, 7, 8, 9, 7 ,7, 5, 8, 6, 5, 9,10, 9, 6, 5, 7, /* S */
 9, 7, 8, 8, 6, 7, 8, 8, 7, 8, 6, 8, 7, 5, 8, 9,11, 3, 5, 8, /* T */
 2,10, 4, 1, 0, 3, 1, 1, 5, 3, 6 ,5, 4, 8, 2, 6, 3,25, 8, 2, /* W */
 5, 4, 6, 4, 8, 4, 4, 3, 8, 7, 7, 4, 6,15, 3, 5, 5, 8,18, 6, /* Y */
 8, 6, 6, 6, 6, 6, 6, 7, 6,12,10, 6,10, 7, 7, 7, 8, 2, 6,12  /* V */
};

/* PAM120 matrix with each cell offset by 8 to make all scores non-negative */
char pam120_matrix[21][21]={
11, 5, 8, 8, 5, 7, 8, 9, 5, 7, 5, 6, 6, 4, 9, 9, 9, 1, 4, 8, 8, /* A */
 5,14, 7, 5, 4, 9, 5, 4, 9, 6, 4,10, 7, 4, 7, 7, 6, 9, 2, 5, 8, /* R */
 8, 7,12,10, 3, 8, 9, 8,10, 6, 4, 9, 5, 4, 6, 9, 8, 3, 6, 5, 8, /* N */
 8, 5,10,13, 1, 9,11, 8, 8, 5, 3, 7, 4, 1, 6, 8, 7, 0, 3, 5, 8, /* D */
 5, 4, 3, 1,17, 1, 1, 3, 4, 5, 1, 1, 2, 2, 5, 7, 5, 0, 7, 6, 8, /* C */
 7, 9, 8, 9, 1,14,10, 5,11, 5, 6, 8, 7, 2, 8, 6, 6, 2, 3, 5, 8, /* Q */
 8, 5, 9,11, 1,10,13, 7, 7, 5, 4, 7, 4, 2, 7, 7, 6, 0, 4, 5, 8, /* E */
 9, 4, 8, 8, 3, 5, 7,13, 4, 4, 3, 5, 4, 3, 6, 9, 7, 0, 2, 6, 8, /* G */
 5, 9,10, 8, 4,11, 7, 4,15, 4, 5, 6, 4, 6, 7, 6, 5, 3, 7, 5, 8, /* H */
 7, 6, 6, 5, 5, 5, 5, 4, 4,14, 9, 6, 9, 8, 5, 6, 8, 1, 6,11, 8, /* I */
 5, 4, 4, 3, 1, 6, 4, 3, 5, 9,13, 4,11, 8, 5, 4, 5, 3, 5, 9, 8, /* L */
 6,10, 9, 7, 1, 8, 7, 5, 6, 6, 4,13, 8, 2, 6, 7, 7, 3, 2, 4, 8, /* K */
 6, 7, 5, 4, 2, 7, 4, 4, 4, 9,11, 8,16, 7, 5 ,6 ,7, 1, 4, 9, 8, /* M */
 4, 4, 4, 1, 2, 2, 2, 3, 6, 8, 8, 2, 7,16, 3, 5 ,4 ,7,12, 5, 8, /* F */
 9, 7, 6, 6, 5, 8, 7, 6, 7, 5, 5, 6, 5, 3,14, 9, 7, 1, 2, 6, 8, /* P */
 9, 7, 9, 8, 7, 6, 7, 9, 6 ,6, 4, 7, 6, 5, 9,11,10, 6, 5, 6, 8, /* S */
 9, 6, 8, 7, 5, 6, 6, 7, 5, 8, 5, 7, 7, 4, 7,10,12, 2, 5, 8, 8, /* T */
 1, 9, 3, 0, 0, 2, 0, 0, 3, 1, 3 ,3, 1, 7, 1, 6, 2,20, 7, 0, 8, /* W */
 4, 2, 6, 3, 7, 3, 4, 2, 7, 6, 5, 2, 4,12, 2, 5, 5, 7,16, 5, 8, /* Y */
 8, 5, 5, 5, 6, 5, 5, 6, 5,11, 9, 4, 9, 5, 6, 6, 8, 0, 5,13, 8, /* V */
 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 0  /* J */
};

/*=======================================================================*/
/* Number to amino acid                  */
char *num_to_aachar(num)
int num;
{
  switch (num) {
    case 0: return("A");
    case 1: return("R");
    case 2: return("N");
    case 3: return("D");
    case 4: return("C");
    case 5: return("Q");
    case 6: return("E");
    case 7: return("G");
    case 8: return("H");
    case 9: return("I");
    case 10: return("L");
    case 11: return("K");
    case 12: return("M");
    case 13: return("F");
    case 14: return("P");
    case 15: return("S");
    case 16: return("T");
    case 17: return("W");
    case 18: return("Y");
    case 19: return("V");
    case 20: return("J");
    case -1: return(".");
    default: return("*");		/* Should never happen */
    }
}
/*======================================================================*/
/* Amino acid to number                  */
int aachar_to_num(ch)
char ch;
{
  switch (ch) {
    case 'A': return(0);
    case 'R': return(1);
    case 'N': return(2);
    case 'D': return(3);
    case 'B': return(3);
    case 'C': return(4);
    case 'Q': return(5);
    case 'E': return(6);
    case 'Z': return(6);
    case 'G': return(7);
    case 'H': return(8);
    case 'I': return(9);
    case 'L': return(10);
    case 'K': return(11);
    case 'M': return(12);
    case 'F': return(13);
    case 'P': return(14);
    case 'S': return(15);
    case 'T': return(16);
    case 'W': return(17);
    case 'Y': return(18);
    case 'V': return(19);
    case 'J': return(20);
    case 'O': return(20);
    case 'X': return(20);
    case '.': return(-1);
    default: return(-1);
    }
}
/*=====================================================================*/
/* Use internal numerical representation to print amino acid: */
void pr_num_to_aa(num)
char num;
{
  switch (num) {
    case 0: printf("A"); break;
    case 1: printf("R"); break;
    case 2: printf("N"); break;
    case 3: printf("D"); break;
    case 4: printf("C"); break;
    case 5: printf("Q"); break;
    case 6: printf("E"); break;
    case 7: printf("G"); break;
    case 8: printf("H"); break;
    case 9: printf("I"); break;
    case 10: printf("L"); break;
    case 11: printf("K"); break;
    case 12: printf("M"); break;
    case 13: printf("F"); break;
    case 14: printf("P"); break;
    case 15: printf("S"); break;
    case 16: printf("T"); break;
    case 17: printf("W"); break;
    case 18: printf("Y"); break;
    case 19: printf("V"); break;
    case 20: printf("."); break;
    case -1: printf("."); break;
    default: printf("*");		/* Should never happen */
    }
}
/*======================================================================*/
void pr_num_to_aa_space(c)
char c;
{
	pr_num_to_aa(c);
	printf(" ");
}
/*=======================================================================
      getscore reads a file containing a scoring matrix and
      loads it into Score[MATSIZE][MATSIZE].  Assumes alphabet for file is
      listed on first non-blank line.
=========================================================================*/
void getscore(matrix)
struct score *matrix;
{
   char filename[FNAMELEN], line[MAXLINE], chigh[6], *ptr;
   FILE *fin, *fstp;
   int alpha[MATSIZE+10], nrows, ncols, row, col, i;

   if ((fstp = fopen("protomat.stp", "rt")) == NULL)
   {
      fin = NULL;
      strcpy(filename, "def");
   }
   else
   {
      line[0] = filename[0] = '\0';
      while(fgets(line, sizeof(line), fstp) != NULL)
      {
	  if (strncmp(line, "SCORE", 5) == 0)
	  {
	     ptr = strtok(line, " ,\t\n");
	     if (ptr != NULL)
	     {
		ptr = strtok(NULL, " ,\t\n");
		if (ptr != NULL) strcpy(filename, ptr);
	     }
	     if ((fin = fopen(filename, "rt")) == NULL)
	     {
		printf("Could not open %s, using default PAM scoring matrix\n",
			    filename);
		strcpy(filename, "def");
	     }
	  }
	  else if (strncmp(line, "HIGH", 4) == 0)
	  {
	     ptr = strtok(line, " ,\t\n");
	     if (ptr != NULL)
	     {
		ptr = strtok(NULL, " ,\t\n");
		if (ptr != NULL) strcpy(chigh, ptr);
		matrix->highpass = atoi(chigh);
	     }
	  }
      }
      fclose(fstp);
   }

/*----------Read file until first non-blank line --------------*/
   if (fin != NULL)
   {
      printf("\nUsing scoring matrix from %s\n", filename);
      line[0] = '\0';
      while (strlen(line) < 1 && fgets(line, sizeof(line), fin) != NULL)
	    ;
/*------See if the first line has characters on it ------------*/
      for (col=0; col < 30; col++) alpha[col] = -1;
      if (strstr(line, "A") != NULL)	/* This line has characters */
      {
	 row = 0;	/* # of alphabetic characters on the line */
	 for (i=0; i<strlen(line); i++)
	 {
	    col = aachar_to_num(line[i]);
	    if (col >= 0)
	    {
	       alpha[row] = col;
	       row++;
	    }
	    else if (isalpha(line[i])) row++; /* skip over other alpha */
	 }
      }
/*-------Get the data values now ------------*/
      for (row=0; row<MATSIZE; row++)
	for (col=0; col<MATSIZE; col++)
	   matrix->scores[row][col] = -99;		/* Null value */
      nrows = 0;
      line[0] = '\0';
      while (fgets(line, sizeof(line), fin) != NULL)
      {
	 if (strlen(line) > 1)
	 {
	    if (alpha[nrows] >= 0 && alpha[nrows] < MATSIZE)
	    {
	       row = alpha[nrows]; ncols = 0;
	       ptr = strtok(line, " ,\t\n");
	       while (ptr != NULL)
	       {
		  if (strspn(ptr, "+-0123456789") == strlen(ptr))
		  {
		     col = alpha[ncols];
		     if (col >= 0 && col < MATSIZE)
			matrix->scores[row][col] = atoi(ptr);
		     ncols++;
		  }
		  ptr = strtok(NULL, " ,\t\n");
	       }
	    }
	    nrows++;
	 }
      }

/*-------If some entries are still missing, assume symmetry ---------*/
      for (row=0; row<MATSIZE; row++)
      {
	for (col=0; col<MATSIZE; col++)
	{
	   if (matrix->scores[row][col] == -99)
		  matrix->scores[row][col] = matrix->scores[col][row];
/*	   printf("%2d ", matrix->scores[row][col]); */
	}
/*	printf("\n"); */
      }
      fclose(fin);
   }
   else    /*   no input file  */
   {
      printf("\nUsing PAM120 scoring matrix.\n");
      matrix->highpass = HIGHPASS;
      for (row=0; row<MATSIZE; row++)
      {
	 for (col=0; col<MATSIZE; col++)
	 {
	    matrix->scores[row][col] = pam120_matrix[row][col];
/*	    printf("%2d ", matrix->scores[row][col]);*/
	 }
/*	 printf("\n");*/
      }
   }
   matrix->highpass *= 100;
/*   printf("HighPass = %d", matrix->highpass);*/
}   /* end of getscore() */
/*======================================================================*/
void init_dbs(dbs)
struct db_info *dbs[];
{
   int i;

   for (i=0; i<MAXDB; i++)
   {
      dbs[i] = (struct db_info *) malloc(sizeof(struct db_info));
      dbs[i]->type = (char *) malloc(10 * sizeof(char));
      dbs[i]->start = (char *) malloc(12 * sizeof(char));
      dbs[i]->desc = (char *) malloc(12 * sizeof(char));
      dbs[i]->seq = (char *)  malloc(12 * sizeof(char));
      dbs[i]->end = (char *)  malloc(6 * sizeof(char));
   }

   dbs[GB]->type = "GENBANK";
   dbs[GB]->start = "LOCUS";
   dbs[GB]->desc = "DEFINITION";
   dbs[GB]->seq = "ORIGIN";
   dbs[GB]->end = "//";
   dbs[GB]->title_offset = 12;
   dbs[GB]->seq_offset = 10;

   dbs[PIR]->type = "PIR";
   dbs[PIR]->start = "ENTRY";
   dbs[PIR]->desc = "TITLE";
   dbs[PIR]->seq = "SEQUENCE";
   dbs[PIR]->end = "//";
   dbs[PIR]->title_offset = 12;
   dbs[PIR]->seq_offset = 10;

   dbs[EMBL]->type = "EMBL";
   dbs[EMBL]->start = "ID";
   dbs[EMBL]->desc = "DE";
   dbs[EMBL]->seq = "SQ";
   dbs[EMBL]->end = "//";
   dbs[EMBL]->title_offset = 5;
   dbs[EMBL]->seq_offset = 5;

   dbs[UNI]->type = "UNI";
   dbs[UNI]->start = ">";
   dbs[UNI]->desc = ">";
   dbs[UNI]->seq = "";
   dbs[UNI]->end = "*";
   dbs[UNI]->title_offset = 1;
   dbs[UNI]->seq_offset = 0;

   dbs[VMS]->type = "VMS";
   dbs[VMS]->start = ">";		/* first line */
   dbs[VMS]->desc = "";			/* second line */
   dbs[VMS]->seq = "";			/* third line */
   dbs[VMS]->end = "*";
   dbs[VMS]->title_offset = 4;		/* first line only */
   dbs[VMS]->seq_offset = 0;

}   /*  end of init_dbs */
/*======================================================================
      type_dbs() determines what type a database is from the allowable
      types
========================================================================*/
int type_dbs(fin, dbs)
FILE *fin;
struct db_info *dbs[];
{
   int db, i;
   char line[MAXLINE];

   db = -1;
/*---------  Figure out what type of input file it is ------------------*/
   while (db < 0 && fgets(line, sizeof(line), fin) != NULL)
      for (i=0; i<MAXDB; i++)
	 if (strncmp(line, dbs[i]->start, strlen(dbs[i]->start)) == 0)
	    db = i;
   if (db == VMS && line[3] != ';')	/* start=='>' is ambiguous */
	db = UNI;
   if (db < 0 || db >= MAXDB) db = -1;	/* can't tell what it is */

   rewind (fin);
   return(db);
}  /* end of type_dbs */
/*====================================================================*/
/*  This is Kernighan & Ritchie's ASCII to integer conversion (p. 58) */
/*====================================================================*/
int kr_atoi(s)
char s[];
{
   int i, n, sign;

   for (i=0; s[i]==' ' || s[i]=='\n' || s[i]=='\t'; i++)
		;
   sign = 1;
   if (s[i] == '+' || s[i] == '-')
      sign = (s[i++]=='+') ? 1 : -1;
   for (n=0; s[i] >= '0' && s[i] <= '9'; i++)
      n = 10 * n + s[i] - '0';
   return(sign * n);
}  /* end of kr_atoi */
/*=====================================================================*/
/*  This is Kernighan & Ritchie's integer to ASCII conversion (p. 60) */
/*====================================================================*/
void kr_itoa(n, s)
int n;
char s[];
{
   int c, i, j, sign;

   sign = n;
   if (sign < 0) n = -n;
   i = 0;
   do {
      s[i++] = n % 10 + '0';
   }  while ( (n /= 10) > 0);
   if (sign < 0) s[i++] = '-';
   s[i] = '\0';
   for (i=0, j=strlen(s)-1; i<j; i++, j--)
   {
      c = s[i];
      s[i] = s[j];
      s[j] = c;
   }
}   /*  end of kr_itoa */
/*=====================================================================
     Locate the directory, file name and extension in a file name
========================================================================*/
struct split_name *split_names(filename)
char *filename;
{
   struct split_name *new;
   int i, ext_len;

   ext_len = 0;
   new = (struct split_name *) malloc(sizeof(struct split_name));
   new->dir_len=new->file_len=new->name_len=0;
   i = strlen(filename);
   /*-------  Read the file name backwards ---------------*/
   while (i>=0 && (!ext_len || !new->dir_len))
   {
      /*---  Last period in string => file extension ----*/
      if (filename[i] == '.') ext_len = strlen(filename)-i;
      /*--- Last slash in string => directory -----------*/
      if (filename[i] == '/' && new->dir_len == 0) new->dir_len = i+1;
      /*--- Last colon and no slash after it => DOS directory -----*/
/*      if (filename[i] == ':' && new->dir_len == 0) new->dir_len = i+1; */
      i--;
   }
   new->file_len = strlen(filename)-new->dir_len;
   new->name_len = new->file_len - ext_len;

   return(new);
}
/*========================================================================
  dir_dos(): DOS code to get name of directory from line and
       create it if necessary
===========================================================================*/
/*char *dir_dos(line)
char *line;
{
   char tname[FNAMELEN], mem[MAXLINE], pros[FNAMELEN], *ptr;
   char filename[FNAMELEN];
   int test;
   FILE *ftmp;

   pros[0] = '\0';
   if (line[0] != '>' &&
	  (strstr(line, "\\") != NULL || strstr(line, ":") != NULL) )
   {
         ptr = strtok(line, "\n\r");
         strcpy(pros, ptr);
	 if (pros[strlen(pros)-1] != '\\') strcat(pros, "\\");
   }
*/
/*-------------------Create the directory ---------------------------------*/
/*
   if (strlen(pros))
   {
      tmpnam(tname);
      strcpy(filename, pros);
      strcat(filename, tname);
      if ( (ftmp=fopen(filename, "w"))== NULL)
      {
	 strcpy(tname, pros);     
	 tname[strlen(pros)-1] = '\0';
	 sprintf(mem, "md %s", tname);
	 test = system(mem);
	 if (test == 0) printf("\nCreated directory %s", tname);
	 else
	 {
	    printf("\nUnable to create directory %s", tname);
	    printf("\nProtein files will be placed in current directory");
	    pros[0] = '\0';
	 }
      }
      else
      {
	 fclose(ftmp);
	 unlink(filename);
      }
   }
   return(pros);
} 
*/
/*=======================================================================
  dir_unix(): UNIX code to get name of directory from line and
	      create it if necessary
==========================================================================*/
char *dir_unix(line)
char *line;
{
   char tname[FNAMELEN], mem[MAXLINE], pros[FNAMELEN], *ptr;
   int test;
   DIR *dp;

   pros[0] = '\0';
   if (line[0] != '>' && (strstr(line, "/") != NULL != NULL) )
   {
         ptr = strtok(line, "\n\r");
         strcpy(pros, ptr);
	 if (pros[strlen(pros)-1] != '/') strcat(pros, "/");
   }
/*-----------------Create the directory ---------------------------------*/
   if (strlen(pros))
   {
      strcpy(tname, pros);
      tname[strlen(pros)-1] = '\0';
      sprintf(mem, "mkdir %s", tname);
      if ((dp=opendir(tname))==NULL)
      {
         test=system(mem);
         if (test == 0) printf("\nCreated directory %s\n", tname);
         else
         {
            printf("\nUnable to create directory %s", tname);
            printf("\nProtein files will be placed in current directory\n");
            pros[0] = '\0';
         }
      }
   }
   return(pros);
}    /* end of dir_unix */
/*======================================================================
     Create & intialize a db_id structure
========================================================================*/
struct db_id *makedbid()
{
   struct db_id *new;

   new = (struct db_id *) malloc (sizeof(struct db_id));
   new->entry[0] = '\0';
   new->ps[0] = '\0';
   new->len = 0;
   new->rank = new->score = 0;
   new->lst = NO;
   new->found = NO;
   new->block = NO;
   new->frag = NO;
   new->search = NO;
   new->pvalue = (double) 0.0; 
   new->next = NULL;
   new->prior = NULL;
   return(new);
}  /*  end of makedbid */

/*======================================================================
     get_ids() reads a .lis or .lst file & inserts the sequences
       found in it into a list sorted by sequence name.
========================================================================*/
int get_ids(flis, ids)
FILE *flis;
struct db_id *ids;
{
   char line[MAXLINE], ctemp[10], *ptr;
   struct db_id *id, *last, *new;
   int len, nids = 0;

   while(!feof(flis) && fgets(line, MAXLINE, flis) != NULL &&
	 strlen(line) > 2)
   {		/* skip over title or directory lines */
      if (line[0] != '>' && strstr(line, "/") == NULL &&
	  strstr(line,"\\") == NULL && strstr(line, ":") == NULL)
      {
	 nids += 1;
/*-----  Copy up to the first space or carriage return ------------*/
	 len = strcspn(line, " \t\r\n");
	 if (len > IDLEN) len = IDLEN;	/* No id should be longer that this*/
	 new = makedbid();
	 strncpy(new->entry, line, len);  new->entry[len] = '\0';
	 last = ids;  id = ids->next;
/*-------- Get any other information on the .lis file line --------*/
	 if (strstr(line, "FRAGMENT") != NULL) new->frag = YES;
	 if (strstr(line, "LST") != NULL) new->lst = YES;
	 ptr = strstr(line, "PS=");
	 if (ptr != NULL)
	    strncpy(new->ps, ptr+3, 1); new->ps[1] = '\0';
	 ptr = strstr(line, "LENGTH=");
	 if (ptr != NULL)
	 {
	    len = strcspn(ptr+7, " \t\r\n");
	    if (len > 0)
	    {
	       strncpy(ctemp, ptr+7, len); ctemp[len] = '\0';
	       new->len = atoi(ctemp);
	    }
	 }
/*------  Insert id into a sorted list ----------------------------*/
	 while (id != NULL && id->entry != NULL &&
	     strcmp(id->entry, new->entry) < 0)
	 {
	    last = id;
	    id = id->next;
	 }
	 new->prior = last;
	 new->next = id;
	 last->next = new;
	 if (id != NULL) id->prior = new;
	 new = NULL;
      }
   }
   return(nids);
}  /*  end of get_ids */
/*======================================================================
     Extracts the sequences from file fin if they appear in the
     sorted list ids.
=======================================================================*/
int extract_seqs(nids, dbs, fin, ids, pros)
int nids;
struct db_info *dbs[MAXDB];
FILE *fin;
struct db_id *ids;
char *pros;
{
   FILE *fout;
   struct db_id *id;
   int nseq, i, db;
   char line[MAXLINE], title[MAXLINE], temp[MAXLINE], *ptr;
   char foutname[FNAMELEN];

   nseq = 0;
   db = type_dbs(fin, dbs);
   if (db < 0 || db >= MAXDB)
   {
      printf("\nCannot determine type of input file");
      return(-1);
   }
   printf("\nProcessing input file as %s", dbs[db]->type);
   if (db==VMS)
   {
      printf("\n WARNING: Titles are sometimes truncated in this format;");
      printf("\n          I may not be able to distinguish fragments.");
   }

/*-------- Modify the entry ids depending on the file type.  If it is
      EMBL, change any % to $, since the ids have $ in the database.
      If it is UNI, change any $ to % since the ids have % in the database */
/*   id = ids->next;
   while(id != NULL)
   {
      for (i=0; i<strlen(id->entry); i++)
      {
	 if      (db == EMBL && id->entry[i] == '%') id->entry[i] = '$';
	 else if (db == UNI  && id->entry[i] == '$') id->entry[i] = '%';
      }
      id = id->next;
   }
*/

/*---  Since the list of requested entries is sorted by id, read through
       the list until past current db entry for each db entry until all
       requested ids have been extracted -------------------------------*/
   while (!feof(fin) && fgets(line, MAXLINE, fin) != NULL &&
	   nseq < nids)
   {
      if (strncmp(line, dbs[db]->start, strlen(dbs[db]->start)) == 0)
      {
	 if (db == VMS )	/* get the 2nd line, too */
	 {
	   ptr=strtok(line, "\n\r");		/* get rid of CRLF */
	   strcat(line, " ");
	   if (fgets(temp, MAXLINE, fin) != NULL) strcat(line, temp);
	 }
/*--------------Check to see if this sequence is in the list -----------*/
	 id = ids->next;
	 while (id != NULL &&
		strncmp(id->entry, &line[dbs[db]->title_offset],
		    strlen(id->entry)) <= 0)
	 {
	    if (strncmp(id->entry, &line[dbs[db]->title_offset],
			strlen(id->entry)) == 0 &&
			id->found == NO)
	    {
	       nseq += 1;  id->found = YES; id->len = 0;
     /*----  Check to see if this sequence is a fragment ----- */
	       if (strstr(line, "FRAGMENT") != NULL ||
		   strstr(line, "fragment") != NULL) id->frag = 1;

/*------  Fix up the file names:  NOTE: since id->entry is 10 chars
    and DOS file names are only 8, DOS file name may be truncated and
    therefore not unique --------*/
	       for (i=0; i<strlen(id->entry); i++)
		  if (db == UNI  &&  id->entry[i] == '%')
			  id->entry[i] = '$';
	       strcpy(foutname, pros);
	       strcpy(temp, id->entry);  temp[SNAMELEN-1] = '\0';
	       strcat(foutname, temp);
	       if (db == GB) strcat(foutname, ".dna");
	       else strcat(foutname, ".pro");
	       printf("\n%d. Entry %s found...Creating %s",
		      nseq, id->entry, foutname);
/*------Open file: should check whether it already exists... */
	       if ( (fout = fopen(foutname, "w+t")) == NULL)
	       {
		  printf("\nCannot open %s\n", foutname);
		  return(-1);
	       }
	       if (db == UNI || db == VMS)
	       {
		   strcpy(temp, &line[dbs[db]->title_offset]);
		   for (i=0; i<strlen(temp); i++)
		   {
		      if (temp[i] == '$')  temp[i] = '%';
		      else if (temp[i] == '\n' || temp[i] == '\r')
			    temp[i] = ' ';
		   }
/*		   fprintf(fout, ">%s$\n", temp); */
                   fprintf(fout, ">%s\n", temp);
	       }
	       else		/*  EMBL, GENBANK, PIR */
	       {
		  ptr = strtok(&line[dbs[db]->title_offset], " ");
		  strcpy(title, ptr);
	       }
/*----- Now write out the sequence, counting line lengths at same time ---*/
	      while (!feof(fin) && fgets(line, MAXLINE, fin) != NULL &&
		 strncmp(line, dbs[db]->end, strlen(dbs[db]->end)) != 0 &&
		 strncmp(line, dbs[db]->start, strlen(dbs[db]->start)) != 0)
	      {
		  if (db == UNI || db == VMS)
		  {
		     fputs(line, fout);
		     for (i=0; i<strlen(line); i++)
		       if (line[i] != ' '  && line[i] != '\n' &&
			   line[i] != '\r' && line[i] != '\t')
			      id->len = id->len + 1;
		  }
		  else
		  {
		     if (strncmp(line,dbs[db]->desc,strlen(dbs[db]->desc))==0 &&
			 strlen(title) + strlen(line) < MAXLINE )
		     {
			strcpy(temp, title); strcat(temp, " ");
			strcat(temp, &line[dbs[db]->title_offset]);
			ptr = strtok(temp, "\n\r");   /* remove CRLF */
			strcpy(title, temp);
		     }
/*------ Process the sequence ------------------------------------------*/
		     if (strncmp(line,dbs[db]->seq,strlen(dbs[db]->seq))==0)
		     {
	      /*--- Replace any $ characters with % in the sequence title.*/
/*			for (i=0; i<strlen(temp); i++)
			   if (temp[i] == '$')  temp[i] = '%';
*/
	      /*----  Check to see if this sequence is a fragment ----- */
			if (strstr(title, "FRAGMENT") != NULL ||
			    strstr(title, "fragment") != NULL) id->frag = 1;
			fprintf(fout, ">%s\n", title);  /*print title*/
			while(!feof(fin) &&
			   fgets(line, MAXLINE, fin) != NULL &&
			  strncmp(line,dbs[db]->end,strlen(dbs[db]->end))!=0)
			{
			    fprintf(fout, &line[dbs[db]->seq_offset]);
			    for (i=0; i<(strlen(line)-dbs[db]->seq_offset);i++)
			     if (line[i] != ' '  && line[i] != '\n' &&
				 line[i] != '\r' && line[i] != '\t')
				    id->len = id->len + 1;
			 }
		      }
		  }   /*  end of if not UNI and not VMS */
	       }  /*  end of records for entry */
/*	       if (db != VMS) fprintf(fout, "*\n");   VMS has * already */
	       if (db == VMS &&
		   strncmp(line, dbs[db]->start, strlen(dbs[db]->start)) == 0)
	       {                           /* get the 2nd line, too */
		   ptr = strtok(line, "\n\r");        /* remove CRLF */
		   strcat(line, " ");
		   if (fgets(temp, MAXLINE, fin) != NULL) strcat(line, temp);
	       }
	       fclose(fout);
	    }  /*  end of found entry */
	    id = id->next;
	 }  /* end of list of requested entries */
      }   /*  end of a db entry */
   }  /*  end of db file */

   return(nseq);
}  /* end of extract_seqs */
