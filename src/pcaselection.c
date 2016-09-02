#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  
#include "globals.h"  

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#define WVERSION   "1000" 

#define MAXFL  50   
#define MAXSTR  512
#define MAXSIZE 8.0e9

typedef enum outputmodetype inputmodetype;

extern int packmode;
extern int malexhet;
extern int verbose;
extern int plotmode;

char *trashdir = "/var/tmp";
int qtmode = NO;

/* major data structures */
Indiv **indivmarkers;
SNP **snpmarkers;
int numsnps, numindivs;

char *genotypename = NULL; /* name of genotype file */
char *snpname = NULL; /* name of SNP file */
char *indivname = NULL; /* name of sample file */
char *pcaname = NULL; /* name of pca file */
char *imode = "eigenstrat"; /* input mode */
char *outputname = NULL; /* name of output file */
int numpc = 10; /* number of principal components
 to correct */

/* 
 If these are to be global, remove them from function parameter lists.
 If they're going to be local, put the rest in (chisq routines)    
 */

int NSAMPLES;
int *outlier;
int L;

inputmodetype inmode;
FILE *fpout; /* output file */

void
readcommands (int argc, char **argv);
void
setinmode (inputmodetype *inmode, char *imode);
int
read_evec (char *filename, double **eval, double **evec, size_t *K, size_t *N);

int
main (int argc, char **argv)
{
  size_t K, N;
  int nignore;
  double rowsum, rowsum1;
  double chisq, Echisq, gamma, denom;

  readcommands (argc, argv);
  if (outputname != NULL)
    openit (outputname, &fpout, "w");
  else
    fpout = stdout;
  fprintf (fpout, "Chisq PCASELECTION\n");

  setinmode (&inmode, imode);
  packmode = YES;

  numsnps = getsnps (snpname, &snpmarkers, 0.0, NULL, &nignore, 1);

  NSAMPLES = getindivs (indivname, &indivmarkers);

  setstatus (indivmarkers, NSAMPLES, "Case");
  setgenotypename (&genotypename, indivname);
  if (genotypename != NULL)
    {
      getgenos (genotypename, snpmarkers, indivmarkers, numsnps, NSAMPLES,
                nignore);
    }

  double *eval, *evec;

  if (read_evec (pcaname, &eval, &evec, &K, &N) < 1)
    {
      printf ("Error reading evec file \"%s\".\n", pcaname);
      return 1;
    }

  if (N != NSAMPLES)
    {
      printf ("Number of samples doesn't match: %d != %d", NSAMPLES, N);
      return 1;
    }

    {
      size_t i, j, k;
      double *vg = (double *) malloc (K * sizeof(double));
      double *v1 = (double *) malloc (K * sizeof(double));
      for (i = 0; i < numsnps; i++)
        {
          N = 0;
          double p = 0;

          SNP *cupt = snpmarkers[i];

          for (k = 0; k < K; k++)
            {
              vg[k] = 0;
              v1[k] = 0;
            }

          for (j = 0; j < NSAMPLES; j++)
            {
              int g = getgtypes (cupt, j);

              if (g >= 0)
                {
                  N++;
                  p += g;
                  for (k = 0; k < K; k++)
                    {
                      vg[k] += evec[j * K + k] * g;
                    }
                }

              for (k = 0; k < K; k++)
                {
                  v1[k] += evec[j * K + k];
                }
            }

          p /= 2*N;

          fprintf(fpout, "%s", cupt->ID);
          for (k = 0; k < K; k++) {
              vg[k] -= 2*p*v1[k];
              vg[k] *= vg[k];
              vg[k] /= 2*p*(1-p);
              vg[k] /= eval[k];
              fprintf(fpout, "\t%g", vg[k]);
          }
          fprintf(fpout, "\n");
        }
    }

  return 0;
}

void
readcommands (int argc, char **argv)
{
  int i;
  char *parname = NULL;
  phandle *ph;

  while ((i = getopt (argc, argv, "p:vV")) != -1)
    {

      switch (i)
        {

        case 'p':
          parname = strdup (optarg);
          break;

        case 'v':
          printf ("version: %s\n", WVERSION);
          break;

        case 'V':
          verbose = YES;
          break;

        case '?':
          printf ("Usage: bad params.... \n");
          fatalx ("bad params\n");
        }
    }

  pcheck (parname, 'p');
  printf ("parameter file: %s\n", parname);
  ph = openpars (parname);
  dostrsub (ph);

  getint (ph, "packmode:", &packmode); // controls internals 
  getstring (ph, "genotypename:", &genotypename);
  getstring (ph, "snpname:", &snpname);
  getstring (ph, "indivname:", &indivname);
  getstring (ph, "pcaname:", &pcaname);
  getstring (ph, "outputname:", &outputname);

  writepars (ph);
  closepars (ph);

}

void
setinmode (inputmodetype *inmode, char *imode)
{
  char *ss = strdup (imode);
  int len = strlen (ss);

  int i;
  for (i = 0; i < len; i++)
    {
      ss[i] = tolower(ss[i]);
    }

  *inmode = EIGENSTRAT; /* default */
  if (strcmp (ss, "eigenstrat") == 0)
    *inmode = EIGENSTRAT;
  if (strcmp (ss, "ped") == 0)
    *inmode = PED;
  if (strcmp (ss, "packedped") == 0)
    *inmode = PACKEDPED;
  if (strcmp (ss, "ancestrymap") == 0)
    *inmode = ANCESTRYMAP;
  if (strcmp (ss, "packedancestrymap") == 0)
    *inmode = PACKEDANCESTRYMAP;

}

int
read_evec (char *filename, double **eval, double **evec, size_t *K, size_t *N)
{
  // OPEN THE FILE
  FILE *fp = fopen (filename, "r");
  if (fp == NULL)
    {
      printf ("File \"%s\" could not be opened.\n", filename);
      return -1;
    }

  // READ THE "#eigvals:" FROM THE FIRST LINE
  char buffer[256];
  fscanf (fp, "%s", buffer);

  if (strcmp (buffer, "#eigvals:") != 0)
    {
      printf ("File \"%s\" doesn't begin with \"#eigvals:\".\n", filename);
      fclose (fp);
      return -1;
    }

  // READ THE REST OF THE FIRST LINE OF THE FILE
  char *line = NULL;
  size_t red, len = 0;

  red = getline (&line, &len, fp);

  if (red < 1)
    {
      printf ("File \"%s\" appears to be empty.\n", filename);
      fclose (fp);
      return -1;
    }

  // COUNT THE EVALS AND PUT THEM INTO **eval
  size_t k = 0;
    {
      double val;

      char *pEnd = line;
      while ((val = strtod (pEnd, &pEnd)) != 0.0)
        {
          k++;
        }

      *eval = (double *) malloc (k * sizeof(double));

      pEnd = line;
      size_t i;
      for (i = 0; i < k; i++)
        {
          (*eval)[i] = strtod (pEnd, &pEnd);
        }
    }
  free (line);

  // COUNT THE LINES
  size_t end_of_first_line = ftell (fp);
  size_t n = 0;
  char c = getc (fp);
  while (c != EOF)
    {
      if (c == '\n')
        {
          n++;
        }
      c = getc (fp);
    }

  *K = k;
  *N = n;

  // READ EVEC
  *evec = (double *) malloc (k * n * sizeof(double));

  fseek (fp, end_of_first_line, 0);

    {
      size_t i;
      for (i = 0; i < n; i++)
        {
          fscanf (fp, "%s", buffer);

          line = NULL;
          len = 0;
          red = getline (&line, &len, fp);

          char *pEnd = line;
          size_t j;
          for (j = 0; j < k; j++)
            {
              (*evec)[k * i + j] = strtod (pEnd, &pEnd);
            }
          free (line);
        }
    }

  fclose (fp);
  return (k);
}

