#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>
#include <globals.h> 

#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  
#include "egsubs.h"  
#include "exclude.h" 

#define WVERSION   "420" 

// badpairsname added 
/** 
 does nothing but read the data 
 and print snps
 sometimes a good place to start !!

 New I/O (mcio.c) added
 New admutils (snpindex hash)
 mcio bug fixed (large files)
 */

#define MAXFL  50   
#define MAXSTR  512

extern int packmode;

char *trashdir = "/var/tmp";
extern int verbose;
int qtmode = NO;
Indiv **indivmarkers;
SNP **snpmarkers;
int numsnps, numindivs;

char *genotypename = NULL;
char *snpname = NULL;
char *genooutfilename = NULL;
char *indoutfilename = NULL;
char *indivname = NULL;
char *badsnpname = NULL;
char *goodsnpname = NULL;
char *badpairsname = NULL;
char *markername = NULL;
char *idname = NULL;

char *outputname = NULL;
FILE *ofile;

double fakespacing = 0.0;

char unknowngender = 'U';

void
readcommands (int argc, char **argv);
void
dophyscheck (SNP **snpm, int numsnps);

int
main (int argc, char **argv)
{

  int i, j, k, g;
  SNP *cupt;
  Indiv *indx;
  int ch1, ch2;

  int numvind, nignore, numrisks = 1;
  int markernum, idnum;

  ofile = stdout;
  packmode = YES;
  readcommands (argc, argv);
  if (indivname == NULL)
    {
      printf ("no indivname\n");
      return 0;
    }
  if (outputname != NULL)
    openit (outputname, &ofile, "w");

  numsnps = getsnps (snpname, &snpmarkers, fakespacing, badsnpname, &nignore,
                     numrisks);

// fakespacing 0.0 (default)

  numindivs = getindivs (indivname, &indivmarkers);
  setstatus (indivmarkers, numindivs, "Case");

  setgenotypename (&genotypename, indivname);

  printf ("genotypename:  %s\n", genotypename);

  if (genotypename != NULL)
    {
      getgenos (genotypename, snpmarkers, indivmarkers, numsnps, numindivs,
                nignore);

      /**
       if (badpairsname != NULL) {
       loadbadpsc(snpmarkers, numsnps, NO, goodsnpname) ;
       dobadpairs(badpairsname, snpmarkers, numsnps) ;
       }
       */
    }
  dophyscheck (snpmarkers, numsnps);

  numvind = numvalidind (indivmarkers, numindivs);
  printf ("\n\n");
  printf ("numindivs: %d valid: %d numsnps: %d nignore: %d\n", numindivs,
          numvind, numsnps, nignore);

  if (verbose)
    {
      for (i = 0; i < numindivs; ++i)
        {
          indx = indivmarkers[i];
          printf ("%20s ", indx->ID);
          for (j = 0; j < numsnps; ++j)
            {
              cupt = snpmarkers[j];
              if (cupt->ignore)
                continue;
              g = getgtypes (cupt, i);
              if (g < 0)
                g = 9;
              printf ("%1d", g);
            }
          printf ("  %20s", indx->egroup);
          printnl ();
        }
    }
// numsnps includes fakes

  if (markername != NULL)
    {
      markernum = snpindex (snpmarkers, numsnps, markername);
      if (markernum < 0)
        fatalx ("markername %s not found\n", markername);
      cupt = snpmarkers[markernum];
      printf ("markername: %s  %d   %9.3f %12.0f\n", cupt->ID, cupt->chrom,
              cupt->genpos, cupt->physpos);
      for (i = 0; i < numindivs; ++i)
        {
          indx = indivmarkers[i];
          g = getgtypes (cupt, i);
          printf ("%20s %20s %2d\n", cupt->ID, indx->ID, g);
        }
    }

  if (idname != NULL)
    {
      idnum = indindex (indivmarkers, numindivs, idname);
      if (idnum < 0)
        fatalx ("idname %s not found\n", idname);
      indx = indivmarkers[idnum];
      printf ("idname: %20s  %c   %20s\n", indx->ID, indx->gender,
              indx->egroup);
      for (j = 0; j < numsnps; ++j)
        {
          cupt = snpmarkers[j];
          if (cupt->ignore)
            continue;
          g = getgtypes (cupt, idnum);
          printf ("%20s %20s %2d", cupt->ID, indx->ID, g);
          printf ("  %3d %12.0f", cupt->chrom, cupt->physpos);
          printf (" %c %c", cupt->alleles[0], cupt->alleles[1]);
          printnl ();
        }
    }

  /**
   if (genotypename != NULL) {
   printdata(genooutfilename, indoutfilename, snpmarkers, indivmarkers, numsnps, numindivs, NO) ;
   }
   */

  printf ("##end of run\n");
  return 0;
}
void
readcommands (int argc, char **argv)

{
  int i, haploid = 0;
  char *parname = NULL;
  phandle *ph;
  char str[5000];
  char *tempname;
  int n;

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

  /**
   DIR2:  /fg/nfiles/admixdata/ms2
   SSSS:  DIR2/outfiles 
   genotypename: DIR2/autos_ccshad_fakes
   eglistname:    DIR2/eurlist  
   output:        eurout
   */
  getint (ph, "packmode:", &packmode); // controls internals 

  getstring (ph, "genotypename:", &genotypename);
  getstring (ph, "genooutfilename:", &genooutfilename);
  getstring (ph, "indoutfilename:", &indoutfilename);
  getstring (ph, "snpname:", &snpname);
  getstring (ph, "indivname:", &indivname);
  getstring (ph, "output:", &outputname);
  getstring (ph, "badsnpname:", &badsnpname);
  getstring (ph, "goodsnpname:", &goodsnpname);
  getstring (ph, "badpairsname:", &badpairsname);
  getstring (ph, "markername:", &markername);
  getstring (ph, "idname:", &idname);
  getdbl (ph, "fakespacing:", &fakespacing);
  getint (ph, "familynames:", &familynames);
  writepars (ph);
  closepars (ph);

}

void
dophyscheck (SNP **snpm, int numsnps)
{
// catch places where physpos genpos are in opposite order
  SNP *cupt, *cuptold;
  int i;

  for (i = 0; i < numsnps; i++)
    {
      cupt = snpm[i];
      if (i == 0)
        cuptold = cupt;
      if (cupt->isfake)
        continue;
      if (cupt->ignore)
        continue;
      if (cupt->chrom == cuptold->chrom)
        {
          if (cupt->physpos < cuptold->physpos)
            {
              printf ("physcheck %20s %15s %12.3f %12.3f %13.0f %13.0f\n",
                      cuptold->ID, cupt->ID, cuptold->genpos, cupt->genpos,
                      cuptold->physpos, cupt->physpos);
            }
        }
      cuptold = cupt;
    }
}

