#include <stdio.h>
#include <string.h>
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

#define WVERSION   "1000" 


#define MAXFL  50   
#define MAXSTR  512
#define MAXSIZE 8.0e9

typedef enum outputmodetype inputmodetype;  

extern int packmode ;
extern int malexhet ;  
extern int verbose ;  
extern int plotmode ;  

char *trashdir = "/var/tmp" ;
int qtmode = NO ;

/* major data structures */
Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 

char  *genotypename = NULL ;		/* name of genotype file */
char  *snpname = NULL ; 		/* name of SNP file */
char  *indivname = NULL ;		/* name of sample file */
char  *pcaname = NULL ;			/* name of pca file */
char *imode = "eigenstrat";		/* input mode */
char *outputname = NULL ;		/* name of output file */
int numpc = 10;				/* number of principal components
						to correct */

/* 
   If these are to be global, remove them from function parameter lists.
   If they're going to be local, put the rest in (chisq routines)    
*/

int NSAMPLES;
int *outlier;
int L;

inputmodetype inmode;
FILE *fpout; 				/* output file */

void readcommands(int argc, char **argv) ;
void setinmode(inputmodetype *inmode, char *imode);
void readpcafile(double **Vp, int **outlierp, int *kp, int L, int NSAMPLES);
void getphenos(int NSAMPLES, double **iscasep, int *outlier, 
  double **iscasecorrp, int L, double *V);
double compute_chisq(double *source, double *target);
double compute_chisqE(double *source, double *target);


int main(int argc, char **argv)
{
  double *V;
  double *xx;
  double *iscase;
  double *iscasecorr;
  int K;
  int k,m,n;
  int nignore;
  double rowsum, rowsum1;
  double chisq, Echisq, gamma, denom;

  readcommands(argc, argv) ;
  if (outputname != NULL) 
    openit(outputname, &fpout, "w") ;
  else 
    fpout = stdout;
  fprintf(fpout, "Chisq EIGENSTRAT\n");

  setinmode(&inmode, imode);
  packmode = YES;

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0,  NULL, &nignore, 1) ;

  NSAMPLES = getindivs(indivname, &indivmarkers) ;

  setstatus(indivmarkers, NSAMPLES, "Case") ;
  setgenotypename(&genotypename, indivname) ;
  if (genotypename != NULL)  {
   getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, NSAMPLES, nignore) ;
  }

  /*******************************************************************/
  /*  Free memory:  Usually this is done in outfiles:                */
  /*                                                                 */
  /*  nind = rmindivs(&snpmarkers, numsnps, &indmarkers, NSAMPLES);  */
  /*                                                                 */
  /*  But where is the snpmarkers array released?                    */
  /*******************************************************************/

  L = numpc;
  readpcafile(&V, &outlier, &K, L, NSAMPLES);
  getphenos(NSAMPLES, &iscase, outlier, &iscasecorr, L, V);

  /* main eigenstrat loop here */

  if ((xx = (double *)malloc(NSAMPLES*sizeof(*xx))) == NULL)
  {  fprintf(stderr,"CM\n");  exit(1);  }

  for(m=0;m<numsnps;m++)  {

    SNP *cupt = snpmarkers[m];
    for(n=0; n<NSAMPLES; n++)
    {
      int j = getgtypes(cupt,n);

      if(j == 0)       { xx[n] = 0.0; }
      else if(j == 1)  { xx[n] = 0.5; }
      else if(j == 2 ) { xx[n] = 1.0; }
      else if(j == -1) { xx[n] = -100.0; }

      if(outlier[n] == 1) xx[n] = -100.0;

    }

    /* mean-adjust xx */
    rowsum = 0.0; rowsum1 = 0.0;
    for(n=0; n<NSAMPLES; n++)
    {
      if(qtmode == NO && ((outlier[n]) || (xx[n] < -99.0))) continue;
      if(qtmode == YES && ((outlier[n]) || (xx[n] == -100.0))) continue;
      rowsum += xx[n];
      rowsum1 += 1.0;
    }
    for(n=0; n<NSAMPLES; n++)
    {  
      if(outlier[n]) continue;
      if(qtmode == NO)  {
        if (xx[n] < -99.0) 
          xx[n] = -100.0; /* still keep track */
        else 
	  xx[n] -= rowsum/rowsum1;
      }
      else  {
        if (xx[n] == -100.0) 
          xx[n] = -100.0; /* still keep track */
        else 
	  xx[n] -= rowsum/rowsum1;
      }
    }

    /* Chisq */
    chisq = compute_chisq(xx,iscase);

    /* EIGENSTRAT */
    for(k=0; k<L; k++)
    {
      gamma = 0.0;
      denom = 0.0;
      for(n=0; n<NSAMPLES; n++) 
      {
        if(qtmode == NO && (outlier[n] || xx[n] < -99.0)) continue;
        if(qtmode == YES && (outlier[n] || xx[n] == -100.0)) continue;
        gamma += xx[n]*V[NSAMPLES*n+k];
        denom += V[NSAMPLES*n+k]*V[NSAMPLES*n+k];
      }
      gamma /= denom;
      for(n=0; n<NSAMPLES; n++) 
      {
        if(qtmode == NO && (outlier[n] || xx[n] < -99.0)) continue;
        if(qtmode == YES && (outlier[n] || xx[n] == -100.0)) continue;
        xx[n] -= gamma*V[NSAMPLES*n+k];
      }
    }
    Echisq = compute_chisqE(xx,iscasecorr);

    if(rowsum1 == 0.0)
    {
      chisq = -1.0; Echisq = -1.0;
    }

    if(chisq >= 0.0) fprintf(fpout,"%.04f",chisq);
    else fprintf(fpout,"NA");
    if(Echisq >= 0.0) fprintf(fpout," %.04f\n",Echisq);
    else fprintf(fpout," NA\n");

    if(NSAMPLES*m > MAXSIZE)
    {
      fprintf(stderr,"OOPS genotype file has > %g genotypes\n",MAXSIZE);
      fprintf(fpout,"OOPS genotype file has > %g genotypes\n",MAXSIZE);
      exit(1);
    }
  }
  return 0;
}


void readcommands(int argc, char **argv) 
{
  int i;
  char *parname = NULL ;
  phandle *ph ;

  while ((i = getopt (argc, argv, "p:vV")) != -1) {

    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'v':
	printf("version: %s\n", WVERSION) ; 
	break; 

      case 'V':
	verbose = YES ;
	break; 

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
   pcheck(parname,'p') ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

   getint(ph, "packmode:", &packmode) ; // controls internals 
   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "pcaname:", &pcaname) ;
   getint(ph, "numpc:", &numpc) ;
   getint(ph, "numeigs:", &numpc) ;
   getint(ph, "qtmode:", &qtmode) ;
   getint(ph, "hashcheck:", &hashcheck) ;

   getstring(ph, "outputname:", &outputname);

   writepars(ph) ;
   closepars(ph) ;

}


void setinmode(inputmodetype *inmode, char *imode)  {
  char *ss = strdup(imode);
  int len = strlen(ss);

  int i;
  for(i=0;i<len;i++)  {
    ss[i] = tolower(ss[i]);
  }

  *inmode = EIGENSTRAT;			/* default */
  if ( strcmp(ss, "eigenstrat") == 0 )  
    *inmode = EIGENSTRAT;
  if ( strcmp(ss, "ped") == 0 )  
    *inmode = PED;
  if ( strcmp(ss, "packedped") == 0 )  
    *inmode = PACKEDPED;
  if ( strcmp(ss, "ancestrymap") == 0 )  
    *inmode = ANCESTRYMAP;
  if ( strcmp(ss, "packedancestrymap") == 0 )  
    *inmode = PACKEDANCESTRYMAP;

}

void readpcafile(double **Vp, int **outlierp, int *kp, int L, 
  int NSAMPLES)  {

  /* Build V[] and K and outlier[] */
  char *PCAFILE = strdup(pcaname);
  double tempdouble;
  int K;
  FILE *fppca;
  double *V;
  int *outlier;
  int noutlier;
  int x,n,k;

  if ((V = (double *) malloc(NSAMPLES*NSAMPLES*sizeof(*V))) == NULL)  
  {  fprintf(stderr, "CM\n");  exit(1);  }
  if ((outlier = (int *)malloc(NSAMPLES*sizeof(*outlier))) == NULL)
  {  fprintf(stderr,"CM\n");  exit(1);  }


  if( (fppca = fopen(PCAFILE, "r")) == NULL)
  {
    fprintf(stderr,"Could not open input file %s\n", PCAFILE);  
    exit(1);
  }

  fscanf(fppca,"%d",&K);
  if(L > K)
  {
    fprintf(stderr,"OOPS l=%d is larger than k=%d in %s\n",L,K,PCAFILE);
    fprintf(fpout,"OOPS l=%d is larger than k=%d in %s\n",L,K,PCAFILE);
    exit(1);
  }
  for(x=0; x<K; x++) fscanf(fppca,"%lf",&tempdouble); /* eigenvalues */
  for(n=0; n<NSAMPLES; n++)
  {
    for(k=0; k<K; k++) fscanf(fppca,"%lf",&V[NSAMPLES*n+k]);
    if(feof(fppca))
    {
      fprintf(stderr,"OOPS: %s contains less than %d times %d entries\n",
        PCAFILE,NSAMPLES,K);
      fprintf(fpout,"OOPS: %s contains less than %d times %d entries\n",
        PCAFILE,NSAMPLES,K);
      exit(1);
    }
    /* check for outliers */
    outlier[n] = 1;
    for(k=0; k<K; k++) { if(V[NSAMPLES*n+k] != 0.0) outlier[n] = 0; }
    if(outlier[n] == 1) noutlier++;
  }
  fscanf(fppca,"%lf",&tempdouble);
  if(!(feof(fppca)))
  {
    fprintf(stderr,"OOPS: %s contains too many entries\n",PCAFILE);
    fprintf(fpout,"OOPS: %s contains too many entries\n",PCAFILE);
    exit(1);
  }
  fclose(fppca);

  *kp = K;
  *Vp = V;
  *outlierp = outlier;
}

void getphenos(int NSAMPLES, double **iscasep, int *outlier, 
  double **iscasecorrp, int L, double *V)  {

  int k,n;
  double *iscase, *iscasecorr;
  double rowsum, rowsum1, gamma, denom ;

  /* allocate iscase */
  if ((iscase = (double *)malloc(NSAMPLES*sizeof(*iscase))) == NULL)  
  {  fprintf(stderr,"CM\n");  exit(1);  }
  if ((iscasecorr = (double *)malloc(NSAMPLES*sizeof(*iscasecorr))) == NULL)  
  {  fprintf(stderr,"CM\n");  exit(1);  }

  /* get phenotypes */
  for(n=0; n<NSAMPLES; n++)
  {

    if ( qtmode == NO )  {

      char *grp = strdup(indivmarkers[n]->egroup);
      int i;
      for(i=0;i<strlen(grp);i++)  
        grp[i] = tolower(grp[i]);

      if ( !strcmp(grp,"case") )  {
        iscase[n] = 1.0;
      }
      else if ( !strcmp(grp,"control") )  {
        iscase[n] = 0.0;
      }
      else if ( indivmarkers[n]->ignore == YES )  {
        iscase[n] = -100.0;
      }
      else  {
        fprintf(stderr,"OOPS bad phenotype %s\n",grp);
        fprintf(fpout,"OOPS bad phenotype %s\n",grp);
        exit(1);
      }
      free(grp);

    }
    else  {
      iscase[n] = indivmarkers[n]->qval;
    }
  }

  /* mean-adjust iscase */
  rowsum = 0.0; rowsum1 = 0.0;
  for(n=0; n<NSAMPLES; n++)
  {
    if((outlier[n] == 1) || (iscase[n] == -100) )
      continue;
    rowsum += iscase[n];
    rowsum1 += 1.0;
  }
  for(n=0; n<NSAMPLES; n++)
  {
    if(outlier[n]) continue;
    if(iscase[n] == -100.0) iscase[n] = -100.0; /* still keep track */
    else iscase[n] -= rowsum/rowsum1;
  }

  /* make iscasecorr */
  for(n=0; n<NSAMPLES; n++) 
  {
    if(outlier[n] == 0) iscasecorr[n] = iscase[n];
  }
  for(k=0; k<L; k++)
  {
    gamma = 0.0;
    denom = 0.0;
    for(n=0; n<NSAMPLES; n++) 
    {
      if(qtmode == NO   && ((outlier[n]) || (iscase[n] < -99.0))) continue;
      if(qtmode == YES  && ((outlier[n]) || (iscase[n] == -100.0))) continue;
      gamma += iscasecorr[n]*V[NSAMPLES*n+k];
      denom += V[NSAMPLES*n+k]*V[NSAMPLES*n+k];
    }
    gamma /= denom;
    for(n=0; n<NSAMPLES; n++) 
    {
      if(qtmode == NO   && ((outlier[n]) || (iscase[n] < -99.0))) continue;
      if(qtmode == YES  && ((outlier[n]) || (iscase[n] == -100.0))) continue;
      iscasecorr[n] -= gamma*V[NSAMPLES*n+k];
    }
  }
  *iscasep = iscase;
  *iscasecorrp = iscasecorr;
}

double compute_chisq(double *source, double *target)
{
  int n;
  double sum1, sumx, sumxx, sumy, sumyy, sumxy, numer, denom1, denom2;
  double corr;

  sum1 = 0.0; sumx = 0.0; sumxx = 0.0; sumy = 0.0; sumyy = 0.0; sumxy = 0.0;
  for(n=0; n<NSAMPLES; n++)
  {
    if(outlier[n]) continue;
    if(qtmode == NO && source[n] < -99.0) continue;
    if(qtmode == YES && source[n] == -100.0) continue;

    if(qtmode == NO && target[n] < -99.0) continue;
    if(qtmode == YES && target[n] == -100.0) continue;

    sumx += source[n];
    sumxx += source[n]*source[n];
    sumy += target[n];
    sumyy += target[n]*target[n];
    sumxy += source[n]*target[n];
    sum1 += 1.0;
  }
  if(sumxx == 0.0) return -1.0;
  if(sumyy == 0.0) return -1.0;
  numer = sumxy/sum1 - (sumx/sum1)*(sumy/sum1);
  denom1 = (sumxx/sum1 - (sumx/sum1)*(sumx/sum1));
  denom2 = (sumyy/sum1 - (sumy/sum1)*(sumy/sum1));
  if(denom1 <= 0.0) return -1.0;
  if(denom2 <= 0.0) return -1.0;

  corr = (numer/sqrt(denom1*denom2));
  return (sum1*corr*corr);
}

double compute_chisqE(double *source, double *target)
{
  int n;
  double sum1, sumx, sumxx, sumy, sumyy, sumxy, numer, denom1, denom2;
  double corr;

  sum1 = 0.0; sumx = 0.0; sumxx = 0.0; sumy = 0.0; sumyy = 0.0; sumxy = 0.0;
  for(n=0; n<NSAMPLES; n++)
  {
    if(outlier[n]) continue;
    if(qtmode == NO && source[n] < -99.0) continue;
    if(qtmode == YES && source[n] == -100.0) continue;

    if(qtmode == NO && target[n] < -99.0) continue;
    if(qtmode == YES && target[n] == -100.0) continue;

    sumx += source[n];
    sumxx += source[n]*source[n];
    sumy += target[n];
    sumyy += target[n]*target[n];
    sumxy += source[n]*target[n];
    sum1 += 1.0;
  }
  if(sumxx == 0.0) return -1.0;
  if(sumyy == 0.0) return -1.0;
  numer = sumxy/sum1 - (sumx/sum1)*(sumy/sum1);
  denom1 = (sumxx/sum1 - (sumx/sum1)*(sumx/sum1));
  denom2 = (sumyy/sum1 - (sumy/sum1)*(sumy/sum1));
  if(denom1 <= 0.0) return -1.0;
  if(denom2 <= 0.0) return -1.0;

  corr = (numer/sqrt(denom1*denom2));
  sum1 = sum1 - ((double)(L+1));
  return (sum1*corr*corr);
}

