#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <pthread.h>
#include <globals.h>
#include <workqueue.h>

#include <nicklib.h>
#include <getpars.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  
#include "eigsubs.h"  
#include "egsubs.h"  
#include "qpsubs.h" 


#define WVERSION   "9003" 
/** 
Simple eigenvector analysis
Options to look at groups (simple ANOVA)
Weights allowed for individuals 
missing mode
dotpops added
recompiled with new twtail. Output form at changed
Cleaned up twestxx
fancynorm mode (divide by sqrt(p*(1-p)) 
poplistname supported.  Eigenanalysis just on individuals in population 
But all individuals figure in eigenvector output
New way of computing effective marker size  (twl2mode)
popdifference implemented
ldregression  ldlimit (genetic distance in Morgans)
nostatslim added  
dotpop has new format if many groups
uses new I/O
Supports packmode
Alkes style outlier removal added
Only half XTX computed
xdata (huge array) removed

fst calculation added 
popsizelimit added
divergence added (not useful?) 

SNPs discarded if no data.
Phylipfile now supported

Preparations for parallelization made
Various fixups for EIGENSTRAT and altnormstyle

output capability added (like convertf)

bug fixed (a last iteration needed for outlier removal)
bug fixed (numindivs unlimited)
output files fixed up (NULL OK)

Many Alkes style options added
Support for outliername added (outlier info)
familyname added (ped files)

bugfix: jackrat dies (outlier removes all of population)
bugfix: See ~ap92/REICH/FALL06/EIGENSOFTCODE/bugfix  bugs 1, 3 fixed

nrows, ncols output added
nrows, ncols set each outlier iteration
indivs with no data removed

writesnpeig added 

bugfix: popsize of 1 no anova done
minallelecnt added
chrom:  added
latest greatest handling of chromosome number added.
bad bugfix: numvalidgtypes 

checksizemode added.  Bug fix to version 7520 (> 2B genotypes fixed)
pubmean added

fst on X
fst std errors now fixed

bad bug fixed (outfiles changed indivmarkers)  ...  

fstdetailname added
fsthiprecision added
bug fixed  (getrawcolx)

version compatible with Mac

*/


#define MAXFL  50   
#define MAXSTR  512
#define MAXPOPS 100

/* Initialize default thread parameters to signify "no threads" */
#define THREADING_DEFAULT 0
#define HAVE_PTHREAD 1

char *parname = NULL ;
char *twxtabname = NULL ;
char *trashdir = "/var/tmp" ;
int verbose = NO ;
int qtmode = NO ;
int outnum = 0 ;         /* only used in qpsubs::setwt, not called in EIG */
Indiv **indivmarkers;
SNP **snpmarkers ;

int numsnps, numindivs ; 
int numeigs = 10 ;  /// default
int markerscore = NO ;
int seed = 0 ;
int chisqmode = NO ;  // approx p-value better to use F-stat
int missingmode = NO ;
int plotmode = NO ;
int fancynorm = YES ;
int dotpopsmode = YES ;
int noxdata = YES ;  /* default as pop structure dubious if Males and females */
//int fastmode = NO ;
int fstonly = NO ; 
int pcorrmode = NO ;
int pcpopsonly = YES ;
int nostatslim = 10 ;
int znval = -1 ;
int popsizelimit = -1 ;
int altnormstyle = YES ;  // affects subtle details in normalization formula
int minallelecnt = 1 ;
int maxmissing = 9999999 ;
int lopos = -999999999, hipos = 999999999 ;  // use with xchrom

/* threading parameters */
int numthreads = THREADING_DEFAULT;
int numxtxblocksperside = THREADING_DEFAULT;
int numsnppartitions = THREADING_DEFAULT;

int packout = -1 ;
extern enum outputmodetype outputmode  ;
extern int checksizemode ;
extern int packmode ;
int ogmode = NO ;
int fsthiprec = NO ;
int inbreed = NO ;  // for fst
int easymode = NO ;

int numoutliter = 5, numoutleigs = 10 ;  
double outlthresh = 6.0 ; 
OUTLINFO  **outinfo ;
char *outinfoname = NULL ;  
char *fstdetailsname = NULL ;  


double plo = .001 ;
double phi = .999 ;
double pvhit = .001 ;
double pvjack = 1.0e-6 ;
double *chitot ;
int    *xpopsize ;

char  *genotypename = NULL ;
char  *snpname = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *deletesnpoutname = NULL ;
char *poplistname = NULL ;
char *xregionname = NULL ;     /* physical positions of SNPs to exclude */
char *outliername = NULL ;
char *phylipname = NULL ;
char *snpeigname = NULL ;

char *indoutfilename = NULL ;
char *snpoutfilename = NULL ;
char  *genooutfilename = NULL ;
char  *omode = "packedancestrymap" ;
double blgsize = 0.05 ;  // block size in Morgans */

double r2thresh = -1.0 ;  
double r2genlim = 0.01 ; // Morgans 
double r2physlim = 5.0e6 ;     
int killr2 = NO ;  
int pubmean = YES ;  // change default

int randomfillin = NO ;
int usepopsformissing = NO ; // if YES popmean is used for missing.  Overall mean if all missing for pop

int xchrom = -1 ;
// list of outliers

int ldregress = 0 ;
double ldlimit = 9999.0 ;  /* default is infinity */
/* we only consider markers as in possible LD if gdis <= ldlimit */

char *outputname = NULL ;
char *outputvname = NULL ;
char *weightname = NULL ;
FILE *ofile, *ovfile ;

double twestxx(double *lam, int m, double *pzn,  double *pzvar)  ;
double twnorm(double lam, double m, double n) ;

void readcommands(int argc, char **argv) ;
int loadindx(Indiv **xindlist, int *xindex, Indiv **indivmarkers, int numindivs) ;
void loadxdataind(double *xrow, SNP **snplist, int ind,  int ncols) ;           
void fixxrow(double *xrow, double *xmean, double *xfancy, int len)  ;
void dofancy(double *cc, int n, double *fancy) ;
int fvadjust(double *rr, int n, double *pmean, double *fancy)  ;
void getcol(double *cc, double *xdata, int col, int nrows, int ncols)  ;
void getcolxf(double *xcol, SNP *cupt, int *xindex, 
  int nrows, int col, double *xmean, double *xfancy)  ;          
int getcolxz(double *xcol, SNP *cupt, const int *xindex, int *xtypes, 
  int nrows, int col, double *xmean, double *xfancy, int *n0, int *n1)  ;          

double doinbxx(double *inbans, double *inbsd, SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg, double blgsize, SNP **snpmarkers, Indiv **indm) ;

void putcol(double *cc, double *xdata, int col, int nrows, int ncols)  ;
void calcpopmean(double *wmean, char **elist, double *vec, 
 char **eglist, int numeg, int *xtypes, int len) ;
double dottest(char *sss, double *vec, char **eglist, int numeg, int *xtypes, int len) ;
double yll(double x1, double x2, double xlen) ;
void calcmean(double *wmean, double *vec, int len, int *xtypes, int numeg)  ;
double anova1(double *vec, int len, int *xtypes, int numeg) ;
double anova(double *vec, int len, int *xtypes, int numeg) ;
void publishit(char *sss, int df, double chi) ;

void setmiss(SNP **snpm, int numsnps)  ;
void setfvecs(double *fvecs, double *evecs, int nrows, int numeigs) ;
void dotpops(double *X, char **eglist, int numeg, int *xtypes, int nrows)  ;
void printxcorr(double *X, int nrows, Indiv **indxx)  ;      

void ldreg(double *ldmat, double *ldmat2, 
  double *vv, double *vv2, double *ldvv, 
  double *ldvv2, int rsize, int n)  ;

void clearld(double *ldmat, double *ldvv, int rsize, int n, int nclear)  ;


void fixrho(double *a, int n) ;
void printdiag(double *a, int n) ;

int
ridoutlier(double *evecs, int n, int neigs, 
 double thresh, int *badlist, OUTLINFO **outinfo) ;

void addoutersym(double *X, double *v, int n)  ;
void symit(double *X, int n)  ;

double fstcol(double *estn, double *estd, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2) ;

double oldfstcol(double *estn, double *estd, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2) ;

void jackrat(double *xmean, double *xsd, double *top, double *bot,  int len)  ;
void writesnpeigs(char *snpeigname, SNP **xsnplist, double *ffvecs, int numeigs, int ncols)  ;
double dofstxx(double *fstans, double *fstsd, SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg, double blgsize, SNP **snpmarkers, Indiv **indm) ; 
void fixwt(SNP **snpm, int nsnp, double val) ;

/* Hanna:  represents one partial computation of a portion of the XTX matrix.
 * Computes the block of the XTX matrix bounded by column indices [XTX_col_start, XTS_col_Stop]
 * and row indices [XTX_row_start,XTX_row_stop].  Will take into account SNP data points in the
 * range [SNP_col_start,SNP_col_stop].
 */

typedef struct block_t  {
  // temporary buffer holding just the data processed by this partial computation.
  double *XTX_block;

  // region from which this block originates
  int XTX_col_start, XTX_col_stop, XTX_row_start, XTX_row_stop;

  // SNP columns in the calculation of this block
  int SNP_col_start, SNP_col_stop;

  // input data source for the block
  SNP **xsnplist;
  const int *xindex;
  int *xtypes;
  int weightmode;
  int nrows;
  int ncols;

  // output data and side effect - computation of xmean and xfancy for each SNP column
  double *xmean, *xfancy;

} block;

void compute_XTX( double *XTX, SNP **xsnplist, const int *xindex, double *xmean, double *xfancy, const int weightmode,
  const int nrows, const int ncols, int *xtypes );
void *compute_block( void *block );
void collect_output( double *XTX, const work_queue *work_queue, const block *blocks, const int num_threads );
void partition( const int size, const int n, const int num_partitions, int *start, int *stop );
void domult( double *cc, const block *block );



double nhwfilter = -1;

int main(int argc, char **argv)
{

  char sss[MAXSTR] ;
  int **snppos ;
  int *snpindx ;
  char **snpnamelist, **indnamelist ;
  char **eglist ;
  int  lsnplist, lindlist, numeg ;
  int i, j, k, k1, k2, pos; 
  int *vv ;
  SNP *cupt, *cupt1, *cupt2, *cupt3 ;
  Indiv *indx ;
  double y1, y2, y ;
  FILE *twxtestfp;

  int ch1, ch2 ;
  int fmnum , lmnum ;
  int num, n0, n1, nkill ;

  int nindiv = 0, e, f, lag=1  ;
  double xc[9], xd[4], xc2[9] ;
  double ychi,  tail, tw ;
  int nignore, numrisks = 1 ;
  double  *xrow, *xpt ; 
  SNP **xsnplist  ;
  Indiv **xindlist ;
  int *xindex, *xtypes = NULL ;
  int nrows, ncols, m, nused ;
  double *XTX, *cc, *evecs, *ww, weight, *qvec, *qcoord ; 
  double *lambda ;
  double *tvecs ;
  double zn, zvar ;
  double *fvecs, *fxvecs, *fxscal ;
  double *ffvecs ;
  int weightmode = NO ;
  double chisq, ynrows ;
  int *numhits, t, g, tt ;  
  double *xmean, *xfancy ;
  double *ldmat, *ldmat2, *ldvv, *ldvv2, *vv2  ;
  double *fstans, *fstsd ; 
  double *inbans, *inbsd ; 

  int chrom,  numclear ;
  double gdis ;
  int outliter, numoutiter, *badlist, nbad ;
  int a, b, kmax, kmin ;
  FILE *outlfile, *phylipfile  ;
  double **eigmoment, **eigindmoment ;
  double *eigkurt, *eigindkurt ;
  double *snpsc ; 
  double *wmean ;
  char **elist ; 
  

  int xblock, blocksize=10000 ;   
  double *tblock ;  

  OUTLINFO *outpt ;

  readcommands(argc, argv) ;
  printf("## smartpca version: %s\n", WVERSION) ;
  packmode = YES ;
  setomode(&outputmode, omode) ;

  if (parname == NULL) return 0 ;
  if (xchrom == (numchrom+1)) noxdata = NO ;

  if (usepopsformissing) easymode = YES ;

  if (deletesnpoutname != NULL)  {   /* remove because snplog opens in append mode */ 
    char buff[256];
    sprintf(buff,"rm -f %s", deletesnpoutname);
    system(buff);
  } 

  if (fstonly) { 
   printf("fstonly\n") ;
   numeigs = 0 ; 
   numoutliter = 0 ;
   numoutiter = 0 ;
   outputname = NULL ;
   snpeigname = NULL ;
  }

  if (fancynorm) printf("norm used\n\n") ;
  else printf("no norm used\n\n") ;

  nostatslim = MAX(nostatslim, 3) ;

  outlfile = ofile = stdout; 

  if (outputname != NULL)  openit(outputname, &ofile, "w") ;
  if (outliername != NULL) openit(outliername, &outlfile, "w") ;
  if (fstdetailsname != NULL) openit(fstdetailsname, &fstdetails, "w") ;

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;


  numindivs = getindivs(indivname, &indivmarkers) ;
  k = getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, numindivs, nignore) ;


  if (poplistname != NULL) 
  { 
    ZALLOC(eglist, numindivs, char *) ; 
    numeg = loadlist(eglist, poplistname) ;
    seteglist(indivmarkers, numindivs, poplistname);
  }
  else
  {
    setstatus(indivmarkers, numindivs, NULL) ;
    ZALLOC(eglist, MAXPOPS, char *) ;
    numeg = makeeglist(eglist, MAXPOPS, indivmarkers, numindivs) ;
  }
  for (i=0; i<numeg; i++) 
  {  
    /* printf("%3d %s\n",i, eglist[i]) ; */
  }

  nindiv=0 ;
  for (i=0; i<numindivs; i++) 
  {
    indx = indivmarkers[i] ;
    if(indx -> affstatus == YES) ++nindiv  ;
  }

  for (i=0; i<numsnps; i++)  
  {  
    cupt = snpmarkers[i] ; 
    chrom = cupt -> chrom ;
    if ((noxdata) && (chrom == (numchrom+1))) {
      cupt-> ignore = YES ;
      logdeletedsnp(cupt->ID,"chrom-X",deletesnpoutname);
    }
    if (chrom == 0) {
      cupt -> ignore = YES ;
      logdeletedsnp(cupt->ID,"chrom-0",deletesnpoutname);
    }
    if (chrom > (numchrom+1))  {
      cupt -> ignore = YES ;
      logdeletedsnp(cupt->ID,"chrom-big",deletesnpoutname);
    }
  }

  for (i=0; i<numsnps; i++)  
  {
    cupt = snpmarkers[i] ; 
    pos = nnint(cupt -> physpos) ;
    if ((xchrom>0) && (cupt -> chrom != xchrom))  {
      cupt -> ignore = YES ;
      logdeletedsnp(cupt->ID,"not-chrom",deletesnpoutname);
    }
    if ((xchrom > 0) && (pos < lopos))   {
      cupt -> ignore = YES ;
      logdeletedsnp(cupt->ID,"lopos",deletesnpoutname);
    }
    if ((xchrom > 0) && (pos > hipos))   {
      cupt -> ignore = YES ;
      logdeletedsnp(cupt->ID,"hipos",deletesnpoutname);
    }
    if (cupt -> ignore) continue ;
    if (numvalidgtx(indivmarkers, cupt, YES) <= 1) 
    { 
      printf("nodata: %20s\n", cupt -> ID) ;
      cupt -> ignore = YES ;
      logdeletedsnp(cupt->ID,"nodata",deletesnpoutname);
    }
  }

  if (killr2) {
   nkill = killhir2(snpmarkers, numsnps, numindivs, r2physlim, r2genlim, r2thresh) ;
   if (nkill>0) printf("killhir2.  number of snps killed: %d\n", nkill) ;
  }

  if ( xregionname )  {
    excluderegions(xregionname, snpmarkers, numsnps, deletesnpoutname);
  }

  if ( nhwfilter > 0 )  {
    hwfilter(snpmarkers, numsnps, numindivs, nhwfilter, deletesnpoutname);
  }

  ZALLOC(vv, numindivs, int) ;
  numvalidgtallind(vv, snpmarkers, numsnps,  numindivs) ; 
  for (i=0; i<numindivs; ++i)  { 
  if (vv[i] == 0) {
    indx = indivmarkers[i] ;
    indx -> ignore = YES ; 
   }
  }
  free(vv) ;

  numsnps = rmsnps(snpmarkers, numsnps, deletesnpoutname) ;  //  rid ignorable snps
   
  if (missingmode) 
  {
    setmiss(snpmarkers, numsnps) ;
    fancynorm = NO ;
  }

  if  (weightname != NULL)   
  {  
    weightmode = YES ;
    getweights(weightname, snpmarkers, numsnps) ;
  }

  ZALLOC(xindex, numindivs, int) ;
  ZALLOC(xindlist, numindivs, Indiv *) ;
  ZALLOC(xsnplist, numsnps, SNP *) ;

  if (popsizelimit > 0) 
  {  
    setplimit(indivmarkers, numindivs, eglist, numeg, popsizelimit) ; 
  }


  /*  Load non-ignored individuals into xindlist,xindex:
   *  xindex[i]   = index into indivmarkers
   *  xindlist[i] = pointer to Indiv struct   */

  /*  Load non-ignored SNPs into xsnplist:
   *  xsnplist[i] = pointer to SNP struct     */



  nrows = loadindx(xindlist, xindex, indivmarkers, numindivs) ;
  ncols = loadsnpx(xsnplist, snpmarkers, numsnps, indivmarkers) ;

  ZALLOC(xtypes, nrows, int) ;

  /** DEBUGGING:
  cupt = xsnplist[0] ;
  for (j=0; j<nrows; ++j) {  
   k = xindex[j] ;
   g = getgtypes(cupt, k) ;
   indx = indivmarkers[k] ;
   t = indxindex(eglist, numeg, indx -> egroup) ;
   printf("yy1 %20s %20s %20s %d %d %d\n", cupt ->ID, indx -> ID, indx -> egroup, j, k, g) ;
  }
  printf("yya: ") ; printimat(xindex, 1, nrows) ;
  printf("zzindxa:  %s\n", indivmarkers[230] -> egroup) ;
  **/


  /* printf("## nrows: %d  ncols  %d\n", nrows, ncols) ; */
  ZALLOC(xmean, ncols, double) ;
  ZALLOC(xfancy, ncols, double) ;

  ZALLOC(XTX, nrows*nrows, double) ;
  ZALLOC(evecs, nrows*nrows, double) ;
  ZALLOC(tvecs, nrows*nrows, double) ;

  ZALLOC(lambda, nrows, double) ;
  ZALLOC(cc, nrows, double) ;
  ZALLOC(tblock, nrows*blocksize, double) ;
  ZALLOC(ww, nrows, double) ;
  ZALLOC(badlist, nrows, int) ;

  blocksize = MIN(blocksize, ncols) ; 

  // xfancy is multiplier for column xmean is mean to take off
  // badlist is list of rows to delete (outlier removal) 

  numoutiter = 1 ;  

  if (numoutliter>=1) 
  {
    numoutiter = numoutliter+1 ;
    ZALLOC(outinfo, nrows,  OUTLINFO *) ;  
    for (k=0; k<nrows; k++) 
    {  
      ZALLOC(outinfo[k], 1, OUTLINFO) ;
    }
    /* fprintf(outlfile, "##%18s %4s %6s %9s\n", "ID", "iter","eigvec", "score") ; */
  }

  nkill = 0 ;

  for (outliter = 1; outliter <= numoutiter ; ++outliter)  {

    if (fstonly) { 
     setidmat(XTX, nrows) ;
     vclear(lambda, 1.0, nrows) ;
     break ;
    }
    if (outliter>1) {
     ncols = loadsnpx(xsnplist, snpmarkers, numsnps, indivmarkers) ;
    }

// -------------------------------------------------------------------------------------
    compute_XTX(XTX, xsnplist, xindex, xmean, xfancy, weightmode, nrows, ncols, xtypes);
// ----------------------------------------------------------------------------------

    y = trace(XTX, nrows) / (double) (nrows-1) ;
    if (isnan(y)) fatalx("bad XTX matrix\n") ;
    printf("trace:  %9.3f\n", y) ;
    if (y<=0.0) fatalx("XTX has zero trace (perhaps no data)\n") ;
    vst(XTX, XTX, 1.0/y, nrows * nrows) ;

    eigvecs(XTX, lambda, evecs, nrows) ;

// eigenvalues are in decreasing order 
//
    if (outliter > numoutliter) break ;  
    // last pass skips outliers 
    numoutleigs = MIN(numoutleigs, nrows-1) ;
    nbad = ridoutlier(evecs, nrows, numoutleigs, outlthresh, badlist, outinfo) ;
    if (nbad == 0) break ; 
    for (i=0; i<nbad; i++) 
    {  
      j = badlist[i] ;
      indx = xindlist[j] ;
      outpt = outinfo[j] ;
      fprintf(outlfile, "REMOVED outlier %s iter %d evec %d sigmage %.3f\n", indx -> ID, outliter, outpt -> vecno, outpt -> score) ;
      indx -> ignore = YES ;
    }
    nrows = loadindx(xindlist, xindex, indivmarkers, numindivs) ;
    printf("number of samples after outlier removal: %d\n", nrows) ;
  }

  if (outliername != NULL) fclose(outlfile) ;

  m = numgtz(lambda, nrows)  ;
  /* printf("matrix rank: %d\n", m) ; */
  if (m==0) fatalx("no data\n") ;

  /* Now, print Tracy-Widom stats, if twtable is valid */
   if (settwxtable(twxtabname)<0) 
  {
    printf("\n## To get Tracy-Widom statistics: recompile smartpca with");
    printf(" TWTAB correctly specified in Makefile, or\n");
    printf("   just run twstats (see README file in POPGEN directory)\n");
  }
  else
  {
    /* *** START of code to print Tracy-Widom statistics */
    printf("\n## Tracy-Widom statistics: rows: %d  cols: %d\n", nrows, ncols);
    y = -1.0 ;
    printf("%4s  %12s", "#N", "eigenvalue") ;
    printf("%12s",  "difference") ;
    printf(" %9s %12s", "twstat", "p-value") ;
    printf(" %9s", "effect. n") ;
    printf("\n") ;

    ynrows = (double) nrows ;  

    for (i=0; i<m; ++i) { 
      if (fstonly) break ;
      zn = znval ;
      if (zn>0) zn = MAX(ynrows, zn) ;
      tail = dotwcalc(lambda+i, m-i, &tw, &zn, &zvar, nostatslim) ;
      printf("%4d  %12.6f", i+1, lambda[i]) ;
      if (i==0) printf( "%12s", "NA") ;
      else printf("%12.6f", lambda[i]-lambda[i-1]) ;
      if (tail>=0.0) printf( " %9.3f %12.6g", tw, tail) ;
      else printf( " %9s %12s", "NA", "NA") ;
      if (zn>0.0) 
      {
        printf( " %9.3f", zn) ;
      }
      else 
      {
        printf( " %9s", "NA") ;
      } 
      printf( "\n") ;
    }
    /* END of code to print Tracy-Widom statistics */
  }


  numeigs = MIN(numeigs, nrows) ;
  numeigs = MIN(numeigs, ncols) ;

  /* fprintf(ofile, "##genotypes: %s\n", genotypename) ; */
  /* fprintf(ofile, "##numrows(indivs):: %d\n", nrows) ; */
  /* fprintf(ofile, "##numcols(snps):: %d\n", ncols) ; */
  /* fprintf(ofile, "##numeigs:: %d\n", numeigs) ; */
  fprintf(ofile, "%20s ", "#eigvals:") ;
  for (j=0; j<numeigs; j++) { 
	  fprintf(ofile, "%9.3f ", lambda[j]) ;
  }
  fprintf(ofile, "\n") ;

  if (outputvname != NULL)  {  
    openit(outputvname, &ovfile, "w") ;
    for (j=0; j<nrows; j++) { 
      fprintf(ovfile, "%12.6f\n", lambda[j]) ;
    }
    fclose(ovfile) ;
  }


  ZALLOC(fvecs, nrows*numeigs, double) ;
  ZALLOC(fxvecs, nrows*numeigs, double) ;
  ZALLOC(fxscal, numeigs, double) ;

  ZALLOC(ffvecs, ncols*numeigs, double) ;
  ZALLOC(xrow, ncols, double) ;
  setfvecs(fvecs, evecs, nrows, numeigs) ;

  if (easymode) { 
   for (j=0; j<numeigs; j++) { 
    xpt = fvecs+j*nrows ;
    y = asum2(xpt, nrows) ;
    vst(xpt, xpt, 1.0/sqrt(y), nrows) ;  // norm 1
   }
   for (i=0; i < nrows ; i++)  { 
          indx = xindlist[i] ;
	  fprintf(ofile, "%20s ", indx -> ID) ;
	  for (j=0; j<numeigs; j++) { 
           xpt = fvecs+j*nrows ;
           y = xpt[i] ;
           fprintf(ofile, "%10.4f  ", y) ;
	  }
	  fprintf(ofile, "%15s\n", indx -> egroup) ;
   }
    if (pubmean) { 

     ZALLOC(wmean, numeg, double) ;
     ZALLOC(elist, numeg, char *) ;

     for (j=0; j<numeigs; j++) { 
        xpt = fvecs+j*nrows ;
        calcpopmean(wmean, elist, xpt,  eglist, numeg, xtypes, nrows) ;
        printf ("eig: %d ", j+1) ;
        printf("min: %s %9.3f  ", elist[0], wmean[0]) ;
        printf("max: %s %9.3f  ", elist[numeg-1], wmean[numeg-1]) ;
        printnl() ;
        for (k=0; k<numeg; ++k) {
         printf("%20s ", elist[k]) ;
         printf(" %9.3f\n", wmean[k]) ;
        }
      }
    }
    
   printf("## easymode set. end of smartpca run\n") ;
   return 0 ; 
  }



  for (i=0; i<ncols; i++)  { 
    cupt = xsnplist[i] ;
    getcolxf(cc, cupt, xindex, nrows, i, NULL, NULL) ;
    for (j=0; j<numeigs; j++)  { 
     for (k=0; k<nrows; k++)  {
      ffvecs[j*ncols+i] += fvecs[j*nrows+k]*cc[k] ;
     }
    }
  }

  eigmoment    = initarray_2Ddouble(numeigs, 5, 0.0) ;
  eigindmoment = initarray_2Ddouble(numeigs, 5, 0.0) ;

  ZALLOC(eigkurt, numeigs, double) ;
  ZALLOC(eigindkurt, numeigs, double) ;

  for (j=0; j<numeigs; ++j) {  
   eigkurt[j] = kurtosis(ffvecs+j*ncols, ncols) ;
   eigindkurt[j] = kurtosis(fvecs+j*nrows, nrows) ;
  }

  for (i=0; i<nrows; i++) {

   indx = xindlist[i] ;
   k = indxindex(eglist, numeg, indx -> egroup) ;
   xtypes[i] = k ;

   loadxdataind(xrow, xsnplist, xindex[i],  ncols) ;            
   fixxrow(xrow, xmean, xfancy, ncols) ;


   for (j=0; j<numeigs; j++) { 

     xpt = ffvecs+j*ncols ;
     y = fxvecs[j*nrows+i] = vdot(xrow, xpt, ncols) ;
     fxscal[j] += y*y ;
   
   }
  }


  for (j=0; j<numeigs; j++) { 
    y = fxscal[j] ;                                            
//  fxscal[j] = 10.0/sqrt(y) ; // eigenvectors have norm 10 (perhaps eccentric)
    fxscal[j] = 1.0/sqrt(y) ;  // standard
  }

  for (i=0; i < numindivs ; i++)  { 
          indx = indivmarkers[i] ;
          if (indx -> ignore) continue ;
          loadxdataind(xrow, xsnplist, i,  ncols) ;            
          fixxrow(xrow, xmean, xfancy, ncols) ;
	  fprintf(ofile, "%20s ", indx -> ID) ;
	  for (j=0; j<numeigs; j++) { 
           y = fxscal[j]*vdot(xrow, ffvecs+j*ncols, ncols) ;


           fprintf(ofile, "%10.4f  ", y) ;
	  }
	  if ( qtmode )  {
	    fprintf(ofile, "%15.6e\n", indx -> qval) ;
	  }
	  else  {
	    fprintf(ofile, "%15s\n", indx -> egroup) ;
	  }
  }

 printf("%12s %4s %9s %9s\n", "kurtosis", "", "snps", "indivs") ;

  for (j=0; j<numeigs; ++j)  {  
   y1 = eigkurt[j] ;              
   y2 = eigindkurt[j] ;              
   printf("%12s %4d %9.3f %9.3f\n", "eigenvector", j+1, y1, y2) ;
  }


// output files
  settersemode(YES) ;

  ZALLOC(xpopsize, numeg, int) ;
  for (i=0; i<nrows; i++) { 
    k = xtypes[i] ;
    ++xpopsize[k] ;
  }

  for (i=0; i<numeg; i++) 
  {  
    printf("population: %3d %20s %4d\n",i, eglist[i], xpopsize[i]) ; 
  }


  if (numeg==1) dotpopsmode = NO ;

  if (dotpopsmode == NO)  {
    writesnpeigs(snpeigname, xsnplist, ffvecs, numeigs, ncols) ;
    printxcorr(XTX, nrows, xindlist) ;
    outfiles(snpoutfilename, indoutfilename, genooutfilename, 
     snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode) ;

    printf("##end of smartpca run\n") ;
    return 0 ;
  }

  ZALLOC(chitot, numeg*numeg, double) ;

  dotpops(XTX, eglist, numeg, xtypes, nrows) ;
  ZALLOC(fstans, numeg*numeg, double) ;
  ZALLOC(fstsd , numeg*numeg, double) ;

  setinbreed(inbreed) ;

  if (inbreed) {  
   ZALLOC(inbans, numeg, double) ;
   ZALLOC(inbsd , numeg, double) ;
   doinbxx(inbans, inbsd, xsnplist, xindex, xtypes, 
     nrows, ncols, numeg, blgsize, snpmarkers, indivmarkers) ;
     printf("## inbreeding coeffs:   inbreed    std error\n");
     for (k1=0; k1<numeg; ++k1) {  
      printf(" %20s %10.4f %10.4f\n", eglist[k1],  
       inbans[k1], inbsd[k1]) ;
     }
     free(inbans) ; 
     free(inbsd) ;
  }

  dofstxx(fstans, fstsd, xsnplist, xindex, xtypes, 
    nrows, ncols, numeg, blgsize, snpmarkers, indivmarkers); 

  if ((phylipname == NULL)  && (numeg>10)){
    printf("## Fst statistics between populations:         fst       std error\n");
    for (k1=0; k1<numeg; ++k1) {  
      for (k2=k1+1; k2<numeg; ++k2) {  
        printf(" %20s %20s %9.3f %10.4f\n", eglist[k1], eglist[k2], 
        fstans[k1*numeg+k2], fstsd[k1*numeg+k2]) ;
      }
    }
    printf("\n");
  }
  if (fstdetailsname != NULL)  {
    printf("## Fst statistics between populations:         fst       std error\n");
    for (k1=0; k1<numeg; ++k1) {  
      for (k2=k1+1; k2<numeg; ++k2) {  
        fprintf(fstdetails, " %20s %20s %12.6.3f %12.6f\n", eglist[k1], eglist[k2], 
        fstans[k1*numeg+k2], fstsd[k1*numeg+k2]) ;
      }
    }
    fprintf(fstdetails, "\n");
  }
  
  if (phylipname != NULL) { 
    openit(phylipname, &phylipfile, "w") ;
    fprintf(phylipfile, "%6d\n",numeg) ; 
    sss[10] = CNULL ;
    for (k1=0; k1<numeg; ++k1) { 
      strncpy(sss, eglist[k1], 10) ;  
      fprintf(phylipfile, "%10s", sss) ;
      for (k2=0; k2<numeg; ++k2) {  
        y1 = fstans[k1*numeg+k2] ; 
        y2 = fstans[k2*numeg+k1] ; 
        fprintf(phylipfile, "%6.3f", (0.5*(y1+y2))) ;
      }
      fprintf(phylipfile, "\n") ;
    }
    fclose(phylipfile) ;
  }

  if ((numeg<=10) || fstonly) {
    if (fsthiprec == NO) {
      printf("fst *1000:") ;
      printnl() ;
      printmatz5(fstans, eglist, numeg) ; 
      printnl() ;
    }
    if (fsthiprec == YES) {
      printf("fst *1000000:") ;
      printnl() ;
      printmatz10(fstans, eglist, numeg) ; 
      printnl() ;
    }
  }
  printf("s.dev * 1000000:\n") ;
  vst(fstsd, fstsd, 1000.0, numeg*numeg) ;
  printmatz5(fstsd, eglist, numeg) ; 
  printnl() ;
  fflush(stdout) ;
  if (fstonly) {
    printf("##end of smartpca run\n") ;
    return 0 ;
  }
  vst(fstsd, fstsd, 1.0/1000.0, numeg*numeg) ;
  
  for (j=0; j< numeigs; j++)  {  
   sprintf(sss, "eigenvector %d", j+1) ;
   y=dottest(sss, evecs+j*nrows, eglist, numeg, xtypes, nrows) ;
  }

  printf("\n## Statistical significance of differences beween populations:\n");
  printf("                                pop1                  pop2      chisq          p-value   |pop1|   |pop2|\n");
  for (k1=0; k1<numeg; ++k1) {  
   if (fstonly) break ;
   for (k2=k1+1; k2<numeg; ++k2) {  
    ychi = chitot[k1*numeg+k2] ;
    tail = rtlchsq(numeigs, ychi) ;
    printf("popdifference:  %20s  %20s  %12.3f  %12.6g", eglist[k1], eglist[k2], ychi, tail) ;
    printf ("   %5d", xpopsize[k1]) ;
    printf ("   %5d", xpopsize[k2]) ;
    printf("\n") ;
   }
  }
  printf("\n");
  for (i=0; i<ncols; i++)  {  
   if (markerscore == NO) break;
   cupt = xsnplist[i] ;
   getcolxf(cc, cupt, xindex, nrows, i, NULL, NULL) ;
   sprintf(sss, "%s raw", cupt -> ID) ;
   dottest(sss, cc, eglist, numeg, xtypes, nrows) ;
   for (j=0; j< numeigs; j++)  {  
    sprintf(sss, "%s subtract sing vec %d", cupt ->ID, j+1) ;
    y = vdot(cc, evecs+j*nrows, nrows) ;
    vst(ww, evecs+j*nrows, y, nrows) ;
    vvm(cc, cc, ww, nrows) ;
    dottest(sss, cc, eglist, numeg, xtypes, nrows) ;
   }
  }

  printxcorr(XTX, nrows, xindlist)  ; 


  writesnpeigs(snpeigname, xsnplist, ffvecs, numeigs, ncols) ;
  outfiles(snpoutfilename, indoutfilename, genooutfilename, 
   snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode) ;

  printf("##end of smartpca run\n") ;
  return 0 ;
}

void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n, t ;

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

         
   if (parname==NULL) { 
    fprintf(stderr, "no parameters\n") ;
    return ;
   }

   pcheck(parname,'p') ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "poplistname:", &poplistname) ;
   getstring(ph, "snpeigname:", &snpeigname) ;
   getstring(ph, "snpweightoutname:", &snpeigname) ; /* changed 09/18/07 */
   getstring(ph, "output:", &outputname) ;
   getstring(ph, "outputvecs:", &outputname) ;
   getstring(ph, "evecoutname:", &outputname) ; /* changed 11/02/06 */
   getstring(ph, "outputvals:", &outputvname) ;
   getstring(ph, "evaloutname:", &outputvname) ; /* changed 11/02/06 */
   getstring(ph, "badsnpname:", &badsnpname) ;
   getstring(ph, "outliername:", &outliername) ;
   getstring(ph, "outlieroutname:", &outliername) ; /* changed 11/02/06 */
   getstring(ph, "phylipname:", &phylipname) ;
   getstring(ph, "phylipoutname:", &phylipname) ; /* changed 11/02/06 */
   getstring(ph, "weightname:", &weightname) ;
   getstring(ph, "fstdetailsname:", &fstdetailsname) ;
   getint(ph, "numeigs:", &numeigs) ;
   getint(ph, "numoutevec:", &numeigs) ; /* changed 11/02/06 */
   getint(ph, "markerscore:", &markerscore) ; 
   getint(ph, "chisqmode:", &chisqmode) ; 
   getint(ph, "missingmode:", &missingmode) ; 
   getint(ph, "fancynorm:", &fancynorm) ; 
   getint(ph, "usenorm:", &fancynorm) ;  /* changed 11/02/06 */
   getint(ph, "dotpopsmode:", &dotpopsmode) ; 
   getint(ph, "pcorrmode:", &pcorrmode) ;  /* print correlations */
   getint(ph, "pcpopsonly:", &pcpopsonly) ;  /* but only within population */
   getint(ph, "altnormstyle:", &altnormstyle) ;  
   getint(ph, "hashcheck:", &hashcheck) ;
   getint(ph, "popgenmode:", &altnormstyle) ;  
   getint(ph, "noxdata:", &noxdata) ; 
   getint(ph, "inbreed:", &inbreed) ; 
   getint(ph, "easymode:", &easymode) ; 
   getint(ph, "usepopsformissing:", &usepopsformissing) ; 

   t = -1 ;        
   getint(ph, "xdata:", &t) ; if (t>=0) noxdata = 1-t ;
   getint(ph, "nostatslim:", &nostatslim) ; 
   getint(ph, "popsizelimit:", &popsizelimit) ; 
   getint(ph, "minallelecnt:", &minallelecnt) ; 
   getint(ph, "chrom:", &xchrom) ;  
   getint(ph, "maxmissing:", &maxmissing) ; 
   getint(ph, "lopos:", &lopos) ;  
   getint(ph, "hipos:", &hipos) ;  
   getint(ph, "checksizemode:", &checksizemode) ;  
   getint(ph, "pubmean:", &pubmean) ;
   getint(ph, "fstonly:", &fstonly) ;
   getint(ph, "fsthiprecision:", &fsthiprec) ;

   getint(ph, "ldregress:", &ldregress) ;
   getint(ph, "nsnpldregress:", &ldregress) ; /* changed 11/02/06 */
   getdbl(ph, "ldlimit:", &ldlimit) ;  /* in morgans */
   getdbl(ph, "maxdistldregress:", &ldlimit) ;  /* in morgans */ /* changed 11/02/06 */
   getint(ph, "minleneig:", &nostatslim) ;
   getint(ph, "malexhet:", &malexhet) ;
   getint(ph, "nomalexhet:", &malexhet) ; /* changed 11/02/06 */
   getint(ph, "familynames:", &familynames) ;
   getint(ph, "qtmode:", &qtmode) ;

   getint(ph, "numoutliter:", &numoutliter) ;
   getint(ph, "numoutlieriter:", &numoutliter) ; /* changed 11/02/06 */
   getint(ph, "numoutleigs", &numoutleigs) ;
   getint(ph, "numoutlierevec:", &numoutleigs) ; /* changed 11/02/06 */
   getdbl(ph, "outlthresh:", &outlthresh) ;
   getdbl(ph, "outliersigmathresh:", &outlthresh) ; /* changed 11/02/06 */
   getdbl(ph, "blgsize:", &blgsize) ; 

   getstring(ph, "indoutfilename:", &indoutfilename) ;
   getstring(ph, "indivoutname:", &indoutfilename) ; /* changed 11/02/06 */
   getstring(ph, "snpoutfilename:", &snpoutfilename) ;
   getstring(ph, "snpoutname:", &snpoutfilename) ; /* changed 11/02/06 */
   getstring(ph, "genooutfilename:", &genooutfilename) ;
   getstring(ph, "genotypeoutname:", &genooutfilename) ; /* changed 11/02/06 */
   getstring(ph, "outputformat:", &omode) ;
   getstring(ph, "outputmode:", &omode) ;
   getint(ph, "outputgroup:", &ogmode) ;
   getint(ph, "packout:", &packout) ; /* now obsolete 11/02/06 */
   getstring(ph, "twxtabname:", &twxtabname) ;

   getdbl(ph, "r2thresh:", &r2thresh) ;
   getdbl(ph, "r2genlim:", &r2genlim) ;
   getdbl(ph, "r2physlim:", &r2physlim) ;
   getint(ph, "killr2:",  &killr2) ;

   // threading
   getint(ph, "numxtxblocksperside:", &numxtxblocksperside);
   getint(ph, "numsnppartitions:", &numsnppartitions);
   getint(ph, "numthreads:", &numthreads);

   getint(ph, "numchrom:",  &numchrom) ;
   getstring(ph, "xregionname:", &xregionname) ;
   getdbl(ph, "hwfilter:", &nhwfilter) ;
   getstring(ph, "deletesnpoutname:", &deletesnpoutname);

   printf("### THE INPUT PARAMETERS\n");
   printf("##PARAMETER NAME: VALUE\n");
   writepars(ph);

}

int fvadjust(double *cc, int n, double *pmean, double *fancy) 
/* take off mean  force missing to zero */
/* set up fancy norming  */
{
 double p, ynum, ysum, y, ymean, yfancy = 1.0 ;
 int i, nmiss=0 ;



 ynum = ysum = 0.0 ;
 for (i=0; i<n; i++) {  
  y = cc[i] ;
  if (y < 0.0) { 
   ++nmiss ;
   continue ;
  }
  ++ynum ;
  ysum += y ;
 }
 if (ynum==0.0) { 
  fatalx("(fvadjust) snp has no data\n") ;
 }
 ymean = ysum/ynum ;
 for (i=0; i<n; i++) {  
  y = cc[i] ;
  if (y < 0.0) cc[i] = 0.0 ;
  else cc[i] -= ymean ;
 }
 if (pmean != NULL) *pmean = ymean ;
  if (fancynorm) {  
   p = 0.5*ymean ;   // autosomes
   if (altnormstyle == NO) p = (ysum+1.0)/(2.0*ynum+2.0) ;
   y  = p * (1.0-p) ;  
   if (y>0.0) yfancy = 1.0/sqrt(y) ;
  }
 if (fancy != NULL) *fancy = yfancy ;
 return nmiss ;
}

double
dottest(char *sss, double *vec, char **eglist, int numeg, int *xtypes, int len) 
// vec will always have mean 0 
// perhaps should rewrite to put xa1 etc in arrays
{
   double *w1 ; 
   int *xt ;
   int i, k1, k2, k, j, n, x1, x2 ;
   double y1, y2, ylike, yl0, yl1, yl2 ;
   double ychi ;
   double *wmean ;
   int imax, imin, *isort ;
   static int ncall = 0 ;

   char ss1[MAXSTR] ;
   char ss2[MAXSTR] ;
   char sshit[4] ;
   double tail, ans, ftail, ftailx, ansx ; 

   ZALLOC(wmean, numeg, double) ;
   ZALLOC(w1, len + numeg, double) ;
   ZALLOC(isort, numeg, int) ;
   ZALLOC(xt, len, int) ;
   strcpy(ss1, "") ;

   calcmean(wmean, vec, len, xtypes, numeg) ;
   if (pubmean) {  
    copyarr(wmean, w1, numeg) ;
    sortit(w1, isort, numeg) ; 
    printf("%s:means\n", sss) ;
    for (i=0; i<numeg; i++) {  
     k = isort[i] ;
     printf("%20s ", eglist[k]) ;
     printf(" %9.3f\n", wmean[k]) ;
    }
   }

   vlmaxmin(wmean, numeg, &imax, &imin) ;  
    if (chisqmode) {
     ylike = anova1(vec, len, xtypes, numeg) ;
     ans = 2.0*ylike ;
    }
    else {
     ans = ftail = anova(vec, len, xtypes, numeg) ;
    }
    ++ncall ;

    
    if (numeg>2) {
     sprintf(ss2, "%s %s ", sss, "overall") ;
     publishit(ss2, numeg-1, ans) ;
     printf(" %20s minv: %9.3f %20s maxv: %9.3f\n", 
     eglist[imin], wmean[imin], eglist[imax], wmean[imax]) ;
    }


    for (k1 = 0; k1<numeg; ++k1) { 
     for (k2 = k1+1; k2<numeg; ++k2) { 
      n = 0 ;
      x1 = x2 = 0 ; 
      for (i=0; i<len ; i++)   {  
        k = xtypes[i] ;
        if (k == k1) {  
         w1[n] = vec[i] ; 
         xt[n] = 0 ; 
         ++n ;
         ++x1 ;
        }
        if (k == k2) {  
         w1[n] = vec[i] ; 
         xt[n] = 1 ; 
         ++n ;
         ++x2 ;
        }
      }

     if (x1 <= 1) continue ;
     if (x2 <= 1) continue ;

     ylike = anova1(w1, n, xt, 2) ;
     ychi  = 2.0*ylike ;
     chitot[k1*numeg + k2]  += ychi ;
     if (chisqmode) {
      ansx = ychi ;
     }
     else {
      ansx = ftailx = anova(w1, n, xt, 2) ;
     }

      sprintf(ss2,"%s %s %s ", sss, eglist[k1], eglist[k2]) ;
      publishit(ss2, 1, ansx) ;
      
     }
    }
    free(w1) ;
    free(xt) ;
    free(wmean) ;
    free(isort) ;
    return ans ;
}
double anova(double *vec, int len, int *xtypes, int numeg)
// anova 1 but f statistic
{
   int i, k ; 
   double y1, y2, ylike, top, bot, ftail ;  
   double *w0, *w1, *popsize, *wmean ;

   static int ncall2  = 0 ;

   if (numeg >= len) {
    printf("*** warning: bad anova popsizes too small\n") ;
    return 0.0 ;
   }

   ZALLOC(w0, len, double) ;
   ZALLOC(w1, len, double) ;
   ZALLOC(wmean, numeg, double) ;
   ZALLOC(popsize, numeg, double) ;

   y1 = asum(vec, len)/ (double) len ;  // mean
   vsp(w0, vec, -y1, len) ;

    for (i=0; i<len; i++)  { 
     k = xtypes[i] ;
     ++popsize[k] ;
     wmean[k] += w0[i] ;
    }

/* debug */
    if (numeg == 2)  {  
     ++ncall2 ;
     for (i=0; i<len; ++i) {  
      if (ncall2<0) break ;
      k = xtypes[i] ;
//    printf("yy %4d %4d %12.6f %12.6f\n", i, k, vec[i], w0[i]) ;
     }
    }

    vsp(popsize, popsize, 1.0e-12, numeg) ;
    vvd(wmean, wmean, popsize, numeg) ;

    vvt(w1, wmean, wmean, numeg) ;
    top = vdot(w1, popsize, numeg) ;
    
    for (i=0; i<len ; i++)   {  
     k = xtypes[i] ;
     w1[i] = w0[i] - wmean[k] ;
    }
    bot = asum2(w1, len) / (double) (len-numeg) ;
    bot *= (double) (numeg-1) ;
    ftail = rtlf(numeg-1, len-numeg, top/bot) ;

    free(w0) ; 
    free(w1) ; 
    free(popsize) ;
    free(wmean) ;

    return ftail ;

}
double anova1(double *vec, int len, int *xtypes, int numeg)
{
   int i, k ; 
   double y1, y2, ylike ;  
   double *w0, *w1, *popsize, *wmean ;

   ZALLOC(w0, len, double) ;
   ZALLOC(w1, len, double) ;
   ZALLOC(wmean, numeg, double) ;
   ZALLOC(popsize, numeg, double) ;

   y1 = asum(vec, len)/ (double) len ;  // mean
   vsp(w0, vec, -y1, len) ;

    for (i=0; i<len; i++)  { 
     k = xtypes[i] ;
     ++popsize[k] ;
     wmean[k] += w0[i] ;
    }

    vsp(popsize, popsize, 1.0e-12, numeg) ;
    vvd(wmean, wmean, popsize, numeg) ;
    
    for (i=0; i<len ; i++)   {  
     k = xtypes[i] ;
     w1[i] = w0[i] - wmean[k] ;
    }

    y1 = asum2(w0, len) / (double) len ;
    y2 = asum2(w1, len) / (double) len ;
    ylike = 0.5*((double) len)*log(y1/y2) ;

    free(w0) ; 
    free(w1) ; 
    free(popsize) ;
    free(wmean) ;

    return ylike ;

}
void publishit(char *sss, int df, double chi) 
{
      double tail ;
      char sshit[4] ;
      char ss2[MAXSTR] ;
      int i, n ;
      char cblank, cunder ;
      static int ncall = 0 ;

      ++ncall ;
      cblank = ' ' ;
      cunder = '_' ;
      n = strlen(sss) ;

      strcpy(ss2, sss) ;  
      for (i=0; i< n; ++i)  { 
       if (ss2[i] == cblank)  ss2[i] = cunder ;
      }

      if (chisqmode) {
       if (ncall==1) printf("## Anova statistics for population differences along each eigenvector:\n");
       if (ncall==1) printf("%40s %6s %9s %12s\n", "", "dof", "chisq", "p-value") ;
       printf("%40s %6d %9.3f",ss2, ss2, df, chi) ;
       tail = rtlchsq(df, chi) ;  
       printf(" %12.6g", tail) ;
      }
      else {  
       if (ncall==1) printf("## Anova statistics for population differences along each eigenvector:\n");
       if (ncall==1) printf("%40s %12s\n", "",  "p-value") ;
       printf("%40s ", ss2) ;
       tail = chi ;  
       printf(" %12.6g", tail) ;
      }
      strcpy(sshit, "") ;  
      if (tail < pvhit) strcpy(sshit, "***") ;
      if (tail < pvjack) strcpy(sshit, "+++") ;
      printf(" %s", sshit) ;
      printf("\n") ;
}

void
dotpops(double *X, char **eglist, int numeg, int *xtypes, int nrows) 
{
      double *pp, *npp, val, yy ;
      int *popsize ;
      int i, j, k1, k2 ;


     if (fstonly) return ;
     ZALLOC(pp, numeg * numeg, double) ;
     ZALLOC(npp, numeg * numeg, double) ;
     popsize = xpopsize; 

     ivzero(popsize, numeg) ;

     for (i=0; i<nrows; i++) { 
      k1 = xtypes[i] ;
      ++popsize[k1] ;
      for (j=i+1; j<nrows; j++) { 
       k2 = xtypes[j] ;
       if (k1 < 0) fatalx("bug\n") ;
       if (k2 < 0) fatalx("bug\n") ;
       if (k1>=numeg) fatalx("bug\n") ;
       if (k2>=numeg) fatalx("bug\n") ;
       val = X[i*nrows+i] + X[j*nrows+j] - 2.0*X[i*nrows+j] ;
       pp[k1*numeg+k2] += val ;
       pp[k2*numeg+k1] += val ;
       ++npp[k1*numeg+k2]  ;
       ++npp[k2*numeg+k1]  ;
      }
     }
     vsp(npp, npp, 1.0e-8, numeg*numeg) ;
     vvd(pp, pp, npp, numeg*numeg) ;
// and normalize so that mean on diagonal is 1 
     yy = trace(pp, numeg) / (double) numeg ;
     vst(pp, pp, 1.0/yy, numeg*numeg) ;
     printf("\n## Average divergence between populations:");
     if (numeg<=10) {
      printf("\n") ;
      printf("%10s", "") ;
      for (k1=0; k1<numeg; ++k1) {  
       printf(" %10s", eglist[k1]) ;
      }
      printf("  %10s", "popsize") ;
      printf("\n") ;
      for (k2=0; k2<numeg; ++k2) {  
       printf("%10s", eglist[k2]) ;
       for (k1=0; k1<numeg; ++k1) {  
        val = pp[k1*numeg+k2] ;
        printf(" %10.3f", val) ;
       };
       printf("  %10d", popsize[k2]) ;
       printf("\n") ;
      }
     }
     else {   // numeg >= 10 
      printf("\n") ;
      for (k2=0; k2<numeg; ++k2) {  
       for (k1=k2; k1<numeg; ++k1) {  
        printf("dotp: %10s", eglist[k2]) ;
        printf(" %10s", eglist[k1]) ;
        val = pp[k1*numeg+k2] ;
        printf(" %10.3f", val) ;
        printf("    %10d", popsize[k2]) ;
        printf(" %10d", popsize[k1]) ;
        printf("\n") ;
       }
      }
     }
    printf("\n") ;
    printf("\n") ;
    fflush(stdout) ;
     

    free(pp) ;
    free(npp) ;

}
void printxcorr(double *X, int nrows, Indiv **indxx) 
{
   int k1, k2, t ; 
   double y1, y2, yy, rho ;
   Indiv *ind1, *ind2 ;

   if (pcorrmode == NO) return ;
   for (k1=0; k1<nrows; ++k1) {  
    for (k2=k1+1; k2<nrows; ++k2) {  

     ind1 = indxx[k1] ;
     ind2 = indxx[k2] ;

     t = strcmp(ind1 -> egroup, ind2 -> egroup) ;
     if (pcpopsonly && (t != 0))  continue ;
     

     y1 = X[k1*nrows+k1] ;
     y2 = X[k2*nrows+k2] ;
     yy = X[k1*nrows+k2] ;

     rho = yy/sqrt(y1*y2+1.0e-20) ;
     printf("corr: %20s %20s %20s %20s %9.3f\n", 
      ind1 -> ID, ind2 -> ID, ind1 -> egroup, ind2 -> egroup, rho) ;

    }
   }
}
void clearld(double *ldmat, double *ldvv, int rsize, int n, int nclear)    
{
  int i, j, lo, hi ;

  if (nclear>rsize) fatalx("bad nclear\n") ;

  lo = rsize-nclear ;  
  hi = rsize-1 ;

  for (i=lo; i<=hi; ++i)  { 
    vzero(ldvv+i*n, n) ;
    for (j=0; j<=hi; ++j)  { 
      ldmat[j*rsize+i] = ldmat[i*rsize+j] = 0.0 ;
    }
// force matrix non-singular 
    ldmat[i*rsize+i] = 1.0e-8 ;  
  }
}

void
ldreg(double *ldmat, double *ldmat2, double *vv, double *vv2, double *ldvv, 
  double *ldvv2, int rsize, int n) 
/** ldmat2 is inner product matrix for last rsize columns on exit */
{
   int i, j, k1, k2 ; 
   double *rr, *ans, *tt ;
   double y ;

   ZALLOC(rr, rsize, double) ;
   ZALLOC(ans, rsize, double) ;
   ZALLOC(tt, n, double) ;

   if (rsize>1) 
    copyarr(ldvv, ldvv2+n, n*(rsize-1)) ;
   for (i=0; i<rsize-1 ; i++)  { 
    for (j=0; j<rsize-1 ; j++)  { 
     k1 = i*rsize+j  ;
     k2 = (i+1)*rsize+j+1 ;
     ldmat2[k2] = ldmat[k1] ;
    }
   }
   copyarr(vv, ldvv2, n) ;
   i = 0 ;
   for (j=0; j<rsize ; j++)  { 
     y = rr[j] = vdot(vv, ldvv+j*n, n) ;
     y =  vdot(vv, ldvv2+j*n, n) ;
     if (j==0) y += 1.0e-6 ;
     ldmat2[i*rsize+j] = ldmat2[j*rsize+i] = y ;
   }
  solvit(ldmat, rr, rsize, ans) ; /* solve normal equations */
  copyarr(vv, vv2, n) ;
  for (i=0; i<rsize; i++) {  
   vst(tt, ldvv+i*n, -ans[i], n) ;
   vvp(vv2, vv2, tt, n) ;
  }
  free(rr) ;
  free(ans) ;
  free(tt) ;
}

double dofstxx(double *fstans, double *fstsd, SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg, double blgsize, SNP **snpmarkers, Indiv **indm) 

{

   int t1, t2 ;
   int nblocks, xnblocks ;  
   double y, sd ; 
   int *blstart, *blsize ;
   double *xfst ;

   if ( qtmode == YES )  {
     return;
   }

   nblocks = numblocks(snpmarkers, numsnps, blgsize) ;
   printf("number of blocks for moving block jackknife: %d\n", nblocks) ;
   if ( nblocks <= 1 )  {
     return;
   }

   ZALLOC(blstart, nblocks, int) ;
   ZALLOC(blsize, nblocks, int) ;
   ZALLOC(xfst, numeg*numeg, double) ;


   setblocks(blstart, blsize, &xnblocks, xsnplist, ncols, blgsize)  ;
   fixwt(xsnplist, ncols, 1.0) ;

   dofstnumx(xfst, fstans, fstsd, xsnplist, xindex, xtypes,  
    nrows, ncols, numeg, nblocks, indm, YES) ;

   free(blstart) ; 
   free(blsize)  ; 
   free(xfst)  ;

}
void fixwt(SNP **snpm, int nsnp, double val) 
{
  int k ; 
  SNP *cupt ;

  for (k=0; k<nsnp; ++k) { 
   cupt = snpm[k] ;
   cupt -> weight = val ;
  }

}

double oldfstcol(double *estn, double *estd, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2) 
{
   int c1[2], c2[2], *cc ;
   int *rawcol ;
   int k, g, i ; 
   double ya, yb, yaa, ybb, p1, p2, en, ed ;
   double z, zz, h1, h2, yt ;
   static int ncall = 0; 


   ++ncall ;
   ZALLOC(rawcol, nrows, int) ;

   getrawcol(rawcol, cupt, xindex, nrows)  ;

   ivzero(c1, 2) ;
   ivzero(c2, 2) ;

   for (i=0; i< nrows; i++)  { 
    k = xtypes[i] ;
    cc = NULL ;
    if (k==type1) cc = c1 ; 
    if (k==type2) cc = c2 ; 
    if (cc == NULL) continue ;
    g = rawcol[i] ;
    if (g<0) continue ;  
    cc[0] += g ; 
    cc[1] += 2-g ;
   }
   if (ncall < 0) {
    printf("qq2\n") ;
    printimat(c1, 1, 2) ;
    printimat(c2, 1, 2) ;
   }

   ya = c1[0] ;
   yb = c1[1] ;
   yaa = c2[0] ;
   ybb = c2[1] ;
   z = ya + yb ;
   zz = yaa+ybb ;
   if ((z<0.1) || (zz<0.1)) { 
    *estn = 0.0 ;
    *estd = -1.0 ;
    free(rawcol) ;
    return 0.0;
   }

   yt = ya+yb ;
   p1 = ya/yt ;       
   h1 = ya*yb/(yt*(yt-1.0)) ;

   yt = yaa+ybb ;
   p2 = yaa/yt ;         
   h2 = yaa*ybb/(yt*(yt-1.0)) ;

   en = (p1-p2)*(p1-p2) ;  
   en -= h1/z ; 
   en -= h2/zz ; 
 
   ed = en ; 
   ed += h1 ; 
   ed += h2 ;

   *estn = en ; 
   *estd = ed ; 
   

   free(rawcol) ;
   return z + zz ;

}


double fstcol(double *estn, double *estd, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2) 
{
   int c1[2], c2[2], *cc ;
   int *rawcol ;
   int k, g, i ; 
   double ya, yb, yaa, ybb, p1, p2, en, ed ;
   double z, zz, h1, h2, yt ;
   int **ccc ;
   static int ncall = 0 ;


   ++ncall ; 
   ccc = initarray_2Dint(nrows, 2, 0) ;
   ZALLOC(rawcol, nrows, int) ;

   getrawcolx(ccc, cupt, xindex, nrows, indivmarkers)  ;
   getrawcol(rawcol, cupt, xindex, nrows) ;                 

   ivzero(c1, 2) ;
   ivzero(c2, 2) ;

   for (i=0; i< nrows; i++)  { 
    k = xtypes[i] ;
    cc = NULL ;
    if (k==type1) cc = c1 ; 
    if (k==type2) cc = c2 ; 
    if (cc == NULL) continue ;
    g = ccc[i][0] ;
    if (ncall < 1000)  { 
//    printf("zz %d  %d %d\n", rawcol[i], ccc[i][0], ccc[i][1]) ;
    }
    
    if (g<0) continue ;  
    ivvp(cc, cc, ccc[i], 2) ;
   }

   if (ncall < 0)  {
    printf("qqq\n") ;
    printimat(c1, 1, 2) ;
    printimat(c2, 1, 2) ;
   }

   ya = c1[0] ;
   yb = c1[1] ;
   yaa = c2[0] ;
   ybb = c2[1] ;
   z = ya + yb ;
   zz = yaa+ybb ;
   if ((z<1.1) || (zz<1.1)) { 
    *estn = 0.0 ;
    *estd = -1.0 ;
    free(rawcol) ;
    free2Dint(&ccc, nrows) ;
    return 0.0;
   }

   yt = ya+yb ;
   p1 = ya/yt ;       
   h1 = ya*yb/(yt*(yt-1.0)) ;

   yt = yaa+ybb ;
   p2 = yaa/yt ;         
   h2 = yaa*ybb/(yt*(yt-1.0)) ;

   en = (p1-p2)*(p1-p2) ;  
   en -= h1/z ; 
   en -= h2/zz ; 
 
   ed = en ; 
   ed += h1 ; 
   ed += h2 ;

   *estn = en ; 
   *estd = ed ; 
   

   free(rawcol) ;
   free2Dint(&ccc, nrows) ;
   return z + zz ;

}

void
writesnpeigs(char *snpeigname, SNP **xsnplist, double *ffvecs, int numeigs, int ncols) 
{
// this is called at end and ffvecs overwritten
  double *xpt, y, yscal, *snpsc ;
  int i, j, k, kmax, kmin ;
  SNP * cupt  ;
  FILE *fff ;
  
  for (j=0; j<numeigs; ++j) {  
   xpt = ffvecs+j*ncols ;  
   y = asum2(xpt, ncols) ;  
   yscal = (double) ncols / y ;
   yscal = sqrt(yscal) ;
   vst(xpt, xpt, yscal, ncols) ;
  }


  ZALLOC(snpsc, ncols, double) ;
  vclear(snpsc, -99999, ncols) ;
  for (j=0; j<numeigs; ++j) {  
   for (i=0; i<ncols; ++i) {  
    cupt = xsnplist[i] ;
    if (cupt -> ignore) continue ;
    y = ffvecs[j*ncols+i] ;
    snpsc[i] = fabs(y) ; 
   }
   for (k=0; k<=10; ++k) { 
    vlmaxmin(snpsc, ncols, &kmax, &kmin) ;
    cupt = xsnplist[kmax] ;
    printf("eigbestsnp %4d %20s %2d %12d %9.3f\n", j+1, cupt -> ID, cupt -> chrom, nnint(cupt -> physpos), snpsc[kmax]) ;
    snpsc[kmax] = -1.0 ;
   }
  }
  free(snpsc) ;


  if (snpeigname == NULL) return ;
  openit (snpeigname, &fff, "w") ;

  for (i=0; i<ncols; ++i) {  
   cupt = xsnplist[i] ;
   if (cupt -> ignore) continue ;

   fprintf(fff, "%20s", cupt -> ID) ;
   fprintf(fff,  " %2d", cupt -> chrom) ;
   fprintf(fff,  " %12d", nnint(cupt -> physpos)) ;

   for (j=0; j<numeigs; ++j) {  
    fprintf(fff, " %9.3f", ffvecs[j*ncols+i]) ;  
   }
   fprintf(fff, "\n") ;
  }

  fclose(fff) ;

}

/*  load genotype data for this SNP into rawcol  (call this g[])
 *  in fvadjust:
 *    ymean := mean over all non-missing g[i]
 *    xcol[i] -= ymean if g[i] is not missing
 *    xcol[i]  = 0.0   if g[i] is missing
 *    if (fancynorm == NO)
 *      yfancy = 1.0
 *    if (fancynorm == YES and altnormstyle == NO)
 *      yfancy = (ymean/2)*(1-(ymean/2))
 *    if (fancynorm == YES and altnormstyle == YES)
 *      yfancy = ( sum(g[i])+1 ) / ( 2*N + 2 )
 *        for (sum,N) only over non-missing g[i]
 *   back in getcolxz:
 *     on exit:  
 *     xmean[ s ] = ymean * yfancy
 *     xfancy[ s ] = yfancy
 *     *n0 = sum( g[i] )    non-missing g[i] only
 *     *n1 = sum( 2-g[i] )  non-missing g[i] only                
 *     g[i] set to zero where missing data
 *     */


int
getcolxz(double *xcol, SNP *cupt, const int *xindex, int *xtypes, int nrows, int col,  
 double *xmean, double *xfancy, int *n0, int *n1)             
// side effect set xmean xfancy and count variant and reference alleles
// returns missings after fill in
{
 Indiv *indx  ;
 int  j,  n, g, t, k, kmax ;
 double y, pmean, p, yfancy ;
 int *rawcol ;
 int **ccc ;
 int c0, c1, nmiss ;
 double *popnum, *popsum ;

  if (usepopsformissing) { 
   ZALLOC(popnum, MAXPOPS+1, double) ;
   ZALLOC(popsum, MAXPOPS+1, double) ;
  }


  c0 = c1 = 0 ;
  ZALLOC(rawcol, nrows, int) ;
  n = cupt -> ngtypes ;
  if (n<nrows) fatalx("bad snp: %s %d\n", cupt -> ID, n) ;
  getrawcol(rawcol, cupt, (int *)xindex, nrows) ;

  t = 0 ;
  nmiss = 0 ;
  for (j=0; j<nrows; ++j) { 
   g = rawcol[j] ;  
   if (g<0) { 
    ++nmiss ; 
    continue ;  
   }
   c0 += g   ;
   c1 += 2-g ;
   if (usepopsformissing) { 
    k = xtypes[j] ;  
    popsum[k] += (double) g ;
    popnum[k] += 1.0 ;
    kmax = MAX(kmax, k) ;
    ++t ;
   }
  }
  floatit(xcol, rawcol, nrows) ;
  if ((usepopsformissing) && (nmiss > 0)) {
   pmean = asum(popsum, kmax+1)/asum(popnum, kmax+1) ;
   nmiss = 0 ;
   for (j=0; j<nrows; ++j) { 
    g = rawcol[j] ;  
    if (g>=0) continue ;  
    k = xtypes[j] ; 
    if (popnum[k] > 0.5)  { 
      y = popsum[k]/popnum[k] ; 
      xcol[j] = y ;
      continue ;
    }
    ++nmiss ;
   }
  }

  fvadjust(xcol, nrows, &pmean, &yfancy) ;

  vst(xcol, xcol, yfancy, nrows) ;
  if (xmean != NULL) {
   xmean[col] = pmean*yfancy ; 
   xfancy[col] = yfancy ;
  }
  free(rawcol) ;
  if (n0 != NULL) {
   *n0 = c0 ; 
   *n1 = c1 ;
  }
  if (usepopsformissing) { 
   free(popnum) ;
   free(popsum) ;
  }
  return nmiss ;
}


void
domult(double *cc, const block *block)
{
  int i ;
  double ycheck ;
  ycheck = asum(cc, block->nrows) ;
  if (fabs(ycheck)>.00001) fatalx("bad ycheck\n") ;

  // Optimization; if processing the entire block, only compute the upper triangle of the matrix.
  // Otherwise, we need to assume that the entire matrix is being calculated.
  if( block->XTX_col_start == 0 && block->XTX_col_stop == block->nrows-1 &&
      block->XTX_row_start == 0 && block->XTX_row_stop == block->nrows-1 )
    addoutersym(block->XTX_block, cc, block->nrows);
  else
    addouterblock(block->XTX_block, cc, block->nrows,
                  block->XTX_row_start, block->XTX_row_stop, block->XTX_col_start, block->XTX_col_stop);
}


void
getcolxf(double *xcol, SNP *cupt, int *xindex, int nrows, int col,
 double *xmean, double *xfancy)
// side effect set xmean xfancy
{
 Indiv *indx  ;
 int  j,  n, g, t ;
 double y, pmean, p, yfancy ;
 int *rawcol ;

  if (xmean != NULL) {
   xmean[col] = xfancy[col] = 0.0 ;
  }

  if (cupt -> ignore) {
   vzero(xcol, nrows) ;
   return ;
  }

  ZALLOC(rawcol, nrows, int) ;
  n = cupt -> ngtypes ;
  if (n<nrows) fatalx("bad snp: %s %d\n", cupt -> ID, n) ;
  getrawcol(rawcol, cupt, xindex, nrows) ;

  floatit(xcol, rawcol, nrows) ;
  fvadjust(xcol, nrows, &pmean, &yfancy) ;
  vst(xcol, xcol, yfancy, nrows) ;
  if (xmean != NULL) {
   xmean[col] = pmean*yfancy ;
   xfancy[col] = yfancy ;
  }
  free(rawcol) ;
}

double doinbxx(double *inbans, double *inbsd, SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg, double blgsize, SNP **snpmarkers, Indiv **indm) 
{

   int t1, t2 ;
   int nblocks, xnblocks ;  
   double y, sd ; 
   int *blstart, *blsize ;
   double *xinb ;

  nblocks = numblocks(snpmarkers, numsnps, blgsize) ;

  ZALLOC(blstart, nblocks, int) ;
  ZALLOC(blsize, nblocks, int) ;
  ZALLOC(xinb, numeg, double) ;


  setblocks(blstart, blsize, &xnblocks, xsnplist, ncols, blgsize)  ;
  fixwt(xsnplist, ncols, 1.0) ;

  doinbreed(xinb, inbans, inbsd, xsnplist, xindex, xtypes,  
   nrows, ncols, numeg, nblocks, indm) ;

  free(blstart) ; 
  free(blsize)  ; 
  free(xinb)  ;

}


void calcpopmean(double *wmean, char **elist, double *vec, 
 char **eglist, int numeg, int *xtypes, int len) 
// extracted from dotttest ;
{
   double *w0, *w1 ; 
   int *isort ;
   int i, k ;

   ZALLOC(w0, len, double) ;
   ZALLOC(w1, len, double) ;
   ZALLOC(isort, len, int) ;

   
    calcmean(w0, vec, len, xtypes, numeg) ;

    copyarr(w0, w1, numeg) ;
    sortit(w1, isort, numeg) ; 

    for (i=0; i<numeg; i++) {  
     k = isort[i] ;
     elist[i] = eglist[k] ;
     wmean[i] = w0[k] ;
    }



   free(w0) ; 
   free(w1) ; 
   free(isort) ;


}



void compute_XTX( double *XTX, SNP **xsnplist, const int *xindex, double *xmean, double *xfancy, const int weightmode,
  const int nrows, const int ncols, int *xtypes )
{
  /*
   * Make a best guess as to how to decompose the problem into thread-sized chunks.
   * Number of threads is determined based on how many processors in the system.
   * Number of XTX blocks per side increases as the order of magnitude of # of matrix elements  increases, starting when the matrix has 1M elements.
   * Number of SNP partitions is same as the number of threads.
   */
  if(numthreads == THREADING_DEFAULT) numthreads = sysconf(_SC_NPROCESSORS_ONLN);
  if(numxtxblocksperside == THREADING_DEFAULT) numxtxblocksperside = nrows > 1000 ? (int)log10(nrows*nrows)-4 : 1;
  if(numsnppartitions == THREADING_DEFAULT) numsnppartitions = numthreads;

  if(numxtxblocksperside <= 0) fatalx("Number of xtx blocks per side must be greater than 0.");
  if(numsnppartitions <= 0) fatalx("Number of SNP partitions must be greater than 0.");
  if(numthreads <= 0) fatalx("Number of threads must be greater than 0.");

#ifdef HAVE_PTHREAD
  printf("Multicore implementation; using %d thread%s (%d matrix block%s per side, %d snp partition%s)\n",
         numthreads, numthreads > 1 ? "s" : "",
         numxtxblocksperside, numxtxblocksperside > 1 ? "s" : "",
         numsnppartitions, numsnppartitions > 1 ? "s" : "");
#else
  printf("Single core implementation; using %d snp partition%s\n",
         numxtxblocksperside, numxtxblocksperside > 1 ? "s" : "");
#endif

  const int xtx_size = nrows;
  block* blocks = NULL;
  int block_row, block_col, block_width, block_height;
  int i, snp_partition, thread = 0;

  ZALLOC(blocks,numthreads,block);
  vzero(XTX, nrows*nrows);

  work_queue* queue = NULL;
  block *block = NULL;

  for( block_col = 0; block_col < numxtxblocksperside; block_col++ )
  {
    for( block_row = 0; block_row < numxtxblocksperside; block_row++ )
    {
      for( snp_partition = 0; snp_partition < numsnppartitions; snp_partition++ )
      {
        if( thread == 0 )
          create_work_queue(&queue);

        // Initialize the block
        block = &blocks[thread];
        block->xsnplist = xsnplist;
        block->xindex = xindex;
        block->xtypes = xtypes;
        block->xmean = xmean;
        block->xfancy = xfancy;
        block->weightmode = weightmode;
        block->nrows = nrows;
        block->ncols = ncols;

        // Determine the column and row bounds for this block
        partition( xtx_size, block_col, numxtxblocksperside, &block->XTX_col_start, &block->XTX_col_stop );
        partition( xtx_size, block_row, numxtxblocksperside, &block->XTX_row_start, &block->XTX_row_stop );

        // Only compute blocks which contribute to the UTM.  If the upper right corner of this
        // block is below the diagonal of the matrix, skip computation of the block.
        if( block->XTX_col_stop < block->XTX_row_start )
          continue;

        // Determine which SNPs the system will take into account in determining this block.
        partition( ncols, snp_partition, numsnppartitions, &block->SNP_col_start, &block->SNP_col_stop );

        // Allocate temporary space for the block
        block_width = block->XTX_col_stop - block->XTX_col_start + 1;
        block_height = block->XTX_row_stop - block->XTX_row_start + 1;
        ZALLOC(block->XTX_block, block_height*block_width,double);
        vzero(block->XTX_block, block_height*block_width);

        // create a task to process this block.
        work_task task;
        task.start_routine = compute_block;
        task.argument = block;

        queue_task(queue,&task);

        // If the max number of threads have been reached, collect the results.
        thread = ++thread % numthreads;
        if( thread == 0 )
        {
          collect_output( XTX, queue, blocks, numthreads );
          destroy_work_queue( &queue );
        }
      }
    }
  }

  // Collect remaining output.
  if( thread > 0 )
  {
    collect_output( XTX, queue, blocks, thread );
    destroy_work_queue( &queue );
  }

  symit(XTX, nrows) ;

  free(blocks);
}



void* compute_block( void* block_holder ) {
  int i, j, k, n0, n1 ;

  SNP *cupt, *cupt2 ;

  double *cc, weight;
  int t;
  int chrom, numclear ;
  double gdis ;
  double *ldmat = NULL, *ldmat2 = NULL, *ldvv = NULL, *ldvv2 = NULL, *vv2 = NULL  ;
  int xblock;

  block *block = block_holder;

  ZALLOC(cc, block->nrows, double) ;

  if (ldregress>0)
  {
    ZALLOC(ldvv,  ldregress*numindivs, double) ;
    ZALLOC(ldvv2,  ldregress*numindivs, double) ;
    ZALLOC(vv2,  numindivs, double) ;
    ZALLOC(ldmat,  ldregress*ldregress, double) ;
    ZALLOC(ldmat2,  ldregress*ldregress, double) ;
    setidmat(ldmat, ldregress) ;
    vst(ldmat, ldmat, 1.0e-6, ldregress*ldregress) ;
  }

  for (i = block->SNP_col_start; i <= block->SNP_col_stop; i++)
    {
      cupt = block->xsnplist[i] ;
      chrom = cupt -> chrom ;
      getcolxz(cc, cupt, block->xindex, block->xtypes, block->nrows, i, block->xmean, block->xfancy, &n0, &n1) ;
      t = MIN(n0, n1) ;

      if (t < minallelecnt)  {
        cupt -> ignore = YES ;
        logdeletedsnp(cupt->ID,"minallelecnt",deletesnpoutname);
        vzero(cc, block->nrows) ;
      }

      if (block->weightmode)
        {
          weight = block->xsnplist[i] -> weight ;
          vst(cc, cc, block->xsnplist[i] -> weight, block->nrows) ;
        }
      if (ldregress>0)
        {
          numclear = 0 ;
          for (k=1; k<= ldregress; ++k)
            {
              j = i-k ;
              if (j<0)
                {
                  numclear = ldregress-k+1 ;
                  break ;
                }
              cupt2 = block->xsnplist[j] ;
              if (cupt2 -> chrom != chrom) gdis = ldlimit + 1.0 ;
              else gdis = cupt -> genpos - cupt2 -> genpos ;
              if (gdis>=ldlimit)
                {
                  numclear = ldregress-k+1 ;
                  break ;
                }
            }
          if (numclear>0) clearld(ldmat, ldvv, ldregress, block->nrows, numclear) ;
          ldreg(ldmat, ldmat2, cc, vv2, ldvv, ldvv2, ldregress, block->nrows) ;
          copyarr(ldmat2, ldmat, ldregress*ldregress) ;
          copyarr(vv2, cc, block->nrows) ;
          copyarr(ldvv2, ldvv, ldregress*block->nrows) ;
        }

      domult(cc, block);
    }

  free(cc);

  if(ldvv != NULL) free(ldvv);
  if(ldvv2 != NULL) free(ldvv2);
  if(vv2 != NULL) free(vv2);
  if(ldmat != NULL) free(ldmat);
  if(ldmat2 != NULL) free(ldmat2);

  return NULL;
}

void collect_output( double *XTX, const work_queue* queue, const block* blocks, const int num_threads )
{
  const block *block = NULL;
  int block_width,block_height,thread = 0;

  wait_for_queue_to_complete(queue);

  for( thread = 0; thread < num_threads; thread++ )
  {
      // Determine sizing for the block.
      block = &blocks[thread];
      block_width = block->XTX_col_stop - block->XTX_col_start + 1;
      block_height = block->XTX_row_stop - block->XTX_row_start + 1;

      // Add the block to the big XTX matrix.
      addblocktomatrix(XTX,block->XTX_col_start,block->XTX_row_start,block->nrows,block->XTX_block,block_width,block_height);

      // Free temporary memory.
      free(block->XTX_block);
  }
}

// Calculate the start and stop point of the nth partition of size, given num_partition
// partitions.  Any excess from rounding error will be added to the num_partitions-1th bucket.
void partition( const int size, const int n, const int num_partitions, int *start, int* stop )
{
  *start = n * (size/num_partitions);
  *stop = (n != (num_partitions-1)) ? (n+1)*(size/num_partitions)-1 : size-1;
}






