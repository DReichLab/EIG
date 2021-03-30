#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdint.h>
#include <inttypes.h>

#include <nicklib.h>
#include <getpars.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"
#include "mcmcpars.h"
#include "eigsubs.h"
#include "gval.h"
#include "egsubs.h"
#include "qpsubs.h"
#include "smartsubs.h"
#include "exclude.h"
#include "globals.h"

/** 
 Most of this code written by Nick Patterson 
 (Broad institute, Harvard Medical, Harvard Evolutionary Biology)
 Improvements and elimination of FORTRAN code by Chris Chang (BGI) 

 Code added to support grm output + improved ld rregression by Alexander Gusev 
*/

#define WVERSION   "18140"

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
s. dev -> std. err. in legend 

bad bug fixed (outfiles changed indivmarkers)  ...  

fstdetailsname added
fsthiprecision added
fstsnpout switch added 

bug fixed  (getrawcolx)

bad bug fix.  xtypes not allocated correctly

version compatible with Mac
XTX.dbg commented out

outliermode added

regmode added
maxpops parametric.  Use easymode if large

id2pops added

Threading added Chris Chang) 
fastmode (Kevin Galinski) 
bugfix to ldregx (Angela Yu)
shrinkmode (see doshrinkp) 

// was smshrink.c
ZALLOC improves error message.  
TW code cleaned up.

nochrom option added (for Jackknife)  
hiprec option added (for eigenvector output) 

fastshrink added (Edgar Dobriban) 

fstz code added (lower triangle is Z scores) 

O2 version.  some extra declarations;  cputime and calcmem called
new estimation of effective number of snps
topright added -- flips top 2 eigenvectors if wanted
dotpopsmode NO default
fstnum added
megaoutname added (mega output)
*/

#include <pthread.h>
#define THREAD_RET_TYPE void*
#define THREAD_RETURN return NULL
#define MAX_THREADS 127
#define MAX_THREADS_P1 128

#define MAXFL  50
#define MAXSTR  512
#define MAXPOPS 1000

char *parname = NULL;
char *twxtabname = NULL;
char *trashdir = "/var/tmp";
int qtmode = NO;
Indiv **indivmarkers;
SNP **snpmarkers;

int shrinkp_done = NO ; 

int numsnps, numindivs;
int numeigs = 10;               /// default
int markerscore = NO;
int maxpops = -1;
long seed = 0;
int chisqmode = NO;             // approx p-value better to use F-stat
int missingmode = NO;
int shrinkmode = NO;
int readparsonly = NO ; 
int newshrink = NO;
int fastshrink = NO ;
int autoshrink = NO  ; 
double autothresh = .01 ;  // TW test threshold for significant eigenvalue
int mpshrink = NO ; 
int dotpopsmode = NO;
int noxdata = YES;              /* default as pop structure dubious if Males and females */
int fstonly = NO;
int fstjack = YES ;
int pcorrmode = NO;
int pcpopsonly = YES;
int nostatslim = 10;
int znval = -1;
int popsizelimit = -1;
int altnormstyle = YES;         // affects subtle details in normalization formula
int minallelecnt = 1;
int maxmissing = 9999999;
int lopos = -999999999, hipos = 999999999;      // use with xchrom
int fstverbose = NO ; 

int packout = -1;
extern enum outputmodetype outputmode;
extern int checksizemode;
extern int packmode;
extern int numchrom;
extern int fancynorm;
extern int verbose;
int ogmode = NO;
int fsthiprec = NO;
int fstz = NO;
int hiprec = NO;
int inbreed = NO;               // for fst
int easymode = NO;
int fastmode = NO;
int fastdim = -1;
int fastiter = -1;
int regmode = YES;
int doproject = YES  ;  // No => poplistname pops only projected

int numoutliter = 5, numoutleigs = 10, outliermode = 0;
double outlthresh = 6.0;
OUTLINFO **outinfo;
char *outinfoname = NULL;
char *fstdetailsname = NULL;
char *elloutname = NULL ; 
char *megaoutname = NULL ; 
double ellconf = -1 ; 
int fstsnpout = NO ; 

double *XTX = NULL ; 
int nrxtx = -1 ; 
int *flip = NULL ; 


double plo = .001;
double phi = .999;
double pvhit = .001;
double pvjack = 1.0e-6;
double *chitot;
int *xpopsize;

char *genotypename = NULL;
char *snpname = NULL;
char *indivname = NULL;
char *badsnpname = NULL;
char *deletesnpoutname = NULL;
char *poplistname = NULL;
char *xregionname = NULL;       /* physical positions of SNPs to exclude */
char *outliername = NULL;
char *phylipname = NULL;
char *snpeigname = NULL;

char *indoutfilename = NULL;
char *snpoutfilename = NULL;
char *genooutfilename = NULL;
char *omode = "packedancestrymap";
char *grmoutname = NULL;
int grmbinary = NO;
double blgsize = 0.05;          // block size in Morgans */
char *id2pops = NULL;
char *elllistname = NULL ; 
int nellindivs = 0 ; 
char **elllist ; 
int *ellindex ; 
Indiv **ellindivs ;  

int nblocks, xnblocks;
int *blstart, *blsize;

int printcover = NO ;

char *topright = NULL ; 
int  toprightindex = -1 ;


double r2thresh = -1.0;
double r2genlim = 0.01;         // Morgans 
double r2physlim = 5.0e6;
int killr2 = NO;
int pubmean = YES;              // change default

double nhwfilter = -1.0;

double *xmean, *xfancy;
double *edgarw = NULL  ;  

int rounakmode = NO ; 

int thread_ct_config = 0;

int randomfillin = NO;
int usepopsformissing = NO;     // if YES popmean is used for missing.  Overall mean if all missing for pop

int xchrom = -1;
int zchrom = - 1 ;
// list of outliers

int ldregress = 0;
double ldlimit = 9999.0;        /* default is infinity */
double ldr2lo = 0.01;
double ldr2hi = 0.95;
int ldposlimit = 1000 * 1000 * 1000;
int ldregx (double *gsource, double *gtarget, double *res, int rsize,
            int n, double r2lo, double r2hi);
void bumpldvv (double *gsource, double *newsource, int *pnumld, int maxld,
               int n, int *ldsnpbuff, int newsnpnum);


char *outputname = NULL;
char *outputvname = NULL;
char *weightname = NULL;
FILE *ofile, *ovfile;

double twestxx (double *lam, int m, double *pzn, double *pzvar);
double twnorm (double lam, double m, double n);
double rhoinv (double x, double gam);

void readcommands (int argc, char **argv);
int loadindx (Indiv ** xindlist, int *xindex, Indiv ** indivmarkers,
              int numindivs);
void loadxdataind (double *xrow, SNP ** snplist, int ind, int ncols);
void fixxrow (double *xrow, double *xmean, double *xfancy, int len);
void dofancy (double *cc, int n, double *fancy);
int fvadjust (double *rr, int n, double *pmean, double *fancy);
void getcol (double *cc, double *xdata, int col, int nrows, int ncols);
void getcolxf (double *xcol, SNP * cupt, int *xindex,
               int nrows, int col, double *xmean, double *xfancy);
int getcolxz (double *xcol, SNP * cupt, int *xindex, int *xtypes,
              int nrows, int col, double *xmean, double *xfancy, int *n0,
              int *n1);
int getcolxz_binary1 (int *rawcol, double *xcol, SNP * cupt, int *xindex,
                      int nrows, int col, double *xmean, double *xfancy,
                      int *n0, int *n1);
void getcolxz_binary2 (int *rawcol, uintptr_t * binary_cols,
                       uintptr_t * binary_mmask, uint32_t xblock,
                       uint32_t nrows);

void doinbxx (double *inbans, double *inbsd, SNP ** xsnplist, int *xindex,
              int *xtypes, int nrows, int ncols, int numeg, double blgsize,
              SNP ** snpmarkers, Indiv ** indm);

void putcol (double *cc, double *xdata, int col, int nrows, int ncols);
void calcpopmean (double *wmean, char **elist, double *vec,
                  char **eglist, int numeg, int *xtypes, int len);
double dottest (char *sss, double *vec, char **eglist, int numeg, int *xtypes,
                int len);
double yll (double x1, double x2, double xlen);
void calcmean (double *wmean, double *vec, int len, int *xtypes, int numeg);
double anova1 (double *vec, int len, int *xtypes, int numeg);
double anova (double *vec, int len, int *xtypes, int numeg);
void publishit (char *sss, int df, double chi);

void setmiss (SNP ** snpm, int numsnps);
void setfvecs (double *fvecs, double *evecs, int nrows, int numeigs);
void dotpops (double *X, char **eglist, int numeg, int *xtypes, int nrows);
void printxcorr (double *X, int nrows, Indiv ** indxx);

void fixrho (double *a, int n);
void printdiag (double *a, int n);
void countsn(int *blcnt, int nblocks, SNP **snplist, int nsnps, int ind)  ;

int
ridoutlier (double *evecs, int n, int neigs,
            double thresh, int *badlist, OUTLINFO ** outinfo);

void addoutersym (double *X, double *v, int n);
void symit (double *X, int n);

double fstcol (double *estn, double *estd, SNP * cupt,
               int *xindex, int *xtypes, int nrows, int type1, int type2);

double oldfstcol (double *estn, double *estd, SNP * cupt,
                  int *xindex, int *xtypes, int nrows, int type1, int type2);

void jackrat (double *xmean, double *xsd, double *top, double *bot, int len);
void lsqproj(int blocknum, SNP **xsnplist, int ncols, Indiv **xindlist, int nind, 
  double *fxscal, double *ffvecs, double * acoeffs, double *bcoeffs, int *xtypes, int numeg) ;
void seteigscale(double *eigscale, double *acoeffs, double *bcoeffs, int *xindex, int nrows, int numeigs)  ;

void domult_increment_lookup (pthread_t * threads, uint32_t thread_ct,
                              double *XTX_lower_tri, double *tblock,
                              uintptr_t * binary_cols,
                              uintptr_t * binary_mmask, uint32_t block_size,
                              uint32_t indiv_ct,
                              double *partial_sum_lookup_buf);
void domult_increment_normal (pthread_t * threads, uint32_t thread_ct,
                              double *XTX_lower_tri, double *tblock,
                              int marker_ct, uint32_t indiv_ct);
void writesnpeigs (char *snpeigname, SNP ** xsnplist, double *ffvecs,
                   int numeigs, int ncols);
void dofstxx (double *fstmean, double *fstans, double *fstsd, int *fstnum, SNP ** xsnplist, int *xindex,
              int *xtypes, int nrows, int ncols, int numeg, double blgsize,
              SNP ** snpmarkers, Indiv ** indm);
void fixwt (SNP ** snpm, int nsnp, double val);
void sqz (double *azq, double *acoeffs, int numeigs, int nrows, int *xindex);
void dumpgrm (double *XTX, int *xindex, int nrows, int numsnps,
              Indiv ** indivmarkers, int numindivs, char *grmoutname);

void printevecs (SNP ** snpmarkers, Indiv ** indivmarkers, Indiv ** xindlist,
                 int numindivs, int ncols, int nrows,
                 int numeigs, double *eigenvecs, double *eigenvals,
                 FILE * ofile);

void doshrinkp (double *mmat, int m, int n, int *xindex, SNP ** xsnplist, double *xcoeffs) ;
void estedgar(double *edgarw, double *lambdav, int lentop, int lenspec, double gamm, double yjfac) ;
int setnstw(double *lambda, int len, int nostatslim) ;
void kjg_fpca (size_t K, size_t L, size_t I, double *eval, double *evec) ;

int mpestimate(double *eigs, int neigs, double *peffect, double *psigma) ; 
void printmega(char *outname, char **eglist, int numeg, double *fstsc) ;

uint32_t
triangle_divide (int64_t cur_prod, int32_t modif)
{
  // return smallest integer vv for which (vv * (vv + modif)) is no smaller
  // than cur_prod, and neither term in the product is negative.  (Note the
  // lack of a divide by two; cur_prod should also be double its "true" value
  // as a result.)
  int64_t vv;
  if (cur_prod == 0) {
    if (modif < 0) {
      return -modif;
    }
    else {
      return 0;
    }
  }
  vv = (int64_t) sqrt ((double) cur_prod);
  while ((vv - 1) * (vv + modif - 1) >= cur_prod) {
    vv--;
  }
  while (vv * (vv + modif) < cur_prod) {
    vv++;
  }
  return vv;
}

void
parallel_bounds (uint32_t ct, int32_t start, uint32_t parallel_idx,
                 uint32_t parallel_tot, int32_t * bound_start_ptr,
                 int32_t * bound_end_ptr)
{
  int32_t modif = 1 - start * 2;
  int64_t ct_tot = ((int64_t) ct) * (ct + modif);
  *bound_start_ptr =
    triangle_divide ((ct_tot * parallel_idx) / parallel_tot, modif);
  *bound_end_ptr =
    triangle_divide ((ct_tot * (parallel_idx + 1)) / parallel_tot, modif);
}

// set align to 1 for no alignment
void
triangle_fill (uint32_t * target_arr, uint32_t ct, uint32_t pieces,
               uint32_t parallel_idx, uint32_t parallel_tot, uint32_t start,
               uint32_t align)
{
  int32_t modif = 1 - start * 2;
  uint32_t cur_piece = 1;
  int64_t ct_tr;
  int64_t cur_prod;
  int32_t lbound;
  int32_t ubound;
  uint32_t uii;
  uint32_t align_m1;
  parallel_bounds (ct, start, parallel_idx, parallel_tot, &lbound, &ubound);
  // x(x+1)/2 is divisible by y iff (x % (2y)) is 0 or (2y - 1).
  align *= 2;
  align_m1 = align - 1;
  target_arr[0] = lbound;
  target_arr[pieces] = ubound;
  cur_prod = ((int64_t) lbound) * (lbound + modif);
  ct_tr = (((int64_t) ubound) * (ubound + modif) - cur_prod) / pieces;
  while (cur_piece < pieces) {
    cur_prod += ct_tr;
    lbound = triangle_divide (cur_prod, modif);
    uii = (lbound - ((int32_t) start)) & align_m1;
    if ((uii) && (uii != align_m1)) {
      lbound = start + ((lbound - ((int32_t) start)) | align_m1);
    }
    // lack of this check caused a nasty bug earlier
    if (((uint32_t) lbound) > ct) {
      lbound = ct;
    }
    target_arr[cur_piece++] = lbound;
  }
}

void
symit2 (double *XTX, uintptr_t nrows)
{
  // unpacks LOWER-triangle-only symmetric matrix representation into regular
  // square matrix.
  uintptr_t row_idx;
  uintptr_t col_idx;
  double *read_col;
  double *write_ptr;
  if (nrows < 3) {
    if (nrows < 2) {
      return;
    }
    // special case, need to avoid overlapping memcpy
    XTX[3] = XTX[2];
    XTX[2] = XTX[1];
    return;
  }
  for (row_idx = nrows - 1; row_idx; row_idx--) {
    memcpy (&(XTX[row_idx * nrows]), &(XTX[(row_idx * (row_idx + 1)) / 2]),
            (row_idx + 1) * sizeof (double));
  }
  for (row_idx = 0; row_idx < nrows; row_idx++) {
    read_col = &(XTX[row_idx]);
    write_ptr = &(XTX[row_idx * nrows + row_idx + 1]);
    for (col_idx = row_idx + 1; col_idx < nrows; col_idx++) {
      *write_ptr++ = read_col[col_idx * nrows];
    }
  }
}

void
copy_transposed (double *orig_matrix, uintptr_t orig_row_ct,
                 uintptr_t orig_col_ct, double *transposed_matrix)
{
  uintptr_t new_row_idx;
  uintptr_t new_col_idx;
  double *orig_col_ptr;
  for (new_row_idx = 0; new_row_idx < orig_col_ct; new_row_idx++) {
    orig_col_ptr = &(orig_matrix[new_row_idx]);
    for (new_col_idx = 0; new_col_idx < orig_row_ct; new_col_idx++) {
      *transposed_matrix++ = orig_col_ptr[new_col_idx * orig_col_ct];
    }
  }
}

// make these file scope so multithreading works
static double *g_XTX_lower_tri;
static double *g_tblock;
static uint32_t g_block_size;
static uintptr_t g_indiv_ct;
static uint32_t g_thread_start[MAX_THREADS_P1];

static double *g_weights;
static uintptr_t *g_binary_cols;
static uintptr_t *g_binary_mmask;

int
main (int argc, char **argv)
{

  char sss[MAXSTR];
  char **eglist;
  int numeg;
  int i, j, jj, k, k1, k2, pos;
  int a, b, bl ;
  int *vv;
  SNP *cupt;
  Indiv *indx;
  double y1 = 0, y2, y2l, y, y3, ymem, yjfac;

  int n0, n1, nkill;

  int nindiv = 0, nindignore ;
  double ychi, tail, tw;
  int nignore, numrisks = 1;
  double *xrow, *xpt;
  double **xrowb ; 
  SNP **xsnplist;
  Indiv **xindlist;
  int *xindex, *xtypes = NULL;
  int nrows, ncols, m, nused;
  double  *cc, *evecs, *ww, *evals;
  double *partial_sum_lookup_buf = NULL;
  double *lambda, *esize;
  int **blcnt ; 
  double zn, zvar, zp, zfn, zfp;
  double *fvecs, *fxvecs, *fxscal;
  double *ffvecs;
  int weightmode = NO;
  double ynrows;
  int t, tt;
  double *ldvv = NULL, ynumsnps = 0;    // for grm
  int *ldsnpbuff = NULL;
  int lastldchrom, numld;
  double *fstans, *fstsd, *fstzz;
  double *inbans, *inbsd;
  int *fstnum ; 

  int chrom;
  int outliter, numoutiter, *badlist, nbad;
  FILE *outlfile, *phylipfile;
  double *eigkurt, *eigindkurt;
  double *wmean, *emean ;
  char **elist;
  double *mmat;

  double *trow = NULL, *rhs = NULL, *emat = NULL, *regans = NULL;
  double *eigscale = NULL ;
  double *eigscmat ;  // scaling matrix for lsqproj

  int kk;
  double *acoeffs, *bcoeffs, *xcoeffs = NULL ;
  double **aellc , **bellc ;  
  double *aco, *bco ;
  double *awk, *bwk ;
  int rngmode = NO;


  int xblock;
  int blocksize = 1024;
  double *tblock = NULL;
  int *binary_rawcol = NULL;
  uintptr_t *binary_cols = NULL;
  uintptr_t *binary_mmask = NULL;

  OUTLINFO *outpt;

  pthread_t threads[MAX_THREADS];
  uint32_t thread_ct;

  int numshrink ; // set by autoshrink 
  double fstzsc, *fstmean ;

  double zeffect, zsigma ; // for MP estimation

  double *jwt, *vest, *var ; 
  double **jmean, **vvest, **vvar ; 
  double zval, za1, za2; 
  double zlambda[2], zvecs[4], zangle ; 
  double zzmean[2], zzvar[4], zzest[2] ; 

  FILE *ellfile ;
  int fnum ;  

  readcommands (argc, argv);

  cputime(0) ;
  calcmem(0) ;


  printf ("## smartpca version: %s\n", WVERSION);

  if (readparsonly) { 
   printf("terminating...\n") ;
   return 0 ;
  }
  if (regmode == NO) printf("lsqproject off:  deprecated\n") ;
  if (usepopsformissing && shrinkmode) { 
    fatalx("usepopsformissing and shrinkmode together not supported\n") ;
  }

  packmode = YES;
  setomode (&outputmode, omode);
  setfstsnpout(fstsnpout) ; 

  if (parname == NULL)
    return 0;
  if (xchrom == (numchrom + 1))
    noxdata = NO;

  if (fastmode) {
    if (fastiter < 0)
      fastiter = numeigs;
    if (fastdim < 0)
      fastdim = 2 * numeigs;
    rngmode = YES;
  }

 if (newshrink) shrinkmode = YES ;

 if (mpshrink) { 
  fastshrink = YES ;
  printf("mpshrink set!\n") ;
 }

  if (rounakmode) { 
   fastshrink = YES ; 
  }

  if (autoshrink) { 
   fastshrink = YES ;
  }

  if (fastshrink) {
   shrinkmode = NO ;
  }

  if (popsizelimit > 0)
    rngmode = YES;

  if (rngmode) {
    if (seed == 0)
      seed = seednum ();
    printf ("seed: %ld\n", seed);
    SRAND (seed);
  }


/**
  if (usepopsformissing) {
    printf ("usepopsformissing => easymode\n");
    easymode = YES;
  }
*/

  if (hiprec) printf("high precision evec output\n") ;

  if (deletesnpoutname != NULL) {       /* remove because snplog opens in append mode */
    char buff[256];
    sprintf (buff, "rm -f %s", deletesnpoutname);
    system (buff);
  }

  if (fstonly) {
    printf ("fstonly\n");
    numeigs = 0;
    numoutliter = 0;
    numoutiter = 0;
    outputname = NULL;
    snpeigname = NULL;
    if (fastmode) printf("*** no fastmode with fstonly!\n") ; 
    fastmode = NO ;
  }

  if (fancynorm)
    printf ("norm used\n\n");
  else
    printf ("no norm used\n\n");
  if (regmode)
    printf ("lsqproject used\n");
  if (shrinkmode) printf("shrinkmode used\n") ; 
  if (rounakmode) printf("rounakmode used\n") ; 

  nostatslim = MAX (nostatslim, 3);

  outlfile = ofile = stdout;

//if (outputname != NULL)  openit(outputname, &ofile, "w") ;
  if (shrinkmode)
    openit ("/dev/null", &ofile, "w");
  if ((outputname != NULL) && (!shrinkmode))
    openit (outputname, &ofile, "w");
  if (outliername != NULL)
    openit (outliername, &outlfile, "w");
  if (fstdetailsname != NULL)
    openit (fstdetailsname, &fstdetails, "w");

  if (shrinkmode) {
    easymode = fastmode = NO;
    altnormstyle = YES ;  // maybe bug if NO 
    printf ("shrinkmode set!\n");
  }

  if ((ldlimit <= 0) || (ldposlimit <= 0))
    ldregress = 0;

  numsnps =
    getsnps (snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks);

  numindivs = getindivs (indivname, &indivmarkers);

  if (id2pops != NULL) {
    setid2pops (id2pops, indivmarkers, numindivs);
  }

  k = getgenos (genotypename, snpmarkers, indivmarkers,
                numsnps, numindivs, nignore);

  if (elllistname != NULL) {
    nellindivs = numlines(elllistname) ; 
    ZALLOC (elllist, nellindivs, char *);
    nellindivs = loadlist(elllist, elllistname) ; 
    ZALLOC(ellindex, nellindivs, int) ; 
    printf("nellindivs: %3d\n", nellindivs) ; 
    printstrings(elllist, nellindivs) ; 
    regmode = YES ;
//  shrinkmode = NO ;
    ZALLOC(ellindivs, nellindivs, Indiv *) ;
  }
  if (maxpops < 0) maxpops = MAXPOPS ; 
  if (maxpops > MAXPOPS) { 
   printf("maxpops: %d is large!\n", maxpops) ;
  }

  if (poplistname != NULL) {
    numeg = numlines(poplistname) ; 
    ZALLOC (eglist, numeg, char *);
    numeg = loadlist (eglist, poplistname);
    seteglist (indivmarkers, numindivs, poplistname);
  }
  else {
    setstatus (indivmarkers, numindivs, NULL);
    ZALLOC (eglist, maxpops, char *);
    numeg = makeeglist (eglist, maxpops, indivmarkers, numindivs);
  }
  for (i = 0; i < numeg; i++) {
    /* printf("%3d %s\n",i, eglist[i]) ; */
  }

  if (topright != NULL) { 
   toprightindex = indxstring(eglist, numeg, topright) ;
  }

  nindiv = 0;
  for (i = 0; i < numindivs; i++) {
    indx = indivmarkers[i];
    indx -> flag = 0 ;
    if (indx->affstatus == YES)
      ++nindiv;
  }

  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    chrom = cupt->chrom;
    if ((noxdata) && (chrom == (numchrom + 1))) {
      cupt->ignore = YES;
      logdeletedsnp (cupt->ID, "chrom-X", deletesnpoutname);
      continue ;
    }
    if (chrom == 0) {
      cupt->ignore = YES;
      logdeletedsnp (cupt->ID, "chrom-0", deletesnpoutname);
      continue ;
    }
    if (chrom > (numchrom + 1)) {
      cupt->ignore = YES;
      logdeletedsnp (cupt->ID, "chrom-big", deletesnpoutname);
      continue ;
    }
    if (chrom == zchrom) { 
      cupt->ignore = YES;
      logdeletedsnp (cupt->ID, "chrom-deltet", deletesnpoutname);
      continue ;
    }
  }

  tt = 0 ;
  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    pos = nnint (cupt->physpos);
    if ((xchrom > 0) && (cupt->chrom != xchrom)) {
      cupt->ignore = YES;
      logdeletedsnp (cupt->ID, "not-chrom", deletesnpoutname);
    }
    if ((xchrom > 0) && (pos < lopos)) {
      cupt->ignore = YES;
      logdeletedsnp (cupt->ID, "lopos", deletesnpoutname);
    }
    if ((xchrom > 0) && (pos > hipos)) {
      cupt->ignore = YES;
      logdeletedsnp (cupt->ID, "hipos", deletesnpoutname);
    }
    if (cupt->ignore)
      continue;
    if (numvalidgtx (indivmarkers, cupt, YES) <= 1) {
      ++tt ; 
      if (verbose) printf ("nodata: %20s\n", cupt->ID);
      cupt->ignore = YES;
      logdeletedsnp (cupt->ID, "nodata", deletesnpoutname);
    }
  }
  printf("snps deleted (nodata): %d.  deletesnpoutname: for details", tt) ;

  if (killr2) {
    nkill =
      killhir2 (snpmarkers, numsnps, numindivs, r2physlim, r2genlim,
                r2thresh);
    if (nkill > 0)
      printf ("killhir2.  number of snps killed: %d\n", nkill);
  }

  if (xregionname) {
    excluderegions (xregionname, snpmarkers, numsnps, deletesnpoutname);
  }

  if (nhwfilter > 0) {
    hwfilter (snpmarkers, numsnps, numindivs, nhwfilter, deletesnpoutname);
  }

  ZALLOC (vv, numindivs, int);
  numvalidgtallind (vv, snpmarkers, numsnps, numindivs);
  for (i = 0; i < numindivs; ++i) {
     indx = indivmarkers[i];
    if (vv[i] == 0) {
      indx->ignore = YES;
    }
    if (doproject == NO) { 
      t = indxindex(eglist, numeg, indx -> egroup) ; 
      if (t<0) indx -> ignore = YES ;
    }
  }
  free (vv);

  numsnps = rmsnps (snpmarkers, numsnps, deletesnpoutname);     //  rid ignorable snps

  
  nindignore = 0 ; 
  for (i=0; i<numindivs; ++i)  { 
   indx = indivmarkers[i] ; 
   t = strcmp(indx -> egroup, "Ignore") ; 
   if (t==0) indx -> ignore = YES ;
   if (indx -> ignore == YES) ++nindignore ;
  }

  // printf("zzz %d %d\n", shrinkmode, nindignore) ; 
  if (shrinkmode && (nindignore > 0)) { 
   printf("removing individuals set Ignore\n") ; 
   numindivs = rmindivs(snpmarkers, numsnps, indivmarkers, numindivs) ;
  }


  if (missingmode) {
    setmiss (snpmarkers, numsnps);
    fancynorm = NO;
  }

  if (weightname != NULL) {
    weightmode = YES;
    getweights (weightname, snpmarkers, numsnps);
  }
  if (ldregress > 0) {
    ZALLOC (ldvv, ldregress * numindivs, double);
    ZALLOC (ldsnpbuff, ldregress, int); // index of snps 
  }

  ZALLOC (xindex, numindivs, int);
  ZALLOC (xindlist, numindivs, Indiv *);
  ZALLOC (xsnplist, numsnps, SNP *);

  if (popsizelimit > 0) {
    setplimit (indivmarkers, numindivs, eglist, numeg, popsizelimit);
  }


  /*  Load non-ignored individuals into xindlist,xindex:
   *  xindex[i]   = index into indivmarkers
   *  xindlist[i] = pointer to Indiv struct   */

  ZALLOC (xtypes, numindivs, int);
  ivclear(xtypes, -2, numindivs) ;



  /*  Load non-ignored SNPs into xsnplist:
   *  xsnplist[i] = pointer to SNP struct     */

  nrows = loadindx (xindlist, xindex, indivmarkers, numindivs);
  ncols = loadsnpx (xsnplist, snpmarkers, numsnps, indivmarkers);

  for (i = 0; i < nrows; i++) {
    indx = xindlist[i];
    k = indxindex (eglist, numeg, indx->egroup);
    xtypes[i] = k;
  }

  printf ("number of samples used: %d number of snps used: %d\n", nrows, ncols) ; 
  printf ("number of pops for axes: %d\n",  numeg) ;

  if (shrinkmode) {
    ZALLOC (mmat, nrows * ncols, double);
    regmode = YES;
  }

  if (numeigs==0) fastmode = NO ;

  if (fastmode) {

    setgval (xsnplist, nrows, indivmarkers, numindivs, xindex, xtypes, ncols);
// side-effect monomorphic snps -> ignore

    ZALLOC (evals, numeigs, double);
    ZALLOC (evecs, numeigs * nrows, double);

    kjg_fpca (numeigs, fastdim, fastiter, evals, evecs);

    if (verbose) {
      printf ("##bug: \n");
      printmat (evals, 1, numeigs);
      printmat (evecs, 1, 20);
    }

    transpose (evecs, evecs, nrows, numeigs);

    printevecs (xsnplist, indivmarkers, xindlist,
                numindivs, ncols, nrows, numeigs, evecs, evals, ofile);


    printf ("end of smartpca(fastmode)\n");
    return 0;

  }


  /* printf("## nrows: %d  ncols  %d\n", nrows, ncols) ; */
  tt = (int) sqrt (BIGINT);
  if (nrows >= tt)
    fatalx
      ("size of GRM matrix too large: nrows: %d\n Run with fastmode: YES perhaps\n",
       nrows);

  ZALLOC (xmean, ncols, double);
  ZALLOC (xfancy, ncols, double);

  ZALLOC (XTX, nrows * nrows, double);
  ZALLOC (evecs, nrows * nrows, double);
  if ((!usepopsformissing) && (ldregress == 0)) {
    // should not use lookup table if
    // - usepopsformissing is set (since each population may have a different
    //   mean), or
    // - ldregress > 0
#ifdef __LP64__
    blocksize = 20;
    ZALLOC (partial_sum_lookup_buf, 131072, double);
#else
    blocksize = 10;
    ZALLOC (partial_sum_lookup_buf, 65536, double);
#endif
    ZALLOC (binary_rawcol, nrows, int);
    ZALLOC (binary_cols, nrows, uintptr_t);
    ZALLOC (binary_mmask, nrows, uintptr_t);
    ZALLOC (tblock, 3 * blocksize, double);
  }
  else {
    ZALLOC (tblock, nrows * blocksize, double);
  }


  ZALLOC (lambda, nrows, double);
  ZALLOC (esize, nrows, double);
  ZALLOC (cc, (nrows > 3) ? nrows : 3, double);
  ZALLOC (ww, nrows +5, double);
  ZALLOC (badlist, nrows, int);

  blocksize = MIN (blocksize, ncols);

  // xfancy is multiplier for column xmean is mean to take off
  // badlist is list of rows to delete (outlier removal) 

  numoutiter = 1;

  if (numoutliter >= 1) {
    numoutiter = numoutliter + 1;
    ZALLOC (outinfo, nrows, OUTLINFO *);
    for (k = 0; k < nrows; k++) {
      ZALLOC (outinfo[k], 1, OUTLINFO);
    }
    /* fprintf(outlfile, "##%18s %4s %6s %9s\n", "ID", "iter","eigvec", "score") ; */
    setoutliermode (outliermode);
  }
  else
    setoutliermode (2);

  // try to autodetect number of (virtual) processors, and use that number to
  // set the thread count.  allow the user to override this in the future
  if (thread_ct_config <= 0) {
    i = sysconf (_SC_NPROCESSORS_ONLN);
    if (i == -1) {
      thread_ct = 1;
    }
    else {
      thread_ct = i;
    }
  }
  else {
    thread_ct = thread_ct_config;
  }
  if (thread_ct > 8) {
    if (thread_ct > MAX_THREADS) {
      thread_ct = MAX_THREADS;
    }
    else {
      thread_ct--;
    }
  }
  if (thread_ct > nrows * 2) {
    thread_ct = nrows / 2;
    if (!thread_ct) {
      thread_ct = 1;
    }
  }
  printf ("Using %u thread%s%s.\n", thread_ct, (thread_ct == 1) ? "" : "s",
          (partial_sum_lookup_buf) ? ", and partial sum lookup algorithm" :
          "");
  triangle_fill (g_thread_start, nrows, thread_ct, 0, 1, 0, 1);

  nkill = 0;

  for (outliter = 1; outliter <= numoutiter; ++outliter) {

    if (fstonly) {
      setidmat (XTX, nrows);
      vclear (lambda, 1.0, nrows);
      break;
    }
    if (outliter > 1) {
      ncols = loadsnpx (xsnplist, snpmarkers, numsnps, indivmarkers);
    }

    vzero (XTX, (nrows * (nrows + 1)) / 2);
    xblock = 0;

    vzero (xmean, ncols);
    vclear (xfancy, 1.0, ncols);

    nused = 0;
    for (i = 0; i < nrows; i++) {
      indx = xindlist[i];
      k = indxindex (eglist, numeg, indx->egroup);
      xtypes[i] = k;
    }

    numld = 0;
    lastldchrom = -1;
    ynumsnps = 0;
    if (partial_sum_lookup_buf) {
      for (i = 0; i < nrows; i++) {
        binary_cols[i] = 0;
      }
      for (i = 0; i < nrows; i++) {
        binary_mmask[i] = 0;
      }
      vzero (tblock, 3 * blocksize);
    }
    else {
      vzero (tblock, nrows * blocksize);
    }
    for (i = 0; i < ncols; i++) {
      cupt = xsnplist[i];
      chrom = cupt->chrom;
      if (!partial_sum_lookup_buf) {
        tt =
          getcolxz (cc, cupt, xindex, xtypes, nrows, i, xmean, xfancy,
                    &n0, &n1);
      }
      else {
        tt =
          getcolxz_binary1 (binary_rawcol, cc, cupt, xindex, nrows, i,
                            xmean, xfancy, &n0, &n1);
      }


      t = MIN (n0, n1);

      if ((t < minallelecnt) || (tt > maxmissing) || (tt < 0) || (t == 0)) {
        t = MAX (t, 0);
        tt = MAX (tt, 0);
        cupt->ignore = YES;
        logdeletedsnp (cupt->ID, "minallelecnt", deletesnpoutname);
        vzero (cc, nrows);
        if (nkill < 10)
          printf (" snp %20s ignored . allelecnt: %5d  missing: %5d\n",
                  cupt->ID, t, tt);
        ++nkill;
        continue;
      }

      if (lastldchrom != chrom)
        numld = 0;

      if (!partial_sum_lookup_buf) {
        if (weightmode) {
          vst (cc, cc, xsnplist[i]->weight, nrows);
        }


        if (ldregress > 0) {

          t = ldregx (ldvv, cc, ww, numld, nrows, ldr2lo, ldr2hi);
           if (t < 2) {
            bumpldvv (ldvv, cc, &numld, ldregress, nrows, ldsnpbuff, i);
            lastldchrom = chrom;
            ynumsnps += asum2 (ww, nrows) / asum2 (cc, nrows);
            // don't need to think hard about how cc is normalizes
          }
          else {
            // Ignore this SNP and exclude from further regressions (*ww is returned as all zeroes)
            bumpldvv (ldvv, ww, &numld, ldregress, nrows, ldsnpbuff, i);
            lastldchrom = chrom;
          }
          copyarr (ww, cc, nrows);
        }
        else
          ++ynumsnps;
        copyarr (cc, tblock + xblock * nrows, nrows);
      }
      else {
        getcolxz_binary2 (binary_rawcol, binary_cols, binary_mmask,
                          xblock, nrows);
        if (weightmode) {
          vst (cc, cc, xsnplist[i]->weight, 3);
        }
        ++ynumsnps;
        copyarr (cc, &(tblock[xblock * 3]), 3);
      }

      ++xblock;
      ++nused;

/** this is the key code to parallelize */
      if (xblock == blocksize) {
        if (partial_sum_lookup_buf) {
          domult_increment_lookup (threads, thread_ct, XTX, tblock,
                                   binary_cols, binary_mmask, xblock,
                                   nrows, partial_sum_lookup_buf);
          for (j = 0; j < nrows; j++) {
            binary_cols[j] = 0;
          }
          for (j = 0; j < nrows; j++) {
            binary_mmask[j] = 0;
          }
          vzero (tblock, 3 * blocksize);
        }
        else {
          domult_increment_normal (threads, thread_ct, XTX, tblock,
                                   xblock, nrows);
          vzero (tblock, nrows * blocksize);
        }
        xblock = 0;
      }
    }

    if (xblock > 0) {
      if (partial_sum_lookup_buf) {
        domult_increment_lookup (threads, thread_ct, XTX, tblock,
                                 binary_cols, binary_mmask, xblock,
                                 nrows, partial_sum_lookup_buf);
      }
      else {
        domult_increment_normal (threads, thread_ct, XTX, tblock,
                                 xblock, nrows);
      }
    }
    symit2 (XTX, nrows);
    printf ("total number of snps killed in pass: %d  used: %d\n", nkill,
            nused);

    if (verbose) {
      printdiag (XTX, nrows);
    }

    y = trace (XTX, nrows) / (double) (nrows - 1);
    if (isnan (y))
      fatalx ("bad XTX matrix\n");
    /* printf("trace:  %9.3f\n", y) ; */
    if (y <= 0.0)
      fatalx ("XTX has zero trace (perhaps no data)\n");
    vst (XTX, XTX, 1.0 / y, nrows * nrows);


    eigvecs (XTX, lambda, evecs, nrows);
    
    if (verbose) printf("zztrace: %9.3f nrows: %3d xmat[0] %12.6f lambda[0] %12.6f\n", y, nrows, XTX[0], lambda[0]) ;  

    nrxtx = nrows ;
// eigenvalues are in decreasing order 

    if (outliter > numoutliter)
      break;
    // last pass skips outliers 
    numoutleigs = MIN (numoutleigs, nrows - 1);
    nbad =
      ridoutlier (evecs, nrows, numoutleigs, outlthresh, badlist, outinfo);
    if (nbad == 0)
      break;
    for (i = 0; i < nbad; i++) {
      j = badlist[i];
      indx = xindlist[j];
      outpt = outinfo[j];
      fprintf (outlfile,
               "REMOVED outlier %s iter %d evec %d sigmage %.3f pop: %s\n",
               indx->ID, outliter, outpt->vecno, outpt->score, indx->egroup);
      indx->ignore = YES;
    }
    nrows = loadindx (xindlist, xindex, indivmarkers, numindivs);
    printf ("number of samples after outlier removal: %d\n", nrows);
  }

 if (toprightindex >=0) {  
  
  ZALLOC(emean, numeg, double) ; 
  ZALLOC(flip, numeigs, int) ; 
  ivclear(flip, -1, numeigs) ; 
  for (j=0; j< MIN(numeigs, 2); ++j) { 
   xpt = evecs + j*nrows ; 
   calcmean (emean, xpt, nrows, xtypes, numeg);
   y = emean[toprightindex] ; 
   if (y<0) { vst(xpt, xpt, -1, nrows) ; 
    printf("eigenvector %d flipped!\n", j+1) ;
    flip[j] = 1 ; 
   }
  }
 }
  


// now load up elllist 
  for (j=0; j<nellindivs; ++j) { 
   k = indindex(indivmarkers, numindivs, elllist[j]) ; 
   if ((k>=0) && (indivmarkers[k] -> ignore)) k = -1 ;
   ellindex[j] = k ; 
   if (k<0) fatalx("*** warning no ellipse for %s\n", elllist[j]) ; 
   indx = ellindivs[j] = indivmarkers[k] ; 
   printf("elllist: %3d %20s\n", indx -> ID) ; fflush(stdout) ;
  }

  if (outliername != NULL)
    fclose (outlfile);
  dumpgrm (XTX, xindex, nrows, ynumsnps, indivmarkers, numindivs, grmoutname);
  if (grmoutname != NULL)
    printf ("grm dumped\n");

  t = nrows - numeigs ;     
  if (t <= 10) {
     if (fastshrink) printf("*** turning off fastshrink: problem too small\n") ;
    fastshrink = NO ;
  }

  m = numgtz (lambda, nrows);
  /* printf("matrix rank: %d\n", m) ; */
  if (m == 0)
    fatalx ("no data\n");

  /* Now, print Tracy-Widom stats, if twtable is valid */
  if (settwxtable (twxtabname) < 0) {
    printf ("\n## To get Tracy-Widom statistics: recompile smartpca with");
    printf (" TWTAB correctly specified in Makefile, or\n");
    printf ("   just run twstats (see README file in POPGEN directory)\n");
  }
  else {
    /* *** START of code to print Tracy-Widom statistics */
    ynrows = (double) nrows;
    mpestimate(lambda, m, &zeffect, &zsigma) ;

    if (mpshrink) printf("zeffect:  %9.3f  zsigma: %9.3f\n", zeffect, zsigma) ;

    if (fstonly == NO) {
      printf ("\n## Tracy-Widom statistics: rows: %d  cols: %d\n", nrows,
              ncols);
      y = -1.0;
      printf ("%4s  %12s", "#N", "eigenvalue");
      printf ("%12s", "difference");
      printf (" %9s %12s", "twstat", "p-value");
      printf (" %9s", "effect. n");
      printf ("\n");
    }


    for (i = 0; i < m; ++i) {
      if (fstonly)
        break;
      zn = znval;
      if (zn > 0)
        zn = MAX (ynrows, zn);
      tail = dotwcalc (lambda + i, m - i, &tw, &zn, &zvar, nostatslim);
//    printf("zzeff: %3d %12.3f %12.3f\n", i, zeffect, zn) ; 
      if (mpshrink) {
       zn = zeffect ; 
      }
      esize[i] = zeffect ;  // not used

      printf ("%4d  %12.6f", i + 1, lambda[i]);
      if (i == 0)
        printf ("%12s", "NA");
      else
        printf ("%12.6f", lambda[i] - lambda[i - 1]);
      if (tail >= 0.0)
        printf (" %9.3f %12.6g", tw, tail);
      else
        printf (" %9s %12s", "NA", "NA");
      if (zn > 0.0) {
        printf (" %9.3f", zn);
      }
      else {
        printf (" %9s", "NA");
      }
      printf ("\n");
    }
    /* END of code to print Tracy-Widom statistics */
  }

  numeigs = MIN (numeigs, nrows);
  numeigs = MIN (numeigs, ncols);

  numshrink = numeigs ;  

  if (autoshrink) { 
   numshrink = setnstw(lambda, m, nostatslim) ;
   if (numshrink<0) numshrink = 0 ;
  }


    if (fastshrink) { 
     zn = znval ;  
     dotwcalc (lambda + numshrink, m - numshrink, &tw, &zn, &zvar, nostatslim);
     if (mpshrink) zn = zeffect ;  
     zfp = nrows ; 
     zfn = zn ; 
     printf("fastshrink:  nsamps: %6.0f  numshrink: %d  neff: %9.3f\n", zfp, numshrink, zfn) ; 
     printf("nrows ncols %d %d\n", nrows, ncols) ;  
     ZALLOC(edgarw, m, double) ;  
     vclear(edgarw, 1.0, m) ;  
/*
     y = 4339.067 ; 
     yjfac = (double) nrows * y - (double) numeigs ; 
     estedgar(edgarw, lambda, numshrink, m-numshrink, y, yjfac) ; 
     printf("shrink factors: ") ;
     printmatl(edgarw, 1, numeigs) ; 
     printnl() ;
     printnl() ;

     y = (double) ncols / (double) nrows ; 
     yjfac = (double) nrows * y - (double) numeigs ; 
     estedgar(edgarw, lambda, numshrink, m-numshrink, y, yjfac) ; 
     printf("shrink factors: ") ;
     printmatl(edgarw, 1, numeigs) ; 
     printnl() ;
     printnl() ;

*/
     y = zfn/zfp ; 
     yjfac = (double) nrows * y - (double) numeigs ; 
     estedgar(edgarw, lambda, numshrink, m-numshrink, y, yjfac) ; 
     printf("shrink factors: ") ;
     printmatl(edgarw, 1, numeigs) ; 
     printnl() ;
     printnl() ;

     for (i=0; i<nrows; ++i)  {
      indx = xindlist[i];
      indx -> flag = 7777 ;
     }
    }


  if (numeigs > 0) {
    if (outputvname != NULL) {
      openit (outputvname, &ovfile, "w");
      for (j = 0; j < nrows; j++) {
        fprintf (ovfile, "%12.6f\n", lambda[j]);
      }
      fclose (ovfile);
    }

    fprintf (ofile, "%20s ", "#eigvals:");
    for (j = 0; j < numeigs; j++) {
      fprintf (ofile, "%9.3f ", lambda[j]);
    }
    fprintf (ofile, "\n");
    ZALLOC (fvecs, nrows * numeigs, double);
    ZALLOC (fxvecs, nrows * numeigs, double);
    ZALLOC (fxscal, numeigs, double);

    ZALLOC (ffvecs, ncols * numeigs, double);
    ZALLOC (xrow, ncols, double);
    setfvecs (fvecs, evecs, nrows, numeigs);

    if (easymode) {
      for (j = 0; j < numeigs; j++) {
        xpt = fvecs + j * nrows;
        y = asum2 (xpt, nrows);
        vst (xpt, xpt, 1.0 / sqrt (y), nrows);  // norm 1
      }
      for (i = 0; i < nrows; i++) {
        indx = xindlist[i];
        fprintf (ofile, "%20s ", indx->ID);
        for (j = 0; j < numeigs; j++) {
          xpt = fvecs + j * nrows;
          y = xpt[i];
          if (indx -> flag == 7777) y *= edgarw[j] ;  
          fprintf (ofile, "%10.4f  ", y);  // easymode
        }
        fprintf (ofile, "%15s\n", indx->egroup);
      }
      if (pubmean) {

        ZALLOC (wmean, numeg, double);
        ZALLOC (elist, numeg, char *);

        for (j = 0; j < numeigs; j++) {
          xpt = fvecs + j * nrows;
          calcpopmean (wmean, elist, xpt, eglist, numeg, xtypes, nrows);
          printf ("eig: %d ", j + 1);
          printf ("min: %s %9.3f  ", elist[0], wmean[0]);
          printf ("max: %s %9.3f  ", elist[numeg - 1], wmean[numeg - 1]);
          printnl ();
          for (k = 0; k < numeg; ++k) {
            printf ("%20s ", elist[k]);
            printf (" %9.3f\n", wmean[k]);
          }
        }
      }

      printf ("## easymode set. end of smartpca run\n");
      return 0;
    }
    for (i = 0; i < ncols; i++) {
      cupt = xsnplist[i];
      getcolxf (cc, cupt, xindex, nrows, i, NULL, NULL);

      for (j = 0; j < numeigs; j++) {
        for (k = 0; k < nrows; k++) {
          ffvecs[j * ncols + i] += fvecs[j * nrows + k] * cc[k];
        }
      }
    }

    ZALLOC (eigkurt, numeigs, double);
    ZALLOC (eigindkurt, numeigs, double);

    for (j = 0; j < numeigs; ++j) {
      eigkurt[j] = kurtosis (ffvecs + j * ncols, ncols);
      eigindkurt[j] = kurtosis (fvecs + j * nrows, nrows);
    }

    for (i = 0; i < nrows; i++) {

      indx = xindlist[i];
      k = indxindex (eglist, numeg, indx->egroup);
      xtypes[i] = k;

      loadxdataind (xrow, xsnplist, xindex[i], ncols);
      fixxrow (xrow, xmean, xfancy, ncols);

      for (j = 0; j < numeigs; j++) {

        xpt = ffvecs + j * ncols;
        y = fxvecs[j * nrows + i] = vdot (xrow, xpt, ncols);
        fxscal[j] += y * y;

      }
    }

    for (j = 0; j < numeigs; j++) {
      y = fxscal[j];
      fxscal[j] = 1.0 / sqrt (y);       // standard
    }


    ZALLOC (acoeffs, numindivs * numeigs, double);
    ZALLOC (bcoeffs, numindivs * numeigs, double);
    ZALLOC (xcoeffs, numindivs*numeigs, double) ; 
    if (partial_sum_lookup_buf) {
      free (partial_sum_lookup_buf);
      free (binary_rawcol);
      free (binary_cols);
      free (binary_mmask);
    }
    free (tblock);



    printf ("%12s %4s %9s %9s\n", "kurtosis", "", "snps", "indivs");

    for (j = 0; j < numeigs; ++j) {
      y1 = eigkurt[j];
      y2 = eigindkurt[j];
      printf ("%12s %4d %9.3f %9.3f\n", "eigenvector", j + 1, y1, y2);
    }

    if (regmode) {

      ZALLOC(eigscmat, numeigs * numeigs, double) ; 
      ZALLOC(eigscale, numeigs, double) ;
      lsqproj(-99, xsnplist, ncols, indivmarkers, numindivs, fxscal, ffvecs, acoeffs, bcoeffs, xtypes, numeg) ; 

      if (shrinkp_done == NO) { 
       seteigscale(eigscale, acoeffs, bcoeffs, xindex, nrows, numeigs) ; 
      }

      else {
       seteigscale(eigscale, xcoeffs, bcoeffs, xindex, nrows, numeigs) ; 
      }
 
      setdiag(eigscmat, eigscale, numeigs) ; 
      mulmat(acoeffs, eigscmat, acoeffs,  numeigs, numeigs, numindivs) ; 

   }

// evec output

    for (i = 0; i < numindivs; i++) {
      indx = indivmarkers[i];
      indx -> idnum = i ; 
    }
    for (i = 0; i < numindivs; i++) {
      if (shrinkmode) break ; 
      indx = indivmarkers[i];
      if (indx->ignore) continue;
      fprintf (ofile, "%20s ", indx->ID);
      for (j = 0; j < numeigs; j++) {
        y = acoeffs[j * numindivs + i];
        if (indx -> flag == 7777) y *= edgarw[j] ;  
        if (hiprec) fprintf (ofile, "%12.6f  ", y);
        else fprintf (ofile, "%10.4f  ", y);   // main output
      }
      if (qtmode) {
        fprintf (ofile, "%15.6e\n", indx->qval);
      }
      else {
        fprintf (ofile, "%15s\n", indx->egroup);
      }
    }

  if (nellindivs>0) { 

   nblocks = numblocks (snpmarkers, numsnps, blgsize);
   if (nblocks <= 1) {
    printf("failed to set blocks...\n") ; 
    return -1 ;  
   }

  ++nblocks ;  // dangerous bend 
  printf ("number of blocks for block jackknife (ell): %d\n", nblocks);

  ZALLOC (blstart, nblocks, int);
  ZALLOC (blsize, nblocks, int);
  blcnt = initarray_2Dint(nellindivs, nblocks, 0) ;

  t = numeigs * nellindivs ; 
  aellc = initarray_2Ddouble(nblocks, t, 0) ; 
  bellc = initarray_2Ddouble(nblocks, t, 0) ; 

  ZALLOC(awk, t, double) ; 
  ZALLOC(bwk, t, double) ; 

  setblocks (blstart, blsize, &xnblocks, xsnplist, ncols, blgsize);
  fixwt (xsnplist, ncols, 1.0);

    for (i = 0; i < ncols; ++i) {
      cupt = xsnplist[i];
      if (cupt -> ignore) continue ; 
      if (cupt -> tagnumber < 0) continue ; 
      ++ cupt -> tagnumber ;
    }

// set up block count for ellidivs[jj] 
   for (jj=0; jj < nellindivs; ++jj) { 
      i = ellindex[jj] ; 
      if (i<0) fatalx("bad bug! bad index: %s\n", ellindivs[jj] -> ID) ;   
      countsn(blcnt[jj], nblocks, xsnplist, ncols, i) ; 
   }

   for (bl=0; bl<nblocks; ++bl) {  
    aco = aellc[bl] ;  
    bco = bellc[bl] ;  
    lsqproj(bl, xsnplist, ncols, ellindivs, nellindivs, fxscal, ffvecs, aco, NULL, xtypes, numeg) ; 
    mulmat(aco, eigscmat, aco,  numeigs, numeigs, nellindivs) ; 
   }
    
 }


 if (shrinkmode) { 

  for (i = 0; i < ncols; ++i) {
    cupt = xsnplist[i];
    getcolxf (cc, cupt, xindex, nrows, i, NULL, NULL);
    for (j = 0; j < nrows; ++j) {
      mmat[j * ncols + i] = cc[j];
    }
  }

  fclose (ofile) ; 

  if (outputname != NULL)
    openit (outputname, &ofile, "w");
  else
    ofile = stdout;
  doshrinkp (mmat, nrows, ncols, xindex, xsnplist, xcoeffs) ;

  

  shrinkp_done = YES ;

 for (j=0; j<numeigs; ++j) { 
   if (nellindivs<1) break ; 
   for (jj=0; jj<nellindivs; ++jj) { 
    indx = ellindivs[jj] ; 
    i = indx -> idnum ; 
    awk[j*nellindivs+jj] = acoeffs[j*numindivs+i] ; 
    bwk[j*nellindivs+jj] = xcoeffs[j*numindivs+i] ; 
   }
   aco = awk + j*nellindivs ; 
   bco = bwk + j*nellindivs ; 
   y = vdot(aco, bco, nellindivs) / vdot(aco, aco, nellindivs) ;
   printf("coeff multiplier (eigenvector %d): %9.3f\n", y, j) ; 
   vst(aco, aco, y, nellindivs) ; 
   for (jj=0; jj<nellindivs; ++jj) { 
    indx = ellindivs[jj] ; 
    i = indx -> idnum ; 
    acoeffs[j*numindivs+i] =  awk[j*nellindivs+jj]; 
   }
   
   for (bl=0; bl<nblocks; ++bl) { 
    aco = aellc[bl] + j*nellindivs ; 
    vst(aco, aco, y, nellindivs) ; 
   }
  }

 }

/**
 aco = aellc[0] ; 
 bco = aellc[1] ;
 for (j=0; j<numeigs; ++j) { 
   for (jj=0; jj<nellindivs; ++jj) { 
    indx = ellindivs[jj] ; 
    i = indx -> idnum ; 
    printf("zzcheck: %s %d ", indx -> ID, j) ; 
    if (shrinkmode) printf("%10.4f ", xcoeffs[j*numindivs+i]) ;
    printf("%10.4f ", acoeffs[j*numindivs+i]) ;
    printf("%10.4f ", aco[j*nellindivs+jj]) ;
    printf("%10.4f ", bco[j*nellindivs+jj]) ;
    printnl() ;
  }
 }
*/

 ZALLOC(jwt, nblocks, double) ; 
 jmean = initarray_2Ddouble(nblocks, numeigs, -999) ; 
 vvest = initarray_2Ddouble(nellindivs, numeigs, -999) ; 
 vvar  = initarray_2Ddouble(nellindivs, numeigs*numeigs, -999) ; 

 fflush(stdout) ; 

 t = 0 ; 
 if (elloutname != NULL) openit (elloutname, &ellfile, "w") ; 
  else ellfile = stdout ;

 if (ellconf > 0) {  fprintf(ellfile, "%9.3f ", ellconf) ; 
  fprintf(ellfile, ":: %s ", outputname) ;
  fprintf(ellfile, "\n") ;

 for (jj=0; jj < nellindivs; ++jj) { 
    indx = ellindivs[jj] ; 
    floatit(jwt, blcnt[jj], nblocks) ; 
    for (bl=0; bl<nblocks; ++bl) { 
     for (j=0; j<numeigs; ++j) { 
      jmean[bl][j] = aellc[bl][j*nellindivs+jj] ; 
     }
   }
  
  if (numeigs<2) continue ; 
  wjackvest(vvest[jj], vvar[jj], numeigs, jmean[0], jmean, jwt, nblocks) ; 
//   printf("zzbugq: ") ;  printmatl(vvar[jj], numeigs, numeigs) ; 
  vest = vvest[jj] ; 
  var = vvar[jj] ; 
  copyarr(vest, zzest, 2) ; 
   zzvar[0] = var[0] ; 
   zzvar[1] = zzvar[2] = var[1] ; 
   zzvar[3] = var[1*numeigs+1] ; 
  
//  printf("zzbugz %d ", jj) ;  printmatl(zzvar, 2, 2) ;
  emean = jmean[0] ; 
  fprintf(ellfile, "sample: %20s\n", indx -> ID) ;
  printmatlfile(emean, 1, 2, ellfile) ; 
  printmatlfile(zzest, 1, 2, ellfile) ; 
  fprintf(ellfile, "\n") ; 
  printmatlfile(zzvar, 2, 2, ellfile) ;
//  printf("zzbugy: ") ;  printmatl(zzvar, 2, 2) ; 
  if (ellconf<0) continue ; 
  eigvecs(zzvar, zlambda, zvecs, 2) ; 
  zval = critchi(2, 1.0 - ellconf) ; 
  za1 = 2.0 * sqrt(zval*zlambda[0]) ;  // diameter 
  za2 = 2.0 * sqrt(zval*zlambda[1]) ;  // diameter 
//  printf("zzbugza %12.6f %12.6f : ", za1, za2) ; printmatl(zlambda, 1, 2) ; 
  printmatl(zvecs, 2, 2) ; 
  printnl() ; 
  fprintf(ellfile, "ellcoords: ") ;
  y1 = zvecs[1] ; 
  y2 = zvecs[0] ; 
  if (fabs(y2) < fabs(y1)*1.0e-20) zangle = 90 ; 
  else zangle = rad2deg(atan(y1/y2)) ; 
  fprintf(ellfile, "%10.4f %10.4f ", emean[0], emean[1]) ;
  ww[0] = za1/2 ;  
  ww[1] = za2/2 ;  
  ww[2] = -deg2rad(zangle) ;
  printmatwlxfile(ww, 1, 3, 3, ellfile) ;
  fprintf(ellfile, " :: %9.3f\n", ellconf) ;
  
  ++t ;  
 }
 

  fclose(ellfile) ;
 }


// numeigs > 0 
}
//  printf("zzzend\n") ;

// output files
  settersemode (YES);

  ZALLOC (xpopsize, numeg, int);
  for (i = 0; i < numeg; i++) {
    xpopsize[i] = 0;
  }
  for (i = 0; i < nrows; i++) {
    k = xtypes[i];
    ++xpopsize[k];
  }

  for (i = 0; i < numeg; i++) {
    printf ("population: %3d %20s %4d", i, eglist[i], xpopsize[i]);
    if (xpopsize[i] == 0)
      printf (" ***");
    printnl ();
  }


  if (numeg == 1)
    dotpopsmode = NO;

/**
  if (dotpopsmode == NO) {
    writesnpeigs (snpeigname, xsnplist, ffvecs, numeigs, ncols);
    printxcorr (XTX, nrows, xindlist);
    if (snpoutfilename != NULL) {
      outfiles (snpoutfilename, indoutfilename, genooutfilename,
                snpmarkers, indivmarkers, numsnps, numindivs, packout,
                ogmode);
    }

    if (!shrinkmode) {
      ymem = calcmem(1)/1.0e6 ;
      printf("##end of smartpca: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
    }

    fclose (ofile);
    return 0;

  }
*/

  ZALLOC (chitot, numeg * numeg, double);

  dotpops (XTX, eglist, numeg, xtypes, nrows);
  ZALLOC (fstans, numeg * numeg, double);
  ZALLOC (fstsd, numeg * numeg, double);
  ZALLOC (fstzz, numeg * numeg, double);
  ZALLOC (fstnum, numeg * numeg, int);

  setinbreed (inbreed);

  if (inbreed && (fstjack == YES)) {
    ZALLOC (inbans, numeg, double);
    ZALLOC (inbsd, numeg, double);
    doinbxx (inbans, inbsd, xsnplist, xindex, xtypes,
             nrows, ncols, numeg, blgsize, snpmarkers, indivmarkers);
    printf ("## inbreeding coeffs:   inbreed    std error\n");
    for (k1 = 0; k1 < numeg; ++k1) {
      printf (" %20s %10.4f %10.4f\n", eglist[k1], inbans[k1], inbsd[k1]);
    }
    free (inbans);
    free (inbsd);
  }

  ZALLOC(fstmean, numeg*numeg, double) ;

  dofstxx (fstmean, fstans, fstsd, fstnum, xsnplist, xindex, xtypes,
           nrows, ncols, numeg, blgsize, snpmarkers, indivmarkers);

   copyarr(fstans, fstzz, numeg*numeg) ; 
    for (k1 = 0; k1 < numeg; ++k1) {
      for (k2 = k1 + 1; k2 < numeg; ++k2) {
        fstzsc =  fstans[k1 * numeg + k2] / (fstsd[k1 * numeg + k2] + 1.0e-10) ; 
        fstzz[k2*numeg + k1] = clip(fstzsc, -9.999999, 9.999999) ; 
      }
     }


  printmega(megaoutname, eglist, numeg, fstans) ;
   
  if (((phylipname == NULL) && (numeg > 10)) || fstverbose)  {
    printf ("## Fst statistics between populations:         fst       std error      Z\n");
    for (k1 = 0; k1 < numeg; ++k1) {
      if (inbreed && (xpopsize[k1] <= 1)) continue ; 
      for (k2 = k1 + 1; k2 < numeg; ++k2) {
        if (inbreed && (xpopsize[k2] <= 1)) continue ; 
        fstzsc =  fstans[k1 * numeg + k2] / (fstsd[k1 * numeg + k2] + 1.0e-10) ; 
        fnum = fstnum[k1*numeg+k2] ; 
        if (fsthiprec == NO) {
          printf (" %20s %20s %9.3f %10.4f", eglist[k1], eglist[k2],
                  fstans[k1 * numeg + k2], fstsd[k1 * numeg + k2]);
        }
        if (fsthiprec == YES) {
          printf (" %20s %20s %12.6f %12.6f", eglist[k1],
                  eglist[k2], fstans[k1 * numeg + k2],
                  fstsd[k1 * numeg + k2]);
        }
           printf("  %9.3f ", fstzsc) ;
           printf(" %12d", fnum) ; 
           printnl() ;
      }
    }
      printnl() ;
  }
  printf("## end of Fst statistics between populations\n") ;

  if (fstdetailsname != NULL) {
    fprintf
      (fstdetails, "## Fst statistics between populations:     fst     fstjack       std error      Z    snpnumber\n");
    for (k1 = 0; k1 < numeg; ++k1) {
      for (k2 = k1 + 1; k2 < numeg; ++k2) {
       if (fsthiprec == NO) { 
        fprintf (fstdetails, "F_st %20s %20s %12.6f %12.6f %12.6f",
                 eglist[k1], eglist[k2], fstmean[k1*numeg+k2], fstans[k1 * numeg + k2],
                 fstsd[k1 * numeg + k2]);
       }
       else { 
        fprintf (fstdetails, "F_st %20s %20s %15.9f %15.9f %15.9f",
                 eglist[k1], eglist[k2], fstmean[k1*numeg+k2], fstans[k1 * numeg + k2],
                 fstsd[k1 * numeg + k2]);
       }
       fstzsc =  fstans[k1 * numeg + k2] / (fstsd[k1 * numeg + k2] + 1.0e-10) ; 
       fnum =  fstnum[k1 * numeg + k2] ; 
       fprintf(fstdetails, " %9.3f", fstzsc) ;
       fprintf(fstdetails, " %6d", fnum) ; 
       fprintf(fstdetails, "\n") ;
      }
    }
  }

  if (fstz)  {
    printnl() ;
    if (fsthiprec==NO) { 
     printf("Fst: Z/F *1000 [clipped to -10,10]\n") ;
     printmatz5(fstzz, eglist, numeg) ; 
    }
    else  { 
     printf("Fst: Z/F *1000000 [clipped to -10,10]\n") ;
     printmatz10(fstzz, eglist, numeg) ; 
    }
  }


  printnl() ;
  if (phylipname != NULL) {
    openit (phylipname, &phylipfile, "w");
    fprintf (phylipfile, "%6d\n", numeg);
    sss[10] = CNULL;
    for (k1 = 0; k1 < numeg; ++k1) {
      strncpy (sss, eglist[k1], 10);
      fprintf (phylipfile, "%10s", sss);
      for (k2 = 0; k2 < numeg; ++k2) {
        y1 = fstans[k1 * numeg + k2];
        y2 = fstans[k2 * numeg + k1];
        fprintf (phylipfile, "%6.3f", (0.5 * (y1 + y2)));
      }
      fprintf (phylipfile, "\n");
    }
    fclose (phylipfile);
  }

  if ((numeg <= 10) || fstonly) {
    if (fsthiprec == NO) {
      printf ("fst *1000:");
      printnl ();
      printmatz5 (fstans, eglist, numeg);
      printnl ();
    }
    if (fsthiprec == YES) {
      printf ("fst *1000000:");
      printnl ();
      printmatz10 (fstans, eglist, numeg);
      printnl ();
    }
  }
  printf ("std. err.  * 1000000:\n");
  vst (fstsd, fstsd, 1000.0, numeg * numeg);
  printmatz5 (fstsd, eglist, numeg);

  printnl ();
  fflush (stdout);
  if (fstonly) {
    ymem = calcmem(1)/1.0e6 ;
    printf("##end of smartpca: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
    return 0;
  }
  vst (fstsd, fstsd, 1.0 / 1000.0, numeg * numeg);

  for (j = 0; j < numeigs; j++) {
    sprintf (sss, "eigenvector %d", j + 1);
    y = dottest (sss, evecs + j * nrows, eglist, numeg, xtypes, nrows);
  }

  printf
    ("\n## Statistical significance of differences beween populations:\n");
  printf
    ("                                pop1                  pop2      chisq          p-value   |pop1|   |pop2|\n");
  for (k1 = 0; k1 < numeg; ++k1) {
    if (fstonly)
      break;
    for (k2 = k1 + 1; k2 < numeg; ++k2) {
      ychi = chitot[k1 * numeg + k2];
      tail = rtlchsq (numeigs, ychi);
      printf ("popdifference:  %20s  %20s  %12.3f  %12.6g", eglist[k1],
              eglist[k2], ychi, tail);
      printf ("   %5d", xpopsize[k1]);
      printf ("   %5d", xpopsize[k2]);
      printf ("\n");
    }
  }
  printf ("\n");
  for (i = 0; i < ncols; i++) {
    if (markerscore == NO)
      break;
    cupt = xsnplist[i];
    getcolxf (cc, cupt, xindex, nrows, i, NULL, NULL);

    sprintf (sss, "%s raw", cupt->ID);
    dottest (sss, cc, eglist, numeg, xtypes, nrows);
    for (j = 0; j < numeigs; j++) {
      sprintf (sss, "%s subtract sing vec %d", cupt->ID, j + 1);
      y = vdot (cc, evecs + j * nrows, nrows);
      vst (ww, evecs + j * nrows, y, nrows);
      vvm (cc, cc, ww, nrows);
      dottest (sss, cc, eglist, numeg, xtypes, nrows);
    }
  }

  printxcorr (XTX, nrows, xindlist);


  writesnpeigs (snpeigname, xsnplist, ffvecs, numeigs, ncols);
  if (snpoutfilename != NULL) {
    outfiles (snpoutfilename, indoutfilename, genooutfilename,
              snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode);
  }

  if (!shrinkmode) {
    ymem = calcmem(1)/1.0e6 ;
    printf("##end of smartpca: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
    return 0;
  }

  printf ("##end of primary run \n");

 if (shrinkmode && (shrinkp_done == NO)) {
  for (i = 0; i < ncols; ++i) {
    cupt = xsnplist[i];
    getcolxf (cc, cupt, xindex, nrows, i, NULL, NULL);
    for (j = 0; j < nrows; ++j) {
      mmat[j * ncols + i] = cc[j];
    }
  }


  fclose (ofile);
  if (outputname != NULL)
    openit (outputname, &ofile, "w");
  else
    ofile = stdout;
  doshrinkp (mmat, nrows, ncols, xindex, xsnplist, xcoeffs);
 }

   ymem = calcmem(1)/1.0e6 ;
   printf("##end of smartpca: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;

  return 0;
}

void
readcommands (int argc, char **argv)
{
  int i;
  phandle *ph;
  int t;

  while ((i = getopt (argc, argv, "p:vV")) != -1) {

    switch (i) {

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


  if (parname == NULL) {
    fprintf (stderr, "no parameters\n");
    return;
  }

  pcheck (parname, 'p');
  printf ("parameter file: %s\n", parname);
  ph = openpars (parname);
  dostrsub (ph);

  getstring (ph, "genotypename:", &genotypename);
  getstring (ph, "snpname:", &snpname);
  getstring (ph, "indivname:", &indivname);
  getstring (ph, "poplistname:", &poplistname);
  getstring (ph, "snpeigname:", &snpeigname);
  getstring (ph, "snpweightoutname:", &snpeigname);     /* changed 09/18/07 */
  getstring (ph, "output:", &outputname);
  getstring (ph, "outputvecs:", &outputname);
  getstring (ph, "evecoutname:", &outputname);  /* changed 11/02/06 */
  getstring (ph, "outputvals:", &outputvname);
  getstring (ph, "evaloutname:", &outputvname); /* changed 11/02/06 */
  getstring (ph, "badsnpname:", &badsnpname);
  getstring (ph, "outliername:", &outliername);
  getstring (ph, "outlieroutname:", &outliername);      /* changed 11/02/06 */
  getstring (ph, "phylipname:", &phylipname);
  getstring (ph, "phylipoutname:", &phylipname);        /* changed 11/02/06 */
  getstring (ph, "weightname:", &weightname);
  getstring (ph, "fstdetailsname:", &fstdetailsname);
  getint(ph, "fstsnpout:", &fstsnpout) ;
  getstring (ph, "deletsnpoutname:", &deletesnpoutname);
  getstring (ph, "topright:", &topright);
  getstring (ph, "elloutname:", &elloutname);
  getstring (ph, "megaoutname:", &megaoutname);
  getdbl (ph, "ellconf:", &ellconf);
  getint (ph, "numeigs:", &numeigs);
  getint (ph, "maxpops:", &maxpops);
  getint (ph, "numoutevec:", &numeigs); /* changed 11/02/06 */
  getint (ph, "markerscore:", &markerscore);
  getint (ph, "chisqmode:", &chisqmode);
  getint (ph, "missingmode:", &missingmode);
  getint (ph, "shrinkmode:", &shrinkmode);
  getint (ph, "newshrink:", &newshrink);
  getint (ph, "fastshrink:", &fastshrink);
  getint (ph, "autoshrink:", &autoshrink);
  getint (ph, "mpshrink:", &mpshrink);
  getdbl (ph, "autothresh:", &autothresh);
  getint (ph, "fancynorm:", &fancynorm);
  getint (ph, "usenorm:", &fancynorm);  /* changed 11/02/06 */
  getint (ph, "dotpopsmode:", &dotpopsmode);
  getint (ph, "pcorrmode:", &pcorrmode);        /* print correlations */
  getint (ph, "pcpopsonly:", &pcpopsonly);      /* but only within population */
  getint (ph, "readparsonly:", &readparsonly);      /* treminate after reading par file */
  getint (ph, "altnormstyle:", &altnormstyle);
  getint (ph, "hashcheck:", &hashcheck);
  getint (ph, "popgenmode:", &altnormstyle);
  getint (ph, "noxdata:", &noxdata);
  getint (ph, "inbreed:", &inbreed);
  getint (ph, "easymode:", &easymode);
  getint (ph, "printcover:", &printcover);
  getint (ph, "seed:", &t);
  seed = (long) t;

  getint (ph, "fastmode:", &fastmode);
  getint (ph, "fastdim:", &fastdim);
  getint (ph, "fastiter:", &fastiter);

  getint (ph, "usepopsformissing:", &usepopsformissing);
  getint (ph, "regmode:", &regmode);
  getint (ph, "lsqproject:", &regmode);
  getint (ph, "doproject:", &doproject);
  getint (ph, "rounakmode:", &rounakmode);
  getint (ph, "fstverbose:", &fstverbose);


  t = -1;
  getint (ph, "xdata:", &t);
  if (t >= 0)
    noxdata = 1 - t;
  getint (ph, "nostatslim:", &nostatslim);
  getint (ph, "popsizelimit:", &popsizelimit);
  getint (ph, "minallelecnt:", &minallelecnt);
  getint (ph, "chrom:", &xchrom);
  getint (ph, "nochrom:", &zchrom);
  getint (ph, "maxmissing:", &maxmissing);
  getint (ph, "lopos:", &lopos);
  getint (ph, "hipos:", &hipos);
  getint (ph, "checksizemode:", &checksizemode);
  getint (ph, "pubmean:", &pubmean);
  getint (ph, "fstonly:", &fstonly);
  getint (ph, "fsthiprecision:", &fsthiprec);
  getint (ph, "fstjack:", &fstjack);

  getint (ph, "hiprecision:", &hiprec);
  getint (ph, "hiprec:", &hiprec);
  getint (ph, "fstz:", &fstz);

  getint (ph, "ldregress:", &ldregress);
  getint (ph, "nsnpldregress:", &ldregress);    /* changed 11/02/06 */
  getdbl (ph, "ldlimit:", &ldlimit);    /* in morgans */
  getint (ph, "ldposlimit:", &ldposlimit);      /* bases */
  getdbl (ph, "ldr2lo:", &ldr2lo);
  getdbl (ph, "ldr2hi:", &ldr2hi);
  getdbl (ph, "maxdistldregress:", &ldlimit);   /* in morgans *//* changed 11/02/06 */
  getint (ph, "minleneig:", &nostatslim);
  getint (ph, "malexhet:", &malexhet);
  getint (ph, "nomalexhet:", &malexhet);        /* changed 11/02/06 */
  getint (ph, "familynames:", &familynames);
  getint (ph, "qtmode:", &qtmode);

  getint (ph, "numoutliter:", &numoutliter);
  getint (ph, "numoutlieriter:", &numoutliter); /* changed 11/02/06 */
  getint (ph, "numoutleigs", &numoutleigs);
  getint (ph, "numoutlierevec:", &numoutleigs); /* changed 11/02/06 */
  getdbl (ph, "outlthresh:", &outlthresh);
  getdbl (ph, "outliersigmathresh:", &outlthresh);      /* changed 11/02/06 */
  getint (ph, "outliermode:", &outliermode);    /* test distribution with sample removed. Makes sense for small samples */
  getdbl (ph, "blgsize:", &blgsize);

  getstring (ph, "indoutfilename:", &indoutfilename);
  getstring (ph, "indivoutname:", &indoutfilename);     /* changed 11/02/06 */
  getstring (ph, "snpoutfilename:", &snpoutfilename);
  getstring (ph, "snpoutname:", &snpoutfilename);       /* changed 11/02/06 */
  getstring (ph, "genooutfilename:", &genooutfilename);
  getstring (ph, "genotypeoutname:", &genooutfilename); /* changed 11/02/06 */
  getstring (ph, "outputformat:", &omode);
  getstring (ph, "outputmode:", &omode);
  getint (ph, "outputgroup:", &ogmode);
  getstring (ph, "grmoutname:", &grmoutname);
  getint (ph, "grmbinary:", &grmbinary);
  getint (ph, "packout:", &packout);    /* now obsolete 11/02/06 */
  getstring (ph, "twxtabname:", &twxtabname);
  getstring (ph, "id2pops:", &id2pops);
  getstring (ph, "elllistname:", &elllistname) ; 

  getdbl (ph, "r2thresh:", &r2thresh);
  getdbl (ph, "r2genlim:", &r2genlim);
  getdbl (ph, "r2physlim:", &r2physlim);
  getint (ph, "killr2:", &killr2);

  getint (ph, "numchrom:", &numchrom);
  getstring (ph, "xregionname:", &xregionname);
  getdbl (ph, "hwfilter:", &nhwfilter);

  getint (ph, "numthreads:", &thread_ct_config);

  printf ("### THE INPUT PARAMETERS\n");
  printf ("##PARAMETER NAME: VALUE\n");
  writepars (ph);

}

int
fvadjust (double *cc, int n, double *pmean, double *fancy)

/* take off mean  force missing to zero */

/* set up fancy norming  */
{
  double p, ynum, ysum, y, ymean, yfancy = 1.0;
  int i, nmiss = 0;

  ynum = ysum = 0.0;
  for (i = 0; i < n; i++) {
    y = cc[i];
    if (y < 0.0) {
      ++nmiss;
      continue;
    }
    ++ynum;
    ysum += y;
  }
  if (ynum == 0.0) {
    return -999;
  }
  ymean = ysum / ynum;
  for (i = 0; i < n; i++) {
    y = cc[i];
    if (y < 0.0)
      cc[i] = 0.0;
    else
      cc[i] -= ymean;
  }
  if (pmean != NULL)
    *pmean = ymean;
  if (fancynorm) {
    p = 0.5 * ymean;            // autosomes
    if (altnormstyle == NO)
      p = (ysum + 1.0) / (2.0 * ynum + 2.0);
    y = p * (1.0 - p);
    if (y > 0.0)
      yfancy = 1.0 / sqrt (y);
  }
  if (fancy != NULL)
    *fancy = yfancy;
  return nmiss;
}

int
fvadjust_binary (int c0, int c1, int nmiss, int n, double *cc, double *pmean,
                 double *fancy)
{
  double p, ynum, ysum, y, ymean, yfancy = 1.0;

  if (n == nmiss) {
    return -999;
  }
  ynum = n - nmiss;
  ysum = c0;
  ymean = ysum / ynum;
  cc[0] = -ymean;
  cc[1] = 1.0 - ymean;
  cc[2] = 2.0 - ymean;
  if (fancynorm) {
    p = 0.5 * ymean;
    if (altnormstyle == NO) {
      p = (ysum + 1.0) / (2.0 * ynum + 2.0);
    }
    y = p * (1.0 - p);
    if (y > 0.0) {
      yfancy = 1.0 / sqrt (y);
    }
  }
  if (pmean) {
    *pmean = ymean;
  }
  if (fancy) {
    *fancy = yfancy;
  }
  return nmiss;
}

double
dottest (char *sss, double *vec, char **eglist, int numeg, int *xtypes,
         int len)
// vec will always have mean 0 
// perhaps should rewrite to put xa1 etc in arrays
{
  double *w1;
  int *xt;
  int i, k1, k2, k, n, x1, x2;
  double ylike;
  double ychi;
  double *wmean;
  int imax, imin, *isort;
  static int ncall = 0;

  char ss1[MAXSTR];
  char ss2[MAXSTR];
  double ans, ftail, ftailx, ansx;

  ZALLOC (wmean, numeg, double);
  ZALLOC (w1, len + numeg, double);
  ZALLOC (isort, numeg, int);
  ZALLOC (xt, len, int);
  strcpy (ss1, "");

  calcmean (wmean, vec, len, xtypes, numeg);
  if (pubmean) {
    copyarr (wmean, w1, numeg);
    sortit (w1, isort, numeg);
    printf ("%s:means\n", sss);
    for (i = 0; i < numeg; i++) {
      k = isort[i];
      printf ("%20s ", eglist[k]);
      printf (" %9.3f\n", wmean[k]);
    }
  }

  vlmaxmin (wmean, numeg, &imax, &imin);
  if (chisqmode) {
    ylike = anova1 (vec, len, xtypes, numeg);
    ans = 2.0 * ylike;
  }
  else {
    ans = ftail = anova (vec, len, xtypes, numeg);
  }
  ++ncall;


  if (numeg > 2) {
    sprintf (ss2, "%s %s ", sss, "overall");
    publishit (ss2, numeg - 1, ans);
    printf (" %20s minv: %9.3f %20s maxv: %9.3f\n",
            eglist[imin], wmean[imin], eglist[imax], wmean[imax]);
  }


  for (k1 = 0; k1 < numeg; ++k1) {
    for (k2 = k1 + 1; k2 < numeg; ++k2) {
      n = 0;
      x1 = x2 = 0;
      for (i = 0; i < len; i++) {
        k = xtypes[i];
        if (k == k1) {
          w1[n] = vec[i];
          xt[n] = 0;
          ++n;
          ++x1;
        }
        if (k == k2) {
          w1[n] = vec[i];
          xt[n] = 1;
          ++n;
          ++x2;
        }
      }

      if (x1 <= 1)
        continue;
      if (x2 <= 1)
        continue;

      ylike = anova1 (w1, n, xt, 2);
      ychi = 2.0 * ylike;
      chitot[k1 * numeg + k2] += ychi;
      if (chisqmode) {
        ansx = ychi;
      }
      else {
        ansx = ftailx = anova (w1, n, xt, 2);
      }

      sprintf (ss2, "%s %s %s ", sss, eglist[k1], eglist[k2]);
      publishit (ss2, 1, ansx);

    }
  }
  free (w1);
  free (xt);
  free (wmean);
  free (isort);
  return ans;
}

double
anova (double *vec, int len, int *xtypes, int numeg)
// anova 1 but f statistic
{
  int i, k;
  double y1, top, bot, ftail;
  double *w0, *w1, *popsize, *wmean;

  static int ncall2 = 0;

  if (numeg >= len) {
    printf ("*** warning: bad anova popsizes too small\n");
    return 0.0;
  }

  ZALLOC (w0, len, double);
  ZALLOC (w1, len, double);
  ZALLOC (wmean, numeg, double);
  ZALLOC (popsize, numeg, double);

  y1 = asum (vec, len) / (double) len;  // mean
  vsp (w0, vec, -y1, len);

  for (i = 0; i < len; i++) {
    k = xtypes[i];
    ++popsize[k];
    wmean[k] += w0[i];
  }

/* debug */
  if (numeg == 2) {
    ++ncall2;
    for (i = 0; i < len; ++i) {
      if (ncall2 < 0)
        break;
      k = xtypes[i];
//    printf("yy %4d %4d %12.6f %12.6f\n", i, k, vec[i], w0[i]) ;
    }
  }

  vsp (popsize, popsize, 1.0e-12, numeg);
  vvd (wmean, wmean, popsize, numeg);

  vvt (w1, wmean, wmean, numeg);
  top = vdot (w1, popsize, numeg);

  for (i = 0; i < len; i++) {
    k = xtypes[i];
    w1[i] = w0[i] - wmean[k];
  }
  bot = asum2 (w1, len) / (double) (len - numeg);
  bot *= (double) (numeg - 1);
  ftail = rtlf (numeg - 1, len - numeg, top / bot);

  free (w0);
  free (w1);
  free (popsize);
  free (wmean);

  return ftail;

}

double
anova1 (double *vec, int len, int *xtypes, int numeg)
{
  int i, k;
  double y1, y2, ylike;
  double *w0, *w1, *popsize, *wmean;

  ZALLOC (w0, len, double);
  ZALLOC (w1, len, double);
  ZALLOC (wmean, numeg, double);
  ZALLOC (popsize, numeg, double);

  y1 = asum (vec, len) / (double) len;  // mean
  vsp (w0, vec, -y1, len);

  for (i = 0; i < len; i++) {
    k = xtypes[i];
    ++popsize[k];
    wmean[k] += w0[i];
  }

  vsp (popsize, popsize, 1.0e-12, numeg);
  vvd (wmean, wmean, popsize, numeg);

  for (i = 0; i < len; i++) {
    k = xtypes[i];
    w1[i] = w0[i] - wmean[k];
  }

  y1 = asum2 (w0, len) / (double) len;
  y2 = asum2 (w1, len) / (double) len;
  ylike = 0.5 * ((double) len) * log (y1 / y2);

  free (w0);
  free (w1);
  free (popsize);
  free (wmean);

  return ylike;

}

void
publishit (char *sss, int df, double chi)
{
  double tail;
  char sshit[4];
  char ss2[MAXSTR];
  int i, n;
  char cblank, cunder;
  static int ncall = 0;

  ++ncall;
  cblank = ' ';
  cunder = '_';
  n = strlen (sss);

  strcpy (ss2, sss);
  for (i = 0; i < n; ++i) {
    if (ss2[i] == cblank)
      ss2[i] = cunder;
  }

  if (chisqmode) {
    if (ncall == 1)
      printf
        ("## Anova statistics for population differences along each eigenvector:\n");
    if (ncall == 1)
      printf ("%40s %6s %9s %12s\n", "", "dof", "chisq", "p-value");
    printf ("%40s %6d %9.3f", ss2, df, chi);
    tail = rtlchsq (df, chi);
    printf (" %12.6g", tail);
  }
  else {
    if (ncall == 1)
      printf
        ("## Anova statistics for population differences along each eigenvector:\n");
    if (ncall == 1)
      printf ("%40s %12s\n", "", "p-value");
    printf ("%40s ", ss2);
    tail = chi;
    printf (" %12.6g", tail);
  }
  strcpy (sshit, "");
  if (tail < pvhit)
    strcpy (sshit, "***");
  if (tail < pvjack)
    strcpy (sshit, "+++");
  printf (" %s", sshit);
  printf ("\n");
}

void
dotpops (double *X, char **eglist, int numeg, int *xtypes, int nrows)
{
  double *pp, *npp, val, yy;
  int *popsize;
  int i, j, k1, k2;


  if (dotpopsmode == NO) return ;
  if (fstonly)
    return;
  ZALLOC (pp, numeg * numeg, double);
  ZALLOC (npp, numeg * numeg, double);
  popsize = xpopsize;

  ivzero (popsize, numeg);

  for (i = 0; i < nrows; i++) {
    k1 = xtypes[i];
    ++popsize[k1];
    for (j = i + 1; j < nrows; j++) {
      k2 = xtypes[j];
      if (k1 < 0)
        fatalx ("bug\n");
      if (k2 < 0)
        fatalx ("bug\n");
      if (k1 >= numeg)
        fatalx ("bug\n");
      if (k2 >= numeg)
        fatalx ("bug\n");
      val = X[i * nrows + i] + X[j * nrows + j] - 2.0 * X[i * nrows + j];
      pp[k1 * numeg + k2] += val;
      pp[k2 * numeg + k1] += val;
      ++npp[k1 * numeg + k2];
      ++npp[k2 * numeg + k1];
    }
  }
  vsp (npp, npp, 1.0e-8, numeg * numeg);
  vvd (pp, pp, npp, numeg * numeg);
// and normalize so that mean on diagonal is 1 
  yy = trace (pp, numeg) / (double) numeg;
  vst (pp, pp, 1.0 / yy, numeg * numeg);
  printf ("\n## Average divergence between populations:");
  if (numeg <= 10) {
    printf ("\n");
    printf ("%10s", "");
    for (k1 = 0; k1 < numeg; ++k1) {
      printf (" %10s", eglist[k1]);
    }
    printf ("  %10s", "popsize");
    printf ("\n");
    for (k2 = 0; k2 < numeg; ++k2) {
      printf ("%10s", eglist[k2]);
      for (k1 = 0; k1 < numeg; ++k1) {
        val = pp[k1 * numeg + k2];
        printf (" %10.3f", val);
      };
      printf ("  %10d", popsize[k2]);
      printf ("\n");
    }
  }
  else {                        // numeg >= 10 
    printf ("\n");
    for (k2 = 0; k2 < numeg; ++k2) {
      for (k1 = k2; k1 < numeg; ++k1) {
        printf ("dotp: %10s", eglist[k2]);
        printf (" %10s", eglist[k1]);
        val = pp[k1 * numeg + k2];
        printf (" %10.3f", val);
        printf ("    %10d", popsize[k2]);
        printf (" %10d", popsize[k1]);
        printf ("\n");
      }
    }
  }
  printf ("\n");
  printf ("\n");
  fflush (stdout);


  free (pp);
  free (npp);

}

void
printxcorr (double *X, int nrows, Indiv ** indxx)
{
  int k1, k2, t;
  double y1, y2, yy, rho;
  Indiv *ind1, *ind2;

  if (pcorrmode == NO)
    return;
  for (k1 = 0; k1 < nrows; ++k1) {
    for (k2 = k1 + 1; k2 < nrows; ++k2) {

      ind1 = indxx[k1];
      ind2 = indxx[k2];

      t = strcmp (ind1->egroup, ind2->egroup);
      if (pcpopsonly && (t != 0))
        continue;


      y1 = X[k1 * nrows + k1];
      y2 = X[k2 * nrows + k2];
      yy = X[k1 * nrows + k2];

      rho = yy / sqrt (y1 * y2 + 1.0e-20);
      printf ("corr: %20s %20s %20s %20s %9.3f\n",
              ind1->ID, ind2->ID, ind1->egroup, ind2->egroup, rho);

    }
  }
}

void
bumpldvv (double *gsource, double *newsource, int *pnumld, int maxld, int n,
          int *ldsnpbuff, int newsnpnum)
{

  int numld;
  SNP *cuptnew, *cuptold;
  int pdiff;
  double gdiff;


  numld = *pnumld;

  cuptnew = snpmarkers[newsnpnum];

  for (;;) {
    if (numld == 0)
      break;
    cuptold = snpmarkers[ldsnpbuff[0]];
    pdiff = nnint (cuptnew->physpos - cuptold->physpos);
    gdiff = cuptnew->genpos - cuptold->genpos;
    if ((pdiff <= ldposlimit) && (gdiff <= ldlimit))
      break;
    copyarr (gsource + n, gsource, (maxld - 1) * n);    // overlapping move but copyarr works left to right 
    copyiarr (ldsnpbuff + 1, ldsnpbuff, (maxld - 1));   // overlapping move but copyiarr works left to right 
    --numld;
  }

  if (numld < maxld) {
    copyarr (newsource, gsource + numld * n, n);
    ldsnpbuff[numld] = newsnpnum;
    ++numld;
    *pnumld = numld;
    return;
  }

  if (maxld == numld) {
    copyarr (gsource + n, gsource, (maxld - 1) * n);    // overlapping move but copyarr works left to right 
    copyiarr (ldsnpbuff + 1, ldsnpbuff, (maxld - 1));   // overlapping move but copyiarr works left to right 
    --numld;
  }
  copyarr (newsource, gsource + numld * n, n);
  ldsnpbuff[numld] = newsnpnum;
  ++numld;

  *pnumld = numld;
  return;
}

int
ldregx (double *gsource, double *gtarget, double *res, int rsize,
        int n, double r2lo, double r2hi)
{

/** 
 gsource: array of (normalized) genotypes 
 rsize rows n long.   
 So row 1 is gsource[0]..gsource[n-1] 
 row 2 gsource[n]...gsource[2*n-1] 
 gtarget n long normalized genotype 
 Routine should return residual (n long) 
  
 return code 
 a) 0  Did nothing 
 b) 1  Ran regression 
 c) 2  Residual set 0 
*/

  if (rsize == 0) {
    copyarr (gtarget, res, n);
    return 0;
  }

  // Allocate space for all genotypes to pass
  double *gsource_pass;
  ZALLOC (gsource_pass, rsize * n, double);

  int i, ii;

  // Compute correlation to previous SNPs
  double sum;
  int rsize_pass = 0;
  for (i = 0; i < rsize; i++) {
    sum = 0;
    for (ii = 0; ii < n; ii++) {
      sum += gtarget[ii] * gsource[i * n + ii];
    }
    // Normalize by (n-1) and square to get cor^2
    sum = pow (sum / (2 * (n - 1)), 2);
    // Check if correlation too high
    if (sum > r2hi) {
      // Clean up and exit
      free (gsource_pass);

      // Residual set to all zero
      for (ii = 0; ii < n; ii++)
        res[ii] = 0;
      return 2;
      // Check if correlation not too low
    }
    else if (sum > r2lo) {
      // Retain this SNP for the regression
      for (ii = 0; ii < n; ii++)
        gsource_pass[rsize_pass * n + ii] = gsource[i * n + ii];
      rsize_pass++;
    }
  }

  // Do the regression if correlated SNPs were found
  if (rsize_pass > 0) {
    double *t_gsource_pass, *regans, *www;
    ZALLOC (regans, rsize, double);
    ZALLOC (www, n, double);
    ZALLOC (t_gsource_pass, rsize * n, double);


    // BUG FIX  old call in EIG5 was wrong:
    transpose (t_gsource_pass, gsource_pass, rsize_pass, n);

    regressit (regans, t_gsource_pass, gtarget, n, rsize_pass); //run regression
    mulmat (www, regans, gsource_pass, 1, rsize_pass, n);       //multiply regans and gsource_pass

    vvm (res, gtarget, www, n);


    free (regans);
    free (www);
    free (t_gsource_pass);
    free (gsource_pass);
    return 1;
  }
  else {
    copyarr (gtarget, res, n);
    free (gsource_pass);
    return 0;
  }
}


void
dofstxx (double *xfst, double *fstans, double *fstsd, int *fstnum, SNP ** xsnplist, int *xindex,
         int *xtypes, int nrows, int ncols, int numeg, double blgsize,
         SNP ** snpmarkers, Indiv ** indm)
{


  int nblocks, xnblocks;
  int *blstart, *blsize;

  if (qtmode == YES) {
    return;
  }

  nblocks = numblocks (snpmarkers, numsnps, blgsize);
  printf ("number of blocks for block jackknife: %d\n", nblocks);

/**
  if (nblocks <= 1) {
    return;
  }
*/

  ZALLOC (blstart, nblocks, int);
  ZALLOC (blsize, nblocks, int);


  setblocks (blstart, blsize, &xnblocks, xsnplist, ncols, blgsize);
  fixwt (xsnplist, ncols, 1.0);

  dofstnumx (xfst, fstans, fstsd, fstnum, xsnplist, xindex, xtypes,
             nrows, ncols, numeg, nblocks, indm, YES);

  free (blstart);
  free (blsize);
 
}

void
fixwt (SNP ** snpm, int nsnp, double val)
{
  int k;
  SNP *cupt;

  for (k = 0; k < nsnp; ++k) {
    cupt = snpm[k];
    cupt->weight = val;
  }

}

double
oldfstcol (double *estn, double *estd, SNP * cupt,
           int *xindex, int *xtypes, int nrows, int type1, int type2)
{
  int c1[2], c2[2], *cc;
  int *rawcol;
  int k, g, i;
  double ya, yb, yaa, ybb, p1, p2, en, ed;
  double z, zz, h1, h2, yt;
  static int ncall = 0;


  ++ncall;
  ZALLOC (rawcol, nrows, int);

  getrawcol (rawcol, cupt, xindex, nrows);

  ivzero (c1, 2);
  ivzero (c2, 2);

  for (i = 0; i < nrows; i++) {
    k = xtypes[i];
    cc = NULL;
    if (k == type1)
      cc = c1;
    if (k == type2)
      cc = c2;
    if (cc == NULL)
      continue;
    g = rawcol[i];
    if (g < 0)
      continue;
    cc[0] += g;
    cc[1] += 2 - g;
  }

  ya = c1[0];
  yb = c1[1];
  yaa = c2[0];
  ybb = c2[1];
  z = ya + yb;
  zz = yaa + ybb;
  if ((z < 0.1) || (zz < 0.1)) {
    *estn = 0.0;
    *estd = -1.0;
    free (rawcol);
    return 0.0;
  }

  yt = ya + yb;
  p1 = ya / yt;
  h1 = ya * yb / (yt * (yt - 1.0));

  yt = yaa + ybb;
  p2 = yaa / yt;
  h2 = yaa * ybb / (yt * (yt - 1.0));

  en = (p1 - p2) * (p1 - p2);
  en -= h1 / z;
  en -= h2 / zz;

  ed = en;
  ed += h1;
  ed += h2;

  *estn = en;
  *estd = ed;


  free (rawcol);
  return z + zz;

}


double
fstcol (double *estn, double *estd, SNP * cupt,
        int *xindex, int *xtypes, int nrows, int type1, int type2)
{
  int c1[2], c2[2], *cc;
  int *rawcol;
  int k, g, i;
  double ya, yb, yaa, ybb, p1, p2, en, ed;
  double z, zz, h1, h2, yt;
  int **ccc;
  static int ncall = 0;


  ++ncall;
  ccc = initarray_2Dint (nrows, 2, 0);
  ZALLOC (rawcol, nrows, int);

  getrawcolx (ccc, cupt, xindex, nrows, indivmarkers);
  getrawcol (rawcol, cupt, xindex, nrows);

  ivzero (c1, 2);
  ivzero (c2, 2);

  for (i = 0; i < nrows; i++) {
    k = xtypes[i];
    cc = NULL;
    if (k == type1)
      cc = c1;
    if (k == type2)
      cc = c2;
    if (cc == NULL)
      continue;
    g = ccc[i][0];
    if (ncall < 1000) {
//    printf("zz %d  %d %d\n", rawcol[i], ccc[i][0], ccc[i][1]) ;
    }

    if (g < 0)
      continue;
    ivvp (cc, cc, ccc[i], 2);
  }

  ya = c1[0];
  yb = c1[1];
  yaa = c2[0];
  ybb = c2[1];
  z = ya + yb;
  zz = yaa + ybb;
  if ((z < 1.1) || (zz < 1.1)) {
    *estn = 0.0;
    *estd = -1.0;
    free (rawcol);
    free2Dint (&ccc, nrows);
    return 0.0;
  }

  yt = ya + yb;
  p1 = ya / yt;
  h1 = ya * yb / (yt * (yt - 1.0));

  yt = yaa + ybb;
  p2 = yaa / yt;
  h2 = yaa * ybb / (yt * (yt - 1.0));

  en = (p1 - p2) * (p1 - p2);
  en -= h1 / z;
  en -= h2 / zz;

  ed = en;
  ed += h1;
  ed += h2;

  *estn = en;
  *estd = ed;


  free (rawcol);
  free2Dint (&ccc, nrows);
  return z + zz;

}

void
writesnpeigs (char *snpeigname, SNP ** xsnplist, double *ffvecs, int numeigs,
              int ncols)
{
// this is called at end and ffvecs overwritten
  double *xpt, y, yscal, *snpsc;
  int i, j, k, kmax, kmin;
  SNP *cupt;
  FILE *fff;

  for (j = 0; j < numeigs; ++j) {
    xpt = ffvecs + j * ncols;
    y = asum2 (xpt, ncols);
    yscal = (double) ncols / y;
    yscal = sqrt (yscal);
    vst (xpt, xpt, yscal, ncols);
  }


  ZALLOC (snpsc, ncols, double);
  vclear (snpsc, -99999, ncols);
  for (j = 0; j < numeigs; ++j) {
    for (i = 0; i < ncols; ++i) {
      cupt = xsnplist[i];
      if (cupt->ignore)
        continue;
      y = ffvecs[j * ncols + i];
      snpsc[i] = fabs (y);
    }
    for (k = 0; k < 10; ++k) {
      if (ncols <= 10)
        break;
// was <= 10 Tiny bug
      vlmaxmin (snpsc, ncols, &kmax, &kmin);
      cupt = xsnplist[kmax];
      if (snpsc[kmax] < 0)
        break;
      printf ("eigbestsnp %4d %20s %2d %12d %9.3f\n", j + 1, cupt->ID,
              cupt->chrom, nnint (cupt->physpos), snpsc[kmax]);
      snpsc[kmax] = -1.0;
    }
  }
  free (snpsc);


  if (snpeigname == NULL)
    return;
  openit (snpeigname, &fff, "w");

  for (i = 0; i < ncols; ++i) {
    cupt = xsnplist[i];
    if (cupt->ignore)
      continue;

    fprintf (fff, "%20s", cupt->ID);
    fprintf (fff, " %2d", cupt->chrom);
    fprintf (fff, " %12d", nnint (cupt->physpos));

    for (j = 0; j < numeigs; ++j) {
      if (hiprec) fprintf (fff, " %12.6f", ffvecs[j * ncols + i]);
      else fprintf (fff, " %9.3f", ffvecs[j * ncols + i]);
    }
    fprintf (fff, "\n");
  }

  fclose (fff);

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
getcolxz (double *xcol, SNP * cupt, int *xindex, int *xtypes, int nrows,
          int col, double *xmean, double *xfancy, int *n0, int *n1)
// side effect set xmean xfancy and count variant and reference alleles
// returns missings after fill in
{
  int j, n, g, t, k, kmax = -1;
  double y, pmean, yfancy;
  int *rawcol;
  int c0, c1, nmiss;
  double *popnum = NULL;
  double *popsum = NULL;

  if (usepopsformissing) {
    ZALLOC (popnum, maxpops + 1, double);
    ZALLOC (popsum, maxpops + 1, double);
  }

  c0 = c1 = 0;
  ZALLOC (rawcol, nrows, int);
  n = cupt->ngtypes;
  if (n < nrows)
    fatalx ("bad snp: %s %d\n", cupt->ID, n);
  getrawcol (rawcol, cupt, xindex, nrows);
  nmiss = 0;
  for (j = 0; j < nrows; ++j) {
    g = rawcol[j];
    if (g < 0) {
      ++nmiss;
      continue;
    }
    c0 += g;
    c1 += 2 - g;
    if (usepopsformissing) {
      k = xtypes[j];
      if (k>=0) { 
       popsum[k] += (double) g;
       popnum[k] += 1.0;
       kmax = MAX (kmax, k);
      } 
    }
  }
  floatit (xcol, rawcol, nrows);
  if ((usepopsformissing) && (nmiss > 0)) {
    pmean = asum (popsum, kmax + 1) / asum (popnum, kmax + 1);
    nmiss = 0;
    for (j = 0; j < nrows; ++j) {
      g = rawcol[j];
      if (g >= 0)
        continue;
      k = xtypes[j];
      if (popnum[k] > 0.5) {
        y = popsum[k] / popnum[k];
        xcol[j] = y;
        continue;
      }
      ++nmiss;
    }
  }
  t = fvadjust (xcol, nrows, &pmean, &yfancy);
  if (t < -99) {
    if (xmean != NULL) {
      xmean[col] = 0.0;
      xfancy[col] = 0.0;
    }
    vzero (xcol, nrows);
    free (rawcol);
    if (n0 != NULL) {
      *n0 = -1;
      *n1 = -1;
    }
    return -1;
  }
  vst (xcol, xcol, yfancy, nrows);
  if (xmean != NULL) {
    xmean[col] = pmean * yfancy;
    xfancy[col] = yfancy;
  }
  free (rawcol);
  if (n0 != NULL) {
    *n0 = c0;
    *n1 = c1;
  }
  if (usepopsformissing) {
    free (popnum);
    free (popsum);
  }
  return nmiss;
}

int
getcolxz_binary1 (int *rawcol, double *xcol, SNP * cupt, int *xindex,
                  int nrows, int col, double *xmean, double *xfancy, int *n0,
                  int *n1)
{
  // Modified getcolxz() which converts to a 3-bit-per-genotype representation
  // compatible with PLINK 1.5's partial sum lookup outer product algorithm.
  // (Well, to be more precise, the conversion occurs in getcolxz_binary2();
  // this function handles the other duties of getcolxz().)  Assumes
  // usepopsformissing is NOT set, and ldregress is zero.
  //
  // Main genotype array:
  //   Homozygous minor -> 0
  //   Heterozygous     -> 2
  //   Homozygous major -> 3
  //   Missing          -> 0
  //
  // Missing mask:
  //   Nonmissing       -> 0
  //   Missing          -> 7
  //
  // Suppose person 1 has genotype g_1 and missing mask m_1, and person 2 has
  // genotype g_2 and missing mask m_2.  Then, the operation
  //
  //   (g_1 + g_2) | m_1 | m_2
  //
  // executes the following mapping:
  //
  //   Both genotypes hom minor -> 0
  //   Hom minor + het          -> 2
  //   Hom minor + hom major    -> 3
  //   Het + het                -> 4
  //   Het + hom major          -> 5
  //   Hom major + hom major    -> 6
  //   Either genotype missing  -> 7
  //
  // Construction of the corresponding lookup table is deferred to
  // domult_increment_lookup().

  int j, n, g, t;
  double pmean, yfancy;
  int c0, c1, nmiss;

  c0 = c1 = 0;
  n = cupt->ngtypes;
  if (n < nrows) {
    fatalx ("bad snp: %s %d\n", cupt->ID, n);
  }
  getrawcol (rawcol, cupt, xindex, nrows);
  nmiss = 0;
  for (j = 0; j < nrows; ++j) {
    g = rawcol[j];
    if (g < 0) {
      ++nmiss;
      continue;
    }
    c0 += g;
    c1 += 2 - g;
  }
  // instead of storing an entire column of floating point values,
  t = fvadjust_binary (c0, c1, nmiss, nrows, xcol, &pmean, &yfancy);
  if (t < -99) {
    if (xmean != NULL) {
      xmean[col] = 0.0;
      xfancy[col] = 0.0;
    }
    vzero (xcol, 3);
    if (n0 != NULL) {
      *n0 = -1;
      *n1 = -1;
    }
    return -1;
  }
  vst (xcol, xcol, yfancy, 3);
  if (xmean != NULL) {
    xmean[col] = pmean * yfancy;
    xfancy[col] = yfancy;
  }
  if (n0 != NULL) {
    *n0 = c0;
    *n1 = c1;
  }
  return nmiss;
}

void
getcolxz_binary2 (int *rawcol, uintptr_t * binary_cols,
                  uintptr_t * binary_mmask, uint32_t xblock, uint32_t nrows)
{
  // slightly better to position at 0-3-6-9-12-16-19... instead of
  // 0-3-6-9-12-15-18...
  uint32_t shift_val = (xblock * 3) + (xblock / 5);

  uintptr_t bitfield_or[3];
  uint32_t row_idx;
  int cur_geno;
  bitfield_or[0] = ((uintptr_t) 7) << shift_val;
  bitfield_or[1] = ((uintptr_t) 2) << shift_val;
  bitfield_or[2] = ((uintptr_t) 3) << shift_val;
  for (row_idx = 0; row_idx < nrows; row_idx++) {
    cur_geno = *rawcol++;
    if (cur_geno) {
      if (cur_geno > 0) {
        binary_cols[row_idx] |= bitfield_or[(uint32_t) cur_geno];
      }
      else {
        binary_mmask[row_idx] |= bitfield_or[0];
      }
    }
  }
}

void
join_threads (pthread_t * threads, uint32_t ctp1)
{
  if (!(--ctp1)) {
    return;
  }
  uint32_t uii;
  for (uii = 0; uii < ctp1; uii++) {
    pthread_join (threads[uii], NULL);
  }
}

int32_t
spawn_threads (pthread_t * threads, void *(*start_routine) (void *),
               uintptr_t ct)
{
  uintptr_t ulii;
  if (ct == 1) {
    return 0;
  }
  for (ulii = 1; ulii < ct; ulii++) {
    if (pthread_create
        (&(threads[ulii - 1]), NULL, start_routine, (void *) ulii)) {
      join_threads (threads, ulii);
      return -1;
    }
  }
  return 0;
}

THREAD_RET_TYPE
block_increment_binary (void *arg)
{
  uintptr_t tidx = (uintptr_t) arg;
  uintptr_t cur_indiv_idx = g_thread_start[tidx];
  uintptr_t end_indiv_idx = g_thread_start[tidx + 1];
  uintptr_t *binary_cols = g_binary_cols;
  uintptr_t *binary_mmask = g_binary_mmask;
  double *write_ptr =
    &(g_XTX_lower_tri[(cur_indiv_idx * (cur_indiv_idx + 1)) / 2]);
  double *weights0 = g_weights;
  double *weights1 = &(g_weights[32768]);
#ifdef __LP64__
  double *weights2 = &(g_weights[65536]);
  double *weights3 = &(g_weights[98304]);
#endif
  uintptr_t *geno_ptr;
  uintptr_t *mmask_ptr;
  uintptr_t base_geno;
  uintptr_t base_mmask;
  uintptr_t final_geno;
  uintptr_t indiv_idx2;
  for (; cur_indiv_idx < end_indiv_idx; cur_indiv_idx++) {
    geno_ptr = binary_cols;
    base_geno = binary_cols[cur_indiv_idx];
    mmask_ptr = binary_mmask;
    base_mmask = binary_mmask[cur_indiv_idx];
    if (!base_mmask) {
      // special case: current individual has no missing genotypes in block
      for (indiv_idx2 = 0; indiv_idx2 <= cur_indiv_idx; indiv_idx2++) {
        final_geno = ((*geno_ptr++) + base_geno) | (*mmask_ptr++);
#ifdef __LP64__
        *write_ptr +=
          weights0[(uint16_t) final_geno] +
          weights1[(uint16_t) (final_geno >> 16)] +
          weights2[(uint16_t) (final_geno >> 32)] +
          weights3[final_geno >> 48];
#else
        *write_ptr +=
          weights0[(uint16_t) final_geno] + weights1[final_geno >> 16];
#endif
        write_ptr++;
      }
    }
    else {
      for (indiv_idx2 = 0; indiv_idx2 <= cur_indiv_idx; indiv_idx2++) {
        final_geno =
          ((*geno_ptr++) + base_geno) | ((*mmask_ptr++) | base_mmask);
#ifdef __LP64__
        *write_ptr +=
          weights0[(uint16_t) final_geno] +
          weights1[(uint16_t) (final_geno >> 16)] +
          weights2[(uint16_t) (final_geno >> 32)] +
          weights3[final_geno >> 48];
#else
        *write_ptr +=
          weights0[(uint16_t) final_geno] + weights1[final_geno >> 16];
#endif
        write_ptr++;
      }
    }
  }
  THREAD_RETURN;
}

void
domult_increment_lookup (pthread_t * threads, uint32_t thread_ct,
                         double *XTX_lower_tri, double *tblock,
                         uintptr_t * binary_cols, uintptr_t * binary_mmask,
                         uint32_t block_size, uint32_t indiv_ct,
                         double *partial_sum_lookup_buf)
{
  // PLINK 1.5 partial sum lookup algorithm
  double increments[40];
  double *dptr;
  double *dptr2;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t uoo;
  double partial_incr1;
  double partial_incr2;
  double partial_incr3;
  double partial_incr4;
  uintptr_t ulii;

  // populate lookup buffer
#ifdef __LP64__
  for (uii = 0; uii < 20; uii += 5)
#else
  for (uii = 0; uii < 10; uii += 5)
#endif
  {
    dptr = increments;
    for (ujj = 0; ujj < 5; ujj++) {
      dptr2 = &(tblock[(uii + ujj) * 3]);
      *dptr++ = dptr2[0] * dptr2[0];
      *dptr++ = 0;
      *dptr++ = dptr2[0] * dptr2[1];
      *dptr++ = dptr2[0] * dptr2[2];
      *dptr++ = dptr2[1] * dptr2[1];
      *dptr++ = dptr2[1] * dptr2[2];
      *dptr++ = dptr2[2] * dptr2[2];
      *dptr++ = 0;
    }
    dptr = &(partial_sum_lookup_buf[(uii / 5) * 32768]);
    for (ujj = 0; ujj < 8; ujj++) {
      partial_incr1 = increments[ujj + 32];
      for (ukk = 0; ukk < 8; ukk++) {
        partial_incr2 = partial_incr1 + increments[ukk + 24];
        for (umm = 0; umm < 8; umm++) {
          partial_incr3 = partial_incr2 + increments[umm + 16];
          for (unn = 0; unn < 8; unn++) {
            partial_incr4 = partial_incr3 + increments[unn + 8];
            for (uoo = 0; uoo < 8; uoo++) {
              *dptr++ = partial_incr4 + increments[uoo];
            }
          }
        }
      }
    }
  }
  g_XTX_lower_tri = XTX_lower_tri;
  g_weights = partial_sum_lookup_buf;
  g_binary_cols = binary_cols;
  g_binary_mmask = binary_mmask;
  if (spawn_threads (threads, block_increment_binary, thread_ct)) {
    fatalx ("Error: Failed to create thread.\n");
    return;
  }
  ulii = 0;
  block_increment_binary ((void *) ulii);
  join_threads (threads, thread_ct);
}

THREAD_RET_TYPE
block_increment_normal (void *arg)
{
  uintptr_t tidx = (uintptr_t) arg;
  uintptr_t start_indiv_idx = g_thread_start[tidx];
  uintptr_t end_indiv_idx = g_thread_start[tidx + 1];
  uintptr_t indiv_ct = g_indiv_ct;
  uint32_t block_size = g_block_size;
  double *write_start_ptr =
    &(g_XTX_lower_tri[(start_indiv_idx * (start_indiv_idx + 1)) / 2]);
  double *write_ptr;
  double *tblock;
  double *tblock_read_ptr;
  double cur_tblock_val;
  uintptr_t cur_indiv_idx;
  uintptr_t indiv_idx2;
  uint32_t bidx;
  for (bidx = 0; bidx < block_size; bidx++) {
    write_ptr = write_start_ptr;
    tblock = &(g_tblock[bidx * indiv_ct]);
    for (cur_indiv_idx = start_indiv_idx; cur_indiv_idx < end_indiv_idx;
         cur_indiv_idx++) {
      cur_tblock_val = tblock[cur_indiv_idx];
      tblock_read_ptr = tblock;
      for (indiv_idx2 = 0; indiv_idx2 <= cur_indiv_idx; indiv_idx2++) {
        *write_ptr += cur_tblock_val * (*tblock_read_ptr++);
        write_ptr++;
      }
    }
  }
  THREAD_RETURN;
}

void
domult_increment_normal (pthread_t * threads, uint32_t thread_ct,
                         double *XTX_lower_tri, double *tblock,
                         int block_size, uint32_t indiv_ct)
{
  // General case: tblock[] can have an arbitrary number of distinct values, so
  // can't use bit hacks.
  //
  // This multithreaded implementation is pretty far from optimal; if more
  // speed is needed, use the DGEMM function from a vendor-optimized BLAS.
  // (Sum of k outer products is just equal to the product of a n*k and a k*n
  // matrix.)
  int ii;
  double ycheck;
  uintptr_t ulii;
  for (ii = 0; ii < block_size; ii++) {
    ycheck = asum (tblock + ii * indiv_ct, indiv_ct);
    if (fabs (ycheck) > .00001)
      fatalx ("bad ycheck\n");
  }
  g_XTX_lower_tri = XTX_lower_tri;
  g_tblock = tblock;
  g_block_size = block_size;
  g_indiv_ct = indiv_ct;
  if (spawn_threads (threads, block_increment_normal, thread_ct)) {
    fatalx ("Error: Failed to create thread.\n");
    return;
  }
  ulii = 0;
  block_increment_normal ((void *) ulii);
  join_threads (threads, thread_ct);
}

void
getcolxf (double *xcol, SNP * cupt, int *xindex, int nrows, int col,
          double *xmean, double *xfancy)
// side effect set xmean xfancy
{
  int n;
  double pmean, yfancy;
  int *rawcol;

  if (xmean != NULL) {
    xmean[col] = xfancy[col] = 0.0;
  }

  if (cupt->ignore) {
    vzero (xcol, nrows);
    return;
  }

  ZALLOC (rawcol, nrows, int);
  n = cupt->ngtypes;
  if (n < nrows)
    fatalx ("bad snp: %s %d\n", cupt->ID, n);
  getrawcol (rawcol, cupt, xindex, nrows);
  floatit (xcol, rawcol, nrows);

  fvadjust (xcol, nrows, &pmean, &yfancy);
  vst (xcol, xcol, yfancy, nrows);
  if (xmean != NULL) {
    xmean[col] = pmean * yfancy;
    xfancy[col] = yfancy;
  }
  free (rawcol);
}

void
doinbxx (double *inbans, double *inbsd, SNP ** xsnplist, int *xindex,
         int *xtypes, int nrows, int ncols, int numeg, double blgsize,
         SNP ** snpmarkers, Indiv ** indm)
{


  int nblocks, xnblocks;
  int *blstart, *blsize;

  double *xinb;

  nblocks = numblocks (snpmarkers, numsnps, blgsize);

  ZALLOC (blstart, nblocks, int);
  ZALLOC (blsize, nblocks, int);
  ZALLOC (xinb, numeg, double);


  setblocks (blstart, blsize, &xnblocks, xsnplist, ncols, blgsize);
  fixwt (xsnplist, ncols, 1.0);

  doinbreed (xinb, inbans, inbsd, xsnplist, xindex, xtypes,
             nrows, ncols, numeg, nblocks, indm);

  free (blstart);
  free (blsize);
  free (xinb);

}


void
calcpopmean (double *wmean, char **elist, double *vec,
             char **eglist, int numeg, int *xtypes, int len)
// extracted from dotttest ;
{
  double *w0, *w1;
  int *isort;
  int i, k;

  ZALLOC (w0, len, double);
  ZALLOC (w1, len, double);
  ZALLOC (isort, len, int);


  calcmean (w0, vec, len, xtypes, numeg);

  copyarr (w0, w1, numeg);
  sortit (w1, isort, numeg);

  for (i = 0; i < numeg; i++) {
    k = isort[i];
    elist[i] = eglist[k];
    wmean[i] = w0[k];
  }



  free (w0);
  free (w1);
  free (isort);


}

void
sqz (double *azq, double *acoeffs, int numeigs, int nrows, int *xindex)
{

  int i, j, k;
  // Indiv *indx ;
  static int ncall = 0;

  ++ncall;

  for (k = 0; k < nrows; ++k) {
    i = xindex[k];
    if (i < 0)
      fatalx ("zzyuk!\n");
    // indx = indivmarkers[i] ; 
//  if (ncall == 1) printf("zz %3d %12s %12s %d %d\n", k, indx -> ID, indx -> egroup, indx -> ignore, indx -> affstatus) ;

    for (j = 0; j < numeigs; ++j) {
      azq[j * nrows + k] = acoeffs[j * numindivs + i];
    }
  }
}

void
dumpgrmid (char *fname, Indiv ** indivmarkers, int *xindex, int numid)
{
  FILE *fff;
  int a, b;
  Indiv *indx;

  openit (fname, &fff, "w");
  for (a = 0; a < numid; ++a) {
    b = xindex[a];
    if ((b < 0) || (b >= numindivs))
      fatalx ("(dumpgrmid) bad index\n");
    indx = indivmarkers[b];
    fprintf (fff, "%s\t%s\n", "NA", indx->ID);
  }
  fclose (fff);
}

void
dumpgrmbin (double *XTX, int nrows, int numsnps, Indiv ** indivmarkers,
            int numindivs, char *grmoutname)
{
  int a, b;
  double y;
  char sss[256];
  char *bb;
  int wout, numout, fdes, ret = 0;
  float yfloat;

  if (sizeof (yfloat) != 4)
    fatalx ("grm binary only supported for 4 byte floats\n");

  sprintf (sss, "%s.N.bin", grmoutname);
  ridfile (sss);
  fdes = open (sss, O_CREAT | O_TRUNC | O_RDWR, 0666);

  if (fdes < 0) {
    perror ("bad dumpgrmbin");
    fatalx ("open failed for %s\n", sss);
  }
  if (verbose)
    printf ("file %s opened\n", sss);

//  numout = numsnps*(numsnps+1)/4 ;
  numout = nrows * (nrows + 1) / 2;
  wout = numsnps;
  bb = (char *) &wout;

  for (a = 0; a < numout; ++a) {
    ret = write (fdes, bb, 4);
  }
  if (ret < 0) {
    perror ("write failure");
    fatalx ("(outpack) bad write");
  }
  close (fdes);

  sprintf (sss, "%s.bin", grmoutname);
  ridfile (sss);
  fdes = open (sss, O_CREAT | O_TRUNC | O_RDWR, 0666);

  if (fdes < 0) {
    perror ("bad dumpgrmbin");
    fatalx ("open failed for %s\n", sss);
  }
  if (verbose)
    printf ("file %s opened\n", sss);

  // Re-adjust values based on diagonal normalization
  double y_norm;
  y_norm = trace (XTX, nrows) / (double) nrows;

  bb = (char *) &yfloat;
  for (a = 0; a < nrows; a++) {
    for (b = 0; b <= a; b++) {
      y = XTX[a * nrows + b] / y_norm;  // bugfix
      yfloat = (float) y;
      ret = write (fdes, bb, 4);
    }
  }
  close (fdes);
}

void
dumpgrm (double *XTX, int *xindex, int nrows, int numsnps,
         Indiv ** indivmarkers, int numindivs, char *grmoutname)
{
  int a, b;
  double y;
  FILE *fff;
  char sss[256];

  if (grmoutname == NULL)
    return;

  sprintf (sss, "%s.id", grmoutname);
  dumpgrmid (sss, indivmarkers, xindex, nrows);

  if (grmbinary) {
    dumpgrmbin (XTX, nrows, numsnps, indivmarkers, numindivs, grmoutname);
    return;
  }

  // Re-adjust values based on diagonal normalization
  double y_norm_recip;
  double *d;
  ZALLOC (d, nrows, double);
  getdiag (d, XTX, nrows);
  y_norm_recip = ((double) nrows) / asum (d, nrows);
  free (d);

  openit (grmoutname, &fff, "w");
  for (a = 0; a < nrows; a++) {
    for (b = 0; b <= a; b++) {
      y = XTX[a * nrows + b];   // bugfix: do NOT want to dereference xindex here
      fprintf (fff, "%d %d ", a + 1, b + 1);
      fprintf (fff, "%d ", numsnps);
      fprintf (fff, "%0.6f\n", y * y_norm_recip);
    }
  }
  fclose (fff);

}

void
printevecs (SNP ** snpmarkers, Indiv ** indivmarkers, Indiv ** xindlist,
            int numindivs, int ncols, int nrows,
            int numeigs, double *eigenvecs, double *eigenvals, FILE * ofile)
{

  double *ffvecs, *fvecs, *cc, *xrow, *bcoeffs, y;
  double *fxscal, *xpt, val;
  int i, j, k, a, b ;
  Indiv *indx;
  static int ncall = 0 ; 

/**
  char *isample = "I3122" ; 
  int qindex ;  

  qindex = indindex(indivmarkers, numindivs, isample) ; 
*/

  ++ncall;  
  if (ncall>1) { 
   printf("repeated calls to printevecs ... all but first ignored\n") ;
   return  ; 
  }



  fprintf (ofile, "%20s ", "#eigvals:");
  for (j = 0; j < numeigs; j++) {
    fprintf (ofile, "%9.3f ", eigenvals[j]);
  }
  fprintf (ofile, "\n");

  if ((easymode) || (nrows == numindivs)) {

// should be separate routine
    printf("writing eigenvecotrs for  %d samples\n", nrows) ; 

    ZALLOC (fvecs, nrows * numeigs, double);
    setfvecs (fvecs, eigenvecs, nrows, numeigs);

    for (j = 0; j < numeigs; j++) {
      xpt = fvecs + j * nrows;
      y = asum2 (xpt, nrows);
      vst (xpt, xpt, 1.0 / sqrt (y), nrows);    // norm 1
    }
    for (i = 0; i < nrows; i++) {
      indx = xindlist[i];
      fprintf (ofile, "%20s ", indx->ID);
      for (j = 0; j < numeigs; j++) {
        xpt = fvecs + j * nrows;
        y = xpt[i];
        if (indx -> flag == 7777) y *= edgarw[j] ;
        if (hiprec) fprintf (ofile, "%12.6f  ", y);
        else fprintf (ofile, "%10.4f  ", y);  // printevecs 
        eigenvecs[j*numindivs + i] = y ; 
      }
      fprintf (ofile, "%15s\n", indx->egroup);
      fflush(ofile) ;
    }
    free (fvecs);
   return;
  }

  printf("zznoteasy\n") ;
  fatalx("... not yet implemented!\n") ;
  ZALLOC (ffvecs, ncols * numeigs, double);
  ZALLOC (fvecs, nrows * numeigs, double);
  ZALLOC (cc, nrows, double);
  ZALLOC (xrow, ncols, double);
  ZALLOC (bcoeffs, numeigs * numindivs, double);
  ZALLOC (fxscal, numeigs, double);

  setfvecs (fvecs, eigenvecs, nrows, numeigs);

  for (i = 0; i < ncols; i++) {
    for (j = 0; j < numeigs; j++) {
      for (k = 0; k < nrows; k++) {
        getgval (k, i, &val);
        ffvecs[j * ncols + i] += fvecs[j * nrows + k] * val;
      }
    }
  }

 
  for (i = 0; i < nrows; i++) {

    for (k = 0; k < ncols; ++k) {
      getgval (i, k, &val);
      xrow[k] = val;
    }

    for (j = 0; j < numeigs; j++) {
      xpt = ffvecs + j * ncols;
      y = vdot (xrow, xpt, ncols);
      fxscal[j] += y * y;
    }
  }

  vsqrt (fxscal, fxscal, numeigs);
  vinvert (fxscal, fxscal, numeigs);

  for (i = 0; i < numindivs; i++) {
    indx = indivmarkers[i];
    if (indx->ignore)
      continue;
    for (k = 0; k < ncols; ++k) {
      getggval (i, k, &val);
      xrow[k] = val;
    }

    for (j = 0; j < numeigs; j++) {
      bcoeffs[j * numindivs + i] = y =
        fxscal[j] * vdot (xrow, ffvecs + j * ncols, ncols);
    }
  }

  for (i = 0; i < numindivs; i++) {
    indx = indivmarkers[i];
    if (indx->ignore)
      continue;
    fprintf (ofile, "%20s ", indx->ID);
    for (j = 0; j < numeigs; j++) {
      y = bcoeffs[j * numindivs + i];
      eigenvecs[j * numindivs + i] = y ; 
      if (hiprec) fprintf (ofile, "%12.6f  ", y);
      else fprintf (ofile, "%10.4f  ", y);
    }
    fprintf (ofile, "%15s\n", indx->egroup);
    fflush(ofile) ; 
  }

  writesnpeigs (snpeigname, snpmarkers, ffvecs, numeigs, ncols);


  free (fvecs);
  free (ffvecs);
  free (cc);
  free (xrow);
  free (bcoeffs);
  free (fxscal);
}

void
norme (double *ee, int n)
{
// ee returned as mean 0 , square norm 1 
  double y;
  y = asum (ee, n) / (double) n;
  vsp (ee, ee, -y, n);
  y = asum2 (ee, n);
  vst (ee, ee, 1.0 / sqrt (y), n);

/**
  y = asum(ee, n) ; 
  if (isnan(y)) fatalx("(norme) isnan!\n") ;
*/
}

void
cpr (double *a, double *b, int n)
{
  double y;
  return ; 

// debug
  printf ("zzcpr\n");
  y = vdot (a, a, n);
  printf ("a a:: %12.3f\n", y);
  y = vdot (b, b, n);
  printf ("b b:: %12.3f\n", y);
  y = vdot (a, b, n);
  printf ("a b:: %12.3f\n", y);
}

void
doproj (double *regans, int isample, SNP ** xsnplist, double *fxvecs,
        int numeigs, int ncols)
// work (numeigs^2 + ncols) * numeigs  --- small
{

  double *xrow, *trow, *rhs, *emat;
  int k, kk, j;

  ZALLOC (xrow, ncols, double);
  ZALLOC (trow, ncols, double);
  ZALLOC (rhs, ncols, double);
  ZALLOC (emat, ncols * numeigs, double);

  loadxdataind (xrow, xsnplist, isample, ncols);
  copyarr (xrow, trow, ncols);
  fixxrow (xrow, xmean, xfancy, ncols);
  kk = 0;
  for (k = 0; k < ncols; ++k) {
    if (trow[k] < 0)
      continue;
    rhs[kk] = xrow[k];
    for (j = 0; j < numeigs; j++) {
      emat[kk * numeigs + j] = fxvecs[j * ncols + k];
    }
    ++kk;
  }
  regressit (regans, emat, rhs, kk, numeigs);

  free (xrow);
  free (trow);
  free (rhs);
  free (emat);

}

void
doshrinkp2 (double *mmat, int m, int n, int *xindex, SNP ** xsnplist, double *xcoeffs)
// lsq project mode
// shrink multiple evecs at once
{
  double *mmatt, *xmat, y, yscale;
  double *evecs, *lambda, *ffvecs, *fxvecs;
  double *dv1, *dd, *d, delta, *evec, *ww, *w2, *ediff, *fxv;
  double eco, lam;
  int i, j, k, l, a, b, t;
  char cc;
  double *sold, *snew, *enew, *ss, *ymul, *xpt;
  int debug;
  double *oldffv, *trow, *ffvec;
  int nrows = m ; 

  printf ("doshrink (newshrink) called\n");
  fflush (stdout);
  cputimes(0, 2) ; 

  ZALLOC (xmat, m * m, double);
  ZALLOC (evecs, m * m, double);
  ZALLOC (lambda, m * m, double);
  ZALLOC (dd, m * m, double);
  ZALLOC (dv1, m, double);
  ZALLOC (d, m, double);
  ZALLOC (ww, m, double);
  ZALLOC (w2, m, double);
  ZALLOC (ymul, numeigs, double);
  ZALLOC (ediff, m, double);
  ZALLOC (enew, m, double);
  ZALLOC (oldffv, n, double);
  ZALLOC (ffvec, n, double);
  ZALLOC (ffvecs, numeigs * n, double);
  ZALLOC (fxvecs, numeigs * n, double);

  ZALLOC (sold, numeigs * m, double);
  ZALLOC (snew, numeigs * m, double);
  ss = xcoeffs ; 

  if ( (XTX != NULL)  && (m == nrxtx)) { 
    xmat = XTX ; 
  }
// can we copy XTX

 else {
  printf(" (doshrinkp2) *** recomputing xmat\n") ;  fflush(stdout) ;
  for (i = 0; i < m; i++) {
    for (j = 0; j < m; j++) {
      for (k = 0; k < n; ++k) {
        xmat[i * m + j] += mmat[i * n + k] * mmat[j * n + k];
      }
    }
  }
 }

  yscale = y = trace (xmat, m) / (double) (m - 1);
  y = fabs(yscale-1.0) ; 
  if (y>0.01) printf ("trace:  %9.3f\n", yscale);

  fflush(stdout) ; 
  if (yscale <= 0.0)
    fatalx ("mmat has zero trace (perhaps no data)\n");
       
   vst (xmat, xmat, 1.0 / yscale, m * m);  // force mean eigenvalue to be 1 

  eigvecs (xmat, lambda, evecs, m);


 if (toprightindex >=0) {  
  for (j=0; j< numeigs; ++j) { 
   if (flip == NULL) break ; 
   if (flip[j] != 1) continue ; 
    xpt = evecs + j*nrows ; 
    vst(xpt, xpt, -1, nrows) ; 
    printf("(shrinkp2) eigenvector %d flipped!\n", j+1) ;
   }
  }
  
  for (i = 0; i < numeigs; i++) {
    evec = evecs + i * m;
    norme (evec, m);
    mulmat (ffvec, evec, mmat, 1, m, n);  // work m*n*numeigs  (not too bad)  ffvecs is numeigs*n  
    y = asum2 (ffvec, n) / (double) n;
    vst (ffvecs + i * n, ffvec, 1.0 / sqrt (y), n);
  }

  for (i = 0; i < numindivs; i++) {
    doproj (ww, i, xsnplist, ffvecs, numeigs, n);
    for (j = 0; j < numeigs; ++j) {
      ss[j * numindivs + i] = ww[j];
    }
  }

// old style projection onto eigenvectors
// normalized eigenvector // don't need to do this for indivs in nrows  
   printf(" old style projection complete\n") ;

  copyarr (ffvecs, fxvecs, numeigs * n);
  for (i = 0; i < numeigs; i++) {

    evec = evecs + i * m;
    lam = lambda[i];

    mulmat (ww, evec, xmat, 1, m, m);
    vst (ww, ww, 1.0 / lam, m);
    norme (ww, m);
    copyarr (ww, sold + i * m, m);   // redundant but uniform with snew 
  }

    for (a = 0; a < m; a++) {
     cputimes(0, 3) ; 
     for (i=0; i< numeigs; ++i) {

      evec = evecs + i*m ; 
      lam = lambda[i];

      vst (dv1, xmat + a * m, -1, m);
      vzero (dd, m * m);
      for (k = 0; k < m; k++) {
        dd[a * m + k] = dd[k * m + a] = dv1[k];
      }
      mulmat (ww, dd, evec, m, m, 1);
      delta = vdot (ww, evec, m);
      vzero (ediff, m);
      for (k = 0; k < m; k++) {
        if (k == i)
          continue;
        eco = vdot (evecs + k * m, ww, m);
        eco /= (lam - lambda[k]);
        vst (w2, evecs + k * m, eco, m);
        vvp (ediff, ediff, w2, m);
      }
      if (lam > -delta) {
        y = ymul[i] = lam / (lam + delta);
      }

      else {
        ymul[i] = 1.0;
        printf ("*** bad delta  -- negative eigenvalue!\n");
        printf ("a: %d  i: %d\n", a, i);
        printf ("lam, delta: %9.3f %9.3f\n", lam, delta);
      }

      vvp (enew, evec, ediff, m);
      enew[a] = 0;
      y = asum (enew, m) / (double) (m - 1);
      vsp (enew, enew, -y, m);
      enew[a] = 0;
      norme (enew, m);


      mulmat (ffvec, enew, mmat, 1, m, n);

      y = asum2 (ffvec, n) / (double) n;
      vst (ffvec, ffvec, 1.0 / sqrt (y), n);
      fxv = fxvecs + i*n ; 
      copyarr (ffvec, fxv, n);
    }
      doproj (ww, xindex[a], xsnplist, fxvecs, numeigs, n);
      vvt(ww, ww, ymul, numeigs) ;
      for (i=0; i<numeigs; i++) {
       snew[i * m + a] = ww[i] ;
      }
     y = cputimes(1, 3) ; 
     printf("sample %d shrunk!  cpu time: %9.3f\n", a, y) ;  fflush(stdout) ; 
   }

  for (j = 0; j < numeigs; j++) {
    for (t=0; t<m; ++t) { 
     i = xindex[t] ; 
     ss[j * numindivs + i] = snew[j * m + t];
     continue;
    }
  }

  if (verbose) printf("zzevecs: p2\n") ;
  printevecs (xsnplist, indivmarkers, indivmarkers,
              numindivs, n, numindivs, numeigs, ss, lambda, ofile);

  free (xmat);
  free (evecs);
  free (lambda);
  free (dd);
  free (dv1);
  free (d);
  free (ww);
  free (w2);
  free (ymul);
  free (sold);
  free (snew);
  free (enew);
  free (oldffv) ;
  free (ffvec);
  free (ffvecs);
  free (fxvecs);
  y = cputimes(1, 2) ;
  printf ("doshrinkp2 exited :: cpu time : %9.3f\n", y);
  fflush (stdout);
}

void
doshrinkp (double *mmat, int m, int n, int *xindex, SNP ** xsnplist, double *xcoeffs)
// lsq project mode
{
  double *mmatt, *xmat, y, yscale;
  double *evecs, *lambda, *ffvecs, *fxvecs;
  double *dv1, *dd, *d, delta, *evec, *ww, *w2, *ediff, *fxv;
  double eco, lam;
  int i, j, k, l, a, b, t;
  char cc;
  double *sold, *snew, *enew, *ss, *xpt;
  int debug;
  double *oldffv, *trow, *ffvec;
  double ymul;
  int nrows = m ; 

/**
  char *isample = "I3122" ; 
  int qindex ;  

  qindex = indindex(indivmarkers, numindivs, isample) ; 
*/

  if (shrinkmode == NO)
    return;

  if (newshrink) { 
   doshrinkp2(mmat, m, n, xindex, xsnplist, xcoeffs) ;
   return ; 
  }
  printf ("doshrink called\n");
  fflush (stdout);
  cputimes(0, 2) ; 

  ZALLOC (xmat, m * m, double);
  ZALLOC (evecs, m * m, double);
  ZALLOC (lambda, m * m, double);
  ZALLOC (dd, m * m, double);
  ZALLOC (dv1, m, double);
  ZALLOC (d, m, double);
  ZALLOC (ww, m, double);
  ZALLOC (w2, m, double);
  ZALLOC (ediff, m, double);
  ZALLOC (enew, m, double);
  ZALLOC (oldffv, n, double);
  ZALLOC (ffvec, n, double);
  ZALLOC (ffvecs, numeigs * n, double);
  ZALLOC (fxvecs, numeigs * n, double);

  ZALLOC (sold, numeigs * m, double);
  ZALLOC (snew, numeigs * m, double);
  ss = xcoeffs ; 

  if ( (XTX != NULL)  && (m == nrxtx)) { 
    xmat = XTX ; 
  }
// can we copy XTX

 else {
  printf(" (doshrinkp) *** recomputing xmat\n") ;
  for (i = 0; i < m; i++) {
    for (j = 0; j < m; j++) {
      for (k = 0; k < n; ++k) {
        xmat[i * m + j] += mmat[i * n + k] * mmat[j * n + k];
      }
    }
  }
 }

  yscale = y = trace (xmat, m) / (double) (m - 1);
  printf ("trace:  %9.3f\n", y);
  fflush(stdout) ; 
  if (y <= 0.0)
    fatalx ("mmat has zero trace (perhaps no data)\n");
   vst (xmat, xmat, 1.0 / y, m * m);  // force mean eigenvalue to be 1 

  eigvecs (xmat, lambda, evecs, m);


 if (toprightindex >=0) {  
  for (j=0; j< numeigs; ++j) { 
   if (flip == NULL) break ; 
   if (flip[j] != 1) continue ; 
    xpt = evecs + j*nrows ; 
    vst(xpt, xpt, -1, nrows) ; 
    printf("(shrinkp) eigenvector %d flipped!\n", j+1) ;
   }
  }
  

  for (i = 0; i < numeigs; i++) {
    evec = evecs + i * m;
    norme (evec, m);

    mulmat (ffvec, evec, mmat, 1, m, n);  // work m*n*numeigs  (not too bad)  ffvecs is numeigs*n  
    y = asum2 (ffvec, n) / (double) n;
    vst (ffvecs + i * n, ffvec, 1.0 / sqrt (y), n);
  }

  for (i = 0; i < numindivs; i++) {
    doproj (ww, i, xsnplist, ffvecs, numeigs, n);
    for (j = 0; j < numeigs; ++j) {
      ss[j * numindivs + i] = ww[j];
    }
  }

  copyarr (ffvecs, fxvecs, numeigs * n);
  for (i = 0; i < numeigs; i++) {
    evec = evecs + i * m;
    lam = lambda[i];
    mulmat (ww, evec, xmat, 1, m, m);
    vst (ww, ww, 1.0 / lam, m);
    norme (ww, m);
    copyarr (ww, sold + i * m, m);   // redundant but uniform with snew 

    for (a = 0; a < m; a++) {
      vst (dv1, xmat + a * m, -1, m);
      vzero (dd, m * m);
      for (k = 0; k < m; k++) {
        dd[a * m + k] = dd[k * m + a] = dv1[k];
      }
      mulmat (ww, dd, evec, m, m, 1);
      delta = vdot (ww, evec, m);
      vzero (ediff, m);
      for (k = 0; k < m; k++) {
        if (k == i)
          continue;
        eco = vdot (evecs + k * m, ww, m);
        eco /= (lam - lambda[k]);
        vst (w2, evecs + k * m, eco, m);
        vvp (ediff, ediff, w2, m);
      }
      if (lam > delta) {
        ymul = lam / (lam + delta);
      }

      else {
        ymul = 1.0;
        printf ("*** bad delta  -- negative eigenvalue!\n");
        printf ("a: %d  i: %d\n", a, i);
        printf ("lam, delta: %9.3f %9.3f\n", lam, delta);
      }

      vvp (enew, evec, ediff, m);
      enew[a] = 0;
      y = asum (enew, m) / (double) (m - 1);
      vsp (enew, enew, -y, m);
      enew[a] = 0;
      norme (enew, m);

      mulmat (ffvec, enew, mmat, 1, m, n);
      y = asum2 (ffvec, n) / (double) n;
      vst (ffvec, ffvec, 1.0 / sqrt (y), n);
      fxv = fxvecs + i*n ; 
      copyarr(fxv, oldffv, n) ;
      copyarr (ffvec, fxv, n);
      doproj (ww, xindex[a], xsnplist, fxvecs, numeigs, n);
      snew[i * m + a] = ww[i] * ymul;
      copyarr(oldffv, fxv, n) ;

// why don't we make a new fxvecs and call doproj just once?  

    }
  }

  for (j = 0; j < numeigs; j++) {
    for (i = 0; i < numindivs; i++) {
      t = findfirst (xindex, m, i);
      if (t >= 0) {
        ss[j * numindivs + i] = snew[j * m + t];
        continue;
      }
    }
  }

  if (verbose) printf("zzevecs p0\n") ;
  printevecs (xsnplist, indivmarkers, indivmarkers,
              numindivs, n, numindivs, numeigs, ss, lambda, ofile);

  free (xmat);
  free (evecs);
  free (lambda);
  free (dd);
  free (dv1);
  free (d);
  free (ww);
  free (w2);
  free (sold);
  free (snew);
  free (enew);
  free (oldffv) ;
  free (ffvec);
  free (ffvecs);
  free (fxvecs);
  y = cputimes(1, 2) ;
  printf ("doshrink exited :: cpu time : %9.3f\n", y);
  fflush (stdout);
}

int setnstw(double *lambda, int len, int nostatslim) 
// number of "significant" eigenvectors 
{
  int i ; 
  double tw, zn, zvar, y ;

  printmat(lambda, 1, len) ;  
  for (i=0; i<len-nostatslim; ++i) { 
   zn = znval ;
   y = dotwcalc(lambda+i, len-i, &tw, &zn, &zvar, nostatslim) ;
   if (y < -0.5) return -1 ;
   if (y > autothresh) return i ;  
  }


  return -1 ;



}
void estrounak(double *edgarw, double *lambdav, int lentop, int lenspec, double gamm, double yjfac)
{

 double *top  ; 
 double *w1, *lam, lambda ; 
 double mphat, vphat, dphat, y, y1, y2, y0, ytemp ;
 int i, j, len, jnum ;

 if (lentop <= 0) return ;
 len = lentop + lenspec ;
 printf("zzz %d %d %9.3f\n", lentop, len, gamm) ; 
 ZALLOC(w1, len, double) ; 
 ZALLOC(lam, len, double) ; 
 copyarr(lambdav, lam, len) ;     
 y = bal1(lam, len) ;    
// scales lam to have sum 1    
//  printf("rounak: sum: %9.3f\n", y)  ; 
/** 
Rounak's R code
temp=temp+samp.eval[j]/(samp.eval[i]-samp.eval[j])
spike.est[i]<-samp.eval[i]/(1+(gamma/(p-m))*temp)
*/ 
 top = lam ;  
  for (i=0; i<lentop; ++i) {  
   ytemp = 0 ; 
   lambda = top[i] ;  
   jnum = 0 ; 
   for (j=lentop; j<len; ++j) { 
     y0  = top[j] / (lambda - top[j]) ; 
//   printf("zzmul: %4d %12.6f %12.6f %12.6f\n", j, lambda, top[j], y0) ; 
     ytemp += y0 ;
     ++jnum ; 
   } 
   printf("zzrounak: %d %12.6f %12.6f %12.6f\n", i, top[i] * (double) len, ytemp, yjfac) ;  
// ytemp is scale invariant 
   y1 = gamm/(double) yjfac ;  
   y1 *= ytemp ;   
   edgarw[i] = y2 = 1.0/(1.0+y1) ; 
   printf("ytemp y1 y2: %9.3f %9.3f %9.3f\n", ytemp, y1, y2) ; 
// y2 should be shrinkage factor
  }

  free(w1) ;
  free(lam) ;

} 

void estedgar(double *edgarw, double *lambdav, int lentop, int lenspec, double gamm, double yjfac)
{

 double *top  ; 
 double *spec ; 
 double *w1, *lam, lambda, mhat, vhat, dhat, ell ;
 double mphat, vphat, dphat, y ;
 int i, len ;

 if (rounakmode) { 
  estrounak(edgarw, lambdav, lentop, lenspec, gamm, yjfac) ; 
  return ; 


 }
 if (lentop <= 0) return ;
 len = lentop + lenspec ;
 ZALLOC(w1, len, double) ; 
 ZALLOC(lam, len, double) ; 
 copyarr(lambdav, lam, len) ;     
 y = bal1(lam, len) ; 
//  printf("edgar: sum: %9.3f\n", y)  ; 
 top = lam ;  
 spec = top + lentop ; 
  for (i=0; i<lentop; ++i) {  
   lambda = top[i] ;  
   vsp(w1, spec, -lambda, lenspec) ; 
   vinvert(w1, w1, lenspec) ; 
   mhat = asum(w1, lenspec)/ (double) (lenspec) ; 
   vhat = gamm*mhat - (1.0-gamm)/lambda ; 
   dhat = lambda*mhat*vhat ; 
   ell = 1.0/dhat ;

   vvt(w1, w1, w1, lenspec) ; 
   mphat = asum(w1, lenspec)/ (double) lenspec ;
   vphat = gamm*mphat + (1.0-gamm)/(lambda*lambda) ; 
   dphat = mhat*vhat  + lambda*mphat*vhat + lambda*mhat*vphat ;   ; 
   edgarw[i] = sqrt(mhat/(dphat*ell)) ; 
  }

  free(w1) ;
  free(lam) ;

} 
int mpestimate(double *eigs, int neigs, double *peffect, double *psigma)
{
// eigs sorted in descending order 
 double qlo, qhi, qmed, qrat, xmed ; 
 double y, yn ; 
 
 int x, t ; 

 static int ntab = 0 ;

 static double **mptable = NULL ; 
 double *gval, *rat, *med ; 

 *peffect = *psigma = 0 ; 

 if (neigs < 10) return -1 ;
// y interpolate here 
 qhi  = quartile(eigs, neigs, 0.75) ;
 qmed = quartile(eigs, neigs, 0.50) ;
 qlo  = quartile(eigs, neigs, 0.25) ;

 qrat = (qhi-qlo)/qmed ; 
 printf("mpstats: ") ; 
 printf("%9.3f", qlo) ;
 printf("%9.3f", qmed) ;
 printf("%9.3f", qhi) ;
 printf("%12.6f", qrat) ;
 printnl() ; 

 if (ntab == 0) ntab = loadmptable(&mptable) ;

 gval = mptable[0] ; 
 rat  = mptable[4] ;
 med =  mptable[2] ;  
 t = qinterp(rat, gval, ntab, qrat, &y) ; 
 if (t<=0) y = 1.0e-6 ; 
 *peffect =  (double) neigs / y ;   

 t = qinterp(gval, med, ntab, y, &xmed) ;
 *psigma = xmed/qmed ; 
// printf("zzmp :: %12.6f %12.6f\n", y, xmed) ; 

 return 1 ; 

}

void
countsn(int *blcnt, int nblocks, SNP **snplist, int nsnps, int ind)  
{
  SNP *cupt ; 
  int i, t, g ; 
  static int ncall = 0 ;

  ++ncall ;

  ivzero(blcnt, nblocks) ;
  for (i=0; i<nsnps; ++i) { 
   cupt = snplist[i] ; 
   if (ncall == -99) printf("snp: %s %d %d %d\n", cupt -> ID, cupt -> chrom, cupt -> tagnumber, cupt -> ignore) ;
   if (cupt -> ignore) continue ; 
   t = cupt -> tagnumber ; 
/**
   if ((i==0) || (i == (nsnps-1))) {
    printf("snp: %s %d %d\n", cupt -> ID, cupt -> chrom, cupt -> tagnumber) ;
   }
*/
   if (t<0) continue ; 
   if (t>=nblocks) fatalx("countsn overflow %d %d\n", t, nblocks) ;
   g = getgtypes(cupt, ind) ; 
   if (g>=0) ++blcnt[t] ; 
  }
//if (ncall == 1) printimat(blcnt, 1, nblocks) ; 
}
 
void lsqproj(int blocknum, SNP **xsnplist, int ncols, Indiv **xindlist, int nind, 
  double *fxscal, double *ffvecs, double * acoeffs, double *bcoeffs,int *xtypes, int numeg) 

{

 int qk, i, k, kk, j ; 
 double *azq, *bzq ; 
 double y ;
 static double *xrow, *trow, *rhs, *emat, *regans ; 
 Indiv *indx ; 
 SNP *cupt ; 
 int t, u, x  ; 
 double *xnum, *xden ; 
 
 static int ncall = 0 ;

 fflush(stdout) ;
 if (!regmode) return ; 

  ++ncall ; 

    if (ncall==1) {
      ZALLOC (trow, ncols, double);
      ZALLOC (xrow, ncols, double);
      ZALLOC (rhs, ncols, double);
      ZALLOC (emat, ncols * numeigs, double);
      ZALLOC (regans, numeigs, double);
      printf("lsqproj called!\n") ;

     for (qk = 0; qk < nind; qk++) {
      indx = xindlist[qk];
      i = indx -> idnum ; 
      if (indx->ignore) continue;
      t = xtypes[i] ; 
      if (t<0) continue ; 
      x = numvalids(indx, xsnplist, 0, ncols-1) ;
      if (printcover) printf("coverage: %8s %20s %7d\n", indx -> ID, indx -> egroup, x) ; 

      if (x <= numeigs) {
        indx->ignore = YES;
        printf ("%s ignored (insufficient data\n", indx->ID);
        continue;
      }
     }
    }

    vzero(trow, ncols) ;
    vzero(xrow, ncols) ;
    vzero(rhs, ncols) ;
    vzero(emat, ncols*numeigs) ;
    vzero(regans, numeigs) ;

    if (usepopsformissing) { 
      ZALLOC(xnum, numeg*ncols, double) ; 
      ZALLOC(xden, numeg*ncols, double) ; 

     for (qk = 0; qk < nind; qk++) {
      indx = xindlist[qk];
      i = indx -> idnum ; 
      if (indx->ignore) continue;
      loadxdataind (xrow, xsnplist, i, ncols);
      t = xtypes[i] ; 
      if (t<0) continue ;
      for (j=0; j<ncols; ++j) { 
       if (xrow[j]>=0) { 
         xnum[t*ncols+j] += xrow[j] ; 
         xden[t*ncols+j] += 1 ; 
       }
      }
    }
   }

    for (qk = 0; qk < nind; qk++) {
      indx = xindlist[qk];
      i = indx -> idnum ; 
      if (indx->ignore) continue;
      loadxdataind (xrow, xsnplist, i, ncols);
      if (usepopsformissing) { 
        t = xtypes[i] ; 
        for (j=0; j<ncols; ++j) { 
         if (t<0) break ; 
         u = t*ncols + j ; 
         if ((xrow[j]<0) && (xden[u] > 0.1)) {
          xrow[j] = xnum[u]/xden[u] ; 
        }
       }
      }

      copyarr (xrow, trow, ncols);
      fixxrow (xrow, xmean, xfancy, ncols);

      kk = 0;
      for (k = 0; k < ncols; ++k) {
        cupt = xsnplist[k] ; 
        if (cupt -> tagnumber < 0) continue ; 

        if ((blocknum >= 0) && (blocknum == cupt -> tagnumber)) continue ; 
        if (trow[k] < 0) continue;
        rhs[kk] = xrow[k];
        for (j = 0; j < numeigs; j++) {
          emat[kk * numeigs + j] = fxscal[j] * ffvecs[j * ncols + k];
        }
        ++kk;
      }
      if (kk <= numeigs) {
        indx->ignore = YES;
        printf ("%s ignored (insufficient data\n", indx->ID);
        continue;
      }
      regressit (regans, emat, rhs, kk, numeigs);
      for (j = 0; j < numeigs; ++j) {
        acoeffs[j * nind + qk] = regans[j];
      }
    }


    for (qk = 0; qk < nind; qk++) {
      if (bcoeffs == NULL) break ; 
      
      indx = xindlist[qk];
      i = indx -> idnum ; 
      if (indx->ignore) continue;
      loadxdataind (xrow, xsnplist, i, ncols);

      if (usepopsformissing) { 
       t = xtypes[i] ; 
       for (j=0; j<ncols; ++j) { 
        if (t<0) break ;
        u = t*ncols + j ; 
        if ((xrow[j]<0) && (xden[u] > 0.1)) {
         xrow[j] = xnum[u]/xden[u] ; 
       }
      }
     }

      fixxrow (xrow, xmean, xfancy, ncols);

      for (j = 0; j < numeigs; j++) {
        y = fxscal[j] * vdot (xrow, ffvecs + j * ncols, ncols);
        bcoeffs[j * numindivs + qk] = y;
      }

    }

    if (usepopsformissing) { 
      free(xnum) ; 
      free(xden) ;
    }

    return ; 

}
void seteigscale(double *eigscale, double *acoeffs, double *bcoeffs, int *xindex, int nrows, int numeigs)  
{
     double *azq, *bzq ; 
     double *apt, *bpt ; 
     int j ; 


     ZALLOC (azq, nrows * numeigs, double);
     ZALLOC (bzq, nrows * numeigs, double);

     sqz (azq, acoeffs, numeigs, nrows, xindex);
     sqz (bzq, bcoeffs, numeigs, nrows, xindex);

    for (j = 0; j < numeigs; ++j) {
// rescale a to match b  

      apt = azq + j * nrows;
      bpt = bzq + j * nrows;
      eigscale[j] = vdot (apt, bpt, nrows) / vdot (apt, apt, nrows);
    }

    free(azq) ; 
    free(bzq) ;  

}

void
printmega(char *outname, char **eglist, int numeg, double *fstsc) 
{
 FILE *fff ;
 int i, j ; 
 double y ; 

 if (outname == NULL) return ;  
 openit(outname, &fff, "w") ; 

 fprintf(fff, "#mega\n") ;
 fprintf(fff, "!TITLE distance for %d samples ;", numeg) ; 
 fprintf(fff, "\n\n") ; 

 for (i=0; i<numeg; i++) { 
  fprintf(fff, "#%s\n", eglist[i]) ; 
 }

 fprintf(fff, "\n\n") ; 
 for (i=0; i<numeg; i++) { 
  for (j=0; j<i; j++) { 
   y = fstsc[i*numeg+j] ;
   y = MAX(y, 0) ;
   fprintf(fff, "  %9.6f", y) ;
  } 
  fprintf(fff, "\n") ;
 }
 fclose(fff) ;
}
