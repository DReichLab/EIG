#include <stdio.h>
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
#include "egsubs.h"  
#include "qpsubs.h" 
#include "smartsubs.h" 
#include "exclude.h" 

/** 
 Most of this code written by Nick Patterson 
 (Broad institute and Harvard Medical)
 Some improvements and elimination of FORTRAN code by Chris Chang (BGI) 

 Code added to support grm output + improved ld rregression by Alexander Gusev 
*/

#define WVERSION   "10210" 
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

fstdetailsname added
fsthiprecision added
bug fixed  (getrawcolx)

bad bug fix.  xtypes not allocated correctly

version compatible with Mac
XTX.dbg commented out

outliermode added

regmode added

*/

#if _WIN32
// just in case we try a Windows port in the future
#include <windows.h>
#include <process.h>
#define pthread_t HANDLE
#define THREAD_RET_TYPE unsigned __stdcall
#define THREAD_RETURN return 0
#define MAX_THREADS 63
#define MAX_THREADS_P1 64
#else
#include <pthread.h>
#define THREAD_RET_TYPE void*
#define THREAD_RETURN return NULL
#define MAX_THREADS 127
#define MAX_THREADS_P1 128
#endif

#define MAXFL  50   
#define MAXSTR  512
#define MAXPOPS 100

char *parname = NULL ;
char *twxtabname = NULL ;
char *trashdir = "/var/tmp" ;
int qtmode = NO ;
Indiv **indivmarkers;
SNP **snpmarkers ;

int numsnps, numindivs ; 
int numeigs = 10 ;  /// default
int markerscore = NO ;
int seed = 0 ;
int chisqmode = NO ;  // approx p-value better to use F-stat
int missingmode = NO ;
int shrinkmode = NO ;
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

int packout = -1 ;
extern enum outputmodetype outputmode  ;
extern int checksizemode ;
extern int packmode ;
extern int numchrom ;
extern int fancynorm ;
extern int verbose ;
int ogmode = NO ;
int fsthiprec = NO ;
int inbreed = NO ;  // for fst
int easymode = NO ;
int regmode = NO ;

int numoutliter = 5, numoutleigs = 10, outliermode = 0 ;  
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
char *grmoutname = NULL ;
int  grmbinary  = NO ;
double blgsize = 0.05 ;  // block size in Morgans */

double r2thresh = -1.0 ;  
double r2genlim = 0.01 ; // Morgans 
double r2physlim = 5.0e6 ;     
int killr2 = NO ;  
int pubmean = YES ;  // change default

double nhwfilter = -1.0;

int thread_ct_config = 0;

int randomfillin = NO ;
int usepopsformissing = NO ; // if YES popmean is used for missing.  Overall mean if all missing for pop

int xchrom = -1 ;
// list of outliers

int ldregress = 0 ;
double ldlimit = 9999.0 ;  /* default is infinity */
double ldr2lo = 0.01 ; 
double ldr2hi = 0.95 ; 
int    ldposlimit = 1000*1000*1000 ;
int ldregx(double *gsource, double *gtarget, double *res, int rsize, 
 int n, double r2lo, double r2hi) ;
void bumpldvv(double *gsource, double *newsource, int *pnumld, int maxld, int n, int *ldsnpbuff, int newsnpnum) ;


char *outputname = NULL ;
char *outputvname = NULL ;
char *weightname = NULL ;
FILE *ofile, *ovfile ;

double twestxx(double *lam, int m, double *pzn,  double *pzvar)  ;
double twnorm(double lam, double m, double n) ;
double rhoinv(double x, double gam) ;

void readcommands(int argc, char **argv) ;
int loadindx(Indiv **xindlist, int *xindex, Indiv **indivmarkers, int numindivs) ;
void loadxdataind(double *xrow, SNP **snplist, int ind,  int ncols) ;           
void fixxrow(double *xrow, double *xmean, double *xfancy, int len)  ;
void dofancy(double *cc, int n, double *fancy) ;
int fvadjust(double *rr, int n, double *pmean, double *fancy)  ;
void getcol(double *cc, double *xdata, int col, int nrows, int ncols)  ;
void getcolxf(double *xcol, SNP *cupt, int *xindex, 
  int nrows, int col, double *xmean, double *xfancy)  ;          
int getcolxz(double *xcol, SNP *cupt, int *xindex, int *xtypes, 
  int nrows, int col, double *xmean, double *xfancy, int *n0, int *n1)  ;
int getcolxz_binary1(int* rawcol, double* xcol, SNP* cupt, int* xindex,
                     int nrows, int col, double* xmean, double* xfancy,
                     int* n0, int* n1);
void getcolxz_binary2(int* rawcol, uintptr_t* binary_cols,
                      uintptr_t* binary_mmask, uint32_t xblock,
                      uint32_t nrows);

void doinbxx(double *inbans, double *inbsd, SNP **xsnplist, int *xindex, int *xtypes, 
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
void domult_increment_lookup(pthread_t* threads, uint32_t thread_ct, double *XTX_lower_tri, double* tblock, uintptr_t* binary_cols, uintptr_t* binary_mmask, uint32_t block_size, uint32_t indiv_ct, double* partial_sum_lookup_buf);
void domult_increment_normal(pthread_t* threads, uint32_t thread_ct, double* XTX_lower_tri, double* tblock, int marker_ct, uint32_t indiv_ct);
void writesnpeigs(char *snpeigname, SNP **xsnplist, double *ffvecs, int numeigs, int ncols)  ;
void dofstxx(double *fstans, double *fstsd, SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg, double blgsize, SNP **snpmarkers, Indiv **indm) ; 
void fixwt(SNP **snpm, int nsnp, double val) ;
void sqz(double *azq, double *acoeffs, int numeigs, int nrows, int *xindex) ;
void dumpgrm(double *XTX, int *xindex, int nrows, int numsnps, Indiv **indivmarkers, int numindivs, char *grmoutname) ;

uint32_t
triangle_divide(int64_t cur_prod, int32_t modif)
{
  // return smallest integer vv for which (vv * (vv + modif)) is no smaller
  // than cur_prod, and neither term in the product is negative.  (Note the
  // lack of a divide by two; cur_prod should also be double its "true" value
  // as a result.)
  int64_t vv;
  if (cur_prod == 0) {
    if (modif < 0) {
      return -modif;
    } else {
      return 0;
    }
  }
  vv = (int64_t)sqrt((double)cur_prod);
  while ((vv - 1) * (vv + modif - 1) >= cur_prod) {
    vv--;
  }
  while (vv * (vv + modif) < cur_prod) {
    vv++;
  }
  return vv;
}

void
parallel_bounds(uint32_t ct, int32_t start, uint32_t parallel_idx, uint32_t parallel_tot, int32_t* bound_start_ptr, int32_t* bound_end_ptr)
{
  int32_t modif = 1 - start * 2;
  int64_t ct_tot = ((int64_t)ct) * (ct + modif);
  *bound_start_ptr = triangle_divide((ct_tot * parallel_idx) / parallel_tot, modif);
  *bound_end_ptr = triangle_divide((ct_tot * (parallel_idx + 1)) / parallel_tot, modif);
}

// set align to 1 for no alignment
void
triangle_fill(uint32_t* target_arr, uint32_t ct, uint32_t pieces, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t start, uint32_t align)
{
  int32_t modif = 1 - start * 2;
  uint32_t cur_piece = 1;
  int64_t ct_tr;
  int64_t cur_prod;
  int32_t lbound;
  int32_t ubound;
  uint32_t uii;
  uint32_t align_m1;
  parallel_bounds(ct, start, parallel_idx, parallel_tot, &lbound, &ubound);
  // x(x+1)/2 is divisible by y iff (x % (2y)) is 0 or (2y - 1).
  align *= 2;
  align_m1 = align - 1;
  target_arr[0] = lbound;
  target_arr[pieces] = ubound;
  cur_prod = ((int64_t)lbound) * (lbound + modif);
  ct_tr = (((int64_t)ubound) * (ubound + modif) - cur_prod) / pieces;
  while (cur_piece < pieces) {
    cur_prod += ct_tr;
    lbound = triangle_divide(cur_prod, modif);
    uii = (lbound - ((int32_t)start)) & align_m1;
    if ((uii) && (uii != align_m1)) {
      lbound = start + ((lbound - ((int32_t)start)) | align_m1);
    }
    // lack of this check caused a nasty bug earlier
    if (((uint32_t)lbound) > ct) {
      lbound = ct;
    }
    target_arr[cur_piece++] = lbound;
  }
}

void
symit2(double* XTX, uintptr_t nrows)
{
  // unpacks LOWER-triangle-only symmetric matrix representation into regular
  // square matrix.
  uintptr_t row_idx;
  uintptr_t col_idx;
  double* read_col;
  double* write_ptr;
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
    memcpy(&(XTX[row_idx * nrows]), &(XTX[(row_idx * (row_idx + 1)) / 2]), (row_idx + 1) * sizeof(double));
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
copy_transposed(double* orig_matrix, uintptr_t orig_row_ct, uintptr_t orig_col_ct, double* transposed_matrix)
{
  uintptr_t new_row_idx;
  uintptr_t new_col_idx;
  double* orig_col_ptr;
  for (new_row_idx = 0; new_row_idx < orig_col_ct; new_row_idx++) {
    orig_col_ptr = &(orig_matrix[new_row_idx]);
    for (new_col_idx = 0; new_col_idx < orig_row_ct; new_col_idx++) {
      *transposed_matrix++ = orig_col_ptr[new_col_idx * orig_col_ct];
    }
  }
}

// make these file scope so multithreading works
static double* g_XTX_lower_tri;
static double* g_tblock;
static uint32_t g_block_size;
static uintptr_t g_indiv_ct;
static uint32_t g_thread_start[MAX_THREADS_P1];

static double* g_weights;
static uintptr_t* g_binary_cols;
static uintptr_t* g_binary_mmask;

int main(int argc, char **argv)
{

  char sss[MAXSTR] ;
  char **eglist ;
  int numeg ;
  int i, j, k, k1, k2, pos; 
  int *vv ;
  SNP *cupt ;
  Indiv *indx ;
  double y1 = 0, y2, y2l, y, y3 ;

  int n0, n1, nkill ;

  int nindiv = 0 ;
  double ychi,  tail, tw ;
  int nignore, numrisks = 1 ;
  double  *xrow, *xpt ; 
  SNP **xsnplist  ;
  Indiv **xindlist ;
  int *xindex, *xtypes = NULL ;
  int nrows, ncols, m, nused ;
  double *XTX, *cc, *evecs, *ww ; 
  double* partial_sum_lookup_buf = NULL;
  double *lambda, *esize ;
  double zn, zvar ;
  double *fvecs, *fxvecs, *fxscal ;
  double *ffvecs ;
  int weightmode = NO ;
  double ynrows ;
  int t, tt ;  
  double *xmean, *xfancy ;
  double *ldvv = NULL , ynumsnps = 0 ; // for grm
  int *ldsnpbuff = NULL ;
  int lastldchrom, numld ;
  double *fstans, *fstsd ; 
  double *inbans, *inbsd ; 

  int chrom ;
  int outliter, numoutiter, *badlist, nbad ;
  FILE *outlfile, *phylipfile  ;
  double *eigkurt, *eigindkurt ;
  double *wmean ;
  char **elist ; 
  double *shrink ; 
  double *trow = NULL, *rhs = NULL, *emat = NULL, *regans = NULL ;
  int kk ;
  double *acoeffs, *bcoeffs, *apt, *bpt, *azq, *bzq ;
  

  int xblock ;
  int blocksize = 1024;
  double *tblock = NULL;
  int* binary_rawcol = NULL;
  uintptr_t* binary_cols = NULL;
  uintptr_t* binary_mmask = NULL;

  OUTLINFO *outpt ;

  pthread_t threads[MAX_THREADS];
  uint32_t thread_ct;

  readcommands(argc, argv) ;
  printf("## smartpca version: %s\n", WVERSION) ;
  packmode = YES ;
  setomode(&outputmode, omode) ;

  if (parname == NULL) return 0 ;
  if (xchrom == (numchrom+1)) noxdata = NO ;

  if (usepopsformissing) {
   printf("usepopsformissing => easymode\n") ;
   easymode = YES ;
  }

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

  if ((ldlimit <= 0) || (ldposlimit<=0)) ldregress = 0 ;

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
    if ((noxdata) && (chrom == (numchrom+1)))   {
      cupt-> ignore = YES ;
      logdeletedsnp(cupt->ID,"chrom-X",deletesnpoutname);
    }
    if (chrom == 0)   {
      cupt -> ignore = YES;   
      logdeletedsnp(cupt->ID,"chrom-0",deletesnpoutname);
    }
    if (chrom > (numchrom+1))   {
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
  if (ldregress>0) 
  {  
    ZALLOC(ldvv,  ldregress*numindivs, double) ;
    ZALLOC(ldsnpbuff, ldregress, int) ; // index of snps 
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

  ZALLOC(xtypes, numindivs, int) ;



  /*  Load non-ignored SNPs into xsnplist:
   *  xsnplist[i] = pointer to SNP struct     */

  nrows = loadindx(xindlist, xindex, indivmarkers, numindivs) ;
  ncols = loadsnpx(xsnplist, snpmarkers, numsnps, indivmarkers) ;

  printf("number of samples used: %d number of snps used: %d\n", nrows, ncols) ;

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
  // if ((!usepopsformissing) && (ldregress == 0)) {
  if (0) {
    // should not use lookup table if
    // - usepopsformissing is set (since each population may have a different
    //   mean), or
    // - ldregress > 0
#ifdef __LP64__
    blocksize = 20;
    ZALLOC(partial_sum_lookup_buf, 131072, double);
#else
    blocksize = 10;
    ZALLOC(partial_sum_lookup_buf, 65536, double);
#endif
    ZALLOC(binary_rawcol, nrows, int);
    ZALLOC(binary_cols, nrows, uintptr_t);
    ZALLOC(binary_mmask, nrows, uintptr_t);
    ZALLOC(tblock, 3 * blocksize, double);
  } else {
    ZALLOC(tblock, nrows*blocksize, double) ;
  }

  ZALLOC(lambda, nrows, double) ;
  ZALLOC(esize, nrows, double) ;
  ZALLOC(cc, (nrows < 3)? nrows : 3, double) ;
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
    setoutliermode(outliermode) ;
  }
  else setoutliermode(2) ;

  // try to autodetect number of (virtual) processors, and use that number to
  // set the thread count.  allow the user to override this in the future
#if _WIN32
  SYSTEM_INFO sysinfo;
  if (thread_ct_config <= 0) {
    GetSystemInfo(&sysinfo);
    thread_ct = sysinfo.dwNumberOfProcessors;
  } else {
    thread_ct = thread_ct_config;
  }
#else
  if (thread_ct_config <= 0) {
    i = sysconf(_SC_NPROCESSORS_ONLN);
    if (i == -1) {
      thread_ct = 1;
    } else {
      thread_ct = i;
    }
  } else {
    thread_ct = thread_ct_config;
  }
#endif
  if (thread_ct > 8) {
    if (thread_ct > MAX_THREADS) {
      thread_ct = MAX_THREADS;
    } else {
      thread_ct--;
    }
  }
  if (thread_ct > nrows * 2) {
    thread_ct = nrows / 2;
    if (!thread_ct) {
      thread_ct = 1;
    }
  }
  printf("Using %u thread%s%s.\n", thread_ct, (thread_ct == 1)? "" : "s", (partial_sum_lookup_buf)? ", and partial sum lookup algorithm" : "");
  triangle_fill(g_thread_start, nrows, thread_ct, 0, 1, 0, 1);

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

    vzero(XTX, (nrows*(nrows+1)) / 2) ;
    xblock = 0 ; 

    vzero(xmean, ncols) ;
    vclear(xfancy, 1.0, ncols) ;

    nused = 0 ;
    for (i=0; i<nrows; i++) {
     indx = xindlist[i] ;
      k= indxindex(eglist, numeg, indx -> egroup) ;
      xtypes[i] = k ;
    }


    numld = 0 ;
    lastldchrom = -1 ;
    ynumsnps = 0 ;
    if (partial_sum_lookup_buf) {
      for (i = 0; i < nrows; i++) {
	binary_cols[i] = 0;
      }
      for (i = 0; i < nrows; i++) {
	binary_mmask[i] = 0;
      }
      vzero(tblock, 3 * blocksize);
    } else {
      vzero(tblock, nrows*blocksize) ;
    }
    for (i=0; i<ncols; i++)  {
      cupt = xsnplist[i] ;
      chrom = cupt -> chrom ;
      if (!partial_sum_lookup_buf) {
        tt = getcolxz(cc, cupt, xindex, xtypes, nrows, i, xmean, xfancy, &n0, &n1) ;
      } else {
        tt = getcolxz_binary1(binary_rawcol, cc, cupt, xindex, nrows, i, xmean, xfancy, &n0, &n1);
      }

      t = MIN(n0, n1) ; 

      if ((t < minallelecnt) || (tt >maxmissing) || (tt<0) || (t==0))  {  
	t =  MAX(t, 0) ;
	tt = MAX(tt, 0) ;
	cupt -> ignore = YES ;
	logdeletedsnp(cupt->ID,"minallelecnt",deletesnpoutname);
	vzero(cc, nrows) ; 
	if (nkill < 10) printf(" snp %20s ignored . allelecnt: %5d  missing: %5d\n", cupt -> ID, t, tt) ;
	++nkill ;
	continue ;
      }

      if (lastldchrom != chrom)  numld = 0 ;

      if (!partial_sum_lookup_buf) {
	if (weightmode)
	{
	  vst(cc, cc, xsnplist[i] -> weight, nrows) ;
	}


	if (ldregress>0) 
	{  

	  t = ldregx(ldvv, cc, ww, numld, nrows, ldr2lo, ldr2hi) ; 
	  if (t<2) {
	    bumpldvv(ldvv, cc, &numld, ldregress, nrows, ldsnpbuff, i) ; 
	    lastldchrom = chrom ;
	    ynumsnps += asum2(ww, nrows)/ asum2(cc, nrows) ; 
  // don't need to think hard about how cc is normalizes
	  } else {
	    // Ignore this SNP and exclude from further regressions (*ww is returned as all zeroes)
	    bumpldvv(ldvv, ww, &numld, ldregress, nrows, ldsnpbuff, i) ; 
	    lastldchrom = chrom ;
	  }
	  copyarr(ww, cc, nrows) ;
	}
	else ++ynumsnps ;
        copyarr(cc, tblock+xblock*nrows, nrows) ;
      } else {
	getcolxz_binary2(binary_rawcol, binary_cols, binary_mmask, xblock, nrows);
	if (weightmode) {
	  vst(cc, cc, xsnplist[i]->weight, 3);
	}
	++ynumsnps;
	copyarr(cc, &(tblock[xblock * 3]), 3);
      }

      ++xblock ; 
      ++nused ;

/** this is the key code to parallelize */
      if (xblock==blocksize) 
      {
	if (partial_sum_lookup_buf) {
	  domult_increment_lookup(threads, thread_ct, XTX, tblock, binary_cols, binary_mmask, xblock, nrows, partial_sum_lookup_buf);
	  for (j = 0; j < nrows; j++) {
	    binary_cols[j] = 0;
	  }
	  for (j = 0; j < nrows; j++) {
	    binary_mmask[j] = 0;
	  }
	  vzero(tblock, 3 * blocksize);
	} else {
          domult_increment_normal(threads, thread_ct, XTX, tblock, xblock, nrows);
          vzero(tblock, nrows*blocksize) ;
	}
        xblock = 0 ;
      }
    }

    if (xblock>0) 
    {
      if (partial_sum_lookup_buf) {
	domult_increment_lookup(threads, thread_ct, XTX, tblock, binary_cols, binary_mmask, xblock, nrows, partial_sum_lookup_buf);
      } else {
        domult_increment_normal(threads, thread_ct, XTX, tblock, xblock, nrows);
      }
    }
    symit2(XTX, nrows) ;
    printf("total number of snps killed in pass: %d  used: %d\n", nkill, nused) ;

    if (verbose) 
    {
      printdiag(XTX, nrows) ;
    }

    y = trace(XTX, nrows) / (double) (nrows-1) ;
    if (isnan(y)) fatalx("bad XTX matrix\n") ;
    /* printf("trace:  %9.3f\n", y) ; */
    if (y<=0.0) fatalx("XTX has zero trace (perhaps no data)\n") ;
    vst(XTX, XTX, 1.0/y, nrows * nrows) ;

    eigvecs(XTX, lambda, evecs, nrows) ;
// eigenvalues are in decreasing order 

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
      fprintf(outlfile, "REMOVED outlier %s iter %d evec %d sigmage %.3f pop: %s\n", 
       indx -> ID, outliter, outpt -> vecno, outpt -> score, indx -> egroup) ;
      indx -> ignore = YES ;
    }
    nrows = loadindx(xindlist, xindex, indivmarkers, numindivs) ;
    printf("number of samples after outlier removal: %d\n", nrows) ;
  }

  if (outliername != NULL) fclose(outlfile) ;
  dumpgrm(XTX, xindex, nrows,  ynumsnps, indivmarkers, numindivs, grmoutname) ; 

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
      esize[i] = zn ;
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

  ZALLOC(shrink, numeigs, double) ;
  vclear(shrink, 1.0, numeigs) ;
  t = nrows - numeigs ; 
  if (t>0) y1 = asum(lambda+numeigs, t)/(double) t ;
  y = (double) nrows / esize[numeigs] ;      
  y = MIN(y, 1.0/y) ; // gamma
  for (j=0; j<numeigs; j++) { 
   if (!shrinkmode) break ; 
   if (t<=0) break ;
   if (esize[j] < 0.1) break ;
   y2 = lambda[j]/y1 ;
// this is d after normalization (Baik Silverman);  now estimate true eigenvalue 
   y2l = rhoinv(y2, y) ;
   if (y2l<0.0) break ;  
   y3 = (y2l-1.0)/(y2l+y-1.0) ;
   y3 = MIN(y3, 1.0) ;
   if (y3<0.0) y3 = 1.0 ;
   shrink[j] = y3 ;
   printf("shrink: %3d %9.3f %9.3f %9.3f\n", j, shrink[j], y2, y2l) ;
  }

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

  
  ZALLOC(acoeffs, numindivs*numeigs, double) ; 
  ZALLOC(bcoeffs, numindivs*numeigs, double) ; 
  if (partial_sum_lookup_buf) {
    free(partial_sum_lookup_buf);
    free(binary_rawcol);
    free(binary_cols);
    free(binary_mmask);
  }
  free(tblock);
  if (regmode) {
   ZALLOC(trow, ncols, double) ;
   ZALLOC(rhs, ncols, double) ; 
   ZALLOC(emat, ncols*numeigs, double) ; 
   ZALLOC(regans, numeigs, double) ; 
/**
   for (j=0; j<numeigs; ++j) { 
     xpt = ffvecs+j*ncols ;
     y = asum2(xpt, ncols) ; 
     fxscal[j] = (double) ncols / sqrt(y*y) ;
   }
*/
  }


  for (i=0; i < numindivs ; i++)  { 
          if (!regmode) break ;
          indx = indivmarkers[i] ;
          if (indx -> ignore) continue ;
          loadxdataind(xrow, xsnplist, i,  ncols) ;            
          copyarr(xrow, trow, ncols) ;  
          fixxrow(xrow, xmean, xfancy, ncols) ;

          kk = 0 ;
          for (k=0; k<ncols; ++k) { 
           if (trow[k]<0) continue ; 
           rhs[kk] = xrow[k] ; 
	   for (j=0; j<numeigs; j++) { 
            emat[kk*numeigs+j] = fxscal[j]*ffvecs[j*ncols+k] ;
           }
           ++kk ;
          }
          if (kk <= numeigs) { 
            indx -> ignore = YES ;
            printf("%s ignored (insufficient data\n", indx -> ID) ;
            continue ;
          }
          regressit(regans, emat, rhs, kk, numeigs) ; 
          for (j=0; j<numeigs; ++j) {
           acoeffs[j*numindivs+i] = regans[j] ;
          }
  }

  for (i=0; i < numindivs ; i++)  { 
          indx = indivmarkers[i] ;
          if (indx -> ignore) continue ;
          loadxdataind(xrow, xsnplist, i,  ncols) ;            
          fixxrow(xrow, xmean, xfancy, ncols) ;

	  for (j=0; j<numeigs; j++) { 
           y = fxscal[j]*vdot(xrow, ffvecs+j*ncols, ncols) ;
           if (shrinkmode && (indx -> affstatus == YES)) y *=shrink[j] ;
           bcoeffs[j*numindivs+i] = y ; 
	  }
  }

  if (!regmode) { 
   free(acoeffs) ;
   acoeffs = bcoeffs ;
  }

  ZALLOC(azq, nrows*numeigs, double) ;
  ZALLOC(bzq, nrows*numeigs, double) ;

  sqz(azq, acoeffs, numeigs, nrows, xindex) ;
  sqz(bzq, bcoeffs, numeigs, nrows, xindex) ;

   for (j=0; j<numeigs; ++j) {  
     if (!regmode) break ;
     apt = azq + j*nrows  ;
     bpt = bzq + j*nrows ;
     y = vdot(apt, bpt, nrows)  / vdot(apt, apt, nrows)  ;
     vst(acoeffs+j*numindivs, acoeffs+j*numindivs, y, numindivs) ;
   }


       for (i=0; i < numindivs ; i++)  { 
          indx = indivmarkers[i] ;
          if (indx -> ignore) continue ;
	  fprintf(ofile, "%20s ", indx -> ID) ;
	  for (j=0; j<numeigs; j++) { 
           y = acoeffs[j*numindivs+i] ; 
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
  for (i = 0; i < numeg; i++) {
    xpopsize[i] = 0;
  }
  printf("nrows: %d\n", nrows);
  for (i=0; i<nrows; i++) { 
    k = xtypes[i] ;
    printf("%d ", k);
    ++xpopsize[k] ;
  }
  printf("\n");

  for (i=0; i<numeg; i++) 
  {  
    printf("population: %3d %20s %4d",i, eglist[i], xpopsize[i]) ; 
    if (xpopsize[i] == 0) printf(" ***") ;
    printnl() ;
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
        fprintf(fstdetails, "F_st %20s %20s %12.6f %12.6f\n", eglist[k1], eglist[k2], 
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
  int i ;
  phandle *ph ;
  int t ;

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
   getstring(ph, "deletsnpoutname:", &deletesnpoutname) ;
   getint(ph, "numeigs:", &numeigs) ;
   getint(ph, "numoutevec:", &numeigs) ; /* changed 11/02/06 */
   getint(ph, "markerscore:", &markerscore) ; 
   getint(ph, "chisqmode:", &chisqmode) ; 
   getint(ph, "missingmode:", &missingmode) ; 
   getint(ph, "shrinkmode:", &shrinkmode) ; 
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
   getint(ph, "regmode:", &regmode) ; 
   getint(ph, "lsqproject:", &regmode) ; 

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
   getint(ph, "ldposlimit:", &ldposlimit) ;  /* bases */
   getdbl(ph, "ldr2lo:", &ldr2lo) ;  
   getdbl(ph, "ldr2hi:", &ldr2hi) ; 
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
   getint(ph, "outliermode:", &outliermode) ; /* test distribution with sample removed. Makes sense for small samples */
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
   getstring(ph, "grmoutname:", &grmoutname) ;
   getint(ph, "grmbinary:", &grmbinary) ;
   getint(ph, "packout:", &packout) ; /* now obsolete 11/02/06 */
   getstring(ph, "twxtabname:", &twxtabname) ;

   getdbl(ph, "r2thresh:", &r2thresh) ;
   getdbl(ph, "r2genlim:", &r2genlim) ;
   getdbl(ph, "r2physlim:", &r2physlim) ;
   getint(ph, "killr2:",  &killr2) ;

   getint(ph, "numchrom:",  &numchrom) ;
   getstring(ph, "xregionname:", &xregionname) ;
   getdbl(ph, "hwfilter:", &nhwfilter) ;

   getint(ph, "numthreads:", &thread_ct_config) ;

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
  return -999 ; 
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

int fvadjust_binary(int c0, int c1, int nmiss, int n, double* cc, double* pmean, double* fancy)
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
    p = 0.5*ymean;
    if (altnormstyle == NO) {
      p = (ysum+1.0)/(2.0*ynum+2.0);
    }
    y = p * (1.0-p);
    if (y>0.0) {
      yfancy = 1.0/sqrt(y);
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
dottest(char *sss, double *vec, char **eglist, int numeg, int *xtypes, int len) 
// vec will always have mean 0 
// perhaps should rewrite to put xa1 etc in arrays
{
   double *w1 ; 
   int *xt ;
   int i, k1, k2, k, n, x1, x2 ;
   double ylike ;
   double ychi ;
   double *wmean ;
   int imax, imin, *isort ;
   static int ncall = 0 ;

   char ss1[MAXSTR] ;
   char ss2[MAXSTR] ;
   double ans, ftail, ftailx, ansx ; 

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
   double y1, top, bot, ftail ;  
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
       printf("%40s %6d %9.3f",ss2, df, chi) ;
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
void bumpldvv(double *gsource, double *newsource, int *pnumld, int maxld, int n, int *ldsnpbuff, int newsnpnum) 
{

   int numld ;
   SNP *cuptnew, *cuptold ;
   int pdiff ; 
   double gdiff ;
   

   numld = *pnumld ; 
  
   cuptnew = snpmarkers[newsnpnum] ;
  
   for (;;) { 
    if (numld==0) break ;
    cuptold = snpmarkers[ldsnpbuff[0]] ;
    pdiff = nnint(cuptnew -> physpos - cuptold -> physpos) ;
    gdiff = cuptnew -> genpos - cuptold -> genpos ;
    if ((pdiff <= ldposlimit) && (gdiff<=ldlimit)) break ;
    copyarr(gsource+n, gsource, (maxld-1)*n) ;     // overlapping move but copyarr works left to right 
    copyiarr(ldsnpbuff+1, ldsnpbuff, (maxld-1)) ;  // overlapping move but copyiarr works left to right 
    --numld ;
   }

   if (numld < maxld) { 
    copyarr(newsource, gsource + numld*n, n) ;
    ldsnpbuff[numld] = newsnpnum ; 
    ++numld ;
    *pnumld = numld ; 
    return ;
   }
   
   if (maxld == numld) {
    copyarr(gsource+n, gsource, (maxld-1)*n) ;  // overlapping move but copyarr works left to right 
    copyiarr(ldsnpbuff+1, ldsnpbuff, (maxld-1)) ;  // overlapping move but copyiarr works left to right 
    --numld ;
   }
   copyarr(newsource, gsource + numld*n, n) ;
   ldsnpbuff[numld] = newsnpnum ; 
   ++numld ;

   *pnumld = numld ; 
   return ;
}

int ldregx(double *gsource, double *gtarget, double *res, int rsize, 
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

  if (rsize==0) { 
   copyarr(gtarget, res, n) ;
   return 0 ;
  }

  // Allocate space for all genotypes to pass
  double *gsource_pass ;
  ZALLOC(gsource_pass , rsize * n , double);

  int i,ii;

  // Compute correlation to previous SNPs
  double sum;
  int rsize_pass = 0 ;
  for ( i = 0 ; i < rsize ; i++ ) {
    sum = 0;
    for ( ii = 0 ; ii < n ; ii++ ) {
      sum += gtarget[ii] * gsource[i*n+ii] ;
    }
    // Normalize by (n-1) and square to get cor^2
    sum = pow(sum / (2*(n-1)),2) ;
    // Check if correlation too high
    if ( sum > r2hi ) {
      // Clean up and exit
      free(gsource_pass);

      // Residual set to all zero
      for ( ii = 0 ; ii < n ; ii++ ) res[ii] = 0;
      return 2;
    // Check if correlation not too low
    } else if ( sum > r2lo ) {
      // Retain this SNP for the regression
     for ( ii = 0 ; ii < n ; ii++ ) gsource_pass[rsize_pass*n+ii] = gsource[i*n+ii] ;
     rsize_pass++;
    }
  }

  // Do the regression if correlated SNPs were found
  if ( rsize_pass > 0 ) {
    double *t_gsource_pass , *regans , *www;
    ZALLOC(regans, rsize, double) ;
    ZALLOC(www, n, double) ;
    ZALLOC(t_gsource_pass , rsize * n , double);

    // Transpose gsource_pass to comply with regressit
    transpose(t_gsource_pass,gsource_pass,rsize,n);

    regressit(regans, t_gsource_pass, gtarget, n, rsize_pass) ; 
    mulmat(www, regans, gsource_pass,  1, rsize_pass, n) ;
    vvm(res, gtarget, www, n) ;

    free(regans) ;
    free(www) ;
    free(t_gsource_pass) ;
    free(gsource_pass);
    return 1;
  } 
  else {
    copyarr(gtarget, res, n) ;
    free(gsource_pass);
    return 0;
  }
}


void dofstxx(double *fstans, double *fstsd, SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg, double blgsize, SNP **snpmarkers, Indiv **indm) 

{

   int nblocks, xnblocks ;  
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
   for (k=0; k<10; ++k) { 
// was <= 10 Tiny bug
    vlmaxmin(snpsc, ncols, &kmax, &kmin) ;
    cupt = xsnplist[kmax] ;
    if (snpsc[kmax]<0) break ;
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
getcolxz(double *xcol, SNP *cupt, int *xindex, int *xtypes, int nrows, int col,  
 double *xmean, double *xfancy, int *n0, int *n1)             
// side effect set xmean xfancy and count variant and reference alleles
// returns missings after fill in
{
 int   j,  n, g, t, k, kmax = -1 ;
 double y, pmean, yfancy ;
 int *rawcol ;
 int c0, c1, nmiss ;
 double* popnum = NULL;
 double* popsum = NULL;

  if (usepopsformissing) { 
   ZALLOC(popnum, MAXPOPS+1, double) ;
   ZALLOC(popsum, MAXPOPS+1, double) ;
  }

  c0 = c1 = 0 ;
  ZALLOC(rawcol, nrows, int) ;
  n = cupt -> ngtypes ;
  if (n<nrows) fatalx("bad snp: %s %d\n", cupt -> ID, n) ;
  getrawcol(rawcol, cupt, xindex, nrows) ;
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
  t = fvadjust(xcol, nrows, &pmean, &yfancy) ;
  if (t < -99) {  
   if (xmean != NULL) {
    xmean[col] = 0.0  ; 
    xfancy[col] = 0.0  ;
   }
   vzero(xcol, nrows) ;
   free(rawcol) ;
   if (n0 != NULL) {
    *n0 = -1 ; 
    *n1 = -1 ;
   }
   return -1;
  }
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

int
getcolxz_binary1(int* rawcol, double* xcol, SNP* cupt, int* xindex, int nrows,
                 int col, double* xmean, double* xfancy, int* n0, int* n1)
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
    fatalx("bad snp: %s %d\n", cupt->ID, n);
  }
  getrawcol(rawcol, cupt, xindex, nrows);
  nmiss = 0;
  for (j=0; j<nrows; ++j) { 
    g = rawcol[j];  
    if (g<0) {
      ++nmiss; 
      continue;  
    }
    c0 += g;
    c1 += 2-g;
  }
  // instead of storing an entire column of floating point values,
  t = fvadjust_binary(c0, c1, nmiss, nrows, xcol, &pmean, &yfancy);
  if (t < -99) {
    if (xmean != NULL) {
      xmean[col] = 0.0; 
      xfancy[col] = 0.0;
    }
    vzero(xcol, 3);
    if (n0 != NULL) {
      *n0 = -1; 
      *n1 = -1;
    }
    return -1;
  }
  vst(xcol, xcol, yfancy, 3);
  if (xmean != NULL) {
    xmean[col] = pmean*yfancy; 
    xfancy[col] = yfancy;
  }
  if (n0 != NULL) {
    *n0 = c0 ; 
    *n1 = c1 ;
  }
  return nmiss ;
}

void
getcolxz_binary2(int* rawcol, uintptr_t* binary_cols, uintptr_t* binary_mmask,
                 uint32_t xblock, uint32_t nrows)
{
  // slightly better to position at 0-3-6-9-12-16-19... instead of
  // 0-3-6-9-12-15-18...
  uint32_t shift_val = (xblock * 3) + (xblock / 5);

  uintptr_t bitfield_or[3];
  uint32_t row_idx;
  int cur_geno;
  bitfield_or[0] = ((uintptr_t)7) << shift_val;
  bitfield_or[1] = ((uintptr_t)2) << shift_val;
  bitfield_or[2] = ((uintptr_t)3) << shift_val;
  for (row_idx = 0; row_idx < nrows; row_idx++) {
    cur_geno = *rawcol++;
    if (cur_geno) {
      if (cur_geno > 0) {
        binary_cols[row_idx] |= bitfield_or[(uint32_t)cur_geno];
      } else {
        binary_mmask[row_idx] |= bitfield_or[0];
      }
    }
  }
}

void
join_threads(pthread_t* threads, uint32_t ctp1)
{
  if (!(--ctp1)) {
    return;
  }
#if _WIN32
  WaitForMultipleObjects(ctp1, threads, 1, INFINITE);
#else
  uint32_t uii;
  for (uii = 0; uii < ctp1; uii++) {
    pthread_join(threads[uii], NULL);
  }
#endif
}

#if _WIN32
int32_t
spawn_threads(pthread_t* threads, unsigned (__stdcall *start_routine)(void*), uintptr_t ct)
#else
int32_t
spawn_threads(pthread_t* threads, void* (*start_routine)(void*), uintptr_t ct)
#endif
{
  uintptr_t ulii;
  if (ct == 1) {
    return 0;
  }
  for (ulii = 1; ulii < ct; ulii++) {
#if _WIN32
    threads[ulii - 1] = (HANDLE)_beginthreadex(NULL, 4096, start_routine, (void*)ulii, 0, NULL);
    if (!threads[ulii - 1]) {
      join_threads(threads, ulii);
      return -1;
    }
#else
    if (pthread_create(&(threads[ulii - 1]), NULL, start_routine, (void*)ulii)) {
      join_threads(threads, ulii);
      return -1;
    }
#endif
  }
  return 0;
}

THREAD_RET_TYPE block_increment_binary(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t cur_indiv_idx = g_thread_start[tidx];
  uintptr_t end_indiv_idx = g_thread_start[tidx + 1];
  uintptr_t* binary_cols = g_binary_cols;
  uintptr_t* binary_mmask = g_binary_mmask;
  double* write_ptr = &(g_XTX_lower_tri[(cur_indiv_idx * (cur_indiv_idx + 1)) / 2]);
  double* weights0 = g_weights;
  double* weights1 = &(g_weights[32768]);
#ifdef __LP64__
  double* weights2 = &(g_weights[65536]);
  double* weights3 = &(g_weights[98304]);
#endif
  uintptr_t* geno_ptr;
  uintptr_t* mmask_ptr;
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
	*write_ptr += weights0[(uint16_t)final_geno] + weights1[(uint16_t)(final_geno >> 16)] + weights2[(uint16_t)(final_geno >> 32)] + weights3[final_geno >> 48];
#else
        *write_ptr += weights0[(uint16_t)final_geno] + weights1[final_geno >> 16];
#endif
        write_ptr++;
      }
    } else {
      for (indiv_idx2 = 0; indiv_idx2 <= cur_indiv_idx; indiv_idx2++) {
	final_geno = ((*geno_ptr++) + base_geno) | ((*mmask_ptr++) | base_mmask);
#ifdef __LP64__
	*write_ptr += weights0[(uint16_t)final_geno] + weights1[(uint16_t)(final_geno >> 16)] + weights2[(uint16_t)(final_geno >> 32)] + weights3[final_geno >> 48];
#else
        *write_ptr += weights0[(uint16_t)final_geno] + weights1[final_geno >> 16];
#endif
	write_ptr++;
      }
    }
  }
  THREAD_RETURN;
}

void
domult_increment_lookup(pthread_t* threads, uint32_t thread_ct, double *XTX_lower_tri, double* tblock, uintptr_t* binary_cols, uintptr_t* binary_mmask, uint32_t block_size, uint32_t indiv_ct, double* partial_sum_lookup_buf)
{
  // PLINK 1.5 partial sum lookup algorithm
  double increments[40];
  double* dptr;
  double* dptr2;
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
  if (spawn_threads(threads, block_increment_binary, thread_ct)) {
    fatalx("Error: Failed to create thread.\n");
    return;
  }
  ulii = 0;
  block_increment_binary((void*)ulii);
  join_threads(threads, thread_ct);
}

THREAD_RET_TYPE block_increment_normal(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t start_indiv_idx = g_thread_start[tidx];
  uintptr_t end_indiv_idx = g_thread_start[tidx + 1];
  uintptr_t indiv_ct = g_indiv_ct;
  uint32_t block_size = g_block_size;
  double* write_start_ptr = &(g_XTX_lower_tri[(start_indiv_idx * (start_indiv_idx + 1)) / 2]);
  double* write_ptr;
  double* tblock;
  double* tblock_read_ptr;
  double cur_tblock_val;
  uintptr_t cur_indiv_idx;
  uintptr_t indiv_idx2;
  uint32_t bidx;
  for (bidx = 0; bidx < block_size; bidx++) {
    write_ptr = write_start_ptr;
    tblock = &(g_tblock[bidx * indiv_ct]);
    for (cur_indiv_idx = start_indiv_idx; cur_indiv_idx < end_indiv_idx; cur_indiv_idx++) {
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
domult_increment_normal(pthread_t* threads, uint32_t thread_ct, double* XTX_lower_tri, double* tblock, int block_size, uint32_t indiv_ct)
{
  // tblock[] can have an arbitrary number of distinct values, so can't use
  // bit hacks
  int ii;
  double ycheck;
  uintptr_t ulii;
  for (ii=0; ii<block_size; ii++) {  
    ycheck = asum(tblock+ii*indiv_ct, indiv_ct) ;
    if (fabs(ycheck)>.00001) fatalx("bad ycheck\n");
  }
  g_XTX_lower_tri = XTX_lower_tri;
  g_tblock = tblock;
  g_block_size = block_size;
  g_indiv_ct = indiv_ct;
  if (spawn_threads(threads, block_increment_normal, thread_ct)) {
    fatalx("Error: Failed to create thread.\n");
    return;
  }
  ulii = 0;
  block_increment_normal((void*)ulii);
  join_threads(threads, thread_ct);
}

void
getcolxf(double *xcol, SNP *cupt, int *xindex, int nrows, int col,
 double *xmean, double *xfancy)
// side effect set xmean xfancy
{
 int n ;
 double pmean, yfancy ;
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

void doinbxx(double *inbans, double *inbsd, SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg, double blgsize, SNP **snpmarkers, Indiv **indm) 
{

   int nblocks, xnblocks ;  
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

void
sqz(double *azq, double *acoeffs, int numeigs, int nrows, int *xindex) 
{

   int i, j, k ;
   // Indiv *indx ;
   static int ncall = 0 ; 

   ++ncall ;

   for (k=0; k<nrows; ++k)  {
    i = xindex[k] ;
    if (i<0) fatalx("zzyuk!\n") ;
    // indx = indivmarkers[i] ; 
//  if (ncall == 1) printf("zz %3d %12s %12s %d %d\n", k, indx -> ID, indx -> egroup, indx -> ignore, indx -> affstatus) ;

      for (j=0; j<numeigs; ++j) { 
       azq[j*nrows+k] = acoeffs[j*numindivs+i] ;
      }
   }
}
void dumpgrmid(char *fname, Indiv **indivmarkers, int *xindex, int numid) 
{
  FILE *fff ; 
  int a, b ;
  Indiv *indx ;

  openit (fname, &fff, "w") ;
  for (a=0; a<numid; ++a) { 
   b = xindex[a] ;
   if ((b<0) || (b>=numindivs)) fatalx("(dumpgrmid) bad index\n") ;
   indx = indivmarkers[b] ;
   fprintf(fff, "%s\t%s\n", "NA", indx -> ID) ;
  }
  fclose(fff) ;
}
void
dumpgrmbin(double *XTX, int *xindex, int nrows, int numsnps, Indiv **indivmarkers, int numindivs, char *grmoutname)  
{
  int a, b, xa, xb ;
  double y ;
  char sss[256] ;
  char *bb ;  
  int wout, numout, fdes, ret = 0 ;
  float yfloat ;
  
  if (sizeof(yfloat) != 4) fatalx("grm binary only supported for 4 byte floats\n") ;
  
  sprintf(sss, "%s.N.bin", grmoutname) ;
  ridfile(sss) ;
  fdes = open(sss, O_CREAT | O_TRUNC | O_RDWR, 0666);

  if (fdes<0) {
    perror("bad dumpgrmbin") ;
    fatalx("open failed for %s\n", sss) ;
  }
  if (verbose)
    printf("file %s opened\n", sss) ;

//  numout = numsnps*(numsnps+1)/4 ;
  numout = nrows*(nrows+1)/2 ;
  wout = numsnps ;
  bb = (char *) &wout ;
  
  for (a=0; a<numout; ++a) {
   ret = write(fdes, bb, 4) ;
  }
  if (ret<0) {
    perror("write failure") ;
    fatalx("(outpack) bad write") ;
  }
  close(fdes) ;

  sprintf(sss, "%s.bin", grmoutname) ;
  ridfile(sss) ;
  fdes = open(sss, O_CREAT | O_TRUNC | O_RDWR, 0666);

  if (fdes<0) {
    perror("bad dumpgrmbin") ;
    fatalx("open failed for %s\n", sss) ;
  }
  if (verbose)
    printf("file %s opened\n", sss) ;

  // Re-adjust values based on diagonal normalization
  double y_norm ;
  y_norm = trace(XTX, nrows) / (double) nrows ;

  bb = (char *) &yfloat ;
  for (a = 0; a < nrows; a++ ){
  for (b = 0; b <= a; b++ ){
    xa = xindex[a] ;
    xb = xindex[b] ;
    y = XTX[xa*nrows+xb] / y_norm;
    yfloat = (float) y ;
    ret = write(fdes, bb, 4) ;
  }
  }
  close(fdes) ;
}
void
dumpgrm(double *XTX, int *xindex, int nrows, int numsnps, Indiv **indivmarkers, int numindivs, char *grmoutname)  
{
  int a, b, xa, xb ;
  double y ;
  FILE *fff ;
  char sss[256] ;
  
  if (grmoutname == NULL) return ;

  sprintf(sss, "%s.id", grmoutname) ;
  dumpgrmid(sss, indivmarkers, xindex, nrows) ;

  if (grmbinary) { 
   dumpgrmbin(XTX, xindex, nrows, numsnps, indivmarkers, numindivs, grmoutname)  ;
   return ;
  }

  // Re-adjust values based on diagonal normalization
  double y_norm ;
  double *d ;
  ZALLOC(d, nrows, double) ;
  getdiag(d, XTX, nrows) ;
  y_norm = asum(d,nrows) / (double) nrows ;
  free(d) ;

  openit(grmoutname, &fff, "w") ;
  for (a = 0; a < nrows; a++ ){
  for (b = 0; b <= a; b++ ){
    xa = xindex[a] ;
    xb = xindex[b] ;
    y = XTX[xa*nrows+xb] ;
    fprintf(fff, "%d %d ", a+1, b+1) ;
    fprintf(fff, "%d ", numsnps) ;
    fprintf(fff, "%0.6f\n", y/y_norm) ;
  }
  }
  fclose(fff) ;

}
