#include "qpsubs.h" 
#include "eigsubs.h" 
#include "smartsubs.h" 
extern int fancynorm, verbose, plotmode, outnum ;
extern FILE *fstdetails ;

// static Indiv **indm ;
// static void wjackestx(double *est, double *sig, double mean, double *jmean, double *jwt, int g)  ;
// static void wjackvestx(double *vest, double *var, int d, double *mean, double **jmean, double *jwt, int g)  ;
static int  outliermode = 0 ; 

void setoutliermode(int mode)  
{
 outliermode = mode ;
}
int
ridoutlier(double *evecs, int n, int neigs, 
 double thresh, int *badlist, OUTLINFO **outinfo) 
{
/* badlist contains list of outliers */
 double *ww, *w2, y1 , y2, yy, zz ; 
 int *vbad ;
 int i, j  ;
 int nbad = 0 ; 
 OUTLINFO *outpt;

 if (outliermode > 1) return 0;
 if (n<3) return 0;
 ZALLOC(ww, n, double) ;
 ZALLOC(vbad, n, int) ;
 for(j=0;j<n;j++)  {
   outpt = outinfo[j];
   outpt->vecno = -1;
 }
 for (i=0; i<neigs; ++i) {  
  copyarr(evecs+i*n, ww, n) ;
  if (outliermode == 0) { 
   y1 = asum(ww, n) / (double) n ;
   vsp(ww, ww, -y1, n) ;
   y2 = asum2(ww, n) / (double)  n ;
   y2 = sqrt(y2) ;
   vst(ww, ww, 1.0/y2, n) ;

   for (j=0; j<n; j++) {  
    if (fabs(ww[j])>thresh) { 
     vbad[j] = 1 ;
     outpt = outinfo[j];
     if (outpt->vecno < 0)  {
       outpt->vecno = i;
       outpt->score = ww[j];
     }
    }
   }
  }
  if (outliermode == 1) { 
   ZALLOC(w2, n, double) ;
   for (j=0; j<n; j++) {  
     yy = ww[j] ;
     ww[j] = 0 ;
     y1 = asum(ww, n) / (double) (n-1) ;
     vsp(w2, ww, -y1, n) ;
     w2[j] = 0 ;
     y2 = asum2(w2, n) / (double)  n ;
     y2 = sqrt(y2) ;
     zz = yy-y1  ;
     zz /= y2 ;
     if (fabs(zz)>thresh) { 
      vbad[j] = 1 ;
      outpt = outinfo[j];
      if (outpt->vecno < 0)  {
       outpt->vecno = i;
       outpt->score = zz ;
      }
    }
    ww[j] = yy ;
   }
   free(w2) ;
  }
 }
 for (j=0; j<n; j++) {  
  if (vbad[j] == 1) { 
   badlist[nbad] = j ;
   ++nbad ;
  }
 }
 free(ww) ; 
 free(vbad) ;
 return nbad ;

}
