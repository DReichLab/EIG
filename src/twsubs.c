#include <stdio.h>

#include <limits.h>
#include <math.h>  
#include <nicklib.h> 
#include "admutils.h"  
#include "eigsubs.h" 

/* ********************************************************************* */

extern int verbose ;

int twl2mode = YES ;
int mval = -1 ;
int nval = -1 ;
int numsamp = 100;
double mul1 = 1.0 ;


double xxlike(int m, double a, double var, double logsum, double lsum)    ;
double xxlikex(int m, double a, double logsum, double lsum)    ;
double xxliked(int m, double a, double logsum, double lsum)    ;
double xxliked2(int m, double a, double logsum, double lsum)    ;
double oldtwestxx(double *lam, int m, double *pzn,  double *pzvar)  ;
double doeig2(double *vals, int m, double *pzn, double *ptw)  ;

double twestxx(double *lam, int m, double *pzn,  double *pzvar) 
{
  double tw, y ;

  if (twl2mode == NO)  return oldtwestxx(lam, m, pzn, pzvar) ;
  (void) doeig2(lam,  m, pzn, &tw)  ;

  y = (*pzn) * (double) m ;
  *pzvar = asum(lam, m) / y ;
  return tw ;

}
double oldtwestxx(double *lam, int m, double *pzn,  double *pzvar) 
{
  double lsum, logsum ;
  double *ww ;
  double  a, p, yn, var ;
  double ylike, ybase, y, ylmax, ynmax, yld, yld2, ainc, ym ;
  int k ;


  ZALLOC(ww, m, double) ;
  copyarr(lam, ww, m) ;
  lsum = asum(ww, m) ;
  vlog(ww, ww, m) ;
  logsum = asum(ww, m) ;
  
  ylmax = -1.0e20 ;
  yn = (double) m ;
  ybase = xxlikex(m, yn, logsum, lsum) ; 

  for (k= 1; k<=100; ++k) {  
   a = yn/2.0 ;
   ylike = xxlikex(m, a, logsum, lsum) ; 
   yld = xxliked(m, a, logsum, lsum) ; 
   ylike -= ybase ;
   if (verbose) 
    printf("ynloop %12.3f %12.3f %12.3f\n", yn / (double) m , ylike, yld) ;
   if (ylike < ylmax) break ;
   ylmax = ylike ;
   ynmax = yn ;
   yn *= 1.1 ;
  }
  a = ynmax/2.0 ;
  for (k= 1; k<=10; ++k) {  
// newton iteration
   ylike = xxlikex(m, a, logsum, lsum) ; 
   yld  =  xxliked(m, a, logsum, lsum) ; 
   yld2 =  xxliked2(m, a, logsum, lsum) ; 
   ylike -= ybase ;
   ainc = -yld/yld2 ;
   a += ainc ;
   if (verbose) 
   printf("newton: %3d  %15.9f  %15.9f  %15.9f\n", k, ylike, yld, ainc)  ;
  }
  fflush(stdout) ;
  yn = 2.0*a ;
  ym = (double) m ;
  var = lsum/ (2.0*a*ym) ;

  *pzn = yn ;
  *pzvar = var ;

  free(ww) ;
  return 0 ;
}
double xxlike(int m, double a, double var, double logsum , double lsum)   
{
  double p , yl = 0.0 ;
  double ym , x ;
  int j ;

  p = 0.5* (double) (m+1) ;
  ym = (double) m ;

  yl = -ym*a*log(2.0) ;
  for (j=1;  j<= m; ++j) {  
   x = a - 0.5* (double) (m-j) ;
   yl -= lgamma(x) ;
  }
// so far this is log (C_L)    normalizing constant
  yl -= ym*a*log(var) ;
  yl += (a-p)*logsum ;
  yl -= lsum/(2.0*var) ;

  return yl ;

}
double xxlikex(int m, double a, double logsum , double lsum)   
{
  double p , yl = 0.0 ;
  double ym , x, var, lco ;
  int j ;

  p = 0.5* (double) (m+1) ;
  ym = (double) m ;
  lco = lsum/(2.0*ym) ;
  var = lco/a ;


  yl = -ym*a*log(2.0) ;
  for (j=1;  j<= m; ++j) {  
   x = a - 0.5* (double) (m-j) ;
   yl -= lgamma(x) ;
  }
// so far this is log (C_L)    normalizing constant
  yl -= ym*a*log(var) ;
  yl += (a-p)*logsum ;
  yl -= ym*a ;  // plugging in var 

  return yl ;

}
double xxliked(int m, double a, double logsum , double lsum)   
// first deriv wrt a
{
  double p , yl = 0.0 ;
  double ym , x, var, vard ;
  int j ;

  p = 0.5* (double) (m+1) ;
  ym = (double) m ;
  var = lsum/ (2.0*a*ym) ;
  vard = -var/a ;


  yl = -ym*log(2.0) ;
  for (j=1;  j<= m; ++j) {  
   x = a - 0.5* (double) (m-j) ;
   if (x<0.0) return 100.0 ;
   yl -= psi(x) ;
  }
// so far this is log (C_L)    normalizing constant
  yl -= ym*log(var) ;
  yl -= (ym*a/var)*vard ;
  yl += logsum ;
  yl -= ym ;  // plugging in var 

  return yl ;

}
double xxliked2(int m, double a, double logsum , double lsum)   
// second deriv wrt a
{
  double p , yl = 0.0 ;
  double ym , x, var, vard, vard2, y ;
  int j ;

  p = 0.5* (double) (m+1) ;
  ym = (double) m ;
  var = lsum/ (2.0*a*ym) ;
  vard = -var/a ;
  vard2 = 2.0*var/(a*a) ;


  yl = 0.0 ; 
  for (j=1;  j<= m; ++j) {  
   x = a - 0.5* (double) (m-j) ;
   if (x<0.0) return 100.0 ;
   yl -= tau(x) ;
  }
// so far this is log (C_L)    normalizing constant
  yl -= 2.0*(ym/var)*vard ;
  yl -= (ym*a/var)*vard2 ;
  y = vard/var ;
  yl += (ym*a)*y*y ;

  return yl ;

}
double doeig2(double *vals, int m, double *pzn, double *ptw) 
{
  static int ncall = 0 ;
  double y, tw, tail ;
  double zn, zvar, top, bot ;
  double *evals ;
 
  ++ncall ;
  ZALLOC(evals, m, double) ;
  copyarr(vals, evals, m) ;
  y = (double) m / asum(evals, m) ;
  vst(evals, evals, y, m) ;      
  top = (double) (m*(m+2)) ;
  bot = asum2(evals, m) - (double) m ;
  zn = top/bot ;
  y = evals[0]*zn ;
  tw = twnorm(y, (double) m, zn) ;
  tail = twtail(tw) ;
  free(evals) ;
  *pzn = zn ;
  *ptw = tw ;  
  return tail ;
}
double rhoinv(double x, double gam) 
// Lee et al. page 5 for \rho^{-1} 
{
  double y1, y2 ;  

  y1 = x + 1.0 - gam ;  
  y2 = y1*y1-4.0*x ; 
  if (y2 <= 0.0) return -1.0 ;

  y1 += sqrt(y2) ; 

  return 0.5*y1 ;

}


