#include<math.h>
#include<stdlib.h>


#include <values.h>

#define BIGINT MAXINT    
#define SRAND  srandom
#define LRAND  random
#define DRAND() ( (double) (random() % BIGINT) / (double) (BIGINT)) 
#define DRAND2() ( drand2() ) 
/* random must return random integer in range 0 to BIGINT-1  */


#define NORMAL gauss



double  gauss()  ;
void    gaussa(double *a, int n)  ;
double  gds(double a)  ;
double  poidev(double mean) ;
double  ranpoiss(double mean) ;
double  ranpoissx(double mean) ;
void  ranperm(int *a, int n) ;

double ranexp( void) ;
double rangam(double a) ;
int randis(double *a, int n) ;
void ransamp(int *samp, int nsamp, double *p, int plen) ;
void pick2(int n, int *k1, int *k2) ;
int ranmod(int n)  ;
double ranbeta(double a, double b) ;
int ranbinom(int n, double p) ;
void ewens(int *a, int n, double theta) ;
void genmultgauss(double *rvec, int num, int n, double *covar) ;
double drand2() ;
void ranmultinom(int *samp, int n, double *p, int len)  ;
double ranchi (int d)  ;
double raninvwis(double *wis, int t, int d, double *s)  ;
double uniform(double lo, double hi) ;
void randirichlet(double *x, double *pp, int n)  ;
void randirmult(double *pp, int *aa, int len, int m) ;
