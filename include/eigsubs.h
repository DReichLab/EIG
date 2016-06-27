#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h> 

void
eigvals (double *mat, double *evals, int n);
void
eigvecs (double *mat, double *evals, double *evecs, int n);
void
eigb (double *lam, double *a, double *b, int n);
void
eigc (double *lam, double *a, double *b, int n);
double
twestxx (double *lam, int m, double *pzn, double *pzvar);

typedef struct
{
  int vecno;
  double score;
} OUTLINFO;
;

