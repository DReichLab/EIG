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

#include "admutils.h"
#include "mcio.h"
#include "gval.h"

static SNP **xxsnps = NULL;
static Indiv **xindivmarkers = NULL;
static int xnrows, xncols;
static int xnumindivs;
static int *xxindex = NULL;
static double *xmean, *xfancy;
static double **gtable;

void
setgval (SNP ** xsnps, int nrows, Indiv ** indivmarkers, int numindivs,
         int *xindex, int *xtypes, int ncols)
{

  double *cc;
  int t, n0, n1, i, k, col;
  SNP *cupt;
  double mean, y;

  unsetgval ();

  xxsnps = xsnps;
  xnrows = nrows;
  xncols = ncols;
  xindivmarkers = indivmarkers;
  xnumindivs = numindivs;
  xxindex = xindex;

  for (i = 1; i < nrows; i++) {
    if (xxindex[i] < xxindex[i - 1]) {
      fprintf (stderr, "xindex not sorted\n");
      exit (1);
    }
  }

  ZALLOC (cc, nrows, double);
  ZALLOC (xmean, ncols, double);
  ZALLOC (xfancy, ncols, double);
  vclear (xfancy, 1.0, ncols);
  gtable = initarray_2Ddouble (ncols, 4, 0);

  for (i = 0; i < ncols; ++i) {
    col = i;
    cupt = xsnps[i];

      /**
       if (i>=0) {
       printf("zz: %d %s\n", cupt -> ID) ;  fflush(stdout) ;
       }
       */
    getcolxz (cc, cupt, xindex, xtypes, nrows, i, xmean, xfancy, &n0, &n1);

    mean = xmean[col] / xfancy[col];
    for (k = 0; k < 3; ++k) {
      y = ((double) k) - mean;
      y *= xfancy[col];
      gtable[col][k] = y / sqrt (2.0);
    }
    gtable[col][3] = 0;

    t = MIN (n0, n1);
    if (t == 0)
      cupt->ignore = YES;       // side-effect
  }

  free (cc);
}

void
unsetgval ()
{
  if (xxsnps == NULL)
    return;

  xxsnps = NULL;
  xindivmarkers = NULL;
  xxindex = NULL;

  free2D (&gtable, xncols);

  gtable = NULL;

  free (xmean);
  free (xfancy);
}

int
getgval (int row, int col, double *val)
{

  /**
   if (row>=xnrows) fatalx("row index overflow\n") ;
   if (col>=xncols) fatalx("col index overflow\n") ;
   */

  return getggval (xxindex[row], col, val);

}

int
getggval (int indindx, int col, double *val)
// indindex is index in full array
{
  SNP *cupt;
  int t, z;
  double y, mean;

  *val = 0;
  if (xindivmarkers[indindx]->ignore)
    return -1;
  cupt = xxsnps[col];
  t = getgtypes (cupt, indindx);
  if (t < 0)
    return t;

  *val = gtable[col][t];
  return t;

}

// Unpack lookup table

// macro to unpack a single byte
#define U(n) { ((n) >> 6) & 3, ((n) >> 4) & 3, ((n) >> 2) & 3, (n) & 3 }

// macros to build the u(n)packi(n)g table
#define U1(n)  U(n),  U((n) +  1),  U((n) +  2),  U((n) +  3)
#define U2(n) U1(n), U1((n) +  4), U1((n) +  8), U1((n) + 12)
#define U3(n) U2(n), U2((n) + 16), U2((n) + 32), U2((n) + 48)

// the unpacking table
static const uint8_t UL[256][4] = { U3 (0), U3 (64), U3 (128), U3 (192) };

size_t
get_nrows ()
{
  return (xnrows);
}

size_t
get_ncols ()
{
  return (xncols);
}

/**
 * Unpacks a SNP column
 * @param snp_index
 * @param *y arrayref to store data
 */
void
kjg_geno_get_normalized_row (const size_t snp_index, double *y)
{
  uint8_t *packed = xxsnps[snp_index]->pbuff;
  double *norm_lookup = gtable[snp_index];

  size_t i = 0, j = xxindex[i];
  while (1) {
    size_t k = j / 4;           // packed location
    size_t jf = (k + 1) * 4;    // last index in packed location

    uint8_t p = packed[k];      // packed data
    const uint8_t *u = UL[p];   // unpacked data

    while (j < jf) {
      size_t o = j % 4;         // offset in packed data
      size_t t = u[o];          // unpacked data
      y[i] = norm_lookup[t];    // normalized data

      if (++i == xnrows)        // move onto next entry
        return;                 // break if we are done with SNP
      j = xxindex[i];           // perform the lookup
    }
  }
}

/**
 * Unpacks several SNP coluns
 * @param snp_index index of the SNP
 * @param *unpacked arrayref to store data
 */

size_t
kjg_geno_get_normalized_rows (const size_t i, const size_t r, double *Y)
{
  size_t j;
  for (j = i; j < i + r && j < xncols; j++) {
    kjg_geno_get_normalized_row (j, Y);
    Y += xnrows;
  }
  return (j - i);
}
