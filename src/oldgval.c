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
static uint8_t *xind_mask;
static size_t xtda;

void setgval (
        SNP ** xsnps,
        int nrows,
        Indiv ** indivmarkers,
        int numindivs,
        int *xindex,
        int *xtypes,
        int ncols) {

    double *cc;
    int t, n0, n1, i, k, col;
    SNP *cupt;
    double mean, y;

    unsetgval();

    xxsnps = xsnps;
    xnrows = nrows;
    xncols = ncols;
    xindivmarkers = indivmarkers;
    xnumindivs = numindivs;
    xxindex = xindex;
    ZALLOC(cc, nrows, double);
    ZALLOC(xmean, ncols, double);
    ZALLOC(xfancy, ncols, double);
    vclear(xfancy, 1.0, ncols);
    gtable = initarray_2Ddouble(ncols, 4, 0);

    xtda = (xnrows + 3) / 4;

    for (i = 0; i < ncols; ++i) {
        col = i;
        cupt = xsnps[i];
        /**
         if (i>=0) {
         printf("zz: %d %s\n", cupt -> ID) ;  fflush(stdout) ;
         }
         */
        getcolxz(cc, cupt, xindex, xtypes, nrows, i, xmean, xfancy, &n0, &n1);

        mean = xmean[col] / xfancy[col];
        for (k = 0; k < 3; ++k) {
            y = ((double) k) - mean;
            y *= xfancy[col];
            gtable[col][k] = y / sqrt(2.0);
        }
        gtable[col][3] = 0;

        t = MIN(n0, n1);
        if (t == 0) cupt->ignore = YES;	// side-effect
    }

    set_ind_mask();

    free(cc);
}

void set_ind_mask () {
    size_t i, j, k;
    xind_mask = calloc(xtda, sizeof(uint8_t));
    for (i = 0; i < xnrows; i++) {
        if (xindivmarkers[i]->ignore) {
            j = i / 4;
            k = (3 - (i % 4)) * 2;
            xind_mask[j] |= 3 << k;
        }
    }

}

void unsetgval () {
    if (xxsnps == NULL) return;

    xxsnps = NULL;
    xindivmarkers = NULL;
    xxindex = NULL;

    free2D(&gtable, xncols);

    gtable = NULL;

    free(xmean);
    free(xfancy);
    free(xind_mask);
}

int getgval (int row, int col, double *val) {

    /**
     if (row>=xnrows) fatalx("row index overflow\n") ;
     if (col>=xncols) fatalx("col index overflow\n") ;
     */

    return getggval(xxindex[row], col, val);

}

int getggval (int indindx, int col, double *val)
// indindex is index in full array
{
    SNP *cupt;
    int t, z;
    double y, mean;

    *val = 0;
    if (xindivmarkers[indindx]->ignore) return -1;
    cupt = xxsnps[col];
    t = getgtypes(cupt, indindx);
    if (t < 0) return t;

    *val = gtable[col][t];
    return t;

// dead code
    y = (double) t;
    mean = xmean[col] / xfancy[col];
    y -= mean;
    y *= xfancy[col];

    /**
     z = ranmod(10000000) ;
     if (z==0) {
     printf("zzcheck: %d %d %12.6f %12.6f   %12.6f\n", indindx, col, xmean[col], xfancy[col], y) ;
     }
     */

    *val = y / sqrt(2.0);
    return t;

}

// Unpack lookup table

// macro to unpack a single byte
#define U0(n) { (n >> 6) & 3, (n >> 4) & 3, (n >> 2) & 3, n & 3 }

// macros to build the unpacking table
#define U1(n) U0(n),        U0(n+1),        U0(n+2),        U0(n+3)
#define U2(n) U1(n), U1((n)+(1<<2)), U1((n)+(2<<2)), U1((n)+(3<<2))
#define U3(n) U2(n), U2((n)+(1<<4)), U2((n)+(2<<4)), U2((n)+(3<<4))

// the unpacking table
static const uint8_t UL[256][4] = { U3(0), U3(1<<6), U3(2<<6), U3(3<<6) };

size_t get_nrows () {
    return (xnrows);
}

size_t get_ncols () {
    return (xncols);
}

/**
 * Unpacks a SNP column
 * @param snp_index
 * @param *y arrayref to store data
 */
void kjg_geno_get_normalized_row (const size_t snp_index, double* y) {
    size_t j;

    // Newer method looking up 4 at once

    size_t t = xtda;
    uint8_t* packed = xxsnps[snp_index]->pbuff;
    uint8_t* ind_mask = xind_mask;
    double* norm_lookup = gtable[snp_index];

    while (--t) {
        const uint8_t* u = UL[*(packed++) | *(ind_mask++)];
        for (j = 0; j < 4; j++)
            *(y++) = norm_lookup[*(u++)];
    }

    const uint8_t* u = UL[*packed | *ind_mask];
    for (j = (xtda - 1) * 4; j < xnrows; j++)
        *(y++) = norm_lookup[*(u++)];

    // using getgval (slower)
/*
    for (j = 0; j < xnrows; j++)
	getgval (j, snp_index, y++);
*/

}

/**
 * Unpacks several SNP coluns
 * @param snp_index index of the SNP
 * @param *unpacked arrayref to store data
 */

size_t kjg_geno_get_normalized_rows (const size_t i, const size_t r, double* Y) {
    size_t j;
    for (j = i; j < i + r && j < xncols; j++) {
        kjg_geno_get_normalized_row(j, Y);
        Y += xnrows;
    }
    return (j - i);
}
