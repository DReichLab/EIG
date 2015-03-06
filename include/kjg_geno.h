/*
 * kjg_geno.h
 *
 *  Created on: Jul 31, 2013
 *      Author: kjg063
 */

#ifndef KJG_GENO_H_
#define KJG_GENO_H_

#include <stddef.h>
#include <stdint.h>

#include "admutils.h"

/**
 * Compute the mean of the genotypes.
 *
 * @param *x array of genotypes
 * @param n number of subjects
 * @return mean of the genotypes whose assays didn't fail
 */
double kjg_geno_mean(uint8_t *x, size_t n);

/**
 * Normalize genotypes.
 *
 * @param *x array of genotypes
 * @param *y array to put scaled genotypes
 * @param n number of subjects
 * @return success or failure (-1)
 */
int kjg_geno_normalize(uint8_t *x, double *y, size_t n);

/**
 * Normalize genotypes when you already have the mean.
 *
 * @param m mean of the genotypes
 * @param *x array of genotypes
 * @param *y array to put scaled genotypes
 * @param n number of subjects
 * @return success or failure (-1)
 */
int kjg_geno_normalize_m(double m, uint8_t* x, double* y,
        size_t n);

/**
 * Remap genotypes values to doubles.
 *
 * @param s array of genotype mappings
 * @param *x array of genotypes
 * @param *y array to put scaled genotypes
 * @param n number of subjects
 * @return success (0) or zero-good genos (1)
 */
void kjg_geno_remap(double s[4], uint8_t* x, double* y,
        size_t n);

/**
 * Compute the normalization lookup array.
 *
 * @param m genotype mean
 * @param s[4] array to store the scale
 * @return success (0) or zero-good genos (1)
 */
int kjg_geno_normalization_lookup(double m, double s[4]);

/**
 * Unack genotype array
 *
 * @param *x where to unpack genotypes
 * @param *y packed genotypes
 * @param n number of genotypes
 */

void kjg_geno_unpack(uint8_t* x,
        SNP *cupt, Indiv **indivmarkers,
        size_t numindivs);

/**
 * Get a row in the geno object
 *
 * @param *x unpacked genotype row
 * @param *g geno object
 * @param i row index
 */

void kjg_geno_get_row(uint8_t* x,
        SNP **snpmarkers, Indiv **indivmarkers,
        size_t numsnps, size_t numindivs,
        size_t i);

/**
 * Get a row and normalize it.
 *
 * @param *x unpacked genotype row
 * @param *y normalized genotype row
 * @param *g geno struct
 * @param *M mean array
 * @param i row index
 */

void kjg_geno_get_normalized_row(uint8_t* x, double* y,
        SNP **snpmarkers, Indiv **indivmarkers,
        size_t numsnps, size_t numindivs,
        double* M, size_t i);

/**
 * Get multiple row and normalized rows
 *
 * @param *x unpacked genotype row (will hold last one)
 * @param *Y normalized genotype rows
 * @param *g geno struct
 * @param *M mean array
 * @param i row index
 * @param r number of rows to get
 * @return number of rows retrieved
 */


size_t kjg_geno_get_normalized_rows(uint8_t* x, double* Y,
        SNP **snpmarkers, Indiv **indivmarkers,
        size_t numsnps, size_t numindivs,
        double* M, size_t i, size_t r);

/**
 * Compute the mean genotype of all SNPs in the geno object
 *
 * @param *g geno object
 * @param *M array to store means
 */

void kjg_geno_row_means(SNP **snpmarkers, Indiv **indivmarkers,
        size_t numsnps, size_t numindivs,
        double* M);

#endif /* KJG_GENO_H_ */
