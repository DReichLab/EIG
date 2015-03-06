/*
 * kjg_gsl.h
 *
 *  Created on: Aug 1, 2013
 *      Author: kjg063
 */

#ifndef KJG_GSL_H_
#define KJG_GSL_H_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

/**
 * Prints the matrix tab-delimited
 *
 * @param *stream file pointer to print output
 * @param *m gsl_matrix to print
 * @param *template character template for fprintf
 */

void kjg_gsl_matrix_fprintf(FILE* stream, gsl_matrix* m, char* template);


/**
 * Print the eigenvalues and then eigenvectors below
 *
 * @param *stream file pointer to print output
 * @param *eval eigenvalues
 * @param *evec eigenvectors
 * @param *template character template for fprintf */

void kjg_gsl_evec_fprintf(FILE* stream, gsl_vector* eval, gsl_matrix* evec, char* template);

void kjg_gsl_matrix_fscanf(FILE* stream, gsl_matrix* m);
int kjg_gsl_evec_fscanf(FILE* stream, gsl_vector* eval, gsl_matrix* evec);

/**
 * Initialize random number generation.
 */

gsl_rng *kjg_gsl_rng_init();

/**
 * Initialize the matrix with random unit gaussians
 *
 * @param *m matrix to be set
 * @param *r random number generator
 */

void kjg_gsl_matrix_set_ran_ugaussian(gsl_matrix* m, gsl_rng* r);

/**
 * Normalize the matrix so the frobenius norm is M*N
 *
 * @param *m matrix to normalize
 * @return if error
 */

int kjg_gsl_matrix_frobenius_normalize(gsl_matrix* m);

/**
 * Calculate the norm of a matrix
 *
 * @param norm type of norm to return, see lapack dlange
 * @param *m matrix to find norm of
 * @return norm
 */

double kjg_gsl_dlange(char norm, gsl_matrix* m);

#endif /* KJG_GSL_H_ */
