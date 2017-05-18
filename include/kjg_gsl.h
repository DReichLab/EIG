/**
 * @file kjg_gsl.h
 * @brief Augment GSL functions
 */

#ifndef KJG_GSL_H_
#define KJG_GSL_H_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

/**
 * Prints the matrix tab-delimited
 * @param *stream output file pointer
 * @param *m gsl_matrix to print
 * @param *template character template for fprintf
 */

void kjg_gsl_matrix_fprintf (FILE * stream, gsl_matrix * m,
			     const char *template);

/**
 * Prints the eigenvalues and then eigenvectors below
 * @param *stream output file pointer
 * @param *eval eigenvalues
 * @param *evec eigenvectors
 * @param *template character template for fprintf */

void kjg_gsl_evec_fprintf (FILE * stream,
			   gsl_vector * eval,
			   gsl_matrix * evec, const char *template);

/**
 * Reads a matrix
 * @param *stream input file pointer
 * @param *m matrix to store
 */

void kjg_gsl_matrix_fscanf (FILE * stream, gsl_matrix * m);

/**
 * Reads an evec
 * @param *stream input file pointer
 * @param *eval eigenvalues vector
 * @param *evec eigenvectors matrix
 */

int kjg_gsl_evec_fscanf (FILE * stream, gsl_vector * eval, gsl_matrix * evec);

/**
 * Initializes random number generation.
 */

gsl_rng *kjg_gsl_rng_init ();

/**
 * Initializes the matrix with random unit gaussians
 * @param *m matrix to be set
 * @param *r random number generator
 */

void kjg_gsl_ran_ugaussian_pair (const gsl_rng * r, double x[2]);

/** Fills a matrix with unit Gaussian random variates
 * @param *r random number generator
 * @param *m matrix to be filled
 */

void kjg_gsl_ran_ugaussian_matrix (const gsl_rng * r, gsl_matrix * m);

/**
 * Normalizes the matrix so the Frobenius norm is M*N
 * @param *m matrix to normalize
 * @return if error
 */

int kjg_gsl_matrix_frobenius_normalize (gsl_matrix * m);

/**
 * Calculates the norm of a matrix
 * @param norm type of norm to return, see lapack dlange
 * @param *m matrix to find norm of
 * @return norm
 */

double kjg_gsl_dlange (const char norm, const gsl_matrix * m);

/**
 * Performs the QR decomposition on the matrix and return Q in the matrix
 * @param *m matrix to orthogonalize
 */

void kjg_gsl_matrix_QR (gsl_matrix * m);

/**
 * Calls LAPACK dgeqrf and return R and compacted Q matrix
 * @param *m input matrix
 * @param *tau see LAPACK documentation
 * @return LAPACK return
 */

int kjg_gsl_dgeqrf (gsl_matrix * m, gsl_vector * tau);

/**
 * Calls LAPACK dorgqr to extract Q matrix
 * @param *m matrix with compacted Q and will store unpacked Q
 * @param *tau see LAPACK documentation
 * @return LAPACK return
 */
int kjg_gsl_dorgqr (gsl_matrix * m, gsl_vector * tau);

/**
 * Calls LAPACK dgesvd, keeping u (in m) and s, discarding v^T
 * @param *m input matrix / where u is stored
 * @param *s entries of the diagonal matrix
 * @return LAPACK return
 */
int kjg_gsl_SVD (gsl_matrix * M, gsl_matrix * V, gsl_vector * S);

#endif /* KJG_GSL_H_ */
