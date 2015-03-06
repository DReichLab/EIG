/** @file kjg_fpca.h
 * @brief Runs fastPCA.
 * This module also has methods to multiply a genotype matrix against the GSL
 * matrices.
 */

#ifndef KJG_FPCA_H_
#define KJG_FPCA_H_

#include <gsl/gsl_matrix.h>

extern size_t KJG_FPCA_ROWS;	// number of rows to process at once

/** Performs a fast PCA
 * @param *eval eigenvalues
 * @param *evec eigenvectors
 * @param K number of eigenvalues/vectors
 * @param L width of projection matrix
 * @param I iterations to do exponentiation
 */

void kjg_fpca (size_t K, size_t L, size_t I, double *eval, double *evec);

/** Multiplies B=X*A1 and A2 = XT*B = XT*X*A1
 * @param *A1 some matrix
 * @param *B intermediate matrix
 * @param *A2 next matrix
 */

void kjg_fpca_XTXA (const gsl_matrix * A1, gsl_matrix * B, gsl_matrix * A2);

/** Multiplies B = X*A
 * @param *A some matrix
 * @param *B another matrix
 */

void kjg_fpca_XA (const gsl_matrix * A, gsl_matrix * B);

/** Multiplies A = XT*B
 * @param *B some matrix
 * @param *A another matrix
 */

void kjg_fpca_XTB (const gsl_matrix * B, gsl_matrix * A);

#endif /* KJG_FPCA_H_ */
