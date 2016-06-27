#include <stdio.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

#include <lapacke.h>
#include "kjg_gsl.h"

void
kjg_gsl_matrix_fprintf (FILE * stream, gsl_matrix * m, const char *template)
{
  size_t i, j;
  for (i = 0; i < m->size1; i++)
    {
      fprintf (stream, template, gsl_matrix_get (m, i, 0));
      for (j = 1; j < m->size2; j++)
        {
          fprintf (stream, "\t");
          fprintf (stream, template, gsl_matrix_get (m, i, j));
        }
      fprintf (stream, "\n");
    }
}

void
kjg_gsl_matrix_fscanf (FILE * stream, gsl_matrix * m)
{
  size_t i, j;
  double x;
  for (i = 0; i < m->size1; i++)
    {
      for (j = 0; j < m->size2; j++)
        {
          fscanf (stream, "%lg", &x);
          gsl_matrix_set (m, i, j, x);
        }
    }
}

void
kjg_gsl_evec_fprintf (FILE * stream, gsl_vector * eval, gsl_matrix * evec,
                      const char *template)
{
  size_t i, j;
  fprintf (stream, "#");
  fprintf (stream, template, gsl_vector_get (eval, 0));
  for (i = 1; i < eval->size; i++)
    {
      fprintf (stream, "\t");
      fprintf (stream, template, gsl_vector_get (eval, i));
    }
  fprintf (stream, "\n");
  kjg_gsl_matrix_fprintf (stream, evec, template);
}

int
kjg_gsl_evec_fscanf (FILE * stream, gsl_vector * eval, gsl_matrix * evec)
{
  size_t i, j;
  int r;
  double x;

  r = fscanf (stream, "#%lg", &x);
  if (r != 1)
    return (r);
  gsl_vector_set (eval, 0, x);

  for (i = 1; i < eval->size; i++)
    {
      r = fscanf (stream, "%lg", &x);
      if (r != 1)
        return (r);
      gsl_vector_set (eval, i, x);
    }

  for (i = 0; i < evec->size1; i++)
    {
      for (j = 0; j < evec->size2; j++)
        {
          r = fscanf (stream, "%lg", &x);
          if (r != 1)
            return (r);
          gsl_matrix_set (evec, i, j, x);
        }
    }

  return (0);
}

gsl_rng *
kjg_gsl_rng_init ()
{
  const gsl_rng_type *T;
  gsl_rng *r;
  extern long seed;

  gsl_rng_env_setup ();

  gsl_rng_default_seed = seed;

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

//  fprintf (stderr, "generator type: %s\n", gsl_rng_name (r));
//  fprintf (stderr, "seed = %lu\n", gsl_rng_default_seed);

  return (r);
}

int
kjg_gsl_matrix_frobenius_normalize (gsl_matrix * m)
{
  double s = kjg_gsl_dlange ('F', m);
  double d = m->size1 * m->size2;
  return (gsl_matrix_scale (m, d / s));
}

double
kjg_gsl_dlange (const char norm, const gsl_matrix * m)
{
  return (LAPACKE_dlange (LAPACK_ROW_MAJOR, norm, m->size1, m->size2, m->data,
                          m->tda));
}

int
kjg_gsl_dgeqrf (gsl_matrix * m, gsl_vector * tau)
{
  return (LAPACKE_dgeqrf (LAPACK_ROW_MAJOR, m->size1, m->size2, m->data, m->tda,
                          tau->data));
}

int
kjg_gsl_dorgqr (gsl_matrix * m, gsl_vector * tau)
{
  return (LAPACKE_dorgqr (LAPACK_ROW_MAJOR, m->size2, m->size2, m->size2,
                          m->data, m->tda, tau->data));
}

void
kjg_gsl_ran_ugaussian_pair (const gsl_rng * r, double x[2])
{
  double r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */
      x[0] = -1 + 2 * gsl_rng_uniform_pos (r);
      x[1] = -1 + 2 * gsl_rng_uniform_pos (r);

      /* see if it is in the unit circle */
      r2 = x[0] * x[0] + x[1] * x[1];
    }
  while (r2 > 1.0 || r2 == 0);

  r2 = sqrt (-2.0 * log (r2) / r2);

  x[0] *= r2;
  x[1] *= r2;
}

void
kjg_gsl_ran_ugaussian_matrix (const gsl_rng * r, gsl_matrix * m)
{
  size_t i, j;
  double *data;
  double x, y, r2;

  for (i = 0; i < m->size1; i++)
    {
      data = gsl_matrix_ptr (m, i, 0);

      for (j = 0; j < m->size2 - 1; j += 2)
        {
          kjg_gsl_ran_ugaussian_pair (r, data);
          data += 2;
        }

      if (m->size2 % 2)
        *data = gsl_rng_uniform_pos (r);
    }
}

void
kjg_gsl_matrix_QR (gsl_matrix * m)
{
  gsl_vector *tau = gsl_vector_alloc (m->size2);
  kjg_gsl_dgeqrf (m, tau);
  kjg_gsl_dorgqr (m, tau);
  gsl_vector_free (tau);
}

int
kjg_gsl_SVD (gsl_matrix * M, gsl_matrix * V, gsl_vector * S)
{
  size_t big_enough = M->size1 + V->size2;
  double *superb = malloc (big_enough * sizeof(double));
  double *U;
  int info = LAPACKE_dgesvd (
      LAPACK_ROW_MAJOR,	// row major
      'O', 'S', M->size1, M->size2, M->data, M->tda, S->data, U, big_enough,
      V->data, V->tda, superb);
  free (superb);
  return (info);
}
