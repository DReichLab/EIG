#include <stdio.h>
#include <nicklib.h>
#include <stdlib.h>
#include "eigsubs.h"
#include <unistd.h>

int
main ()
{
  int NSAMPLES, n, k;
  double *eval, *evec, *XTX;

  NSAMPLES = 2;

  /* malloc */
  if ((eval = (double *) malloc (NSAMPLES * sizeof(*eval))) == NULL)
    {
      fprintf (stderr, "CM\n");
      exit (1);
    }
  if ((evec = (double *) malloc (NSAMPLES * NSAMPLES * sizeof(*evec))) == NULL)
    {
      fprintf (stderr, "CM\n");
      exit (1);
    }
  if ((XTX = (double *) malloc (NSAMPLES * NSAMPLES * sizeof(*XTX))) == NULL)
    {
      fprintf (stderr, "CM\n");
      exit (1);
    }

  XTX[0] = 1;
  XTX[1] = 0;
  XTX[2] = 0;
  XTX[3] = 1; /* 2x2 identity matrix */

  eigvecs (XTX, eval, evec, NSAMPLES); /* eigenvector k is evec[k*NSAMPLES+n] */

  /* print eval and evec */
  printf ("The eigenvectors of the 2x2 identity matrix are:\n");
  for (n = 0; n < NSAMPLES; n++)
    {
      for (k = 0; k < NSAMPLES; k++)
        {
          printf (" ");
          printf ("%.02f", evec[k * NSAMPLES + n]);
        }
      printf ("\n");
    }
  return 0;
}
