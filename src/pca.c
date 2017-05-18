#include <stdio.h>
#include <nicklib.h>
#include <stdlib.h>
#include "eigsubs.h"
#include <unistd.h>

int MAXITER, K, TOPK, NSAMPLES, nSNP;
double SIGMATHRESH;

int
main (int argc, char **argv)
{
  int k, n, m, nn, rowvalid, *outlier, i;
  int iter, nonewoutliers, nflags;
  char Xchar;
  double *X, *XTX, rowsum, rowmean = 0, rowmeanbayes = 0, *sigmaoutlier;
  double *eval, *evec, sum, summ, sum1, mean, sdev, sigma;
  FILE *fp, *fpout, *fplog, *fpeval;
  char *INFILE = NULL;
  char *OUTFILE = NULL;
  char *EVALFILE = NULL;
  char *LOGFILE = NULL;

  /* set default values */
  K = 10;
  MAXITER = 5;
  TOPK = 10;
  SIGMATHRESH = 6.0;

  /* process flags */
  nflags = 0;
  while ((i = getopt (argc, argv, "i:k:o:e:l:m:t:s:")) != -1) {
    switch (i) {
    case 'i':                  /* input file */
      INFILE = (char *) strdup (optarg);
      nflags++;
      break;
    case 'k':
      K = atoi (optarg);        /* number of principal components to output */
      break;
    case 'o':                  /* output file */
      OUTFILE = (char *) strdup (optarg);
      nflags++;
      break;
    case 'e':                  /* output eval file */
      EVALFILE = (char *) strdup (optarg);
      nflags++;
      break;
    case 'l':                  /* log file */
      LOGFILE = (char *) strdup (optarg);
      nflags++;
      break;
    case 'm':
      MAXITER = atoi (optarg);  /* max # of outlier removal iterations */
      break;
    case 't':
      TOPK = atoi (optarg);     /* # of PCs along which to remove outliers */
      break;
    case 's':
      SIGMATHRESH = atof (optarg);      /* # sdev to declare as outlier */
      break;
    }
  }
  if (nflags != 4) {
    fprintf (stderr, "Usage: -i -o -e -l flags must all be specified\n");
    exit (1);
  }

  /* open output files */
  if ((fpout = fopen (OUTFILE, "w")) == NULL) {
    fprintf (stderr, "Could not open output file %s\n", OUTFILE);
    exit (1);
  }
  if ((fpeval = fopen (EVALFILE, "w")) == NULL) {
    fprintf (stderr, "Could not open output file %s\n", OUTFILE);
    exit (1);
  }
  if ((fplog = fopen (LOGFILE, "w")) == NULL) {
    fprintf (stderr, "Could not open input file %s\n", LOGFILE);
    exit (1);
  }

  /* print parameters */
  fprintf (fplog, "pca program run using parameters\n");
  fprintf (fplog, " -i %s\n", INFILE);
  fprintf (fplog, " -k %d\n", K);
  fprintf (fplog, " -o %s\n", OUTFILE);
  fprintf (fplog, " -e %s\n", EVALFILE);
  fprintf (fplog, " -l %s\n", LOGFILE);
  fprintf (fplog, " -m %d\n", MAXITER);
  fprintf (fplog, " -t %d\n", TOPK);
  fprintf (fplog, " -s %.03f\n", SIGMATHRESH);
  fprintf (fplog, "\n");

  /* Determine NSAMPLES */
  if ((fp = fopen (INFILE, "r")) == NULL) {
    fprintf (stderr, "Could not open input file %s\n", INFILE);
    exit (1);
  }
  n = 0;
  while (1) {
    fscanf (fp, "%c", &Xchar);
    if (Xchar == '\n')
      break;
    n++;
  }
  NSAMPLES = n;
  fclose (fp);
  if (K > NSAMPLES - 1) {
    fprintf (stderr, "OOPS k=%d is too large for only %d samples\n", K,
             NSAMPLES);
    fprintf (fplog, "OOPS k=%d is too large for only %d samples\n", K,
             NSAMPLES);
    exit (1);
  }

  /* malloc */
  if ((eval = (double *) malloc (NSAMPLES * sizeof (*eval))) == NULL) {
    fprintf (stderr, "CM\n");
    exit (1);
  }
  if ((evec =
       (double *) malloc (NSAMPLES * NSAMPLES * sizeof (*evec))) == NULL) {
    fprintf (stderr, "CM\n");
    exit (1);
  }
  if ((outlier = (int *) malloc (NSAMPLES * sizeof (*outlier))) == NULL) {
    fprintf (stderr, "CM\n");
    exit (1);
  }
  if ((sigmaoutlier =
       (double *) malloc (NSAMPLES * sizeof (*sigmaoutlier))) == NULL) {
    fprintf (stderr, "CM\n");
    exit (1);
  }
  if ((X = (double *) malloc (NSAMPLES * sizeof (*X))) == NULL) {
    fprintf (stderr, "CM\n");
    exit (1);
  }
  if ((XTX = (double *) malloc (NSAMPLES * NSAMPLES * sizeof (*X))) == NULL) {
    fprintf (stderr, "CM\n");
    exit (1);
  }

  nonewoutliers = 0;
  for (n = 0; n < NSAMPLES; n++)
    outlier[n] = 0;
  iter = 0;
  while (nonewoutliers == 0) {
    for (n = 0; n < NSAMPLES; n++)
      sigmaoutlier[n] = 0.0;
    /* initialize XTX */
    for (n = 0; n < NSAMPLES; n++) {
      for (nn = 0; nn < NSAMPLES; nn++)
        XTX[NSAMPLES * n + nn] = 0.0;
    }

    nonewoutliers = 1;
    /* get data */
    if ((fp = fopen (INFILE, "r")) == NULL) {
      fprintf (stderr, "Could not open input file %s\n", INFILE);
      exit (1);
    }
    m = 0;
    while (1) {                 /* do EVERYTHING for SNP m */
      for (n = 0; n < NSAMPLES; n++) {
        fscanf (fp, "%c", &Xchar);
        if (Xchar == '0') {
          X[n] = 0.0;
        }
        else if (Xchar == '1') {
          X[n] = 0.5;
        }
        else if (Xchar == '2') {
          X[n] = 1.0;
        }
        else if (Xchar == '9') {
          X[n] = -100.0;
        }
        else if (!(feof (fp))) {
          fprintf (stderr, "OOPS bad char %c at m=%d n=%d\n", Xchar, m, n);
          fprintf (fplog, "OOPS bad char %c at m=%d n=%d\n", Xchar, m, n);
          exit (1);
        }
        if (outlier[n] == 1)
          X[n] = -100.0;
      }
      if (feof (fp))
        break;
      fscanf (fp, "%c", &Xchar);        /* should be \n character */

      /* mean-adjust this SNP */
      rowvalid = 0;
      rowsum = 0.0;
      for (n = 0; n < NSAMPLES; n++) {
        if (X[n] >= -99.0) {
          rowvalid++;
          rowsum += X[n];
        }
      }
      if (rowvalid > 0) {
        rowmean = (rowsum) / ((double) (rowvalid));
        rowmeanbayes = (rowsum + 0.5) / ((double) (1 + rowvalid));
      }
      for (n = 0; n < NSAMPLES; n++) {
        if (X[n] >= -99.0) {
          X[n] -= rowmean;
          X[n] /= sqrt (rowmeanbayes * (1.0 - rowmeanbayes));
        }
        else
          X[n] = 0.0;
      }

      /* update XTX */
      for (n = 0; n < NSAMPLES; n++) {
        for (nn = n; nn < NSAMPLES; nn++)
          XTX[NSAMPLES * n + nn] += X[n] * X[nn];
      }
      m++;
    }
    nSNP = m;
    if (K > nSNP - 1) {
      fprintf (stderr, "OOPS k=%d is too large for only %d SNPs\n", K, nSNP);
      fprintf (fplog, "OOPS k=%d is too large for only %d SNPs\n", K, nSNP);
      exit (1);
    }
    if (iter == 0) {
      fprintf (fplog, "nSNP=%d NSAMPLES=%d\n", nSNP, NSAMPLES);
    }

    /* complete XTX */
    for (n = 0; n < NSAMPLES; n++) {
      for (nn = n; nn < NSAMPLES; nn++)
        XTX[NSAMPLES * n + nn] /= ((double) nSNP);
    }
    for (n = 0; n < NSAMPLES; n++) {
      for (nn = 0; nn < n; nn++)
        XTX[NSAMPLES * n + nn] = XTX[NSAMPLES * nn + n];
    }

    /* do eigenanalysis */
    eigvecs (XTX, eval, evec, NSAMPLES);        /* eigenvector k is evec[k*NSAMPLES+n] */

    if (iter == MAXITER)
      break;                    /* no need to look for outliers */

    /* find outliers */
    for (k = 0; k < TOPK; k++) {
      sum = 0.0;
      summ = 0.0;
      sum1 = 0.0;
      for (n = 0; n < NSAMPLES; n++) {
        if (outlier[n] == 1)
          continue;
        sum += evec[k * NSAMPLES + n];
        summ += evec[k * NSAMPLES + n] * evec[k * NSAMPLES + n];
        sum1 += 1.0;
      }
      mean = sum / sum1;
      sdev = sqrt (summ / sum1 - mean * mean);
      for (n = 0; n < NSAMPLES; n++) {
        if (outlier[n] == 1)
          continue;
        sigma = (evec[k * NSAMPLES + n] - mean) / sdev;
        if (sigma < 0)
          sigma = -sigma;
        if (sigma > SIGMATHRESH) {
          if (sigma > sigmaoutlier[n])
            sigmaoutlier[n] = sigma;
          nonewoutliers = 0;
        }
      }
    }
    fprintf (fplog, "Outlier removal iteration %d:\n", iter);
    if (nonewoutliers)
      fprintf (fplog, "  no outliers detected\n");
    for (n = 0; n < NSAMPLES; n++) {
      if (sigmaoutlier[n] > 0.0) {
        fprintf (fplog,
                 "  removed outlier individual %d (%.02f sigma)\n", n,
                 sigmaoutlier[n]);
        outlier[n] = 1;
      }
    }
    iter++;
    fclose (fp);
  }

  /* print eval and evec */
  for (k = 0; k < NSAMPLES; k++)
    fprintf (fpeval, "%.06f\n", eval[k]);
  fprintf (fpout, "%d\n", K);
  for (k = 0; k < K; k++)
    fprintf (fpout, "%.04f\n", eval[k]);
  for (n = 0; n < NSAMPLES; n++) {
    for (k = 0; k < K; k++) {
      fprintf (fpout, " ");
      if (evec[k * NSAMPLES + n] > 0)
        fprintf (fpout, " ");
      fprintf (fpout, "%.04f", evec[k * NSAMPLES + n]);
    }
    fprintf (fpout, "\n");
  }
  return 0;
}
