#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define MAXSIZE 2000000000
int *outlier, noutlier;
int K, L, NSAMPLES, nSNP;

int
main (int argc, char **argv)
{
  int k, x, m, n, i, nflags;
  double *xx, *iscase;
  double *iscasecorr, gamma, rowsum, rowsum1, *V, denom;
  double tempdouble;
  double chisq, Echisq;
  char Xchar;
  FILE *fp, *fppheno, *fppca, *fpout;
  char *INFILE = NULL;
  char *OUTFILE = NULL;
  char *PCAFILE = NULL;
  char *PHENOFILE = NULL;

  double compute_chisq (double *source, double *target);
  double compute_chisqE (double *source, double *target);

  /* set default values */
  L = 10;

  /* process flags */
  nflags = 0;
  while ((i = getopt (argc, argv, "i:j:p:l:o:")) != -1) {
    switch (i) {
    case 'i':			/* input genotype file */
      INFILE = (char *) strdup (optarg);
      nflags++;
      break;
    case 'j':			/* input phenotype file */
      PHENOFILE = (char *) strdup (optarg);
      nflags++;
      break;
    case 'p':			/* input phenotype file */
      PCAFILE = (char *) strdup (optarg);
      nflags++;
      break;
    case 'l':
      L = atoi (optarg);	/* number of principal components to correct */
      break;
    case 'o':			/* output file */
      OUTFILE = (char *) strdup (optarg);
      nflags++;
      break;
    }
  }

  if (nflags != 4) {
    fprintf (stderr, "Usage: -i -j -p -o flags must all be specified\n");
    exit (1);
  }

  /* open output file */
  if ((fpout = fopen (OUTFILE, "w")) == NULL) {
    fprintf (stderr, "Could not open output file %s\n", OUTFILE);
    exit (1);
  }

  /* print parameters */
  fprintf (fpout, "eigenstratQTL program run using parameters\n");
  fprintf (fpout, " -i %s\n", INFILE);
  fprintf (fpout, " -j %s\n", PHENOFILE);
  fprintf (fpout, " -p %s\n", PCAFILE);
  fprintf (fpout, " -l %d\n", L);
  fprintf (fpout, " -o %s\n", OUTFILE);
  fprintf (fpout, "\n");
  fprintf (fpout, "Chisq EIGENSTRAT\n");

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

  /* malloc */
  if ((V = (double *) malloc (NSAMPLES * NSAMPLES * sizeof (*V))) == NULL) {
    fprintf (stderr, "CM\n");
    exit (1);
  }
  if ((outlier = (int *) malloc (NSAMPLES * sizeof (*outlier))) == NULL) {
    fprintf (stderr, "CM\n");
    exit (1);
  }
  if ((iscase = (double *) malloc (NSAMPLES * sizeof (*iscase))) == NULL) {
    fprintf (stderr, "CM\n");
    exit (1);
  }
  if ((iscasecorr =
       (double *) malloc (NSAMPLES * sizeof (*iscasecorr))) == NULL) {
    fprintf (stderr, "CM\n");
    exit (1);
  }
  if ((xx = (double *) malloc (NSAMPLES * sizeof (*xx))) == NULL) {
    fprintf (stderr, "CM\n");
    exit (1);
  }

  /* FIRST, do PCAFILE */
  /* Build V[] and K and outlier[] */
  if ((fppca = fopen (PCAFILE, "r")) == NULL) {
    fprintf (stderr, "Could not open input file %s\n", PCAFILE);
    exit (1);
  }
  fscanf (fppca, "%d", &K);
  if (L > K) {
    fprintf (stderr, "OOPS l=%d is larger than k=%d in %s\n", L, K, PCAFILE);
    fprintf (fpout, "OOPS l=%d is larger than k=%d in %s\n", L, K, PCAFILE);
    exit (1);
  }
  for (x = 0; x < K; x++)
    fscanf (fppca, "%lf", &tempdouble);	/* eigenvalues */
  for (n = 0; n < NSAMPLES; n++) {
    for (k = 0; k < K; k++)
      fscanf (fppca, "%lf", &V[NSAMPLES * n + k]);
    if (feof (fppca)) {
      fprintf (stderr,
	       "OOPS: %s contains less than %d times %d entries\n",
	       PCAFILE, NSAMPLES, K);
      fprintf (fpout, "OOPS: %s contains less than %d times %d entries\n",
	       PCAFILE, NSAMPLES, K);
      exit (1);
    }
    /* check for outliers */
    outlier[n] = 1;
    for (k = 0; k < K; k++) {
      if (V[NSAMPLES * n + k] != 0.0)
	outlier[n] = 0;
    }
    if (outlier[n] == 1)
      noutlier++;
  }
  fscanf (fppca, "%lf", &tempdouble);
  if (!(feof (fppca))) {
    fprintf (stderr, "OOPS: %s contains too many entries\n", PCAFILE);
    fprintf (fpout, "OOPS: %s contains too many entries\n", PCAFILE);
    exit (1);
  }

  /* SECOND, do PHENOFILE */
  /* get phenotypes */
  if ((fppheno = fopen (PHENOFILE, "r")) == NULL) {
    fprintf (stderr, "Could not open input file %s\n", PHENOFILE);
    exit (1);
  }
  for (n = 0; n < NSAMPLES; n++) {
    if (feof (fppheno)) {
      fprintf (stderr, "OOPS: %s contains less than %d entries\n",
	       PHENOFILE, NSAMPLES);
      fprintf (fpout, "OOPS: %s contains less than %d entries\n",
	       PHENOFILE, NSAMPLES);
      exit (1);
    }
    fscanf (fppheno, "%lf", &iscase[n]);	/* QTL phenotype */
  }
  fscanf (fppheno, "%lf", &tempdouble);	/* to check # entries */
  if (!(feof (fppheno))) {
    fprintf (stderr, "OOPS: %s contains too many entries\n", PHENOFILE);
    fprintf (fpout, "OOPS: %s contains too many entries\n", PHENOFILE);
    exit (1);
  }
  /* mean-adjust iscase */
  rowsum = 0.0;
  rowsum1 = 0.0;
  for (n = 0; n < NSAMPLES; n++) {
    if ((outlier[n] == 1) || (iscase[n] == -100.0))
      continue;
    rowsum += iscase[n];
    rowsum1 += 1.0;
  }
  for (n = 0; n < NSAMPLES; n++) {
    if (outlier[n])
      continue;
    if (iscase[n] == -100.0)
      iscase[n] = -100.0;	/* still keep track */
    else
      iscase[n] -= rowsum / rowsum1;
  }
  /* make iscasecorr */
  for (n = 0; n < NSAMPLES; n++) {
    if (outlier[n] == 0)
      iscasecorr[n] = iscase[n];
  }
  for (k = 0; k < L; k++) {
    gamma = 0.0;
    denom = 0.0;
    for (n = 0; n < NSAMPLES; n++) {
      if ((outlier[n]) || (iscase[n] == -100.0))
	continue;
      gamma += iscasecorr[n] * V[NSAMPLES * n + k];
      denom += V[NSAMPLES * n + k] * V[NSAMPLES * n + k];
    }
    gamma /= denom;
    for (n = 0; n < NSAMPLES; n++) {
      if ((outlier[n]) || (iscase[n] == -100.0))
	continue;
      iscasecorr[n] -= gamma * V[NSAMPLES * n + k];
    }
  }

  /* THIRD, do INFILE */
  if ((fp = fopen (INFILE, "r")) == NULL) {
    fprintf (stderr, "Could not open input file %s\n", INFILE);
    exit (1);
  }
  m = 0;
  while (1) {			/* do EVERYTHING for SNP m */
    for (n = 0; n < NSAMPLES; n++) {
      fscanf (fp, "%c", &Xchar);
      if (Xchar == '0') {
	xx[n] = 0.0;
      }
      else if (Xchar == '1') {
	xx[n] = 0.5;
      }
      else if (Xchar == '2') {
	xx[n] = 1.0;
      }
      else if (Xchar == '9') {
	xx[n] = -100.0;
      }
      else if (!(feof (fp))) {
	fprintf (stderr, "OOPS bad char %c at m=%d n=%d\n", Xchar, m, n);
	fprintf (fpout, "OOPS bad char %c at m=%d n=%d\n", Xchar, m, n);
	exit (1);
      }
      if (outlier[n] == 1)
	xx[n] = -100.0;
    }
    if (feof (fp))
      break;
    fscanf (fp, "%c", &Xchar);	/* should be \n character */

    /* mean-adjust xx */
    rowsum = 0.0;
    rowsum1 = 0.0;
    for (n = 0; n < NSAMPLES; n++) {
      if ((outlier[n]) || (xx[n] == -100.0))
	continue;
      rowsum += xx[n];
      rowsum1 += 1.0;
    }
    for (n = 0; n < NSAMPLES; n++) {
      if (outlier[n])
	continue;
      if (xx[n] == -100.0)
	xx[n] = -100.0;		/* still keep track */
      else
	xx[n] -= rowsum / rowsum1;
    }

    /* Chisq */
    chisq = compute_chisq (xx, iscase);

    /* EIGENSTRAT */
    for (k = 0; k < L; k++) {
      gamma = 0.0;
      denom = 0.0;
      for (n = 0; n < NSAMPLES; n++) {
	if ((outlier[n]) || (xx[n] == -100.0))
	  continue;
	gamma += xx[n] * V[NSAMPLES * n + k];
	denom += V[NSAMPLES * n + k] * V[NSAMPLES * n + k];
      }
      gamma /= denom;
      for (n = 0; n < NSAMPLES; n++) {
	if ((outlier[n]) || (xx[n] == -100.0))
	  continue;
	xx[n] -= gamma * V[NSAMPLES * n + k];
      }
    }
    Echisq = compute_chisqE (xx, iscasecorr);

    if (rowsum1 == 0.0) {
      chisq = -1.0;
      Echisq = -1.0;
    }

    if (chisq >= 0.0)
      fprintf (fpout, "%.04f", chisq);
    else
      fprintf (fpout, "NA");
    if (Echisq >= 0.0)
      fprintf (fpout, " %.04f\n", Echisq);
    else
      fprintf (fpout, " NA\n");

    m++;
    if (NSAMPLES * m > MAXSIZE) {
      fprintf (stderr, "OOPS genotype file has > %d genotypes\n", MAXSIZE);
      fprintf (fpout, "OOPS genotype file has > %d genotypes\n", MAXSIZE);
      exit (1);
    }
  }
  fclose (fp);
  return 0;
}

double
compute_chisq (double *source, double *target)
{
  int n;
  double sum1, sumx, sumxx, sumy, sumyy, sumxy, numer, denom1, denom2;
  double corr;

  sum1 = 0.0;
  sumx = 0.0;
  sumxx = 0.0;
  sumy = 0.0;
  sumyy = 0.0;
  sumxy = 0.0;
  for (n = 0; n < NSAMPLES; n++) {
    if (outlier[n])
      continue;
    if (source[n] == -100.0)
      continue;
    if (target[n] == -100.0)
      continue;

    sumx += source[n];
    sumxx += source[n] * source[n];
    sumy += target[n];
    sumyy += target[n] * target[n];
    sumxy += source[n] * target[n];
    sum1 += 1.0;
  }
  if (sumxx == 0.0)
    return -1.0;
  if (sumyy == 0.0)
    return -1.0;
  numer = sumxy / sum1 - (sumx / sum1) * (sumy / sum1);
  denom1 = (sumxx / sum1 - (sumx / sum1) * (sumx / sum1));
  denom2 = (sumyy / sum1 - (sumy / sum1) * (sumy / sum1));
  if (denom1 <= 0.0)
    return -1.0;
  if (denom2 <= 0.0)
    return -1.0;

  corr = (numer / sqrt (denom1 * denom2));
  return (sum1 * corr * corr);
}

double
compute_chisqE (double *source, double *target)
{
  int n;
  double sum1, sumx, sumxx, sumy, sumyy, sumxy, numer, denom1, denom2;
  double corr;

  sum1 = 0.0;
  sumx = 0.0;
  sumxx = 0.0;
  sumy = 0.0;
  sumyy = 0.0;
  sumxy = 0.0;
  for (n = 0; n < NSAMPLES; n++) {
    if (outlier[n])
      continue;
    if (source[n] == -100.0)
      continue;
    if (target[n] == -100.0)
      continue;

    sumx += source[n];
    sumxx += source[n] * source[n];
    sumy += target[n];
    sumyy += target[n] * target[n];
    sumxy += source[n] * target[n];
    sum1 += 1.0;
  }
  if (sumxx == 0.0)
    return -1.0;
  if (sumyy == 0.0)
    return -1.0;
  numer = sumxy / sum1 - (sumx / sum1) * (sumy / sum1);
  denom1 = (sumxx / sum1 - (sumx / sum1) * (sumx / sum1));
  denom2 = (sumyy / sum1 - (sumy / sum1) * (sumy / sum1));
  if (denom1 <= 0.0)
    return -1.0;
  if (denom2 <= 0.0)
    return -1.0;

  corr = (numer / sqrt (denom1 * denom2));
  sum1 = sum1 - ((double) (L + 1));
  return (sum1 * corr * corr);
}
