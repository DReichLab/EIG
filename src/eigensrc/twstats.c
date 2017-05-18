#include  <stdio.h>
#include <string.h>
#include <math.h>

#include <nicklib.h>
#include <getpars.h>

int verbose = NO;
double nval = -1;

char *iname = NULL;
char *parname = NULL;
char *oname = NULL;

char *twxtab = NULL;

void readcommands (int argc, char **argv);

#define VERSION    "1000"

int minleneig = 10;


int
main (int argc, char **argv)
{
  FILE *ofile;
  int nlambda = 0;
  int i, m;
  double zn, zvar, tw, tail;
  double *xx[0], *lambda;

  readcommands (argc, argv);
  settwxtable (twxtab);

  if (oname == NULL)
    ofile = stdout;
  else
    openit (oname, &ofile, "w");

  if (iname == NULL)
    fatalx ("i paraameter compulsory\n");
  nlambda = numlines (iname);
  ZALLOC (lambda, nlambda, double);
  xx[0] = lambda;
  nlambda = getxx (xx, nlambda, 1, iname);
  vst (lambda, lambda, -1.0, nlambda);
  sortit (lambda, NULL, nlambda);
  vst (lambda, lambda, -1.0, nlambda);
  m = numgtz (lambda, nlambda);

  fprintf (ofile, "%4s  %12s", "#N", "eigenvalue");
  fprintf (ofile, "%12s", "difference");
  fprintf (ofile, " %9s %12s", "twstat", "p-value");
  fprintf (ofile, " %9s", "effect. n");
  fprintf (ofile, "\n");

  for (i = 0; i < m; ++i) {

    zn = nval;
    tail = dotwcalc (lambda + i, m - i, &tw, &zn, &zvar, minleneig);
    fprintf (ofile, "%4d  %12.6f", i + 1, lambda[i]);
    if (i == 0)
      fprintf (ofile, "%12s", "NA");
    else
      fprintf (ofile, "%12.6f", lambda[i] - lambda[i - 1]);
    if (tail >= 0.0)
      fprintf (ofile, " %9.3f %12.6g", tw, tail);
    else
      fprintf (ofile, " %9s %12s", "NA", "NA");
    if (zn > 0.0) {
      fprintf (ofile, " %9.3f", zn);
    }
    else {
      fprintf (ofile, " %9s", "NA");
    }
    fprintf (ofile, "\n");
  }
  return 0;
}

void
readcommands (int argc, char **argv)
{

  int i;
  char *parname = NULL;
  phandle *ph;

  while ((i = getopt (argc, argv, "i:o:p:n:m:t:V")) != -1) {

    switch (i) {

    case 'i':
      iname = strdup (optarg);
      break;

    case 'o':
      oname = strdup (optarg);
      break;

    case 't':
      twxtab = strdup (optarg);
      break;

    case 'n':
      nval = atof (optarg);
      break;

    case 'm':
      minleneig = atoi (optarg);
      break;

    case 'p':
      parname = strdup (optarg);
      break;

    case 'V':
      verbose = YES;
      break;

    case '?':
      printf ("Usage: bad params.... \n");
      fatalx ("bad params\n");
    }
  }

  if (parname == NULL)
    return;

  printf ("parameter file: %s\n", parname);
  ph = openpars (parname);

  getstring (ph, "input:", &iname);
  getstring (ph, "output:", &oname);
  getdbl (ph, "nval:", &nval);
  getint (ph, "minleneig:", &minleneig);

  writepars (ph);
  closepars (ph);

}
