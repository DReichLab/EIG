#include "qpsubs.h"
#include "eigsubs.h"
#include "smartsubs.h"
extern int fancynorm, verbose, plotmode, outnum;

extern FILE *fstdetails ; 

// static Indiv **indm ;
// static void wjackestx(double *est, double *sig, double mean, double *jmean, double *jwt, int g)  ;
// static void wjackvestx(double *vest, double *var, int d, double *mean, double **jmean, double *jwt, int g)  ;
static int outliermode = 0;

void
setoutliermode (int mode)
{
  outliermode = mode;
}

int
ridoutlierx (double *evecs, int n, int neigs,
	    double thresh, int *badlist, OUTLINFO ** outinfo)
{
/* badlist contains list of outliers */
  double *ww, *w2, y1, y2, yy, zz;
  int *vbad;
  int i, j;
  int nbad = 0;
  OUTLINFO *outpt;

  if (outliermode > 1)
    return 0;
  if (n < 3)
    return 0;
  ZALLOC (ww, n, double);
  ZALLOC (vbad, n, int);
  for (j = 0; j < n; j++) {
    outpt = outinfo[j];
    outpt->vecno = -1;
  }
  for (i = 0; i < neigs; ++i) {
    copyarr (evecs + i * n, ww, n);
    if (outliermode == 0) {
      y1 = asum (ww, n) / (double) n;
      vsp (ww, ww, -y1, n);
      y2 = asum2 (ww, n) / (double) n;
      y2 = sqrt (y2);
      vst (ww, ww, 1.0 / y2, n);

      for (j = 0; j < n; j++) {
	if (fabs (ww[j]) > thresh) {
	  vbad[j] = 1;
	  outpt = outinfo[j];
	  if (outpt->vecno < 0) {
	    outpt->vecno = i;
	    outpt->score = ww[j];
	  }
	}
      }
    }
    if (outliermode == 1) {
      ZALLOC (w2, n, double);
      for (j = 0; j < n; j++) {
	yy = ww[j];
	ww[j] = 0;
	y1 = asum (ww, n) / (double) (n - 1);
	vsp (w2, ww, -y1, n);
	w2[j] = 0;
	y2 = asum2 (w2, n) / (double) n;
	y2 = sqrt (y2);
	zz = yy - y1;
	zz /= y2;
	if (fabs (zz) > thresh) {
	  vbad[j] = 1;
	  outpt = outinfo[j];
	  if (outpt->vecno < 0) {
	    outpt->vecno = i;
	    outpt->score = zz;
	  }
	}
	ww[j] = yy;
      }
      free (w2);
    }
  }
  for (j = 0; j < n; j++) {
    if (vbad[j] == 1) {
      badlist[nbad] = j;
      ++nbad;
    }
  }
  free (ww);
  free (vbad);
  return nbad;

}

double
dofstnumx_eig (double *fst, double *fstest, double *fstsig, int *qfstnum, SNP ** xsnplist,
           int *xindex, int *xtypes, int nrows, int ncols, int numeg,
           int nblocks, Indiv ** indivmarkers, int fstmode)
// fstmode is classic mode (smartpca)
// fstmode 2  is fstdmode
{

  int t1, t2;
  int a, b;
  int c1[2], c2[2], *cc;
  int *rawcol, *popall, *pop0, *pop1;
  int t, k, g, i, col, j;
  double ya, yb, y, jest, jsig, mean;
  SNP *cupt;
  double *top, *bot, *djack, *wjack, *gtop, *gbot, *wbot, *wtop;
  double **btop, **bbot, wt;
  double *w1, *w2, *w3;
  double ytop, ybot;
  double y1, y2, yscal;
  int bnum;
  int nloop = 0, fstdnum = 0;
  double *ztop, *zbot, qtop, qbot;
  int *fstnum ; 
  char **eglist;
  Indiv **indm ;

  indm = indivmarkers;
  if ((fstdetails != NULL) && (indm == NULL))
    fatalx ("bug in dofstnumx\n");

  ZALLOC (eglist, numeg, char *);
  for (k = 0; k < nrows; ++k) {
    if (indm == NULL)
      break;
    j = xtypes[k];
    if (j < 0)
      continue;
    if (j >= numeg)
      continue;
    t = xindex[k];
    eglist[j] = indm[t]->egroup;
  }

  ZALLOC (w1, numeg * numeg, double);
  ZALLOC (w2, numeg * numeg, double);
  ZALLOC (w3, numeg * numeg, double);
  ZALLOC (gtop, numeg * numeg, double);
  ZALLOC (gbot, numeg * numeg, double);
  ZALLOC (wtop, numeg * numeg, double);
  ZALLOC (wbot, numeg * numeg, double);
  ZALLOC (djack, nblocks, double);
  ZALLOC (wjack, nblocks, double);
  ZALLOC (ztop, numeg * numeg, double);
  ZALLOC (zbot, numeg * numeg, double);
  ZALLOC (fstnum, numeg * numeg, int);
  btop = initarray_2Ddouble (nblocks, numeg * numeg, 0.0);
  bbot = initarray_2Ddouble (nblocks, numeg * numeg, 0.0);

  if (nblocks == 1)
    printf ("number of blocks 1: no standard error\n");

  vzero (fst, numeg * numeg);
  vzero (fstest, numeg * numeg);
  vzero (fstsig, numeg * numeg);
  ivzero (fstnum, numeg * numeg);


  for (col = 0; col < ncols; ++col) {
    cupt = xsnplist[col];
    if (cupt->ignore)
      continue;
    wt = cupt->weight;
    if (wt <= 0.0)
      continue;
    bnum = cupt->tagnumber;
    if (bnum < 0)
      continue;
    ++wjack[bnum];
    top = btop[bnum];
    bot = bbot[bnum];

    fstcolyy (ztop, zbot, cupt, xindex, xtypes, nrows, numeg);

    for (a = 0; a < numeg; a++) {
      for (b = a + 1; b < numeg; b++) {
        k = a * numeg + b;
        ytop = ztop[k];
        ybot = zbot[k];
        if (fstdetails != NULL) {
          if (fstdnum == 0) {
            fprintf (fstdetails, "%20s ", "## pop 1");
            fprintf (fstdetails, "%20s ", "pop 2");
            fprintf (fstdetails, "%20s ", "snpname");
            fprintf (fstdetails, "%20s ", "N");
            fprintf (fstdetails, "%20s ", "D");
            fprintf (fstdetails, "%12s ", "Ratio");
            fprintf (fstdetails, "\n");
          }
          fprintf (fstdetails, "%20s ", eglist[a]);
          fprintf (fstdetails, "%20s ", eglist[b]);
          fprintf (fstdetails, "%20s ", cupt->ID);
          fprintf (fstdetails, "%12.6f ", ytop);
          fprintf (fstdetails, "%12.6f ", ybot);
          if (ybot > 0.0)
            fprintf (fstdetails, "%12.6f", ytop / ybot);
          else
            fprintf (fstdetails, "%12s", "-");
          fprintf (fstdetails, "\n");
          ++fstdnum;
        }


        if (ybot < 0.0)
          continue;

        ++fstnum[k] ;

        if (fstmode == NO) {
          top[k] += wt * ytop;
          bot[k] += 1.0;
        }

        if (fstmode == YES) {
          top[k] += ytop;
          bot[k] += ybot;
        }

        if (fstmode == 2) {
          top[k] += ytop;
          bot[k] += 1.0 / wt;
        }

        w1[k] += ytop;
        w2[k] += ybot;
// classic fst estimate

      }
    }
  }
// symmetrize
  for (a = 0; a < numeg; a++) {
    for (b = a + 1; b < numeg; b++) {
      top[b * numeg + a] = top[a * numeg + b];
      bot[b * numeg + a] = bot[a * numeg + b];
      w1[b * numeg + a] = w1[a * numeg + b];
      w2[b * numeg + a] = w2[a * numeg + b];
      fstnum[b*numeg+a] = fstnum[a*numeg+b] ;
    }
  }


  vsp (w2, w2, 1.0e-10, numeg * numeg);
  vvd (fst, w1, w2, numeg * numeg);


  for (k = 0; k < nblocks; k++) {
    top = btop[k];
    bot = bbot[k];
    vvp (gtop, gtop, top, numeg * numeg);
    vvp (gbot, gbot, bot, numeg * numeg);
  }

  for (k = 0; k < nblocks; k++) {
    top = btop[k];
    bot = bbot[k];
    vvm (wtop, gtop, top, numeg * numeg);
    vvm (wbot, gbot, bot, numeg * numeg);
    vsp (wbot, wbot, 1.0e-10, numeg * numeg);
    vvd (top, wtop, wbot, numeg * numeg);       // delete-block estimate
  }
  vsp (gbot, gbot, 1.0e-10, numeg * numeg);
  vvd (gtop, gtop, gbot, numeg * numeg);

  copyarr(gtop, fst, numeg*numeg) ;

  for (i = 0; i < numeg; i++) {
    for (j = i + 1; j < numeg; j++) {
      for (k = 0; k < nblocks; k++) {
        top = btop[k];
        djack[k] = top[i * numeg + j];
      }

      ++nloop;
      mean = gtop[i * numeg + j];
      jest = mean;
      jsig = 0;
      if (nblocks > 1) {
        wjackest (&jest, &jsig, mean, djack, wjack, nblocks);
      }
      fstest[i * numeg + j] = fstest[j * numeg + i] = jest;
      fstsig[i * numeg + j] = fstsig[j * numeg + i] = jsig;

      if (nloop == -1) {
        printf ("fstest\n");
        printf ("mean: %9.3f\n", mean);
        printmat (djack, 1, nblocks);
        printnl ();
        printmat (wjack, 1, nblocks);
        printf ("%9.3f %9.3f\n", jest, jsig);
      }
    }
  }


/**
   printf("fst:\n") ;
   printmat(fst, numeg, numeg) ;
   printnl() ;
   printmat(fstest, numeg, numeg) ;
   printnl() ;
   printmat(fstsig, numeg, numeg) ;
   printnl() ;
*/

  yscal = 1.0;
  if (fstmode != YES) {
    copyarr (fstsig, w3, numeg * numeg);
    vsp (w3, w3, 1.0e-10, numeg * numeg);
    vvd (w1, fst, w3, numeg * numeg);
    vvd (w2, fstest, w3, numeg * numeg);
//   now do regression  w1 = yscal * w2
    y1 = vdot (w1, w2, numeg * numeg);
    y2 = vdot (w2, w2, numeg * numeg);
    yscal = y1 / y2;
    vst (fstest, fstest, yscal, numeg * numeg);
    vst (fstsig, fstsig, yscal, numeg * numeg);
  }

  if (qfstnum != NULL) copyiarr(fstnum, qfstnum, numeg*numeg) ;
 
  free(fstnum) ;
  free (eglist);
  free (w1);
  free (w2);
  free (w3);

  free (gbot);
  free (wtop);
  free (wbot);
  free (ztop);
  free (zbot);
  free (djack);
  free (wjack);

  free2D (&btop, nblocks);
  free2D (&bbot, nblocks);

  return yscal;

}

double
dofstnum_eig (double *fst, double *fstest, double *fstsig, SNP ** xsnplist,
          int *xindex, int *xtypes, int nrows, int ncols, int numeg,
          int nblocks)
{

  return dofstnumx_eig (fst, fstest, fstsig, NULL, xsnplist, xindex, xtypes, nrows, ncols,
             numeg, nblocks, NULL, NO);

}

