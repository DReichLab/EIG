#include "regsubs.h" 

extern int verbose;
void
squishx (double *xmat, double *mat, int nrow, int oldc, int *cols, int newc);

double
regressit (double *ans, double *eq, double *rhs, int m, int n)
{
  double *co, *rr, *ww, *w2;
  double vres, vbase, ynum, y, trace;
  double *traneq;
  int i, j, k, ret;

  ZALLOC(co, n*n, double);
  ZALLOC(rr, n, double);
  ZALLOC(ww, m, double);
  ZALLOC(w2, m, double);

  for (i = 0; i < m; i++)
    {
      for (j = 0; j < n; j++)
        {
          rr[j] += eq[i * n + j] * rhs[i];
          for (k = j; k < n; k++)
            {
              co[j * n + k] = co[k * n + j] += eq[i * n + j] * eq[i * n + k];
            }
        }
    }

  /**
   y = 1.00001 ;
   for (j=0; j<n ; j++) {  
   co[j*n+j] *= y ;
   }
   */

  if (verbose)
    {
      printf ("coeffs:\n");
      printmat (co, n, n);
      printf ("\n\n");
      printmat (rr, n, 1);

      for (i = 0; i < n; i++)
        {
          printf ("diag: %3d %9.3f\n", i, co[i * n + i]);
        }

      fflush (stdout);
    }

  ret = solvit (co, rr, n, ans);
  if (ret < 0)
    return -1000.0;
  for (i = 0; i < m; i++)
    {
      ww[i] = rhs[i] - vdot (ans, eq + i * n, n);
    }

  ynum = (double) m;
  vres = asum2 (ww, m) / ynum;
  vbase = asum2 (rhs, m) / ynum;

  /**
   printf("zzreg  %15.9f  %15.9f\n", log(vbase), log(vres)) ;
   printmat(rr, 1, n) ;
   printmat(co, n, n) ;
   printf("\n") ;
   */

  free (co);
  free (rr);
  free (ww);
  free (w2);
  return ynum * log (vbase / vres);
}

void
regressitall (char **vname, double *eq, double *rhs, int m, int n)
{
  double *ans;
  int i, j, k, wt;
  int npow;
  int **tab, *tweight, *cols;
  double *teq;
  double yscore;

  npow = (int) pow (2.0, (double) n);
  ZALLOC(tab, npow, int *);
  ZALLOC(tweight, npow, int);
  ZALLOC(cols, n, int);

  for (k = 0; k < npow; k++)
    {
      ZALLOC(tab[k], n, int);
    }
  for (k = 1; k < npow; k++)
    {
      add1 (tab[k], tab[k - 1], n);
      tweight[k] = intsum (tab[k], n);
    }
  ZALLOC(ans, n, double);
  ZALLOC(teq, m*n, double);

  for (wt = 1; wt <= n; ++wt)
    {
      printf ("weight: %d\n", wt);
      for (k = 0; k < npow; ++k)
        {
          if (tweight[k] != wt)
            continue;
          for (i = 0, j = 0; i < n; i++)
            {
              if (tab[k][i] == 0)
                continue;
              cols[j] = i;
              ++j;
            }
          squishx (teq, eq, m, n, cols, wt);
          yscore = regressit (ans, teq, rhs, m, wt);
          printf ("chisq: %9.3f\n", yscore);
          for (i = 0, j = 0; i < n; i++)
            {
              if (tab[k][i] == 0)
                continue;
              printf ("%15s %9.3f\n", vname[i], ans[j]);
              ++j;
            }
          printf ("\n");
        }
    }

  free (ans);
  free (teq);

  for (k = 0; k < npow; k++)
    {
      free (tab[k]);
    }
  free (tab);
  free (tweight);
  free (cols);
}

void
add1 (int *a, int *b, int n)
// b is 0, 1 vector as base 2 integer.  a = b + 1
{
  if (n == 0)
    return;
  copyiarr (b, a, n);
  a[n - 1] = b[n - 1] + 1;
  if (a[n - 1] == 2)
    {
      a[n - 1] = 0;
      add1 (a, b, n - 1);
    }
}

// now logistic regression stuff

double
logregressit (double *ans, double *eq, double **rhs, int neq, int nv)
// return log likelihood NOT chi-sq
{
  double *p, *z, *q;
  double *n0, *n1, *tans;
  double *grad, *hess, rr[2];
  double y0, y1, y, ylike, ybase, yold;
  int i, j;
  int iter, numiter = 10;
  int ret;

  ZALLOC(p, neq, double);
  ZALLOC(q, neq, double);
  ZALLOC(z, neq, double);
  ZALLOC(n0, neq, double);
  ZALLOC(n1, neq, double);
  ZALLOC(tans, neq, double);
  ZALLOC(grad, nv, double);
  ZALLOC(hess, nv*nv, double);

  for (i = 0; i < neq; i++)
    {
      y0 = n0[i] = rhs[i][0];
      y1 = n1[i] = rhs[i][1];
      y0 += 1.0;
      y1 += 1.0;
      y = y1 / (y0 + y1);
      y = MIN(y, 0.75);
      y = MAX(y, 0.25);
// may need changing for some problems
      p[i] = y;
    }
  y0 = asum (n0, neq);
  y1 = asum (n1, neq);
  y = y1 / (y0 + y1);
  for (i = 0; i < neq; i++)
    {
      p[i] = (p[i] + y) / 2.0;
      if (p[i] < 0.0)
        fatalx ("bugbug\n");
      if (p[i] > 1.0)
        fatalx ("bugbug\n");
    }

  if (verbose)
    {
      vzero (rr, 2);
      for (j = 0; j < neq; j++)
        {
          vvp (grad, grad, eq + j * nv, nv);
          vvp (rr, rr, rhs[j], 2);
          addouter (hess, eq + j * nv, nv);
        }
      y = 1.0 / (double) neq;
      vst (grad, grad, y, nv);
      vst (rr, rr, y, 2);
      vst (hess, hess, y, nv * nv);
      printf ("## averages\n");
      printmat (grad, 1, nv);
      printmat (rr, 1, 2);
      printmat (hess, nv, nv);
    }

  ptoz (p, z, neq);
  regressit (ans, eq, z, neq, nv);
  for (j = 0; j < neq; j++)
    {
      z[j] = vdot (eq + j * nv, ans, nv);
    }
  ybase = zlike (eq, n0, n1, ans, neq, nv);
  ztop (p, z, neq);

  calcgh (grad, hess, eq, z, n0, n1, neq, nv);
  y = .001;
  for (i = 0; i < nv; i++)
    {
      if (!verbose)
        break;
      copyarr (ans, tans, nv);
      tans[i] += y;
      ylike = zlike (eq, n0, n1, tans, neq, nv);
      printf ("zzgrad %3d %12.6f %12.6f\n", i, ylike - ybase, grad[i] * y);
    }

  for (iter = 1; iter <= numiter; ++iter)
    {
      calcgh (grad, hess, eq, z, n0, n1, neq, nv);
      ret = solvit (hess, grad, nv, tans);
      if (ret < 0)
        return -1000.0;
      if (verbose)
        {
          printf ("zzzz\n");
          printmat (ans, 1, nv);
          printmat (grad, 1, nv);
          printmat (hess, nv, nv);
          printmat (tans, 1, nv);
          printf ("\n\n");
        }
      vvp (ans, ans, tans, nv);
      for (j = 0; j < neq; j++)
        {
          z[j] = vdot (eq + j * nv, ans, nv);
        }
      ylike = zlike (eq, n0, n1, ans, neq, nv);
      /**
       if (verbose)  {
       printf("iter: %3d  llike: %15.9f incr: %15.9f\n", iter, ylike, ylike-ybase) ;
       printmat(ans, 1, nv) ;
       }
       */
      if ((iter > 1) && (ylike < (yold + .0001)))
        break;
      yold = ylike;
    }

  free (p);
  free (q);
  free (z);
  free (n0);
  free (n1);
  free (tans);
  free (grad);
  free (hess);

  return ylike;
}

void
ptoz (double *p, double *z, int n)
{
  double *w1, *w2;

  ZALLOC(w1, n, double);
  ZALLOC(w2, n, double);

  vst (w2, p, -1.0, n);
  vsp (w2, w2, 1.0, n);  // q 
  vvd (w1, p, w2, n);
  vlog (z, w1, n);
  free (w1);
  free (w2);
}
void
ztop (double *p, double *z, int n)
{
  double *ww, *w1;

  ZALLOC(ww, n, double);
  ZALLOC(w1, n, double);

  vexp (ww, z, n);
  vsp (w1, ww, 1.0, n); // 1 + e^z 
  vvd (p, ww, w1, n);  // p 

  free (ww);
  free (w1);
}

void
calcgh (double *grad, double *hess, double *eq, double *z, double *n0,
        double *n1, int neq, int nv)

{

  double *ww, *w1, *w2, *x0, *x1;
  int j;

  ZALLOC(ww, neq, double);
  ZALLOC(w1, neq, double);
  ZALLOC(w2, neq, double);
  ZALLOC(x0, neq, double);
  ZALLOC(x1, neq, double);

  vexp (ww, z, neq);
  vsp (w1, ww, 1.0, neq);
  vvt (w2, w1, w1, neq);
  vvt (x0, n0, ww, neq);
  vvm (x0, x0, n1, neq);
  vvd (x0, x0, w1, neq);

  vvp (x1, n0, n1, neq);
  vvt (x1, x1, ww, neq);
  vvd (x1, x1, w2, neq);

  vzero (grad, nv);
  vzero (hess, nv * nv);

  for (j = 0; j < neq; j++)
    {
      vst (ww, eq + j * nv, x0[j], nv);
      vvm (grad, grad, ww, nv);
      vst (ww, eq + j * nv, sqrt (x1[j]), nv);
      addouter (hess, ww, nv); // actually -hess 
    }
  free (ww);
  free (w1);
  free (w2);
  free (x0);
  free (x1);
}
double
zlike (double *eq, double *n0, double *n1, double *ans, int neq, int nv)
{
  double *z, *p, *q;
  double ylike, pprob, qprob, y0, y1, ybase;
  int j;

  ZALLOC(z, neq, double);
  ZALLOC(p, neq, double);
  ZALLOC(q, neq, double);

  y0 = asum (n0, neq);
  y1 = asum (n1, neq);

  y0 += 1.0e-10;
  y1 += 1.0e-10;

  pprob = y1 / (y0 + y1);
  qprob = y0 / (y0 + y1);

  ybase = y1 * log (pprob) + y0 * log (qprob);

  for (j = 0; j < neq; j++)
    {
      z[j] = vdot (eq + j * nv, ans, nv);
    }
  ztop (p, z, neq);
  vst (q, p, -1.0, neq);
  vsp (q, q, 1.0, neq);
  ylike = vldot (n1, p, neq) + vldot (n0, q, neq);
  ylike -= ybase;
  free (z);
  free (p);
  free (q);
  return ylike;
}
double
logrscore (double *eq, double **rhs, int neq, int nv)
// test significance of last regressor  
{
  double *teq, *ans;
  double y1, y2, ychi;
  int i;

  ZALLOC(teq, neq*nv, double);
  ZALLOC(ans, nv, double);

  squish (teq, eq, neq, nv, nv - 1);

  y1 = logregressit (ans, teq, rhs, neq, nv - 1);
  y2 = logregressit (ans, eq, rhs, neq, nv);

  ychi = 2.0 * (y2 - y1);

  free (teq);
  free (ans);

  return ychi;

}

void
squish (double *xmat, double *mat, int nrow, int oldc, int newc)
// in place legal !
{
  int i;
  double *ww;

  ZALLOC(ww, nrow*newc, double);

  for (i = 0; i < nrow; i++)
    {
      copyarr (mat + i * oldc, ww + i * newc, newc);
    }

  copyarr (ww, xmat, nrow * newc);
  free (ww);

}
void
squishx (double *xmat, double *mat, int nrow, int oldc, int *cols, int newc)
// copy cols of mat to xmat 
{
  int i, j, k;
  for (i = 0; i < nrow; i++)
    {
      for (j = 0; j < newc; ++j)
        {
          k = cols[j];
          xmat[i * newc + j] = mat[i * oldc + k];
        }
    }
}
void
calcres (double *res, double *ans, double *eq, double *rhs, int neq, int nv)
/**   
 calculate residual
 */
{

  int i;
  for (i = 0; i < neq; i++)
    {
      res[i] = rhs[i] - vdot (eq + i * nv, ans, nv);
    }
}

