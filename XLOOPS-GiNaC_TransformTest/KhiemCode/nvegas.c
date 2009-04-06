/*
 * NAME
 *   nvegas.c
 *   Implementation of G.P.Lepage's VEGAS-algorithm.
 *
 * SYNOPSIS
 *   void vegas(double regn[], int ndim, void (*fxn)(double x[], double f[]),
 *              int init, unsigned long ncall, int itmx, int nprn,
 *              int fcns, int pdim, int wrks,
 *              double tgral[], double sd[], double chi2a[]);
 *
 *     regn[]: array specifying the region to be integrated, 2*ndim entries
 *     ndim: dimensionality of space
 *     (*fxn)(x[],f[]): pointer to function to be evaluated (must be MT-safe!)
 *     init: initialization level (start with 0, then 1, later 2)
 *     ncall: number of samples points per iteration
 *     itmx: number of iterations in one call
 *     nprn: bit field, see constants NPRN_* below
 *     fcns: actual number of integrands (<=FNMX),
 *           if > 1 additional function accumulators are used
 *     pdim: dimension of parallel space, 0==autosense, higher=="manual tuning"
 *     wrks: number of parallel working units
 *     tgral[]: pointer to estimate of result (maybe array)
 *     sd[]: pointer to estimate of standard deviation (maybe array)
 *     chi2a[]: pointer to chi-squared over ndim (maybe array)
 *     NOTE: pdim and wrks are unused here, they are provided only for
 *           object-compatibility with the vegas function from pvegas.c.
 *
 * DESCRIPTION
 *   pvegas is a parallel farmer-worker implementation of the popular
 *   VEGAS-algorithm. It splits up some dimensions (the first
 *   ndim_par ones) into separate chunks and then lets each
 *   worker-thread evaluate one chunk of these dimensions (and all the
 *   remaining dimensions). The random-numbers from gfsr are
 *   parallelized as well.
 *
 *   This is the reference-version without any parallelism. It is
 *   provided in order to make runs on one processor and subsequently
 *   compare the numerical results with the parallel versions pvegas.c
 *   and pvegas_mpi.c.  Please feel free to contact
 *   <Richard.Kreckel@Uni-Mainz.DE> if you encounter any
 *   implementation-specific problems.
 *
 *   No external random number generator (RNG) needs to be supplied. A
 *   shift register-generator (SR) is implemented which is initialized
 *   with the system-time. You can force reproducible results by
 *   defining REPRO to be some positive integer different from zero
 *   and sticking to that value. All versions of vegas provided by the
 *   author should return the same numerical values if the same values
 *   for REPRO are used.
 *
 *   Note that the RNG is guaranteed to work properly if your
 *   architecture adheres to the LP64-model, i.e. ints are 32 bit
 *   long. Another, more crucial, assumption is that chars are 8 bit
 *   long. If you are on some strange hardware you are well advised to
 *   check the random numbers manually by consulting the supplied
 *   sample program vegastest.c!
 *
 *   This version may differ considerably from other implementations
 *   found on the net with regards to the arguments passed to it. The
 *   design goal was to have a uniform interface across all versions
 *   of vegas supplied by the author (nvegas.c, pvegas.c,
 *   pvegas_mpi.c) and to make optimal use of parallel facilities.
 *   Consult vegas.h and the samples to make proper use of the
 *   interface.
 *
 * AUTHOR
 * Richard Kreckel, ThEP, Univ. Mainz, October 1996 - October 2000 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "vegas.h"

#define TINY 1.0e-68         /* small, since we are in double-precision      */
#define REPRO 1              /* 0 = default, others used for comparison      */

unsigned int gfsr_m[SR_P];   /* status-field for the GFSR-generator          */
int gfsr_k;                  /* pointer into status-field                    */
static int gfsr_not_initialized = 1;  /* flag: need to initialize GFSR-field */
unsigned int rdum;           /* linear congruential counter in kw_rand()     */
double gfsr_norm;            /* will be set such that gfsr is normalized     */
int functions;               /* copy of (*ctl).fcns                          */

/*
 * This routine is used by vegas. It rebins a vector of densities xi into
 * new bins defined by a vector r.
 */
void rebin(double rc, int nd, double r[], double xin[], double xi[])
{
  int i;
  int k = 0;
  double dr = 0.0;
  double xn = 0.0;
  double xo = 0.0;
  
  for (i=0; i<nd-1; i++) {
    while (rc > dr) {
      dr += r[k++];
    }
    if (k > 1) xo = xi[k-2];
    xn = xi[k-1];
    dr -= rc;
    xin[i] = xn-(xn-xo)*dr/r[k-1];
  }
  for (i=0; i<nd-1; i++) xi[i] = xin[i];
  xi[nd-1] = 1.0;
}

/*
 * gfsr produces the random numbers once the starting values have been
 * initialized using gfsr_init.
 */
double gfsr_rand(unsigned int w[], int *k)
{
  int j;
  
  (*k)++;
  if (*k >= SR_P) *k = 0;
  j = *k + SR_Q;
  if (j >= SR_P) j -= SR_P;
  w[*k] = w[*k] ^ w[j];
  return((double)w[*k] * gfsr_norm);
}

/*
 * Simple linear congruential generator used to initialize SR.
 * The multiplier and increment were provided by KOMA, Univ.Mainz.
 * (kw stands for Kalos & Whitlock)
 */
unsigned int kw_rand(void)
{
  rdum = (1812433253*rdum + 314159265);
  return (rdum);
}

/*
 * gfsr_init initializes the sequences using values from kw_rand.
 */
void gfsr_init(long seed)
{
  int i, j;
  
  printf("Initializing SR-sequences with seed %ld\n",seed);
  gfsr_norm = (double)(1/(pow(2.0,(long double)(8*sizeof(int))) - 1));
  rdum = (unsigned int)seed;
  for (j=0; j<SR_P; j++)
    gfsr_m[j] = kw_rand();      /* initialize starting values */
  gfsr_k = -1;                  /*         initialize pointer */
  gfsr_not_initialized = 0;
}

void vegas(double regn[], int ndim, void (*fxn)(double x[], double f[]),
           int init, unsigned long ncall, int itmx, int nprn,
           int fcns, int pdim, int wrks,
           double tgral[], double sd[], double chi2a[])
{
  static int ndo;            /*                                    (ndo)     */
  int it;                    /* iteration counter                  (it)      */
  static int ittot;          /* iteration counter across init>1              */
  int i, j, k;               /* counters                           (i, j, k) */
  int nd;                    /* slices in grid (c.f. NDMX)         (nd)      */
  int ng;                    /*                                    (ng)      */
  unsigned int npg;          /* number of calls within bin         (npg)     */
  static int mds;            /* ==1: statified smpl.               (mds)     */
  int ia[MXDIM];             /*                                    (ia[])    */
  int kg[MXDIM];             /*                                    (kg[])    */
  double calls;              /* real total number of calls to fxn  (calls)   */
  double dv2g;               /*                                    (dv2g)    */
  double dxg;                /*                                    (dxg)     */
  double rc;                 /*                                    (rc)      */
  double wgt;                /* weight                             (wgt)     */
  double xn;                 /*                                    (xn)      */
  double xnd;                /*                                    (xnd)     */
  double xo;                 /*                                    (xo)      */
  double xJac;               /* Jacobian of integration            (xjac)    */
  typedef struct {
    double Wgt;              /* weight                             (wgt)     */
    double sWgt;             /* cumulative sum for weights         (swgt)    */
    double sChi;             /* cumulative sum for chi^2           (schi)    */
    double sInt;             /* cumulative sum for integral        (si)      */
  } iterAccu;                /* accumulator for stuff at end of iteration... */
  static iterAccu Ai[FNMX];  /* ...one for each integrand                    */
  typedef struct {
    double ti;               /* sum for f over bins                (ti)      */
    double tsi;              /* sum for variances over bins        (tsi)     */
  } binAccu;                 /* accumulator over bins / hypercubes...        */
  binAccu Ab[FNMX];          /* ...one for each integrand                    */
  typedef struct {
    double f2;               /* f squared                          (f2)      */
    double fb;               /* sum for f within bin               (fb)      */
    double f2b;              /* sum for f2 within bin              (f2b)     */
    unsigned long npg;       /* number of calls within bin f != 0            */
  } pointAccu;               /* accumulator over points x within bins...     */
  pointAccu Ax[FNMX];        /* ...one for each integrand                    */
  double f[FNMX];            /* array passed into fxn for evaluation at x    */
  double x[MXDIM];           /* evaluation point                   (x[])     */
  double d[NDMX][MXDIM];     /*                                    (d[][])   */
  double di[NDMX][MXDIM];    /* delta i                            (di[][])  */
  double dt[MXDIM];          /*                                    (dt[])    */
  double r[NDMX];            /*                                    (r[])     */
  static double xi[MXDIM][NDMX];  /*                               (xi[][])  */
  double xin[NDMX];          /* aux. variable for rebinning        (xin[])   */
  double dx[MXDIM];          /* width of integration region        (dx[])    */
  double xrand;              /* random sample 0.0 <= xrand < 1.0   (xrand)   */

  functions = (fcns<FNMX) ? (fcns) : (FNMX);
#if (REPRO != 0)
  if (gfsr_not_initialized) gfsr_init(REPRO);
#else
  if (gfsr_not_initialized) gfsr_init((long)time(NULL));
#endif
  if (init <= 0) {    /* entry for cold start        */
    mds = ndo = 1;
    for (j=0; j<ndim; j++) xi[j][0] = 1.0;
  }
  if (init <= 1) {    /* inherit the previous grid   */
    for (j=0; j<functions; j++) {
      Ai[j].sInt = 0.0;
      Ai[j].sWgt = 0.0;
      Ai[j].sChi = 0.0;
    }
    ittot = 1;
  }
  if (init <= 2) {    /* inherit grid and results    */
    nd = NDMX;
    ng = 1;
    if (mds) {
      ng = (int)pow(ncall/2.0+0.25,1.0/ndim);
      mds = 1;
      if ((2*ng-NDMX) >= 0) {
        mds = -1;
        npg = ng/NDMX+1;
        nd = ng/npg;
        ng = npg*nd;
      }
    }
    for (k=1,i=0; i<ndim; i++) k *= ng;
    npg = (ncall/k>2) ? (ncall/k) : (2);
    calls = (double)npg * (double)k;
    dxg = 1.0/ng;
    for (dv2g=1,i=0; i<ndim; i++) dv2g *= dxg;
    dv2g = calls*calls*dv2g*dv2g/npg/npg/(npg-1.0);
    xnd = nd;
    dxg *= xnd;
    xJac = 1.0/calls;
    for (j=0; j<ndim; j++) {
      dx[j] = regn[j+ndim]-regn[j];
      xJac *= dx[j];
    }
    if (nd != ndo) {
      for (i=0; i<(nd>ndo?nd:ndo); i++) r[i] = 1.0;
      for (j=0; j<ndim; j++) rebin(ndo/xnd,nd,r,xin,xi[j]);
      ndo = nd;
    }
    if (nprn & NPRN_INPUT) {
      printf("%s:  ndim= %3d  ncall= %8.0f\n",
             " Input parameters for vegas",ndim,calls);
      printf("%28s  ittot=%5d  itmx=%5d\n"," ",ittot,itmx);
      printf("%28s  nprn=0x%04x  ALPH=%5.2f\n"," ",nprn,ALPH);
      printf("%28s  mds=%3d  nd=%4d%15s npg=%d\n"," ",mds,nd," ",npg);
      for (j=0; j<ndim; j++) {
        printf("%30s xl[%2d]= %11.4g xu[%2d]= %11.4g\n",
               " ",j,regn[j],j,regn[j+ndim]);
      }
    }
  }
  for (it=ittot; it<=itmx+ittot-1; it++) {
    for (j=0; j<functions; j++) Ab[j].ti = Ab[j].tsi = 0.0;
    for (j=0; j<ndim; j++) {
      kg[j] = 1;
      for (i=0; i<nd; i++) d[i][j] = di[i][j] = 0.0;
    }
    for (;;) {
      for (j=0; j<functions; j++) {
        Ax[j].fb = 0.0;
        Ax[j].f2b = 0.0;
        Ax[j].npg = 0;
      }
      for (k=0; k<npg; k++) {
        wgt = xJac;
        for (j=0; j<ndim; j++) {
          xrand = gfsr_rand(gfsr_m,&gfsr_k);
          xn = (kg[j]-xrand)*dxg+1.0;
          ia[j] = ((int)xn<NDMX) ? ((int)xn) : (NDMX);
          ia[j] = (ia[j]>1) ? (ia[j]) : (1);
          if (ia[j] > 1) {
            xo = xi[j][ia[j]-1]-xi[j][ia[j]-2];
            rc = xi[j][ia[j]-2]+(xn-ia[j])*xo;
          } else {
            xo = xi[j][ia[j]-1];
            rc = (xn-ia[j])*xo;
          }
          x[j] = regn[j]+rc*dx[j];
          wgt *= xo*xnd;
        }
        fxn(x,f);   /* call integrand at point x */
        for (j=0; j<functions; j++) {
          if (f[j] != 0.0) ++Ax[j].npg;
          f[j] *= wgt;
          Ax[j].f2 = f[j]*f[j];
          Ax[j].fb += f[j];
          Ax[j].f2b += Ax[j].f2;
        }
        for (j=0; j<ndim; j++) {
          di[ia[j]-1][j] += f[0];
          if (mds >= 0) d[ia[j]-1][j] += Ax[0].f2;
        }
      }  /* end of loop within hypercube */
      for (j=0; j<functions; j++) {
        Ax[j].f2b = sqrt(Ax[j].f2b*Ax[j].npg);
        Ax[j].f2b = (Ax[j].f2b-Ax[j].fb)*(Ax[j].f2b+Ax[j].fb);
        if (Ax[j].f2b <= 0.0) Ax[j].f2b = TINY;
        Ab[j].ti += Ax[j].fb;
        Ab[j].tsi += Ax[j].f2b;
      }
      if (mds < 0) {
        for (j=0; j<ndim; j++) d[ia[j]-1][j] += Ax[0].f2b;
      }
      for (k=ndim-1; k>=0; k--) {
        kg[k] %= ng;
        if (++kg[k] != 1) break;
      }
      if (k < 0) break;
    }  /* end of loop over hypercubes */
    for (j=0; j<functions; j++) {
      Ab[j].tsi *= dv2g;
      Ai[j].Wgt = 1.0/Ab[j].tsi;
      Ai[j].sInt += Ai[j].Wgt*Ab[j].ti;
      Ai[j].sChi += Ai[j].Wgt*Ab[j].ti*Ab[j].ti;
      Ai[j].sWgt += Ai[j].Wgt;
      tgral[j] = Ai[j].sInt/Ai[j].sWgt;
      chi2a[j] = (Ai[j].sChi-Ai[j].sInt*tgral[j])/(it-0.9999);
      if (chi2a[j] < 0.0) chi2a[j] = 0.0;
      sd[j] = sqrt(1.0/Ai[j].sWgt); 
      Ab[j].tsi = sqrt(Ab[j].tsi);
    }
    if (nprn & NPRN_RESULT) {
      printf("%s %3d : integral = %14.7g +/-  %9.2g\n",
             " iteration no.",it,Ab[0].ti,Ab[0].tsi);
      printf("%s integral =%14.7g+/-%9.2g  chi^2/IT n = %9.2g\n",
             " all iterations:  ",tgral[0],sd[0],chi2a[0]);
    }
    if (nprn & NPRN_SECRES) {
      for (i=1; i<functions; i++) {
        printf("   %4d%s%14.7g+/-%9.2g  chi^2/IT n = %9.2g\n",
               i,".additional integral= ",tgral[i],sd[i],chi2a[i]);
      }
    }
    if (nprn & (NPRN_GRID | NPRN_GRID_2 | NPRN_GRID_4 | NPRN_GRID_8)) {
      for (j=0; j<ndim; j++) {
        printf(" data for axis  %2d\n",j);
        printf("%6s%13s%11s%13s%11s%13s\n", 
               "X","delta i","X","delta i","X","delta i");
        for (i=0; i<nd; i += 3) {
          for (k=0; k<3 && i+k<nd; k++) {
            printf("%8.5f%12.4g    ",xi[j][i+k],di[i+k][j]);
          }
          printf("\n");
          if (nprn & NPRN_GRID_8) k = 3*(8-1);
          if (nprn & NPRN_GRID_4) k = 3*(4-1);
          if (nprn & NPRN_GRID_2) k = 3*(2-1);
          if (nprn & NPRN_GRID) k = 3*(1-1);
          i += k;
        }
      }
    }
    if (nprn) fflush(NULL);
    for (j=0; j<ndim; j++) {
      xo = d[0][j];
      xn = d[1][j];
      d[0][j] = (xo+xn)/2.0;
      dt[j] = d[0][j];
      for (i=1; i<nd-1; i++) {
        rc = xo+xn;
        xo = xn;
        xn = d[i+1][j];
        d[i][j] = (rc+xn)/3.0;
        dt[j] += d[i][j];
      }
      d[nd-1][j] = (xo+xn)/2.0;
      dt[j] += d[nd-1][j];
    }
    for (j=0; j<ndim; j++) {
      rc = 0.0;
      for (i=0; i<nd; i++) {
        if (d[i][j] < TINY) d[i][j] = TINY;
        r[i] = pow((1.0-d[i][j]/dt[j])/
                   (log(dt[j])-log(d[i][j])),ALPH);
        rc += r[i];
      }
      rebin(rc/xnd,nd,r,xin,xi[j]);
    }
  }
  ittot += itmx;
}

#undef SR_P
#undef SR_Q
#undef REPRO
#undef TINY
