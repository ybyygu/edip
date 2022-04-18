/* forces-edip.c
   -------------

   Version 1.0c

   Force and Energy Calculation with the
   Environment-Dependent Interatomic Potential

   written by Martin Z. Bazant,
   Department of Physics, Harvard University
   April - October 1997
   (based on forces.c, June 1994)

   Current address (2000):
   Professor Martin Z. Bazant
   Department of Mathematics 2-363B
   Massachusetts Institute of Technology
   Cambridge, MA 02139-4307

   E-mail:
   bazant@math.mit.edu


   COPYRIGHT NOTICE
   ----------------

   forces-edip, copyright 1997 by Martin Z. Bazant and Harvard University.
   Permission is granted to use forces-edip.c for academic use only, at
   no cost. Unauthorized sale or commerical use of this software
   is prohibited by United States copyright law. Any publication describing
   research involving this software should contain the following citations,
   at once and in this order, to give proper credit to the theoretical work and
   fitting that produced EDIP and this subroutine:

     1.  M. Z. Bazant and E. Kaxiras, Phys. Rev. Lett. 77, 4370 (1996).
     2.  M. Z. Bazant, E. Kaxiras, J. F. Justo, Phys. Rev. B 56, 8542 (1997).
     3.  J. F. Justo, M. Z. Bazant, E. Kaxiras, V. V. Bulatov, and S. Yip,
           Phys. Rev. B 58, 2539 (1998).

   This software has been extensively tested for molecular dynamics simulations
   on Sun, SGI and IBM architectures, but no guarantees are made.


   WEBSITE
   -------

   Updated versions of this software are available at the EDIP distribution
   site, http://pelion.eas.harvard.edu/software/EDIP/ Postscript files of
   related papers are available at the Kaxiras group web site in the Department
   of Physics at Harvard University, http://pelion.eas.harvard.edu, under
   'Empirical Methods'. A description of the algorithm used in this subroutine
   can be found in the Ph.D. Thesis of M. Z. Bazant (1997), chapter 6, on the
   web at

  HISTORY
  -------
  2022-04-18: Updated for easy to debug and port to Rust (ybyygu@gmail.com)
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* DEBUGGING FLAGS - systematically turn on various pieces of the potential. */
#define V2_on 1
#define V2Z_on 1
#define V3_on 1
#define V3g_on 1
#define V3h_on 1
#define V3Z_on 1
#define Zfast 1

/* Total number of particles */
static int N_own = 5;
static int MAX_NBRS = 30;
static int MAX_PART = 4097;
static int MAX_NBRS_1 = 4096;

typedef double SCALAR_t;
typedef struct {
  SCALAR_t x, y, z;
} VECTOR_t;
typedef VECTOR_t export_t; /* r */
typedef VECTOR_t import_t; /* f */

/* Verlet neighbor list with bonds double-counted */
static int neighbors[100];
static int p_nbrs[100];

/* EDIP Si PARAMETERS Justo et al., Phys. Rev. B 58, 2539 (1998).

     5.6714030     2.0002804     1.2085196     3.1213820     0.5774108
     1.4533108     1.1247945     3.1213820     2.5609104    78.7590539
     0.6966326   312.1341346     1.4074424     0.0070975     3.1083847

connection between these parameters and Justo et al., Phys. Rev. B 58, 2539
(1998):

A((B/r)**rh-palp*exp(-bet*Z*Z)) = A'((B'/r)**rh-exp(-bet*Z*Z))

so in the paper (')
A' = A*palp
B' = B * palp**(-1/rh)
eta = detla/Qo

*/
static double A, B, rh, sig, a;
static double lam, gam, b, c;
static double mu, Qo, bet, alp;
static double bg;   /* cutoff for g(r) */
static double palp; /* justo prefactor for bond order - delete later */
static double delta, eta, zet;

/* tau(Z) (Ismail & Kaxiras, 1993) */
static const double u1 = -0.165799;
static const double u2 = 32.557;
static const double u3 = 0.286198;
static const double u4 = 0.66;

/* measurements in force routine */
static double virial, v2sum, E_potential, coord_total, V2, V3;

/********************************************/

void init_EDIP() {
  A = 5.6714030;
  B = 2.0002804;
  rh = 1.2085196;
  a = 3.1213820;
  sig = 0.5774108;
  lam = 1.4533108;
  gam = 1.1247945;
  b = 3.1213820;
  c = 2.5609104;
  delta = 78.7590539;
  mu = 0.6966326;
  Qo = 312.1341346;
  palp = 1.4074424;
  bet = 0.0070975;
  alp = 3.1083847;
  printf("\n EDIP-Si Parameters: \n");
  printf("%lf %lf %lf %lf %lf \n", A, B, rh, a, sig);
  printf("%lf %lf %lf %lf %lf \n", lam, gam, b, c, delta);
  printf("%lf %lf %lf %lf %lf \n", mu, Qo, palp, bet, alp);
  bg = a;
  eta = delta / Qo;
}

/********************************************/

void compute_forces_EDIP(measure, version, pos, f) int measure;
int version;
export_t pos[MAX_PART];
import_t f[MAX_PART];
{
  /*------------------------- VARIABLE DECLARATIONS -------------------------*/

  int i, j, k, l, n;
  double dx, dy, dz, r, rsqr, asqr;
  double rinv, rmainv, xinv, xinv3, den, Z, fZ;
  double dV2j, dV2ijx, dV2ijy, dV2ijz, pZ, dp;
  double temp0, temp1, temp2;
  double Qort, muhalf, u5;
  double rmbinv, winv, dwinv, tau, dtau, lcos, x, H, dHdx, dhdl;
  double dV3rij, dV3rijx, dV3rijy, dV3rijz;
  double dV3rik, dV3rikx, dV3riky, dV3rikz;
  double dV3l, dV3ljx, dV3ljy, dV3ljz, dV3lkx, dV3lky, dV3lkz;
  double dV2dZ, dxdZ, dV3dZ;
  double dEdrl, dEdrlx, dEdrly, dEdrlz;
  double bmc, cmbinv;
  double fjx, fjy, fjz, fkx, fky, fkz;

  typedef struct {
    double t0, t1, t2, t3; /* various V2 functions and derivatives */
    double dx, dy, dz;     /* unit separation vector */
    double r;              /* bond length (only needed for virial) */
  } store2;
  store2 s2[MAX_NBRS_1]; /* two-body interaction quantities, r<a */
  int n2;                /* size of s2[] */
  int num2[MAX_NBRS_1];  /* atom ID numbers for s2[] */

  typedef struct {
    double g, dg;      /* 3-body radial function and its derivative */
    double rinv;       /* 1/r */
    double dx, dy, dz; /* unit separation vector */
    double r;          /* bond length (only needed for virial) */
  } store3;
  store3 s3[MAX_NBRS_1]; /* three-body interaction quantities, r<b */
  int n3;                /* size of s3[] */
  int num3[MAX_NBRS_1];  /* atom ID numbers for s3[] */

  typedef struct {
    double df;         /* derivative of neighbor function f'(r) */
    double sum;        /* array to accumulate coordination force prefactors */
    double dx, dy, dz; /* unit separation vector */
    double r;          /* bond length (only needed for virial) */
  } storez;
  storez sz[MAX_NBRS_1]; /* coordination number stuff, c<r<b */
  int nz;                /* size of sz[] */
  int numz[MAX_NBRS_1];  /* atom ID numbers for sz[] */

  int nj, nk, nl; /* indices for the store arrays */

  /*---------------------------------------------------------------------------*/

  /* INITIALIZE FORCES AND GLOBAL SUMS */

  for (i = 0; i < N_own; i++) {
    f[i].x = 0.0;
    f[i].y = 0.0;
    f[i].z = 0.0;
  }
  if (measure) {
    coord_total = 0.0;
    virial = 0.0;
    V2 = 0.0; /* E_potential = V2 + V3 */
    V3 = 0.0;
  }

  /* COMBINE COEFFICIENTS */

  asqr = a * a;
  Qort = sqrt(Qo);
  muhalf = mu * 0.5;
  u5 = u2 * u4;
  bmc = b - c;
  cmbinv = 1.0 / (c - b);

  /*--- LEVEL 1: OUTER LOOP OVER ATOMS ---*/

  for (i = 0; i < N_own; i++) {

    /* RESET COORDINATION AND NEIGHBOR NUMBERS */

    Z = 0.0;
    n2 = 0;
    n3 = 0;
    nz = 0;

    /*--- LEVEL 2: LOOP PREPASS OVER PAIRS ---*/

    for (n = p_nbrs[i]; n < p_nbrs[i + 1]; n++) {
      j = neighbors[n];

      /* TEST IF WITHIN OUTER CUTOFF */

      dx = pos[j].x - pos[i].x;
      if (fabs(dx) < a) {
        dy = pos[j].y - pos[i].y;
        if (fabs(dy) < a) {
          dz = pos[j].z - pos[i].z;
          if (fabs(dz) < a) {
            rsqr = dx * dx + dy * dy + dz * dz;
            if (rsqr < asqr) {
              r = sqrt(rsqr);

              /* PARTS OF TWO-BODY INTERACTION r<a */

              num2[n2] = j;
              rinv = 1.0 / r;
              dx *= rinv;
              dy *= rinv;
              dz *= rinv;
              rmainv = 1.0 / (r - a);
              s2[n2].t0 = A * exp(sig * rmainv);
              s2[n2].t1 = pow(B * rinv, rh);
              s2[n2].t2 = rh * rinv;
              s2[n2].t3 = sig * rmainv * rmainv;
              s2[n2].dx = dx;
              s2[n2].dy = dy;
              s2[n2].dz = dz;
              s2[n2].r = r;
              n2++;

              /* RADIAL PARTS OF THREE-BODY INTERACTION r<b */

              if (r < bg) {

                num3[n3] = j;
                rmbinv = 1.0 / (r - bg);
                temp1 = gam * rmbinv;
                temp0 = exp(temp1);
#if V3g_on
                s3[n3].g = temp0;
                s3[n3].dg = -rmbinv * temp1 * temp0;
#else
                s3[n3].g = 1;
                s3[n3].dg = 0;
#endif
                s3[n3].dx = dx;
                s3[n3].dy = dy;
                s3[n3].dz = dz;
                s3[n3].rinv = rinv;
                s3[n3].r = r;
                n3++;

                /* COORDINATION AND NEIGHBOR FUNCTION c<r<b */

                if (r < b) {
                  if (r < c)
                    Z += 1.0;
                  else {
                    xinv = bmc / (r - c);
                    xinv3 = xinv * xinv * xinv;
                    den = 1.0 / (1 - xinv3);
                    temp1 = alp * den;
                    fZ = exp(temp1);
                    Z += fZ;
                    numz[nz] = j;
                    sz[nz].df = fZ * temp1 * den * 3.0 * xinv3 * xinv *
                                cmbinv; /* df/dr */
                    sz[nz].dx = dx;
                    sz[nz].dy = dy;
                    sz[nz].dz = dz;
                    sz[nz].r = r;
                    nz++;
                  }
                }
              }
            }
          }
        }
      }
    }

    if (measure)
      coord_total += Z;

    /* ZERO ACCUMULATION ARRAY FOR ENVIRONMENT FORCES */

    for (nl = 0; nl < nz; nl++)
      sz[nl].sum = 0.0;

      /* ENVIRONMENT-DEPENDENCE OF PAIR INTERACTION */

#if V2_on

#if V2Z_on
    temp0 = bet * Z;
    pZ = palp * exp(-temp0 * Z); /* bond order */
    dp = -2.0 * temp0 * pZ;      /* derivative of bond order */
#else
    pZ = palp * exp(-bet * 16);
    dp = 0.0;
#endif

    /*--- LEVEL 2: LOOP FOR PAIR INTERACTIONS ---*/

    for (nj = 0; nj < n2; nj++) {

      temp0 = s2[nj].t1 - pZ;

      /* two-body energy V2(rij,Z) */

      if (measure)
        V2 += temp0 * s2[nj].t0;

      /* two-body forces */

      dV2j = -(s2[nj].t0) *
             ((s2[nj].t1) * (s2[nj].t2) + temp0 * (s2[nj].t3)); /* dV2/dr */
      dV2ijx = dV2j * s2[nj].dx;
      dV2ijy = dV2j * s2[nj].dy;
      dV2ijz = dV2j * s2[nj].dz;
      f[i].x += dV2ijx;
      f[i].y += dV2ijy;
      f[i].z += dV2ijz;
      j = num2[nj];
      f[j].x -= dV2ijx;
      f[j].y -= dV2ijy;
      f[j].z -= dV2ijz;

      /* dV2/dr contribution to virial */

      if (measure)
        virial -= s2[nj].r * (dV2ijx * s2[nj].dx + dV2ijy * s2[nj].dy +
                              dV2ijz * s2[nj].dz);

      /*--- LEVEL 3: LOOP FOR PAIR COORDINATION FORCES ---*/

      dV2dZ = -dp * s2[nj].t0;
      for (nl = 0; nl < nz; nl++)
        sz[nl].sum += dV2dZ;
    }
#endif

    /* COORDINATION-DEPENDENCE OF THREE-BODY INTERACTION */

#if V3_on

#if V3Z_on
    winv = Qort * exp(-muhalf * Z); /* inverse width of angular function */
    dwinv = -muhalf * winv;         /* its derivative */
    temp0 = exp(-u4 * Z);
    tau = u1 + u2 * temp0 * (u3 - temp0); /* -cosine of angular minimum */
    dtau = u5 * temp0 * (2 * temp0 - u3); /* its derivative */
#else
    winv = Qort * exp(-muhalf * 4);
    dwinv = 0.0;
    tau = 1.0 / 3.0;
    dtau = 0.0;
#endif

    /*--- LEVEL 2: FIRST LOOP FOR THREE-BODY INTERACTIONS ---*/

    for (nj = 0; nj < (n3 - 1); nj++) {

      j = num3[nj];

      /*--- LEVEL 3: SECOND LOOP FOR THREE-BODY INTERACTIONS ---*/

      for (nk = nj + 1; nk < n3; nk++) {

        k = num3[nk];

        /* angular function h(l,Z) */

        lcos = s3[nj].dx * s3[nk].dx + s3[nj].dy * s3[nk].dy +
               s3[nj].dz * s3[nk].dz;
        x = (lcos + tau) * winv;
        temp0 = exp(-x * x);
#if V3h_on
        H = lam * (1 - temp0 + eta * x * x);
        dHdx = 2 * lam * x * (temp0 + eta);
        dhdl = dHdx * winv;
#else
        H = 1.0;
        dhdl = 0.0;
#endif

        /* three-body energy */

        temp1 = s3[nj].g * s3[nk].g;
        if (measure) {
          V3 += temp1 * H;
        }

        /* (-) radial force on atom j */

        dV3rij = s3[nj].dg * s3[nk].g * H;
        dV3rijx = dV3rij * s3[nj].dx;
        dV3rijy = dV3rij * s3[nj].dy;
        dV3rijz = dV3rij * s3[nj].dz;
        fjx = dV3rijx;
        fjy = dV3rijy;
        fjz = dV3rijz;

        /* (-) radial force on atom k */

        dV3rik = s3[nj].g * s3[nk].dg * H;
        dV3rikx = dV3rik * s3[nk].dx;
        dV3riky = dV3rik * s3[nk].dy;
        dV3rikz = dV3rik * s3[nk].dz;
        fkx = dV3rikx;
        fky = dV3riky;
        fkz = dV3rikz;

        /* (-) angular force on j */

        dV3l = temp1 * dhdl;
        dV3ljx = dV3l * (s3[nk].dx - lcos * s3[nj].dx) * s3[nj].rinv;
        dV3ljy = dV3l * (s3[nk].dy - lcos * s3[nj].dy) * s3[nj].rinv;
        dV3ljz = dV3l * (s3[nk].dz - lcos * s3[nj].dz) * s3[nj].rinv;
        fjx += dV3ljx;
        fjy += dV3ljy;
        fjz += dV3ljz;

        /* (-) angular force on k */

        dV3lkx = dV3l * (s3[nj].dx - lcos * s3[nk].dx) * s3[nk].rinv;
        dV3lky = dV3l * (s3[nj].dy - lcos * s3[nk].dy) * s3[nk].rinv;
        dV3lkz = dV3l * (s3[nj].dz - lcos * s3[nk].dz) * s3[nk].rinv;
        fkx += dV3lkx;
        fky += dV3lky;
        fkz += dV3lkz;

        /* apply radial + angular forces to i, j, k */

        f[j].x -= fjx;
        f[j].y -= fjy;
        f[j].z -= fjz;
        f[k].x -= fkx;
        f[k].y -= fky;
        f[k].z -= fkz;
        f[i].x += fjx + fkx;
        f[i].y += fjy + fky;
        f[i].z += fjz + fkz;

        /* dV3/dR contributions to virial */

        if (measure) {
          virial -=
              s3[nj].r * (fjx * s3[nj].dx + fjy * s3[nj].dy + fjz * s3[nj].dz);
          virial -=
              s3[nk].r * (fkx * s3[nk].dx + fky * s3[nk].dy + fkz * s3[nk].dz);
        }

        /* prefactor for 4-body forces from coordination */
#if V3Z_on
        dxdZ = dwinv * (lcos + tau) + winv * dtau;
        dV3dZ = temp1 * dHdx * dxdZ;

        /*--- LEVEL 4: LOOP FOR THREE-BODY COORDINATION FORCES ---*/

        for (nl = 0; nl < nz; nl++)
          sz[nl].sum += dV3dZ;
#endif
      }
    }
#endif

    /*--- LEVEL 2: LOOP TO APPLY COORDINATION FORCES ---*/

    for (nl = 0; nl < nz; nl++) {

      dEdrl = sz[nl].sum * sz[nl].df;
      dEdrlx = dEdrl * sz[nl].dx;
      dEdrly = dEdrl * sz[nl].dy;
      dEdrlz = dEdrl * sz[nl].dz;
      f[i].x += dEdrlx;
      f[i].y += dEdrly;
      f[i].z += dEdrlz;
      l = numz[nl];
      f[l].x -= dEdrlx;
      f[l].y -= dEdrly;
      f[l].z -= dEdrlz;

      /* dE/dZ*dZ/dr contribution to virial */

      if (measure)
        virial -= sz[nl].r * (dEdrlx * sz[nl].dx + dEdrly * sz[nl].dy +
                              dEdrlz * sz[nl].dz);
    }
  }

  E_potential = V2 + V3;
  virial /= 3.0;
}

void main() {
  VECTOR_t pos[5];
  VECTOR_t f[5];
  pos[0].x = 0.715452;
  pos[0].y = 0.944157;
  pos[0].z = -0.270162;

  pos[1].x = 0.822738;
  pos[1].y = -2.210230;
  pos[1].z = -1.453488;

  pos[2].x = -1.171290;
  pos[2].y = -1.456535;
  pos[2].z = 1.156440;

  pos[3].x = 1.251249;
  pos[3].y = -1.227517;
  pos[3].z = 0.765830;

  pos[4].x = -1.006124;
  pos[4].y = -0.587667;
  pos[4].z = -1.143735;

  f[0].x = 0.0;
  f[0].y = 0.0;
  f[0].z = 0.0;

  f[1].x = 0.0;
  f[1].y = 0.0;
  f[1].z = 0.0;

  f[2].x = 0.0;
  f[2].y = 0.0;
  f[2].z = 0.0;

  f[3].x = 0.0;
  f[3].y = 0.0;
  f[3].z = 0.0;

  f[4].x = 0.0;
  f[4].y = 0.0;
  f[4].z = 0.0;

  p_nbrs[0] = 0;
  p_nbrs[1] = 4;
  p_nbrs[2] = 8;
  p_nbrs[3] = 12;
  p_nbrs[4] = 16;
  p_nbrs[5] = 20;

  neighbors[0] = 1;
  neighbors[1] = 2;
  neighbors[2] = 3;
  neighbors[3] = 4;
  /* 0 */
  neighbors[4] = 0;
  neighbors[5] = 2;
  neighbors[6] = 3;
  neighbors[7] = 4;
  /* 1 */
  neighbors[8] = 0;
  neighbors[9] = 1;
  neighbors[10] = 3;
  neighbors[11] = 4;
  /* 2 */
  neighbors[12] = 0;
  neighbors[13] = 1;
  neighbors[14] = 2;
  neighbors[15] = 4;
  /* 3 */
  neighbors[16] = 0;
  neighbors[17] = 1;
  neighbors[18] = 2;
  neighbors[19] = 3;
  /* 4 */

  init_EDIP();
  compute_forces_EDIP(1, 0, pos, f);

  printf("energy: %lf\n", E_potential);
  printf("virial: %lf\n", virial);
  printf("forces\n");
  for (int i = 0; i < 5; i++) {
    printf("%i.x: %lf\n", i, f[i].x);
    printf("%i.y: %lf\n", i, f[i].y);
    printf("%i.z: %lf\n", i, f[i].z);
  }
}
