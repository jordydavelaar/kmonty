

/*
   HARM model specWification routines
 */
#include "decs.h"
#define global
#include "BHAC3D_model.h"
#undef global

struct of_spectrum ***spect;
extern struct of_spectrum ***shared_spect;
double ***Xi_spec;
extern double ***shared_Xi_spec;
double ***ispec;
extern double ***shared_ispec;

#pragma omp threadprivate(spect, Xi_spec, ispec)

void init_model(char *args[]) {
  set_units(args[3], args[4]);

  init_bhac3d_data(args[2]);

  {
    int i, j, k;
    spect = (struct of_spectrum ***)malloc(N_THBINS *
                                           sizeof(struct of_spectrum **));
    for (i = 0; i < N_THBINS; i++) {
      spect[i] = (struct of_spectrum **)malloc(N_PHIBINS *
                                               sizeof(struct of_spectrum *));
      for (j = 0; j < N_PHIBINS; j++) {
        spect[i][j] =
            (struct of_spectrum *)malloc(N_EBINS * sizeof(struct of_spectrum));
      }
    }
    Xi_spec = (double ***)malloc(N1 * sizeof(double **));
    ispec = (double ***)malloc(N1 * sizeof(double **));
    for (i = 0; i < N1; i++) {
      Xi_spec[i] = (double **)malloc(N2 * sizeof(double *));
      ispec[i] = (double **)malloc(N2 * sizeof(double *));
      for (j = 0; j < N2; j++) {
        Xi_spec[i][j] = (double *)malloc(N3 * sizeof(double));
        ispec[i][j] = (double *)malloc(N3 * sizeof(double));
        for (k = 0; k < N3; k++) {
          Xi_spec[i][j][k] = 0;
          ispec[i][j][k] = 0;
        }
      }
    }
  }
  //	get_profiles();
  Rh = 1 + sqrt(1. - a * a);

  /* make look-up table for hot cross sections */
  init_hotcross();

  /* make table for solid angle integrated emissivity and K2 */
  init_emiss_tables();

  /* make table for superphoton weights */
  init_weight_table();

  /* make table for quick evaluation of ns_zone */
  init_nint_table();
}

/*
   make super photon

 */

int n2gen = -1;
double dnmax;
int zone_i, zone_j, zone_k;

void make_super_photon(struct of_photon *ph, int *quit_flag) {

  while (n2gen <= 0) {
    /*dnmax=number of superphotons needed to be crated in a given zone*/
    n2gen = get_zone(&zone_i, &zone_j, &zone_k, &dnmax);
    // fprintf(stdout,"i=%d j=%d k=%d n2gen=%d
    // dnmax=%g\n",zone_i,zone_j,zone_k,n2gen,dnmax);
  }

  n2gen--;

  if (zone_i == N1)
    *quit_flag = 1;
  else
    *quit_flag = 0;

  if (*quit_flag != 1) {
    /* Initialize the superphoton energy, direction, weight, etc. */
    sample_zone_photon(zone_i, zone_j, zone_k, dnmax, ph);
    //	  fprintf(stdout,"nu=%g [Hz]\n",ph->E*ME*CL*CL/HPL);
  }
  return;
}

double bias_func(double Te, double w) {

  double bias;

  bias = 100. * Te * Te / bias_norm;
#pragma omp atomic
  bias /= max_tau_scatt;
  if (bias < 1)
    bias = 1.;
  if (bias > 3000.0)
    bias = 3000.0;
  return 0; // bias;
  /* double bias,max,c,avg_num_scatt;
     avg_num_scatt = 2.+N_scatt / (1. * N_superph_recorded + 1.);

     avg_num_scatt = N_scatt / (1. * N_superph_recorded + 1.);
     bias =  Te * Te/(bias_norm  *(avg_num_scatt + 2));
     #pragma omp atomic
     bias/= max_tau_scatt;
     if (bias < 1)
     bias = 1.;
     //c=0.1*bias_norm;
     //max = 0.5 * w / WEIGHT_MIN;
     // bias =  Te * Te/(bias_norm*max_tau_scatt);
     //bias=0;
     // if(bias>5.)bias=5.;
     // if(bias<0.1){ bias=0.1;}
     // bias /= (max_tau_scatt*avg_num_scatt);
     // if(bias<1.0) bias=1./(max_tau_scatt*avg_num_scatt);
     // fprintf(stderr,"bias %g %g\n", bias, max);
     return bias;*/
}

/*

   these supply basic model data to grmonty

 */

// get zonal variables without any interpolations, this should not be changed
void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B,
                    double Ucon[NDIM], double Bcon[NDIM], int *ACCZONE) {
  int l, m;
  double Ucov[NDIM], Bcov[NDIM];
  double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
  double bsq, beta, beta_trans, b2, trat, two_temp_gam, Th_unit, Be;
  double X[NDIM];

  coord(i, j, k, X);
  // no interpolations and such , this is needed to genereate distribution
  // function
  *Ne = p[KRHO][i][j][k] * Ne_unit;

  Bp[1] = p[B1][i][j][k];
  Bp[2] = p[B2][i][j][k];
  Bp[3] = p[B3][i][j][k];

  Vcon[1] = p[U1][i][j][k];
  Vcon[2] = p[U2][i][j][k];
  Vcon[3] = p[U3][i][j][k];

  /* Get Ucov */
  VdotV = 0.;
  for (l = 1; l < NDIM; l++)
    for (m = 1; m < NDIM; m++)
      VdotV += geom[i][j].gcov[l][m] * Vcon[l] * Vcon[m];
  Vfac = sqrt(-1. / geom[i][j].gcon[0][0] * (1. + fabs(VdotV)));
  Ucon[0] = -Vfac * geom[i][j].gcon[0][0];
  for (l = 1; l < NDIM; l++)
    Ucon[l] = Vcon[l] - Vfac * geom[i][j].gcon[0][l];
  lower(Ucon, geom[i][j].gcov, Ucov);

  /* Get B and Bcov */
  UdotBp = 0.;
  for (l = 1; l < NDIM; l++)
    UdotBp += Ucov[l] * Bp[l];
  Bcon[0] = UdotBp;
  for (l = 1; l < NDIM; l++)
    Bcon[l] = (Bp[l] + Ucon[l] * UdotBp) / Ucon[0];
  lower(Bcon, geom[i][j].gcov, Bcov);

  bsq = Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] + Bcon[2] * Bcov[2] +
        Bcon[3] * Bcov[3] + 1e-40;
  *B = sqrt(bsq) * B_unit;

  gam = 4. / 3.;
  beta = p[UU][i][j][k] * (gam - 1.) / 0.5 / bsq; // beta plasma
  // beta=p[KRHO][i][j][k]/bsq/0.5; // beta plasma
  /*electron temperature depending on the plasma magnetization*/

  Be = (-(1. + gam * p[UU][i][j][k] / p[KRHO][i][j][k]) * Ucov[0]);
  beta_trans = 1.;
  b2 = pow(beta / beta_trans, 2);
  trat = trat_d * b2 / (1. + b2) + trat_j / (1. + b2);
  // two_temp_gam = 0.5 * ((1. + 2. / 3. * (trat + 1.) / (trat + 2.)) + gam);
  two_temp_gam = 4. / 3.;
  Th_unit = (two_temp_gam - 1.) * (MP / ME) / (1. + trat);
  *Thetae = (p[UU][i][j][k] / p[KRHO][i][j][k]) * Th_unit;
  Be = (-(1. + two_temp_gam * p[UU][i][j][k] / p[KRHO][i][j][k] +
          bsq / 2. / p[KRHO][i][j][k]) *
        Ucov[0]);

  if (bsq / p[KRHO][i][j][k] > 1. || exp(X[1]) > 200.) {
    *Ne = 0;
  }

  if (Be > 1.02 && bsq / p[KRHO][i][j][k] < 1.0) {
    *ACCZONE = 1;
  }
}

void get_fluid_params(double X[NDIM], double gcov[NDIM][NDIM], double *Ne,
                      double *Thetae, double *B, double Ucon[NDIM],
                      double Ucov[NDIM], double Bcon[NDIM], double Bcov[NDIM],
                      int *ACCZONE) {
  int i, j, k;
  double del[NDIM];
  double rho, uu;
  double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
  double gcon[NDIM][NDIM], coeff[4];
  double interp_scalar(double ***var, int i, int j, int k, double del[4]);
  double bsq, beta, beta_trans, b2, trat, two_temp_gam, Th_unit, Be;

  if (X[1] < startx[1] || X[1] > stopx[1] || X[2] < startx[2] ||
      X[2] > stopx[2]) {
    *Ne = 0.;
    return;
  }

  Xtoijk(X, &i, &j, &k, del);

  coeff[1] = del[1];
  coeff[2] = del[2];
  coeff[3] = del[3];

  // this is for transfer so interpolations needed
  // interpolate (cubiclinear interp.) like in mibothros
  rho = interp_scalar(p[KRHO], i, j, k, coeff);
  uu = interp_scalar(p[UU], i, j, k, coeff);
  *Ne = rho * Ne_unit;

  // here unlike in mibothros it was interpolating scalars and
  // reconstructing velocity and magnetic field based on interpolated
  // coefficients
  Bp[1] = interp_scalar(p[B1], i, j, k, coeff);
  Bp[2] = interp_scalar(p[B2], i, j, k, coeff);
  Bp[3] = interp_scalar(p[B3], i, j, k, coeff);

  Vcon[1] = interp_scalar(p[U1], i, j, k, coeff);
  Vcon[2] = interp_scalar(p[U2], i, j, k, coeff);
  Vcon[3] = interp_scalar(p[U3], i, j, k, coeff);

  gcon_func(X, gcon);

  /* Get Ucov based on zonal or reconstructed values*/
  /*        VdotV = 0.;
          for (i = 1; i < NDIM; i++)
                  for (j = 1; j < NDIM; j++)
                          VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
          Vfac = sqrt(-1. / gcon[0][0] * (1. + fabs(VdotV)));
          Ucon[0] = -Vfac * gcon[0][0];
          for (i = 1; i < NDIM; i++)
                  Ucon[i] = Vcon[i] - Vfac * gcon[0][i];
          lower(Ucon, gcov, Ucov);

          UdotBp = 0.;
          for (i = 1; i < NDIM; i++)
                  UdotBp += Ucov[i] * Bp[i];
          Bcon[0] = UdotBp;
          for (i = 1; i < NDIM; i++)
                  Bcon[i] = (Bp[i] + Ucon[i] * UdotBp) / Ucon[0];
   */
  // From gammav^i to v^i
  VdotV = 0.;
  for (int i = 1; i < NDIM; i++)
    for (int j = 1; j < NDIM; j++)
      VdotV += gcov[i][j] * Vcon[i] * Vcon[j]; // ww
  for (int i = 1; i < NDIM; i++)
    Vcon[i] = Vcon[i] / sqrt(1.0 + VdotV);

  VdotV = 0;
  for (int i = 1; i < NDIM; i++)
    for (int j = 1; j < NDIM; j++)
      VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
  Ucon[0] = 1. / (sqrt(-1 / gcon[0][0] * (1 - VdotV)));
  for (int i = 1; i < NDIM; i++)
    Ucon[i] = Ucon[0] * (Vcon[i] / sqrt(-gcon[0][0]) + gcon[0][i] / gcon[0][0]);

  // lower(Ucon, gcov, Ucov); // Gammie's lowering function
  lower(Ucon, gcov, Ucov);

  Bcon[0] = 0;
  for (i = 1; i < NDIM; i++) {
    for (int l = 0; l < NDIM; l++) {
      Bcon[0] += Bp[i] * (gcov[i][l] * Ucon[l]) * sqrt(-gcon[0][0]);
    }
  }

  for (int i = 1; i < NDIM; i++)
    Bcon[i] =
        (Bp[i] / Ucon[0]) * sqrt(-gcon[0][0]) + (Bcon[0] * Ucon[i] / Ucon[0]);

  lower(Bcon, gcov, Bcov);

  bsq = Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] + Bcon[2] * Bcov[2] +
        Bcon[3] * Bcov[3] + 1e-40;

  *B = sqrt(bsq) * B_unit;
  gam = 4. / 3.;

  /*electron temperature depending on the plasma magnetization*/
  beta = uu * (gam - 1.) / 0.5 / bsq;
  Be = (-(1. + gam * uu / (rho)) * Ucov[0]);
  beta_trans = 1.;
  b2 = pow(beta / beta_trans, 2);
  trat = trat_d * b2 / (1. + b2) + trat_j / (1. + b2);

  two_temp_gam =
      4. / 3.; // 0.5 * ((1. + 2. / 3. * (trat + 1.) / (trat + 2.)) + gam);
  Th_unit = (two_temp_gam - 1.) * (MP / ME) / (1. + trat);

  *Thetae = (uu / rho) * Th_unit;

  Be = (-(1. + two_temp_gam * uu / rho + bsq / 2. / rho) * Ucov[0]);

  if (bsq / rho > 1. || exp(X[1]) > 200.) {
    *Ne = 0;
  }
}

/*
   Current metric: modified Kerr-Schild, squashed in theta
   to give higher resolution at the equator
 */

/* mnemonics for dimensional indices */
#define TT 0
#define RR 1
#define TH 2
#define PH 3

void gcon_func(double *X, double gcon[][NDIM]) {

  int k, l;
  double sth, cth, irho2;
  double r, th;
  double hfac;
  /* required by broken math.h */
  void sincos(double in, double *sth, double *cth);

  DLOOP gcon[k][l] = 0.;

  bl_coord(X, &r, &th);

  // sincos(th, &sth, &cth);
  // sth = fabs(sth) + SMALL;
  cth = cos(th);
  sth = fabs(sin(th));
  if (sth < SMALL)
    sth = SMALL;

  irho2 = 1. / (r * r + a * a * cth * cth);
  // transformation for Kerr-Schild -> modified Kerr-Schild
  hfac = M_PI + (hslope)*M_PI * cos(2. * M_PI * X[2]);
  // hfac = M_PI;

  gcon[TT][TT] = -1. - 2. * r * irho2;
  gcon[TT][1] = 2. * irho2;

  gcon[1][TT] = gcon[TT][1];
  gcon[1][1] = irho2 * (r * (r - 2.) + a * a) / (r * r);
  gcon[1][3] = a * irho2 / r;

  gcon[2][2] = irho2 / (hfac * hfac);

  gcon[3][1] = gcon[1][3];
  gcon[3][3] = irho2 / (sth * sth);
}

void gcov_func(double *X, double gcov[][NDIM]) {
  int k, l;
  double sth, cth, s2, rho2;
  double r, th;
  double tfac, rfac, hfac, pfac;

  DLOOP gcov[k][l] = 0.;

  bl_coord(X, &r, &th);

  cth = cos(th);
  sth = fabs(sin(th));
  if (sth < SMALL)
    sth = SMALL;
  s2 = sth * sth;
  rho2 = r * r + a * a * cth * cth;

  /* transformation for Kerr-Schild -> modified Kerr-Schild */
  tfac = 1.;
  rfac = r - R0;
  hfac = M_PI + (hslope)*M_PI * cos(2. * M_PI * X[2]);
  //  hfac = M_PI;
  pfac = 1.;

  gcov[TT][TT] = (-1. + 2. * r / rho2) * tfac * tfac;
  gcov[TT][1] = (2. * r / rho2) * tfac * rfac;
  gcov[TT][3] = (-2. * a * r * s2 / rho2) * tfac * pfac;

  gcov[1][TT] = gcov[TT][1];
  gcov[1][1] = (1. + 2. * r / rho2) * rfac * rfac;
  gcov[1][3] = (-a * s2 * (1. + 2. * r / rho2)) * rfac * pfac;

  gcov[2][2] = rho2 * hfac * hfac;

  gcov[3][TT] = gcov[TT][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][3] = s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2)) * pfac * pfac;
}

#undef TT
#undef RR
#undef TH
#undef PH

/*

   connection calculated analytically for modified Kerr-Schild
   coordinates

   this gives the connection coefficient
   \Gamma^{i}_{j,k} = conn[..][i][j][k]
   where i = {1,2,3,4} corresponds to, e.g., {t,ln(r),theta,phi}
 */

void get_connection(double X[4], double lconn[4][4][4]) {

  double r1, r2, r3, r4, sx, cx;
  double th, dthdx2, dthdx22, d2thdx22, sth, cth, sth2, cth2, sth4, cth4, s2th,
      c2th;
  double a2, a3, a4, rho2, irho2, rho22, irho22, rho23, irho23, irho23_dthdx2;
  double fac1, fac1_rho23, fac2, fac3, a2cth2, a2sth2, r1sth2, a4cth4;
  /* required by broken math.h */
  void sincos(double th, double *sth, double *cth);

  r1 = exp(X[1]);
  r2 = r1 * r1;
  r3 = r2 * r1;
  r4 = r3 * r1;

  sincos(2. * M_PI * X[2], &sx, &cx);

  /* HARM-2D MKS */
  // th = M_PI*X[2] + 0.5*(1-hslope)*sx;
  // dthdx2 = M_PI*(1.+(1-hslope)*cx);
  // d2thdx22 = -2.*M_PI*M_PI*(1-hslope)*sx;

  /* HARM-3D MKS */
  // th = th_beg + th_len*X[2] + hslope*sx;
  // dthdx2 = th_len + 2.*M_PI*hslope*cx;
  // d2thdx22 = -4.*M_PI*M_PI*hslope*sx;	/* d^2(th)/dx2^2 */

  /*new for Hotaka data*/
  th = startx[2] + M_PI * X[2];
  //  dthdx2 = M_PI;
  dthdx2 = M_PI + (hslope)*M_PI * cos(2. * M_PI * X[2]);

  d2thdx22 = 0.0; /* d^2(th)/dx2^2 */

  dthdx22 = dthdx2 * dthdx2;

  sincos(th, &sth, &cth);
  sth2 = sth * sth;
  r1sth2 = r1 * sth2;
  sth4 = sth2 * sth2;
  cth2 = cth * cth;
  cth4 = cth2 * cth2;
  s2th = 2. * sth * cth;
  c2th = 2 * cth2 - 1.;

  a2 = a * a;
  a2sth2 = a2 * sth2;
  a2cth2 = a2 * cth2;
  a3 = a2 * a;
  a4 = a3 * a;
  a4cth4 = a4 * cth4;

  rho2 = r2 + a2cth2;
  rho22 = rho2 * rho2;
  rho23 = rho22 * rho2;
  irho2 = 1. / rho2;
  irho22 = irho2 * irho2;
  irho23 = irho22 * irho2;
  irho23_dthdx2 = irho23 / dthdx2;

  fac1 = r2 - a2cth2;
  fac1_rho23 = fac1 * irho23;
  fac2 = a2 + 2 * r2 + a2 * c2th;
  fac3 = a2 + r1 * (-2. + r1);

  lconn[0][0][0] = 2. * r1 * fac1_rho23;
  lconn[0][0][1] = r1 * (2. * r1 + rho2) * fac1_rho23;
  lconn[0][0][2] = -a2 * r1 * s2th * dthdx2 * irho22;
  lconn[0][0][3] = -2. * a * r1sth2 * fac1_rho23;

  // lconn[0][1][0] = lconn[0][0][1];
  lconn[0][1][1] = 2. * r2 * (r4 + r1 * fac1 - a4cth4) * irho23;
  lconn[0][1][2] = -a2 * r2 * s2th * dthdx2 * irho22;
  lconn[0][1][3] = a * r1 * (-r1 * (r3 + 2 * fac1) + a4cth4) * sth2 * irho23;

  // lconn[0][2][0] = lconn[0][0][2];
  // lconn[0][2][1] = lconn[0][1][2];
  lconn[0][2][2] = -2. * r2 * dthdx22 * irho2;
  lconn[0][2][3] = a3 * r1sth2 * s2th * dthdx2 * irho22;

  // lconn[0][3][0] = lconn[0][0][3];
  // lconn[0][3][1] = lconn[0][1][3];
  // lconn[0][3][2] = lconn[0][2][3];
  lconn[0][3][3] = 2. * r1sth2 * (-r1 * rho22 + a2sth2 * fac1) * irho23;

  lconn[1][0][0] = fac3 * fac1 / (r1 * rho23);
  lconn[1][0][1] = fac1 * (-2. * r1 + a2sth2) * irho23;
  lconn[1][0][2] = 0.;
  lconn[1][0][3] = -a * sth2 * fac3 * fac1 / (r1 * rho23);
  // lconn[1][1][0] = lconn[1][0][1];
  lconn[1][1][1] =
      (r4 * (-2. + r1) * (1. + r1) +
       a2 * (a2 * r1 * (1. + 3. * r1) * cth4 + a4cth4 * cth2 + r3 * sth2 +
             r1 * cth2 * (2. * r1 + 3. * r3 - a2sth2))) *
      irho23;
  lconn[1][1][2] = -a2 * dthdx2 * s2th / fac2;
  lconn[1][1][3] = a * sth2 * (a4 * r1 * cth4 + r2 * (2 * r1 + r3 - a2sth2) +
                               a2cth2 * (2. * r1 * (-1. + r2) + a2sth2)) *
                   irho23;

  // lconn[1][2][0] = lconn[1][0][2];
  // lconn[1][2][1] = lconn[1][1][2];
  lconn[1][2][2] = -fac3 * dthdx22 * irho2;
  lconn[1][2][3] = 0.;

  // lconn[1][3][0] = lconn[1][0][3];
  // lconn[1][3][1] = lconn[1][1][3];
  // lconn[1][3][2] = lconn[1][2][3];
  lconn[1][3][3] =
      -fac3 * sth2 * (r1 * rho22 - a2 * fac1 * sth2) / (r1 * rho23);

  lconn[2][0][0] = -a2 * r1 * s2th * irho23_dthdx2;
  lconn[2][0][1] = r1 * lconn[2][0][0];
  lconn[2][0][2] = 0.;
  lconn[2][0][3] = a * r1 * (a2 + r2) * s2th * irho23_dthdx2;

  // lconn[2][1][0] = lconn[2][0][1];
  lconn[2][1][1] = r2 * lconn[2][0][0];
  lconn[2][1][2] = r2 * irho2;
  lconn[2][1][3] =
      (a * r1 * cth * sth *
       (r3 * (2. + r1) +
        a2 * (2. * r1 * (1. + r1) * cth2 + a2 * cth4 + 2 * r1sth2))) *
      irho23_dthdx2;

  // lconn[2][2][0] = lconn[2][0][2];
  // lconn[2][2][1] = lconn[2][1][2];
  lconn[2][2][2] = -a2 * cth * sth * dthdx2 * irho2 + d2thdx22 / dthdx2;
  lconn[2][2][3] = 0.;

  // lconn[2][3][0] = lconn[2][0][3];
  // lconn[2][3][1] = lconn[2][1][3];
  // lconn[2][3][2] = lconn[2][2][3];
  lconn[2][3][3] = -cth * sth *
                   (rho23 + a2sth2 * rho2 * (r1 * (4. + r1) + a2cth2) +
                    2. * r1 * a4 * sth4) *
                   irho23_dthdx2;

  lconn[3][0][0] = a * fac1_rho23;
  lconn[3][0][1] = r1 * lconn[3][0][0];
  lconn[3][0][2] = -2. * a * r1 * cth * dthdx2 / (sth * rho22);
  lconn[3][0][3] = -a2sth2 * fac1_rho23;

  // lconn[3][1][0] = lconn[3][0][1];
  lconn[3][1][1] = a * r2 * fac1_rho23;
  lconn[3][1][2] = -2 * a * r1 * (a2 + 2 * r1 * (2. + r1) + a2 * c2th) * cth *
                   dthdx2 / (sth * fac2 * fac2);
  lconn[3][1][3] = r1 * (r1 * rho22 - a2sth2 * fac1) * irho23;

  // lconn[3][2][0] = lconn[3][0][2];
  // lconn[3][2][1] = lconn[3][1][2];
  lconn[3][2][2] = -a * r1 * dthdx22 * irho2;
  lconn[3][2][3] =
      dthdx2 * (0.25 * fac2 * fac2 * cth / sth + a2 * r1 * s2th) * irho22;

  // lconn[3][3][0] = lconn[3][0][3];
  // lconn[3][3][1] = lconn[3][1][3];
  // lconn[3][3][2] = lconn[3][2][3];
  lconn[3][3][3] = (-a * r1sth2 * rho22 + a3 * sth4 * fac1) * irho23;
}

/* stopping criterion for geodesic integrator */
/* K not referenced intentionally */

#define RMAX 1e5
#define ROULETTE 1.e3
int stop_criterion(struct of_photon *ph) {
  double wmin, X1min, X1max;

  wmin = WEIGHT_MIN; /* stop if weight is below minimum weight */

  X1min = log(Rh);         /* this is coordinate-specific; stop
                              at event horizon */
  X1max = log(RMAX * 1.1); /* this is coordinate and simulation
                              specific: stop at large distance */

  if (ph->X[1] < X1min)
    return 1;

  if (ph->X[1] > X1max) {
    if (ph->w < wmin) {
      if (monty_rand() <= 1. / ROULETTE) {
        ph->w *= ROULETTE;
      } else
        ph->w = 0.;
    }
    return 1;
  }

  if (ph->w < wmin) {
    if (monty_rand() <= 1. / ROULETTE) {
      ph->w *= ROULETTE;
    } else {
      ph->w = 0.;
      return 1;
    }
  }

  if (isnan(ph->K[0])) {
    fprintf(stderr, "isnan K[0]\n");
    return 1;
  }
  return (0);
}

/* criterion for recording photon */
int record_criterion(struct of_photon *ph) {
  //  const double X1max = log(RMAX);
  const double X1max = log(RMAX * 1.1);

  if (ph->X[1] > X1max)
    return (1);
  else
    return (0);
}

/* EPS really ought to be related to the number of
   zones in the simulation. */
#define EPS 0.08
//#define EPS   0.01

double stepsize(double X[NDIM], double K[NDIM], double tau) {

  double dl, dlx1, dlx2, dlx3;
  double idlx1, idlx2, idlx3;

  dlx1 = EPS / (fabs(K[1]) + SMALL * SMALL);
  dlx2 = EPS * GSL_MIN(X[2], 1 - X[2]) / (fabs(K[2]) + SMALL * SMALL);
  dlx3 = EPS / (fabs(K[3]) + SMALL * SMALL);

  idlx1 = 1. / (fabs(dlx1) + SMALL * SMALL);
  idlx2 = 1. / (fabs(dlx2) + SMALL * SMALL);
  idlx3 = 1. / (fabs(dlx3) + SMALL * SMALL);

  dl = 1. / (idlx1 + idlx2 + idlx3);

  return (dl);
}

/*
   record contribution of super photon to spectrum.

   This routine should make minimal assumptions about the
   coordinate system.

 */
void record_super_photon(struct of_photon *ph) {
  int iE, ix2, ix3, ii, jj, kk;
  double lE, dx2, dx3;

  if (isnan(ph->w) || isnan(ph->E)) {
    fprintf(stderr, "record isnan: %g %g\n", ph->w, ph->E);
    return;
  }
#pragma omp critical(MAXTAU)
  {
    if (ph->tau_scatt > max_tau_scatt)
      max_tau_scatt = ph->tau_scatt;
  }
  /* currently, bin in x2 coordinate */

  /* get theta bin, while folding around equator */
  /*
     dx2 = (stopx[2] - startx[2]) / (2. * N_THBINS);
     if (ph->X[2] < 0.5 * (startx[2] + stopx[2]))
     ix2 = (int) (ph->X[2] / dx2);
     else
     ix2 = (int) ((stopx[2] - ph->X[2]) / dx2);
   */

  double phi = fmod(ph->X[3], 2 * M_PI);
  if (phi < 0)
    phi += 2 * M_PI;

  /***************************************************/
  /*  get theta bin, without folding around equator   */
  dx2 = (stopx[2] - startx[2]) / (N_THBINS);
  ix2 = (int)(ph->X[2] / dx2);
  dx3 = (stopx[3] - startx[3]) / (N_PHIBINS);
  ix3 = (int)(phi / dx3);
  /***************************************************/

  /* check limits */
  if (ix2 < 0 || ix2 >= N_THBINS)
    return;
  if (ix3 < 0 || ix3 >= N_PHIBINS) {
    fprintf(stderr, "outside phi bin %e %d\n", ph->X[3], ix3);
    return;
  }

  double Xi[4], del[4], gcov[NDIM][NDIM], g;
  Xi[0] = 0;
  Xi[1] = ph->X1i;
  Xi[2] = ph->X2i;
  Xi[3] = ph->X3i;
  // Xi[3] = fmod(Xi[3],stopx[3]);
  // if(Xi[3] < 0.) Xi[3] += stopx[3];

  gcov_func(Xi, gcov);
  g = gdet_func(gcov);

  double dix1, dix2, dix3;
  /*dix1=(stopx[1]-startx[1])/N1;
     dix2=(stopx[2]-startx[2])/N2;
     dix3=(stopx[3]-startx[3])/N3;

     ii= (int)((Xi[1]-log(Rin)) / dix1);
     jj= (int)(Xi[2] / dix2);
     kk= (int)(Xi[3] /dix3);
   */
  Xtoijk(Xi, &ii, &jj, &kk, del);
  /* get energy bin */
  lE = log(ph->E);
  iE = (int)((lE - lE0) / dlE + 2.5) - 2; /* bin is centered on iE*dlE + lE0 */

  /* check limits */
  if (iE < 0 || iE >= N_EBINS)
    return;
  if (ii >= N1 || jj >= N2 || kk >= N3) {
    fprintf(stderr, "%d %d %d %d %d %d\n", ii, jj, kk, N1, N2, N3);
    fprintf(stderr, "%f %f %f\n", Xi[1], Xi[2], Xi[3]);
    return;
  }

  //#pragma omp atomic
  N_superph_recorded++;
  //#pragma omp atomic
  N_scatt += ph->nscatt;

  /* sum in photon */
  spect[ix2][ix3][iE].dNdlE += ph->w;
  spect[ix2][ix3][iE].dEdlE += ph->w * ph->E;
  spect[ix2][ix3][iE].tau_abs += ph->w * ph->tau_abs;
  spect[ix2][ix3][iE].tau_scatt += ph->w * ph->tau_scatt;
  spect[ix2][ix3][iE].X1iav += ph->w * ph->X1i;
  spect[ix2][ix3][iE].X2isq += ph->w * (ph->X2i); // * ph->X2i);
  spect[ix2][ix3][iE].X3fsq += ph->w * (ph->X3i); // ph->X3i);
  spect[ix2][ix3][iE].ne0 += ph->w * (ph->ne0);
  spect[ix2][ix3][iE].b0 += ph->w * (ph->b0);
  spect[ix2][ix3][iE].thetae0 += ph->w * (ph->thetae0);
  spect[ix2][ix3][iE].nscatt += ph->nscatt;
  spect[ix2][ix3][iE].nph += 1.;
  if (exp(lE) * ME * CL * CL / HPL < pow(10, 19) &&
      exp(lE) * ME * CL * CL / HPL > pow(10, 9)) {
    Xi_spec[ii][jj][kk] +=
        ph->w / (pow(L_unit, 3)); //*dix1*dix2*dix3*g);//WEIGHT_MIN/WEIGHT_MIN;
    ispec[ii][jj][kk] += ph->w * ph->E / (exp(lE) * ME * CL * CL / HPL);
  }
  // dx1dx2dx3sqrt(-g)
}

void omp_reduce_spect() {
  /* Combine partial spectra from each OpenMP process		*
   * Inefficient, but only called once so doesn't matter	*/

  int i, j, k;

#pragma omp critical(UPDATE_SPECT)
  {
    for (i = 0; i < N_THBINS; i++) {
      for (k = 0; k < N_PHIBINS; k++) {
        for (j = 0; j < N_EBINS; j++) {
          shared_spect[i][k][j].dNdlE += spect[i][k][j].dNdlE;
          shared_spect[i][k][j].dEdlE += spect[i][k][j].dEdlE;
          shared_spect[i][k][j].tau_abs += spect[i][k][j].tau_abs;
          shared_spect[i][k][j].tau_scatt += spect[i][k][j].tau_scatt;
          shared_spect[i][k][j].X1iav += spect[i][k][j].X1iav;
          shared_spect[i][k][j].X2isq += spect[i][k][j].X2isq;
          shared_spect[i][k][j].X3fsq += spect[i][k][j].X3fsq;
          shared_spect[i][k][j].ne0 += spect[i][k][j].ne0;
          shared_spect[i][k][j].b0 += spect[i][k][j].b0;
          shared_spect[i][k][j].thetae0 += spect[i][k][j].thetae0;
          shared_spect[i][k][j].nscatt += spect[i][k][j].nscatt;
          shared_spect[i][k][j].nph += spect[i][k][j].nph;
        }
      }
    }
    for (i = 0; i < N1; i++) {
      for (j = 0; j < N2; j++) {
        for (k = 0; k < N3; k++) {
          shared_Xi_spec[i][j][k] += Xi_spec[i][j][k];
          shared_ispec[i][j][k] += ispec[i][j][k];
        }
      }
    }
  }
#pragma omp barrier
#pragma omp master
  {
    for (i = 0; i < N_THBINS; i++) {
      for (k = 0; k < N_PHIBINS; k++) {
        for (j = 0; j < N_EBINS; j++) {
          spect[i][k][j].dNdlE = shared_spect[i][k][j].dNdlE;
          spect[i][k][j].dEdlE = shared_spect[i][k][j].dEdlE;
          spect[i][k][j].tau_abs = shared_spect[i][k][j].tau_abs;
          spect[i][k][j].tau_scatt = shared_spect[i][k][j].tau_scatt;
          spect[i][k][j].X1iav = shared_spect[i][k][j].X1iav;
          spect[i][k][j].X2isq = shared_spect[i][k][j].X2isq;
          spect[i][k][j].X3fsq = shared_spect[i][k][j].X3fsq;
          spect[i][k][j].ne0 = shared_spect[i][k][j].ne0;
          spect[i][k][j].b0 = shared_spect[i][k][j].b0;
          spect[i][k][j].thetae0 = shared_spect[i][k][j].thetae0;
          spect[i][k][j].nscatt = shared_spect[i][k][j].nscatt;
          spect[i][k][j].nph = shared_spect[i][k][j].nph;
        }
      }
    }
    double norm_ispec = 0;
    for (i = 0; i < N1; i++) {
      for (j = 0; j < N2; j++) {
        for (k = 0; k < N3; k++) {
          shared_ispec[i][j][k] =
              (ME * CL * CL) * shared_ispec[i][j][k] * (1. / dlE);
          norm_ispec += shared_ispec[i][j][k];
        }
      }
    }
    fprintf(stderr, "norm ispec = %g\n", norm_ispec);
    for (i = 0; i < N1; i++) {
      for (j = 0; j < N2; j++) {
        for (k = 0; k < N3; k++) {
          Xi_spec[i][j][k] = shared_Xi_spec[i][j][k];
          ispec[i][j][k] = shared_ispec[i][j][k] / norm_ispec;
          // if(Xi_spec[i][j][k]>0)
          //    fprintf(stderr, "%10.5g\n", Xi_spec[i][j][k]);
        }
      }
    }
  }
}

void sum_struct_ts(void *in, void *inout, int *len, MPI_Datatype *type) {
  /* ignore type, just trust that it's our struct type */

  struct of_spectrum *invals = in;
  struct of_spectrum *inoutvals = inout;
  int i, j, k;

  for (int i = 0; i < *len; i++) {

    inoutvals[i].dNdlE += invals[i].dNdlE;
    inoutvals[i].dEdlE += invals[i].dEdlE;
    inoutvals[i].tau_abs += invals[i].tau_abs;
    inoutvals[i].tau_scatt += invals[i].tau_scatt;
    inoutvals[i].X1iav += invals[i].X1iav;
    inoutvals[i].X2isq += invals[i].X2isq;
    inoutvals[i].X3fsq += invals[i].X3fsq;
    inoutvals[i].ne0 += invals[i].ne0;
    inoutvals[i].b0 += invals[i].b0;
    inoutvals[i].thetae0 += invals[i].thetae0;
    inoutvals[i].nscatt += invals[i].nscatt;
    inoutvals[i].nph += invals[i].nph;
  }
  return;
}

void mpi_reduce_spect() {
  /* Combine partial spectra from each OpenMP process		*
   * Inefficient, but only called once so doesn't matter	*/
  void sum_struct_ts(void *in, void *inout, int *len, MPI_Datatype *type);
  MPI_Datatype tstype;
  MPI_Op sumstruct;

  const int count = 13;
  int blocklens[count];
  MPI_Datatype types[count];
  MPI_Aint disps[count];

  for (int i = 0; i < count; i++) {
    types[i] = MPI_DOUBLE;
    blocklens[i] = 1;
  }

  disps[0] = offsetof(type_spectr, dNdlE);
  disps[1] = offsetof(type_spectr, dEdlE);
  disps[2] = offsetof(type_spectr, nph);
  disps[3] = offsetof(type_spectr, nscatt);
  disps[4] = offsetof(type_spectr, X1iav);
  disps[5] = offsetof(type_spectr, X2isq);
  disps[6] = offsetof(type_spectr, X3fsq);
  disps[7] = offsetof(type_spectr, tau_abs);
  disps[8] = offsetof(type_spectr, tau_scatt);
  disps[9] = offsetof(type_spectr, ne0);
  disps[10] = offsetof(type_spectr, thetae0);
  disps[11] = offsetof(type_spectr, b0);
  disps[12] = offsetof(type_spectr, E0);

  MPI_Type_create_struct(count, blocklens, disps, types, &tstype);
  MPI_Type_commit(&tstype);

  MPI_Op_create(sum_struct_ts, 1, &sumstruct);

  //      MPI_Reduce(&(spect[0][0][0]), &(shared_spect[0][0][0]),
  //      N_THBINS*N_PHIBINS*N_EBINS, tstype, MPI_SUM, 0,
  //               MPI_COMM_WORLD);
  for (int i = 0; i < N_THBINS; i++) {
    for (int j = 0; j < N_PHIBINS; j++) {
      for (int k = 0; k < N_EBINS; k++) {
        MPI_Reduce(&(spect[i][j][k]), &(shared_spect[i][j][k]), 1, tstype,
                   sumstruct, 0, MPI_COMM_WORLD);
      }
    }
  }
}

/*

   output spectrum to file

 */

#define SPECTRUM_FILE_NAME "grmonty_kappa4.spec"

void report_spectrum(double N_superph_made) {
  double JY_FACTOR = 1.e23;
  int i, j, k;
  double dx2, dx3, dOmega, nuLnu, tau_scatt, L, Inu, d;
  FILE *fp;
  //       strcat (args[2],".overall.spec");
  //     fprintf(stderr, "%s", args[2]);
  fp = fopen(SPECTRUM_FILE_NAME, "w");
  if (fp == NULL) {
    fprintf(stderr, "trouble opening spectrum file\n");
    exit(0);
  }
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  /* output */
  max_tau_scatt = 0.;
  L = 0.;
  for (i = 0; i < N_EBINS; i++) {

    /* output log_10(photon energy/(me c^2)) */
    //		fprintf(fp, "%10.5g ", (i * dlE + lE0) / M_LN10);
    fprintf(fp, "%10.5g ",
            exp(i * dlE + lE0) * ME * CL * CL / HPL); // nu in [Hz]

    //	  fprintf(fp, "%10.5g ", exp(i*dlE + lE0)*ME*CL*CL/HPL/2.417898940e17);
    ////nu in [keV]

    for (k = 0; k < N_PHIBINS; k++) {
      for (j = 0; j < N_THBINS; j++) {

        /* convert accumulated photon number in each bin
           to \nu L_\nu, in units of Lsun */
        // to fold the SED
        // dx2 = (stopx[2] - startx[2]) / (2. * N_THBINS);

        dx2 = (stopx[2] - startx[2]) / (1. * N_THBINS);
        dx3 = (stopx[3] - startx[3]) / (1. * N_PHIBINS);

        /* factor of 2 accounts for folding around equator */
        // to fold the sed
        // dOmega = 2. * dOmega_func(j * dx2, (j + 1) * dx2);

        dOmega = dOmega_func(j * dx2, (j + 1) * dx2, k * dx3, (k + 1) * dx3);
        d = 2.6228263e22; // SgrA*//RMAX*L_unit;

        Inu = 1 / (4. * M_PI * d * d) * JY_FACTOR;
        nuLnu = (ME * CL * CL) * (4. * M_PI / dOmega) * (1. / dlE);
        nuLnu *= shared_spect[j][k][i].dEdlE / ((double)world_size);
        // nuLnu /= LSUN;

        tau_scatt =
            shared_spect[j][k][i].tau_scatt / ((double)world_size) /
            (shared_spect[j][k][i].dNdlE / ((double)world_size) + SMALL);
        Inu *= nuLnu / (exp(i * dlE + lE0) * ME * CL * CL / HPL);

        // for time variability//

        /*                               fprintf(fp,"%10.5g %10.5g %10.5g %10.5g
           %10.5g %10.5g %10.5g ",
                                               nuLnu,
                                               spect[j][k][i].tau_abs /
           (spect[j][k][i].dNdlE + SMALL),
                                               tau_scatt,
                                               spect[j][k][i].X1iav /
           (spect[j][k][i].dNdlE + SMALL),
                                               ( (spect[j][k][i].X2isq /
           (spect[j][k][i].dNdlE + SMALL))),
                                               (spect[j][k][i].X3fsq /
           (spect[j][k][i].dNdlE + SMALL)),
                                               (spect[j][k][i].thetae0 /
           (spect[j][k][i].dNdlE + SMALL))
                                               );
         */
        // for quick seds
        fprintf(fp, "%10.5g %10.5g ", nuLnu, Inu);
        if (tau_scatt > max_tau_scatt)
          max_tau_scatt = tau_scatt;
        L += nuLnu * dOmega * dlE;
      }
    }
    fprintf(fp, "\n");

    //		exit(1);
  }
  fprintf(
      stderr,
      //		"luminosity %g, dMact %g [Msun/yr], efficiency %g,
      //L/Ladv %g, max_tau_scatt %g\n",
      "luminosity %g, dMact %g [Msun/yr], mdot[mdot edd]= %g,efficiency %g, "
      "L/Ladv %g, max_tau_scatt %g\n",
      L, dMact * M_unit / T_unit / (MSUN / YEAR),
      dMact * M_unit / T_unit / (MSUN / YEAR) / 3.3,
      L / (dMact * M_unit * CL * CL / T_unit),
      L * LSUN / (Ladv * M_unit * CL * CL / T_unit), max_tau_scatt);
  fprintf(stderr, "\n");
  fprintf(stderr, "N_superph_made: %g\n", N_superph_made);
  fprintf(stderr, "N_superph_recorded: %g\n", N_superph_recorded_total);

  fclose(fp);
}
/*
   void report_xi_backup(double N_superph_made, char *args[]){
    int i, j, k;
    double X[4],th,phi,r;
    FILE *fp;
    strcat (args[2],"_xi.txt");
    fprintf(stderr, "%s", args[2]);
    fp = fopen(args[2], "w");
    if (fp == NULL) {
        fprintf(stderr, "trouble opening xi file\n");
        exit(0);
    }

    for(i = 0; i < N1; i++){
        X[1] = startx[1] + (i+0.5)*dx[1];
        for(j = 0; j < N2; j++){
            X[2] = startx[2] + (j+0.5)*dx[2];
            for(k = 0; k < N3; k++){
                X[3] = startx[3] + (k+0.5)*dx[3];
                //delV=dx[1]*dx[2]*dx[3]*g[i][j];
                r  = exp(X[1]) + R0;
                //back to the cureent set up, th_len=pi, adn hslope=0 in our
   case
                th = M_PI*X[2]; //th_len*X[2] + th_beg + hslope*sin(2. * M_PI *
   X[2]);
                phi = fmod(X[3],stopx[3]);
                if(phi < 0.) phi = stopx[3]+phi;
                if(log(Xi_spec[i][j][k])>5. && (k==0 || k==0+N3/2))
                    fprintf(fp,"%10.5g %10.5g %10.5g %10.5g\n",
   r,th,phi,log(Xi_spec[i][j][k]+1));
            }
        }
    }
    fclose(fp);
   }*/

#define N 3
#define NN 1
void report_xi(double N_superph_made, char *args[]) {

  FILE *fp, *fp2;
  int i, j, k;
  double dum_r, dum_r1, dum_r2, dum_r3;
  double X[4], th, phi, r;
  float dum_rho[NN];
  float dum_rho1[NN];
  float dum_rho2[NN];
  float dum_rho3[NN];
  float dum_coord[N];
  strcat(args[2], "_xi.vtk");
  /* Opening the outputfile*/
  fp = fopen(args[2], "w");

  /* There are five basic parts to the VTK file format
     1. Write file version and identifier (Header)*/
  fprintf(fp, "# vtk DataFile Version 2.0\n");

  /* 2. Tilte */
  fprintf(fp, "GRMHD Simulation Result\n");

  /* 3. Data type (ASCII or BINARY)*/
  fprintf(fp, "ASCII\n");

  /* 4.  Dataset structure */

  fprintf(fp, "DATASET STRUCTURED_GRID\n");
  fprintf(fp, "DIMENSIONS %d %d %d\n", N1, N2, N3);
  fprintf(fp, "POINTS %d float\n", N1 * (N2)*N3);

  for (k = 0; k < N3; k++) {
    for (j = 0; j < N2; j++) {
      for (i = 0; i < N1; i++) {
        X[1] = startx[1] + (i + 0.5) * dx[1];
        X[2] = startx[2] + (j + 0.5) * dx[2];
        X[3] = startx[3] + (k + 0.5) * dx[3];
        // delV=dx[1]*dx[2]*dx[3]*g[i][j];
        r = exp(X[1]) + R0;
        // back to the cureent set up, th_len=pi, adn hslope=0 in our case
        th = M_PI * X[2]; // th_len*X[2] + th_beg + hslope*sin(2. * M_PI *
                          // X[2]);
        phi = fmod(X[3], stopx[3]);
        if (phi < 0.)
          phi = stopx[3] + phi;

        if (k == 0 || k == N3 - 1) {
          fprintf(fp, "%g %g %g\n", r * sin(th) * cos(phi), r * cos(th), 0.0);
          Xi_spec[i][j][k] = (Xi_spec[i][j][N3 - 1] + Xi_spec[i][j][0]) / 2.;
        } else
          fprintf(fp, "%g %g %g\n", r * sin(th) * cos(phi), r * cos(th),
                  r * sin(th) * sin(phi));
      }
    }
  }
  /*now dumps data*/

  /* 5 . Data */
  // You need to add newline character for safty after writing data in binary
  fprintf(fp, "\nPOINT_DATA %d\n", N1 * (N2)*N3);

  /********* Write log of density density ********/
  fprintf(fp, "SCALARS log10w float\n");
  fprintf(fp, "LOOKUP_TABLE default\n");

  for (k = 0; k < N3; k++) {
    for (j = 0; j < N2; j++) {
      for (i = 0; i < N1; i++) {
        dum_r = log10(Xi_spec[i][j][k] + 1);
        fprintf(fp, "%g\n", dum_r);
      }
    }
  }

  /********* Write log of density density ********/
  fprintf(fp, "SCALARS ispec float\n");
  fprintf(fp, "LOOKUP_TABLE default\n");

  for (k = 0; k < N3; k++) {
    for (j = 0; j < N2; j++) {
      for (i = 0; i < N1; i++) {
        dum_r = ispec[i][j][k];
        fprintf(fp, "%g\n", dum_r);
      }
    }
  }

  fclose(fp);

  return;
}
