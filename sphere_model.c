/*

   sphere model specification routines
   //model from B. Ryan and C. Gammie

 */

#include "decs.h"
#define global
#include "sphere_model.h"
#undef global

struct of_spectrum ***spect;
extern struct of_spectrum ***shared_spect;

#pragma omp threadprivate(spect)

/*

   encapsulates all initialization routines

 */

void init_model(char *args[]) {
  /* find dimensional quantities from black hole
     mass and its accretion rate */
  set_units(args[3]);

  fprintf(stderr, "getting simulation data...\n");
  init_sphere_data(args[2]); /* read in HARM simulation data */
  #pragma omp parallel
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
     }


  init_geometry();
  fprintf(stderr, "init geom done\n");
  /* make look-up table for hot cross sections */
  init_hotcross();
  fprintf(stderr, "hot cross done\n");
  /* make table for solid angle integrated emissivity and K2 */
  init_emiss_tables();
  fprintf(stderr, "emiss tables done\n");

  /* make table for superphoton weights */
  init_weight_table();
  fprintf(stderr, "weight tables done\n");

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
   // fprintf(stderr,"making some superphotons...\n");
   while (n2gen <= 0) {
     // fprintf(stderr,"n2gen is %d %d\n",zone_i,n2gen);
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

/*

   produces a bias (> 1) for probability of Compton scattering
   as a function of local temperature

 */

double bias_func(double Te, double w) {
  double bias, max;

  max = 0.5 * w / WEIGHT_MIN;

  bias = Te * Te / (5. * max_tau_scatt);
  // bias = 100. * Te * Te / (bias_norm * max_tau_scatt);

  if (bias > max)
    bias = max;

  return  bias ;//fmax(bias*0.15,1.); /// TP_OVER_TE;
}

/*

   these supply basic model data to grmonty

 */

void get_fluid_zone_new(int i, int j, int k, double *Ne, double *Thetae, double *B,
       	double Ucon[NDIM], double Bcon[NDIM],  int *ACCZONE)
{
  // get grid zone center
  double X[4] = { 0. };
  X[0] = 0.;
  X[1] = (i+0.5) * dx[1] + startx[1];
  X[2] = (j+0.5) * dx[2] + startx[2];
  X[3] = (k+0.5) * dx[3] + startx[3];

  // fill gcov
  double gcov[4][4];
  gcov_func(X, gcov);

  // create dummy vectors
  double Ucov[4] = { 0. };
  double Bcov[4] = { 0. };

  double sigma,beta,dx_local;
  int igrid_c;

  // call into analytic model
  get_fluid_params(X, gcov, Ne, Thetae, B,&sigma,&beta, Ucon, Ucov, Bcon, Bcov, ACCZONE,&dx_local,&igrid_c);
}

void get_fluid_zone(int i, int j, int k,double *Ne, double *Thetae, double *B,
                    double Ucon[NDIM], double Bcon[NDIM], int *ACCZONE) {
  double Ucov[NDIM], Bcov[NDIM];
  double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
  double sig;
  double X[NDIM];
  coord(i, j, X);
  *Ne = p[KRHO][i][j][k] * Ne_unit;
  *Thetae = (p[UU][i][j][k] / p[KRHO][i][j][k]) * Thetae_unit;
  Bp[1] = p[B1][i][j][k];
  Bp[2] = p[B2][i][j][k];
  Bp[3] = p[B3][i][j][k];

  Vcon[1] = p[U1][i][j][k];
  Vcon[2] = p[U2][i][j][k];
  Vcon[3] = p[U3][i][j][k];

  // Get Ucov
  VdotV = 0.;
  for (int l = 1; l < NDIM; l++)
    for (int m = 1; m < NDIM; m++)
      VdotV += geom[i][j].gcov[l][m] * Vcon[l] * Vcon[m];
  Vfac = sqrt(-1. / geom[i][j].gcon[0][0] * (1. + fabs(VdotV)));
  Ucon[0] = -Vfac * geom[i][j].gcon[0][0];
  for (int l = 1; l < NDIM; l++)
    Ucon[l] = Vcon[l] - Vfac * geom[i][j].gcon[0][l];
  lower(Ucon, geom[i][j].gcov, Ucov);
  // Get B and Bcov
  UdotBp = 0.;
  for (int l = 1; l < NDIM; l++)
    UdotBp += Ucov[l] * Bp[l];
  Bcon[0] = UdotBp;
  for (int l = 1; l < NDIM; l++)
    Bcon[l] = (Bp[l] + Ucon[l] * UdotBp) / Ucon[0];

  lower(Bcon, geom[i][j].gcov, Bcov);

  *B = sqrt(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] + Bcon[2] * Bcov[2] +
            Bcon[3] * Bcov[3]) *
       B_unit +1e-30;

  sig = pow(*B / B_unit, 2) / (*Ne / Ne_unit);
  if (sig > 1.) {
    *Thetae = SMALL; // *Ne = 1.e-10*Ne_unit;
    *Ne=0;
   // fprintf(stderr, "yup its high!\n");
  }
}

void check_acczone(int *ACCZONE, double var1, double var2) {

  if (var1 > 1.02 && var2 < 1.) {
    *ACCZONE = 1;
  }
}

int get_fluid_params_new(double X[NDIM], double gcov[NDIM][NDIM], double *Ne, double *Thetae, double *B, double *sigma, double *beta,
                     double Ucon[NDIM], double Ucov[NDIM], double Bcon[NDIM], double Bcov[NDIM], int *ACCZONE, double *dx_local, int *igrid_c){

double gcon[NDIM][NDIM];
  // model parameters // TODO, maybe load these from model parameters
double  MODEL_R_0 = 100.;
double  MODEL_TAU_0 = 1.e-5;
double  MODEL_BETA_0 = 20.;
double  MODEL_THETAE_0 = 10.;
double  MODEL_TP_OVER_TE = 3.;
double  MODEL_GAM = 13./9.;
double   MODEL_MBH = 4.1e6;

  // derive model Ne (in cgs)
double  model_Ne0 = MODEL_TAU_0 / SIGMA_THOMSON / MODEL_R_0 / L_unit;

  // since B = B(pressure), we need to specify the thermodynamics to
  // find pressure = pressure(Thetae)
  double gam = MODEL_GAM;
  double game = 4./3;
  double gamp = 5./3;

  // as implemented in RAPTOR + kmonty
double  THETAE_UNIT = MP/ME * (gam-1.) / (1. + MODEL_TP_OVER_TE);

  // now we can find B (again, in gauss)
  double model_B0 = CL * sqrt(8 * M_PI * (gam-1.) * (MP+ME) / MODEL_BETA_0) * sqrt( model_Ne0 * MODEL_THETAE_0 ) / sqrt( THETAE_UNIT );

  if (X[1] > 100) {
    *Ne = 0.;
    *Thetae = 0;
    *B = 0;
    return;
  }

  *Ne = model_Ne0;
  *Thetae = MODEL_THETAE_0;
  *B = model_B0;

  double r = X[1];
  double h = X[2];

  Ucon[0] = 1;
  Ucon[1] = 0.;
  Ucon[2] = 0.;
  Ucon[3] = 0.;

  Bcon[0] = 0.;
  Bcon[1] = model_B0 * cos(h);
  Bcon[2] = - model_B0 * sin(h) / (r + 1.e-8);
  Bcon[3] = 0.;

  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  lower(Ucon, gcov, Ucov);
  lower(Bcon, gcov, Bcov);

}

int get_fluid_params(double X[NDIM], double gcov[NDIM][NDIM], double *Ne, double *Thetae, double *B, double *sigma, double *beta,
                     double Ucon[NDIM], double Ucov[NDIM], double Bcon[NDIM], double Bcov[NDIM], int *ACCZONE, double *dx_local, int *igrid_c){
  int i, j,k;
  double rho, uu;
  double del[NDIM];
  double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
  double gcon[NDIM][NDIM], coeff[4];
  double interp_scalar(double ***var, int i, int j, double del[4]);
  double sig;

  if (X[1] < startx[1] || X[1] > stopx[1] || X[2] < startx[2] ||
      X[2] > stopx[2]) {

    *Ne = 0;

    return 0;
  }

  Xtoijk(X, &i, &j,&k, del);

  coeff[0] = (1. - del[1]) * (1. - del[2]);
  coeff[1] = (1. - del[1]) * del[2];
  coeff[2] = del[1] * (1. - del[2]);
  coeff[3] = del[1] * del[2];

  rho = interp_scalar(p[KRHO],i,j,coeff);//p[KRHO][i][j][k];
  uu = interp_scalar(p[UU],i,j,coeff);;//p[UU][i][j][k];
  *Ne = rho * Ne_unit;
  *Thetae = (uu) / (rho)*Thetae_unit;

  Bp[1] = interp_scalar(p[B1],i,j,coeff);//p[B1][i][j][k];
  Bp[2] = interp_scalar(p[B2],i,j,coeff);//p[B2][i][j][k];
  Bp[3] = interp_scalar(p[B3],i,j,coeff);//p[B3][i][j][k];

  Vcon[1] = interp_scalar(p[U1],i,j,coeff);//p[U1][i][j][k];
  Vcon[2] = interp_scalar(p[U2],i,j,coeff);//p[U2][i][j][k];
  Vcon[3] = interp_scalar(p[U3],i,j,coeff);//p[U3][i][j][k];
  
  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  // Get Ucov
  VdotV = 0.;
  for (int i = 1; i < NDIM; i++)
    for (int j = 1; j < NDIM; j++)
      VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
  Vfac = sqrt(-1. / gcon[0][0] * (1. + fabs(VdotV)));
  Ucon[0] = -Vfac * gcon[0][0];
  for (int i = 1; i < NDIM; i++)
    Ucon[i] = Vcon[i] - Vfac * gcon[0][i];
  lower(Ucon, gcov, Ucov);

  // Get B and Bcov
  UdotBp = 0.;
  for (int i = 1; i < NDIM; i++)
    UdotBp += Ucov[i] * Bp[i];
  Bcon[0] = UdotBp;
  for (int i = 1; i < NDIM; i++)
    Bcon[i] = (Bp[i] + Ucon[i] * UdotBp) / Ucon[0];

  lower(Bcon, gcov, Bcov);

  *B = sqrt(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] + Bcon[2] * Bcov[2] +
            Bcon[3] * Bcov[3]) *
       B_unit +1e-30;


  sig = pow(*B / B_unit, 2) / (*Ne / Ne_unit);
  if (sig > 1.){
    *Thetae = SMALL; // *Ne = 1.e-10*Ne_unit;
   *Ne=0;
   }  
//
  return 1;
}

/*
   Current metric: miknowski spherical coordinates
 */

void gcov_func(double X[NDIM], double gcov[NDIM][NDIM]) {
  MUNULOOP gcov[mu][nu] = 0.;

  gcov[0][0] = -1.;

  gcov[1][1] = 1.;

  gcov[2][2] = pow(X[1], 2);

  gcov[3][3] = pow(X[1] * sin(X[2]), 2);
}


void gcon_func(double gcov[][NDIM], double gcon[][NDIM]) {
  MUNULOOP gcon[mu][nu] = 0.;
	gcon[0][0]=1/gcov[0][0];
	gcon[1][1]=1/gcov[1][1];
	gcon[2][2]=1/(gcov[2][2]+1e-10);
	gcon[3][3]=1/(gcov[3][3]+1e-10);
  // done!
}

/*

   connection calculated analytically for modified Kerr-Schild
    coordinates


   this gives the connection coefficient
   \Gamma^{i}_{j,k} = conn[..][i][j][k]
   where i = {1,2,3,4} corresponds to, e.g., {t,ln(r),theta,phi}
 */

#define DEL (1.e-7)
void get_connection(double X[NDIM], double conn[NDIM][NDIM][NDIM]) {
  double tmp[NDIM][NDIM][NDIM];
  double Xh[NDIM], Xl[NDIM];
  double gcon[NDIM][NDIM];
  double gcov[NDIM][NDIM];
  double gh[NDIM][NDIM];
  double gl[NDIM][NDIM];

  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  for (int k = 0; k < NDIM; k++) {
    for (int l = 0; l < NDIM; l++)
      Xh[l] = X[l];
    for (int l = 0; l < NDIM; l++)
      Xl[l] = X[l];
    Xh[k] += DEL;
    Xl[k] -= DEL;
    gcov_func(Xh, gh);
    gcov_func(Xl, gl);

    for (int i = 0; i < NDIM; i++) {
      for (int j = 0; j < NDIM; j++) {
        conn[i][j][k] = (gh[i][j] - gl[i][j]) / (Xh[k] - Xl[k]);
      }
    }
  }

  // Rearrange to find \Gamma_{ijk}
  for (int i = 0; i < NDIM; i++)
    for (int j = 0; j < NDIM; j++)
      for (int k = 0; k < NDIM; k++)
        tmp[i][j][k] = 0.5 * (conn[j][i][k] + conn[k][i][j] - conn[k][j][i]);

  // G_{ijk} -> G^i_{jk}
  for (int i = 0; i < NDIM; i++) {
    for (int j = 0; j < NDIM; j++) {
      for (int k = 0; k < NDIM; k++) {
        conn[i][j][k] = 0.;
        for (int l = 0; l < NDIM; l++)
          conn[i][j][k] += gcon[i][l] * tmp[l][j][k];
      }
    }
  }
}
#undef DEL

/* stopping criterion for geodesic integrator */
/* K not referenced intentionally */

#define RMAX 10000.
#define ROULETTE 1.e4
int stop_criterion(struct of_photon *ph) {
  double wmin, X1min, X1max;

  wmin = WEIGHT_MIN; /* stop if weight is below minimum weight */

  // Stop at event horizon
  X1min = 0.;

  // Stop at large distance
  X1max = 1.1*RMAX;

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

  return (0);
}

/* criterion for recording photon */

int record_criterion(struct of_photon *ph) {
  const double X1max = 1.1*RMAX;
  /* this is coordinate and simulation
     specific: stop at large distance */

  if (ph->X[1] > X1max)
    return (1);

  else
    return (0);
}

/* EPS really ought to be related to the number of
   zones in the simulation. */
#define EPS 0.05
//#define EPS   0.01

double stepsize(double X[NDIM], double K[NDIM], double tau) {
  double dl, dlx1,dlxx1, dlx2, dlx3;
  double idlx1,idlxx1, idlx2, idlx3;

  dlx1 = EPS * X[1] / (fabs(K[1]) + SMALL); // *fmin(fmax(exp(-tau),1e-6),1);
  dlxx1 = EPS/ (fabs(K[1]) + SMALL); // *fmin(fmax(exp(-tau),1e-6),1);
  dlx2 = EPS * GSL_MIN(X[2], stopx[2] - X[2]) / (fabs(K[2]) + SMALL);
  dlx3 = EPS / (fabs(K[3]) + SMALL);
  idlx1 = 1. / (fabs(dlx1) + SMALL);
  idlxx1 = 1. / (fabs(dlxx1) + SMALL);
  idlx2 = 1. / (fabs(dlx2) + SMALL);
  idlx3 = 1. / (fabs(dlx3) + SMALL);

  dl = 1. / (idlx1 + idlx2 + idlx3);

  if(X[1]>stopx[1]*1.05 || (X[1] < stopx[1] * 0.9))
      dl = 0.2/ (idlxx1 + idlx2 + idlx3);

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
#if(FOLDING)
     dx2 = (stopx[2] - startx[2]) / (2. * N_THBINS);
     if (ph->X[2] < 0.5 * (startx[2] + stopx[2]))
     ix2 = (int) (ph->X[2] / dx2);
     else
     ix2 = (int) ((stopx[2] - ph->X[2]) / dx2);
#else
  dx2 = (stopx[2] - startx[2]) / (N_THBINS);
  ix2 = (int)(ph->X[2] / dx2);
#endif

  double phi = fmod(ph->X[3], 2 * M_PI);
  if (phi < 0)
    phi += 2 * M_PI;

  /***************************************************/
  /*  get theta bin, without folding around equator   */
  dx3 = (stopx[3] - startx[3]) / (N_PHIBINS);
  ix3 = (int)(phi / dx3 );
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
  g = gdet_func(gcov,Xi);

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

  }
}

/*

   output spectrum to file

 */

 #define SPECTRUM_FILE_NAME "grmonty_pwl.spec"

 void report_spectrum(double N_superph_made) {
   double JY_FACTOR = 1.e23;
   int i, j, k;
   double dx2, dx3, dOmega, nuLnu, tau_scatt, L, Inu, d;
   FILE *fp;
   char filename[100];
   //       strcat (args[2],".overall.spec");
   //     fprintf(stderr, "%s", args[2]);
   sprintf(filename,"spec_sphere_%d.dat",index);
   fp = fopen(filename, "w");
   if (fp == NULL) {
     fprintf(stderr, "trouble opening spectrum file\n");
     exit(0);
   }
 #if MPI
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
 #endif
   dx2 = (stopx[2] - startx[2]) / (N_THBINS); 
   //fprintf(stderr, "world size is %d\n", world_size);
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
#if(FOLDING)
         dx2 = (stopx[2] - startx[2]) / (2. * N_THBINS*M_PI);
#else
         dx2 = (1.) / (1. * N_THBINS);
#endif
         dx3 = (2 * M_PI) / (1. * N_PHIBINS);

         /* factor of 2 accounts for folding around equator */
         // to fold the sed
#if(FOLDING)
         dOmega = 2. * dOmega_func(j * dx2, (j + 1) * dx2,k * dx3, (k + 1) * dx3);
#else
         dOmega = dOmega_func(j * dx2, (j + 1) * dx2, k * dx3, (k + 1) * dx3);
#endif
         d = 2.6228263e22;
     //5.061e25; M87
         //                                d= 2.6228263e22; //SgrA*//RMAX*L_unit;

         Inu = 1. / (4. * M_PI * d * d) * JY_FACTOR;
         nuLnu = (ME * CL * CL) * (4. * M_PI / dOmega) * (1. / dlE);
 #if MPI
         nuLnu *= shared_spect[j][k][i].dEdlE / ((double)world_size);
 #endif

 #if OPENMP
         nuLnu *= shared_spect[j][k][i].dEdlE;
 #endif
         // nuLnu /= LSUN;
 #if MPI
         tau_scatt =
             shared_spect[j][k][i].tau_scatt / ((double)world_size) /
             (shared_spect[j][k][i].dNdlE / ((double)world_size) + SMALL);
 #endif

 #if OPENMP
         tau_scatt =
             shared_spect[j][k][i].tau_scatt /
             (shared_spect[j][k][i].dNdlE + SMALL);
 #endif
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
         L += nuLnu * dOmega * dlE/(4.*M_PI);
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
   fprintf(stderr, "N_superph_recorded: %g\n", N_superph_recorded);

   fclose(fp);
 }
