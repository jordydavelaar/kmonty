#include "decs.h"
#include <stdio.h>
#include <stdlib.h>
/* HDF5 v1.6 API */
//#include <H5LT.h>

/* HDF5 v1.8 API */

/* Harm3d globals */

extern double ****bcon;
extern double ****bcov;
extern double ****ucon;
extern double ****ucov;
extern double ****p;
extern double ***ne;
extern double ***thetae;
extern double ***b;

void init_bhac3d_data(char *fname) {
  int i, j, k, l, m;
  double X[NDIM], UdotU, ufac;
  double gcov[NDIM][NDIM], gcon[NDIM][NDIM], g;
  double Ucon[NDIM], Ucov[NDIM];
  double dV, Thetae, V, Be;
  double Th_unit, two_temp_gam;
  double UdotBp, Bcon[NDIM], Bcov[NDIM], Bp[NDIM];
  double bsq, beta, beta_trans, b2;
  double trat;
  double th_end, th_cutout;
  FILE *fp;
  double x[4];
  int prec;

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_rank == 0) {
    fprintf(stderr, "\nGRMHD FILE\n");
    fprintf(stderr, "File reading: %s\n", fname);
  }

  /* header variables not used except locally */
  double t, rho_floor, p_floor;
  char *buffer;
  int nstep, metric_type, code_type, Ndim;

  fp = fopen(fname, "r");

  if (fp == NULL) {
    fprintf(stderr, "can't open sim data file %s\n", fname);
    exit(1);
  }

  size_t buffersize = 50;
  buffer = (char *)malloc(buffersize * sizeof(char));

  double buffer2;
  fread(&N1, sizeof(int), 1, fp);
  fread(&N2, sizeof(int), 1, fp);
  fread(&N3, sizeof(int), 1, fp);
  fread(&buffer2, sizeof(double), 1, fp);
  t = (double)buffer2;
  fread(&buffer2, sizeof(double), 1, fp);
  a = (double)buffer2;
  fread(&buffer2, sizeof(double), 1, fp);
  startx[1] = (double)buffer2;
  fread(&buffer2, sizeof(double), 1, fp);
  startx[2] = (double)buffer2;
  fread(&buffer2, sizeof(double), 1, fp);
  startx[3] = (double)buffer2;
  fread(&buffer2, sizeof(double), 1, fp);
  stopx[1] = (double)buffer2;
  fread(&buffer2, sizeof(double), 1, fp);
  stopx[2] = (double)buffer2;
  fread(&buffer2, sizeof(double), 1, fp);
  stopx[3] = (double)buffer2;
  fread(&buffer2, sizeof(double), 1, fp);
  gam = (double)buffer2;
  fread(&buffer2, sizeof(double), 1, fp);
  hslope = (double)buffer2;
  fread(&metric_type, sizeof(int), 1, fp);
  fread(&code_type, sizeof(int), 1, fp);
  fread(&Ndim, sizeof(int), 1, fp);
  fread(&prec, sizeof(int), 1, fp);

  if (world_rank == 0) {
    fprintf(stderr, "\nHEADER\n");
    fprintf(stderr, "time %e ", t);
    fprintf(stderr, "spin %e ", a);
    fprintf(stderr, "gamma %e ", gam);
    fprintf(stderr, "hslope %e ", hslope);
    fprintf(stderr, "N1 %d N2 %d N3 %d\n", N1, N2, N3);
    fprintf(stderr, "Sim range x1, x2:  %e %e, %e %e\n", startx[1], stopx[1],
            startx[2], stopx[2]);
  }

  //	hslope = hslope;
  n_within_horizon = 0;
  init_storage();
  th_cutout = 0.;
  Rin = 0;
  Rout = 1000.;

  startx[2] /= M_PI;
  stopx[2] /= M_PI;
  th_beg = startx[2];

  dx[1] = (stopx[1] - startx[1]) / N1;
  dx[2] = (stopx[2] - startx[2]) / (N2);
  dx[3] = (stopx[3] - startx[3]) / N3;

  //#if(metric==MMKS)
  //      hslope=1.;
  //#endif
  R0 = 0;
  // printf("done reading header, %lf\n",hslope);
  //    hslope = 1 - hslope;
  // nominal non-zero values for axisymmetric simulations
  stopx[0] = 1.;

  dx[0] = 1.;
  dV = dx[1] * dx[2] * dx[3];

  dMact = Ladv = 0.;
  bias_norm = 0.;
  V = 0.;
  init_geometry();
  two_temp_gam =
      0.5 * ((1. + 2. / 3. * (TP_OVER_TE + 1.) / (TP_OVER_TE + 2.)) + gam);
  Thetae_unit = (two_temp_gam - 1.) * (MP / ME) / (1. + TP_OVER_TE);

  double sizebit;
  double bufferd3;
  float bufferf3;
  if (prec == 0) {
    sizebit = sizeof(double);
  } else {
    sizebit = sizeof(float);
  }
  fseek(fp, -NPRIM * N1 * N2 * N3 * sizebit, SEEK_END);
  int count = 0;
  for (int k = 0; k < N3; k++) {
    for (int j = 0; j < N2; j++) {
      for (int i = 0; i < N1; i++) {
        for (int l = 0; l < NPRIM; l++) {
          if (prec == 0) {
            fread(&bufferd3, sizebit, 1, fp);
            p[l][i][j][k] = (double)bufferd3;
          } else {
            fread(&bufferf3, sizebit, 1, fp);
            p[l][i][j][k] = (double)bufferf3;
          }
        }
        p[UU][i][j][k] = p[UU][i][j][k] / (4. / 3. - 1);
        p[U2][i][j][k] = p[U2][i][j][k] / M_PI;
        //				fprintf(stderr,"rho
        //%e\n",p[KRHO][i][j][k]);
      }
    }
  }
  for (int k = 0; k < NPRIM; k++) {
    for (int j = 0; j < N2; j++) {
      for (int i = 0; i < N1; i++) {
        double p0 = p[k][i][j][0];
        double pN3 = p[k][i][j][N3 - 1];
        p[k][i][j][0] = (p0 + pN3) / 2.;
        p[k][i][j][N3 - 1] = (p0 + pN3) / 2.;
      }
    }
  }

  X[0] = 0.;
  X[3] = 0.;

  dMact = Ladv = 0.;
  bias_norm = 0.;
  V = 0.;

  for (i = 0; i < N1; i++) {
    X[1] = startx[1] + (i + 0.5) * dx[1];
    for (j = 0; j < N2; j++) {
      X[2] = startx[2] + (j + 0.5) * dx[2];

      gcov_func(X, gcov);
      gcon_func(X, gcon);
      g = gdet_func(gcov);
      for (k = 0; k < N3; k++) {
        UdotU = 0.;
        for (l = 1; l < NDIM; l++)
          for (m = 1; m < NDIM; m++)
            UdotU +=
                gcov[l][m] * p[U1 + l - 1][i][j][k] * p[U1 + m - 1][i][j][k];
        ufac = sqrt(-1. / gcon[0][0] * (1 + fabs(UdotU)));
        Ucon[0] = -ufac * gcon[0][0];
        for (l = 1; l < NDIM; l++)
          Ucon[l] = p[U1 + l - 1][i][j][k] - ufac * gcon[0][l];
        lower(Ucon, gcov, Ucov);
        /*calculate local electron temperature depending on the geometry*/
        Be = -(1. + p[UU][i][j][k] / p[KRHO][i][j][k] * gam) * Ucov[0];
        // depending on the magnetic properties
        Bp[1] = p[B1][i][j][k];
        Bp[2] = p[B2][i][j][k];
        Bp[3] = p[B3][i][j][k];
        UdotBp = 0.;
        for (l = 1; l < NDIM; l++)
          UdotBp += Ucov[l] * Bp[l];
        Bcon[0] = UdotBp;
        for (l = 1; l < NDIM; l++)
          Bcon[l] = (Bp[l] + Ucon[l] * UdotBp) / Ucon[0];
        lower(Bcon, geom[i][j].gcov, Bcov);
        bsq = Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] + Bcon[2] * Bcov[2] +
              Bcon[3] * Bcov[3] + 1e-20;
        beta = p[UU][i][j][k] * (gam - 1.) / 0.5 / bsq;
        // beta=p[UU][i][j][k]/(bsq +1e-20); // beta plasma
        beta_trans = 1.;
        b2 = pow(beta / beta_trans, 2);
        trat = trat_d * b2 / (1. + b2) + trat_j / (1. + b2);
        two_temp_gam =
            4. /
            3.; // 0.5 * ((1. + 2. / 3. * (trat + 1.) / (trat + 2.)) + gam);
        Th_unit = (two_temp_gam - 1.) * (MP / ME) / (1. + trat);
        Thetae = (p[UU][i][j][k] / p[KRHO][i][j][k]) * Th_unit;
        if (bsq / p[KRHO][i][j][k] > 1) {
          Thetae = 0;
        }

        // if( Be > 1.){// (Be > 1.0) && {((j<=20 || j>=70)||i<=30)){
        //    Thetae=Theta_e_jet;
        // }
        bias_norm += g * dV * Thetae * Thetae;
        V += g * dV;
        if (isnan(bias_norm) || isnan(V)) {
          fprintf(stderr, "%d %d %d %e %e %e %e %e %e\n", i, j, k, g, dV,
                  Thetae, beta, bsq, V);
          exit(1);
        }

        if (i <= 20)
          dMact += g * p[KRHO][i][j][k] * Ucon[1];
        if (i >= 20 && i < 40)
          Ladv += g * p[UU][i][j][k] * Ucon[1] * Ucov[0];
      }
    }
  }

  bias_norm /= V;
  dMact *= dx[3] * dx[2];
  dMact /= 21.;
  Ladv *= dx[3] * dx[2];
  Ladv /= 21.;
  if (world_rank == 0) {
    fprintf(stderr, "\nGLOBALS\n");
    fprintf(stderr, "dMact: %g, Ladv: %g\n", dMact, Ladv);
  }
}
