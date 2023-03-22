

/*
   HARM model specWification routines
 */
#include "decs.h"
#include "BHAC_model.h"

struct of_spectrum ***spect;
extern struct of_spectrum ***shared_spect;
double ***Xi_spec;
extern double ***shared_Xi_spec;
double ***ispec;
extern double ***shared_ispec;

int LFAC, XI;

int N1, N2, N3;


#pragma omp threadprivate(spect, Xi_spec, ispec)

void init_model(char *args[]) {
    set_units(args[3], args[4]);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0)
        fprintf(stderr, "READING DATA\n");

    init_bhac_amr_data(args[2]);

#pragma omp parallel
    {
        int i, j, k;
        spect = (struct of_spectrum ***)malloc(N_THBINS *
                                               sizeof(struct of_spectrum **));
        for (i = 0; i < N_THBINS; i++) {
            spect[i] = (struct of_spectrum **)malloc(
                N_PHIBINS * sizeof(struct of_spectrum *));
            for (j = 0; j < N_PHIBINS; j++) {
                spect[i][j] = (struct of_spectrum *)malloc(
                    N_EBINS * sizeof(struct of_spectrum));
            }
        }
        /*
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
        */
    }
    //	get_profiles();
    Rh = 1 + sqrt(1. - a * a);

    init_geometry();
    if (world_rank == 0)
        fprintf(stderr, "init geom done\n");
    /* make look-up table for hot cross sections */
    init_hotcross();
    if (world_rank == 0)
        fprintf(stderr, "hot cross done\n");

    /* make table for solid angle integrated emissivity and K2 */
    init_emiss_tables();
    if (world_rank == 0)
        fprintf(stderr, "emiss tables done\n");

    /* make table for superphoton weights */
    init_weight_table();

    if (world_rank == 0)
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

double bias_func(double Te, double w) {

    double bias;

    bias = 100. * Te * Te / bias_norm;
#pragma omp atomic
    bias /= max_tau_scatt;
      if (bias < 1)
        bias = 1.;
    //  if (bias > 3000.0)
    //    bias = 3000.0;
    return bias;
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

int find_igrid(double x[4], struct block *block_info, int igrid_cur) {
    int igrid = 0;
    double small = 1e-9;

    for (int igrid = 0; igrid < nleafs; igrid++) {

        if (x[1] + small >= block_info[igrid].lb[0] &&
            x[1] + small <
                block_info[igrid].lb[0] + (block_info[igrid].size[0]) *
                                              block_info[igrid].dxc_block[0] &&
            x[2] + small >= block_info[igrid].lb[1] &&
            x[2] + small <
                block_info[igrid].lb[1] + (block_info[igrid].size[1]) *
                                              block_info[igrid].dxc_block[1] &&
            x[3] + small >= block_info[igrid].lb[2] &&
            x[3] + small <
                block_info[igrid].lb[2] + (block_info[igrid].size[2]) *
                                              block_info[igrid].dxc_block[2]) {
            return igrid;
        }
    }
    return -1;
}

int find_cell(double x[4], struct block *block_info, int igrid) {

    int i = (int)((x[1] - block_info[igrid].lb[0]) /
                  block_info[igrid].dxc_block[0]);
    int j = (int)((x[2] - block_info[igrid].lb[1]) /
                  block_info[igrid].dxc_block[1]);
    int k = (int)((x[3] - block_info[igrid].lb[2]) /
                  block_info[igrid].dxc_block[2]);

    int cell = i + j * block_info[igrid].size[0] +
               k * block_info[igrid].size[0] * block_info[igrid].size[1];

    return cell;
}

// get zonal variables without any interpolations, this should not be changed
void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B,
                    double *sigma, double *beta, double Ucon[NDIM],
                    double Bcon[NDIM], int *ACCZONE, double *dx_local,
                    int *igrid_c) {
    int igrid = *igrid_c;
    int c;
    double del[NDIM];
    double rho, uu;
    double Bp[NDIM], gVcon[NDIM], Vcon[NDIM], Vfac, gVdotgV, UdotBp;
    double gcon[NDIM][NDIM], gcov[NDIM][NDIM];
    double Ucov[NDIM], Bcov[NDIM];

    // *IN_VOLUME = 1;

    double smalll = 1.e-50;
    double small = 0;

    igrid = i;
    c = j;
    k = 0;

    double Xcent[4], X[4];
    /*
      Xcent[0] = 0.;
      Xcent[1] = Xgrid[igrid][c][0]; // + block_info[igrid].dxc_block[0] * 0.5;
      Xcent[2] = Xgrid[igrid][c][1]; // + block_info[igrid].dxc_block[1] * 0.5;
      Xcent[3] = Xgrid[igrid][c][2]; // + block_info[igrid].dxc_block[2] * 0.5;
                                     // //cell centered coordinate
    */
    // coord(i, j, k, Xcent);
    // fprintf(stderr,"Xgrid %e %e %e\n",Xcent[0],Xcent[1],Xcent[2]);
    calc_coord(c, nx, ndimini, block_info[igrid].lb,
               block_info[igrid].dxc_block, Xcent);
    // fprintf(stderr,"Xcent %e %e %e\n",Xcent[0],Xcent[1],Xcent[2]);

    *dx_local = block_info[igrid].dxc_block[0];

   double r  = get_r(Xcent);

    if (r < 1 + sqrt(1 - a * a) || r > 40. ||
        r < 2) { 
        *Ne = 0;
        *B = 0;
        *Thetae = 0;
        *dx_local = 1e100;
        return 1;
    }
    //  gcon_func(Xcent,gcon); // cell centered, nearest neighbour so need
    //  metric at cell position gcov_func(Xcent, gcov);

    // inteprolatie van je primitieve variabelen
    rho = p[KRHO][igrid][c][0];
    uu = p[UU][igrid][c][0];
    // bepalen van de plasma number density en electron temperatuur
    *Ne = rho * Ne_unit + smalll;
    *Thetae = uu / rho * Thetae_unit;

    Bp[1] = p[B1][igrid][c][0]; // interp_scalar_3d(p[B1], i, j,k,del);
    Bp[2] = p[B2][igrid][c][0]; // interp_scalar_3d(p[B2], i, j,k, del);
    Bp[3] = p[B3][igrid][c][0]; // interp_scalar_3d(p[B3], i, j,k, del);

    gVcon[1] = p[U1][igrid][c][0];
    gVcon[2] = p[U2][igrid][c][0];
    gVcon[3] = p[U3][igrid][c][0];

    double gamma_dd[4][4];
    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
            gamma_dd[i][j] =
                geom[igrid][c]
                    .gcov[i][j]; // + gcon[0][i]*gcon[0][j]/(-gcon[0][0]);
        }
    }
    double shift[4];
    for (int j = 1; j < 4; j++) {
        shift[j] = geom[igrid][c].gcon[0][j] / (-geom[igrid][c].gcon[0][0]);
    }
    double alpha = 1 / sqrt(-geom[igrid][c].gcon[0][0]);
    gVdotgV = 0.;
    Ucon[0] = 0.;
    for (int i = 1; i < NDIM; i++) {
        for (int j = 1; j < NDIM; j++) {
            gVdotgV += gamma_dd[i][j] * gVcon[i] * gVcon[j];
        }
    }

    double lfac = sqrt(gVdotgV + 1.);

    Vcon[1] = gVcon[1] / lfac;
    Vcon[2] = gVcon[2] / lfac;
    Vcon[3] = gVcon[3] / lfac;

    Ucon[0] = lfac / alpha;

    for (int i = 1; i < NDIM; i++) {
        Ucon[i] = 0;
        Ucon[i] = Vcon[i] * lfac - shift[i] * lfac / alpha;
    }

    lower(Ucon, geom[igrid][c].gcov, Ucov);
    /* double sum=0;
     for(int i; i < NDIM; i++){
           sum+=Ucon[i]*Ucov[i];
     }
     fprintf(stderr,"sum should be -1: %e\n",sum);
   */
    Bcon[0] = 0;
    for (i = 1; i < NDIM; i++) {
        for (int l = 1; l < NDIM; l++) {
            Bcon[0] += lfac * Bp[i] * (gamma_dd[i][l] * Vcon[l]) / alpha;
        }
    }

    for (int i = 1; i < NDIM; i++) {
        Bcon[i] = 0;
        Bcon[i] = (Bp[i] + alpha * Bcon[0] * Ucon[i]) / lfac;
    }
    lower(Bcon, geom[igrid][c].gcov, Bcov);

    // sterkte van het magneetveld
    *B = sqrt(fabs(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] + Bcon[2] * Bcov[2] +
                   Bcon[3] * Bcov[3])) *
             B_unit +
         smalll;
    double Bsq = fabs(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] +
                      Bcon[2] * Bcov[2] + Bcon[3] * Bcov[3]);
#if (DEBUG)
    if (isnan(Bsq) || isnan(*B)) {
        double R2 =
            Xcent[0] * Xcent[0] + Xcent[1] * Xcent[1] + Xcent[2] * Xcent[2];
        double a2 = a * a;
        double r2 =
            (R2 - a2 +
             sqrt((R2 - a2) * (R2 - a2) + 4. * a2 * Xcent[2] * Xcent[2])) *
            0.5;
        fprintf(stderr, "B isnan r %e rmin %e\n", r, cutoff_inner);
        fprintf(stderr, "B isnan X %e %e %e %e\n", X[0], X[1], X[2], X[3]);
        fprintf(stderr, "B isnan Xgrid %e %e %e\n", Xcent[0], Xcent[1],
                Xcent[2]);
        fprintf(stderr, "B isnan Bsq %e Bcon %e %e %e %e Ucon %e %e %e %e\n",
                Bsq, Bcon[0], Bcon[1], Bcon[2], Bcon[3], Ucon[0], Ucon[1],
                Ucon[2], Ucon[3]);
        fprintf(stderr, "B isnan Bp %e %e %e\n", Bp[1], Bp[2], Bp[3]);
        fprintf(stderr, "B isnan gVcon %e %e %e\n", gVcon[1], gVcon[2],
                gVcon[3]);
        fprintf(stderr, "B isnan Vcon %e %e %e\n", Vcon[1], Vcon[2], Vcon[3]);
        fprintf(stderr, "B isnan VdotV %e\n", VdotV);
        fprintf(stderr, "B isnan lapse %e\n", sqrt(-geom[igrid][c].gcon[0][0]));
        fprintf(stderr, "B isnan shift %e %e %e\n", geom[igrid][c].gcon[0][1],
                geom[igrid][c].gcon[0][2], geom[igrid][c].gcon[0][3]);
        exit(1);
    }
#endif

    double gam = neqpar[0];
    double beta_trans = 1.0;
    *beta = uu * (gam - 1.) / (0.5 * (Bsq + smalll) * beta_trans);

    double b2 = pow(uu * (gam - 1.) / (0.5 * (Bsq + smalll) * beta_trans), 2.);

    *sigma = Bsq / (rho); // *(1.+ uu/rho*gam));

    double beta_max = 1 / (4. * (*sigma));
    double trat = 3.; // trat_d * b2 / (1. + b2) + trat_j / (1. + b2);

    double two_temp_gam =
        0.5 * ((1. + 2. / 3. * (trat + 1.) / (trat + 2.)) + gam);

    Thetae_unit = (gam - 1.) * ((MP / ME)) / (trat + 1);

    *Thetae = (uu / rho) * Thetae_unit;

    if (*sigma > 1.0 || r < 2.0 ||
        r > 40) { 
        *Ne = 0;
        *B = 0;
        *Thetae = 0;
        *dx_local = 1e100;
    }

    return 1;
}

int get_fluid_params(double X[NDIM], double gcov[NDIM][NDIM], double *Ne,
                     double *Thetae, double *B, double *sigma, double *beta,
                     double Ucon[NDIM], double Ucov[NDIM], double Bcon[NDIM],
                     double Bcov[NDIM], int *ACCZONE, double *dx_local,
                     int *igrid_c) {
    int igrid = *igrid_c;
    int i, j, k, c;
    double del[NDIM];
    double rho, uu;
    double Bp[NDIM], gVcon[NDIM], Vcon[NDIM], Vfac, gVdotgV, UdotBp;
    double gcon[NDIM][NDIM];
    //  *IN_VOLUME = 1;

    double smalll = 1.e-50;
    double small = 0;

    if (X[1] > stopx[1] || X[1] < startx[1] || X[2] < startx[2] ||
        X[2] > stopx[2] || X[3] < startx[3] || X[3] > stopx[3]) {
        *Ne = 0;
        *dx_local = 1e100;

        return 0;
    }

    if (igrid == -1 || X[1] + small < block_info[igrid].lb[0] ||
        X[1] + small >
            block_info[igrid].lb[0] +
                (block_info[igrid].size[0]) * block_info[igrid].dxc_block[0] ||
        X[2] + small < block_info[igrid].lb[1] ||
        X[2] + small >
            block_info[igrid].lb[1] +
                (block_info[igrid].size[1]) * block_info[igrid].dxc_block[1] ||
        X[3] + small < block_info[igrid].lb[2] ||
        X[3] + small >
            block_info[igrid].lb[2] +
                (block_info[igrid].size[2]) * block_info[igrid].dxc_block[2]) {
        *igrid_c = find_igrid(X, block_info, igrid);
        igrid = *igrid_c;
    }
    if (igrid == -1) {
        fprintf(stderr, "issues with finding igrid, exiting... %e %e %e\n",
                X[1], X[2], X[3]);
        return 0;
    }
    *dx_local = block_info[igrid].dxc_block[0];

    c = find_cell(X, block_info, igrid);

    double Xcent[4];

    calc_coord(c, nx, ndimini, block_info[igrid].lb,
               block_info[igrid].dxc_block, Xcent);

    double r = get_r(X);

    if (r < Rh || r > 75) {
        *Ne = 0;
        *B = 0;
        *Thetae = 0;
        *dx_local = 1e100;
        return 1;
    }
    gcon_func(Xcent, gcon);

    gcov_func(Xcent, gcov);

    rho = p[KRHO][igrid][c][0];
    uu = p[UU][igrid][c][0];

    *Ne = rho * Ne_unit + smalll;
    *Thetae = uu / rho * Thetae_unit;

    Bp[1] = p[B1][igrid][c][0];
    Bp[2] = p[B2][igrid][c][0];
    Bp[3] = p[B3][igrid][c][0];

    gVcon[1] = p[U1][igrid][c][0];
    gVcon[2] = p[U2][igrid][c][0];
    gVcon[3] = p[U3][igrid][c][0];

    double gamma_dd[4][4];
    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
            gamma_dd[i][j] = gcov[i][j];
        }
    }
    double shift[4];
    for (int j = 1; j < 4; j++) {
        shift[j] = gcon[0][j] / (-gcon[0][0]);
    }
    double alpha = 1 / sqrt(-gcon[0][0]);
    gVdotgV = 0.;
    Ucon[0] = 0.;
    for (int i = 1; i < NDIM; i++) {
        for (int j = 1; j < NDIM; j++) {
            gVdotgV += gamma_dd[i][j] * Vcon[i] * Vcon[j];
        }
    }


    double lfac = sqrt(gVdotgV + 1.);

    Vcon[1] = gVcon[1] / lfac;
    Vcon[2] = gVcon[2] / lfac;
    Vcon[3] = gVcon[3] / lfac;

    Ucon[0] = lfac / alpha;

    for (int i = 1; i < NDIM; i++) {
        Ucon[i] = 0;
        Ucon[i] = Vcon[i] * lfac - shift[i] * lfac / alpha;
    }

    lower(Ucon, gcov, Ucov);

    Bcon[0] = 0;
    for (i = 1; i < NDIM; i++) {
        for (int l = 1; l < NDIM; l++) {
            Bcon[0] += lfac * Bp[i] * (gamma_dd[i][l] * Vcon[l]) / alpha;
        }
    }

    for (int i = 1; i < NDIM; i++) {
        Bcon[i] = 0;
        Bcon[i] = (Bp[i] + alpha * Bcon[0] * Ucon[i]) / lfac;
    }

    lower(Bcon, gcov, Bcov);

    *B = sqrt(fabs(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] + Bcon[2] * Bcov[2] +
                   Bcon[3] * Bcov[3])) *
             B_unit +
         smalll;
    double Bsq = fabs(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] +
                      Bcon[2] * Bcov[2] + Bcon[3] * Bcov[3]);
#if (DEBUG)
    if (isnan(Bsq) || isnan(*B)) {
        r = get_r(X);
        fprintf(stderr, "B isnan r %e rmin %e\n", r, cutoff_inner);
        fprintf(stderr, "B isnan X %e %e %e %e\n", X[0], X[1], X[2], X[3]);
        fprintf(stderr, "B isnan Xgrid %e %e %e\n", Xcent[0], Xcent[1],
                Xcent[2]);
        fprintf(stderr, "B isnan Bsq %e Bcon %e %e %e %e Ucon %e %e %e %e\n",
                Bsq, Bcon[0], Bcon[1], Bcon[2], Bcon[3], Ucon[0], Ucon[1],
                Ucon[2], Ucon[3]);
        fprintf(stderr, "B isnan Bp %e %e %e\n", Bp[1], Bp[2], Bp[3]);
        fprintf(stderr, "B isnan gVcon %e %e %e\n", gVcon[1], gVcon[2],
                gVcon[3]);
        fprintf(stderr, "B isnan Vcon %e %e %e\n", Vcon[1], Vcon[2], Vcon[3]);
        fprintf(stderr, "B isnan VdotV %e\n", VdotV);
        fprintf(stderr, "B isnan lapse %e\n", sqrt(-gcon[0][0]));
        fprintf(stderr, "B isnan shift %e %e %e\n", gcon[0][1], gcon[0][2],
                gcon[0][3]);
        exit(1);
    }
#endif

    double gam = neqpar[0];

    double beta_trans = 1.0;
    *beta = uu * (gam - 1.) / (0.5 * (Bsq + smalll) * beta_trans);

    double b2 = pow(uu * (gam - 1.) / (0.5 * (Bsq + smalll) * beta_trans), 2.);

    *sigma = Bsq / (rho); 

    double beta_max = 1 / (4 * (*sigma));
    double trat = trat_d * b2 / (1. + b2) + trat_j / (1. + b2);
        
    double two_temp_gam =
        0.5 * ((1. + 2. / 3. * (trat + 1.) / (trat + 2.)) + gam);

    Thetae_unit = (gam - 1.) * ((MP / ME)) / (trat + 1);

    *Thetae = (uu / rho) * Thetae_unit;

    if (*sigma > 1.0 || r < 2.0 ||
        r > 40) { 
        *Ne = 0;
        *Thetae = 0;
        *B = 0;
        *dx_local = 1e100;
    }
    return 1;
}

double get_r(double X[4]) {
    double r;

#if (CKS)
    double R2 = X[1] * X[1] + X[2] * X[2] + X[3] * X[3];
    double a2 = a * a;
    double r2 =
        0.5 * (R2 - a2 + sqrt((R2 - a2) * (R2 - a2) + 4. * a2 * X[3] * X[3]));
    r = sqrt(r2);
#elif (MKS)
    r = exp(X[1]);
#endif

    return r;
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

#if (MKS)
    for(int i=0;i<4;i++)
       	for(int j=0;j<4;j++)
               	gcon[i][j]=0;

    double r = exp(X[1]);
    double theta = X[2] + 0.5 * hslope * sin(2. * X[2]);

    double sinth = sin(theta);
    double sin2th = sinth * sinth;
    double costh = cos(theta);

    double irho2 = 1. / (r * r + a * a * costh * costh);

    double hfac = 1. + hslope * cos(2. * X[2]);

    gcon[0][0] = -1. - 2. * r * irho2;
    gcon[0][1] = 2. * irho2;

    gcon[1][0] = gcon[0][1];
    gcon[1][1] = irho2 * (r * (r - 2.) + a * a) / (r * r);
    gcon[1][3] = a * irho2 / r;

    gcon[2][2] = irho2 / (hfac * hfac);

    gcon[3][1] = gcon[1][3];
    gcon[3][3] = irho2 / (sin2th);

#elif (CKS)
    int i, j, k;
    double R2 = X[1] * X[1] + X[2] * X[2] + X[3] * X[3];
    double a2 = a * a;
    double r2 =
        (R2 - a2 + sqrt((R2 - a2) * (R2 - a2) + 4. * a2 * X[3] * X[3])) * 0.5;
    double r = sqrt(r2);

    double sig = (r2 * r2 + a2 * X[3] * X[3]);
    double del = (r2 + a2);
    double isig, idel;

    if (r < 1)
        isig = 1 / (sig + 1e-3);
    else
        isig = 1 / sig;
    idel = 1 / (del);

    double f = 2. * r2 * r * isig; 

    double l[4];

    l[0] = -1.;
    l[1] = (r * X[1] + a * X[2]) * idel; 
    l[2] = (r * X[2] - a * X[1]) * idel; 
    l[3] = (X[3]) / (r);

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            if (i == j && i == 0)
                gcon[i][j] = -1. - f * l[i] * l[j];
            else if (i == j && i != 0)
                gcon[i][j] = 1. - f * l[i] * l[j];
            else
                gcon[i][j] = -f * l[i] * l[j];
        }
    }
#endif
}

// dd
void gcov_func(double *X, double gcov[][NDIM]) {

#if (MKS)
    for(int i=0;i<4;i++)
	for(int j=0;j<4;j++)
		gcov[i][j]=0;

    double r = exp(X[1]);
    double theta = X[2] + 0.5 * hslope * sin(2. * X[2]);

    double sinth = sin(theta);
    double sin2th = sinth * sinth;
    double costh = cos(theta);
    double tfac, rfac, hfac, pfac;
    double rho2 = r * r + a * a * costh * costh;

    tfac = 1.;
    rfac = r;
    hfac = 1. + hslope * cos(2. * X[2]);
    pfac = 1.;

    gcov[0][0] = (-1. + 2. * r / rho2) * tfac * tfac;
    gcov[0][1] = (2. * r / rho2) * tfac * rfac;
    gcov[0][3] = (-2. * a * r * sin2th / rho2) * tfac * pfac;

    gcov[1][0] = gcov[0][1];
    gcov[1][1] = (1. + 2. * r / rho2) * rfac * rfac;
    gcov[1][3] = (-a * sin2th * (1. + 2. * r / rho2)) * rfac * pfac;

    gcov[2][2] = rho2 * hfac * hfac;

    gcov[3][0] = gcov[0][3];
    gcov[3][1] = gcov[1][3];
    gcov[3][3] =
        sin2th * (rho2 + a * a * sin2th * (1. + 2. * r / rho2)) * pfac * pfac;


#elif (CKS)
    int i, j, k;
    double R2 = X[1] * X[1] + X[2] * X[2] + X[3] * X[3];
    double a2 = a * a;
    double r2 =
        (R2 - a2 + sqrt((R2 - a2) * (R2 - a2) + 4. * a2 * X[3] * X[3])) * 0.5;
    double r = sqrt(r2);

    double sig = (r2 * r2 + a2 * X[3] * X[3]);
    double del = (r2 + a2);
    double isig, idel;

    if (r < 1)
        isig = 1 / (sig + 1e-3);
    else
        isig = 1 / sig;
    idel = 1 / (del);

    double f = 2. * r2 * r * isig;     double l[4];

    l[0] = 1;
    l[1] = (r * X[1] + a * X[2]) * idel; 
    l[2] = (r * X[2] - a * X[1]) * idel; 
    l[3] = (X[3]) / (r);

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            if (i == j && i == 0)
                gcov[i][j] = -1. + f * l[i] * l[j];
            else if (i == j && i != 0)
                gcov[i][j] = 1. + f * l[i] * l[j];
            else
                gcov[i][j] = f * l[i] * l[j];
        }
    }
#endif
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
    int i, j, k;
    double temp[4][4][4];
    double dg[4][4][4];

    LOOP_ijk { lconn[i][j][k] = 0; }

    // Obtain metric at current position (contravariant form)
    double g_uu[4][4], g_dd[4][4], g_dd_m[4][4], g_dd_p[4][4];
    double X_u_temp[4];
    gcon_func(X, g_uu);
    gcov_func(X, g_dd);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            X_u_temp[j] = X[j];
        }

        X_u_temp[i] += delta_num;

        gcov_func(X_u_temp, g_dd_p);

        X_u_temp[i] -= 2. * delta_num;
        gcov_func(X_u_temp, g_dd_m);

        dg[i][0][0] = (g_dd_p[0][0] - g_dd_m[0][0]) / (2. * delta_num);
        dg[i][1][0] = (g_dd_p[1][0] - g_dd_m[1][0]) / (2. * delta_num);
        dg[i][2][0] = (g_dd_p[2][0] - g_dd_m[2][0]) / (2. * delta_num);
        dg[i][3][0] = (g_dd_p[3][0] - g_dd_m[3][0]) / (2. * delta_num);

        dg[i][0][1] = (g_dd_p[0][1] - g_dd_m[0][1]) / (2. * delta_num);
        dg[i][1][1] = (g_dd_p[1][1] - g_dd_m[1][1]) / (2. * delta_num);
        dg[i][2][1] = (g_dd_p[2][1] - g_dd_m[2][1]) / (2. * delta_num);
        dg[i][3][1] = (g_dd_p[3][1] - g_dd_m[3][1]) / (2. * delta_num);

        dg[i][0][2] = (g_dd_p[0][2] - g_dd_m[0][2]) / (2. * delta_num);
        dg[i][1][2] = (g_dd_p[1][2] - g_dd_m[1][2]) / (2. * delta_num);
        dg[i][2][2] = (g_dd_p[2][2] - g_dd_m[2][2]) / (2. * delta_num);
        dg[i][3][2] = (g_dd_p[3][2] - g_dd_m[3][2]) / (2. * delta_num);

        dg[i][0][3] = (g_dd_p[0][3] - g_dd_m[0][3]) / (2. * delta_num);
        dg[i][1][3] = (g_dd_p[1][3] - g_dd_m[1][3]) / (2. * delta_num);
        dg[i][2][3] = (g_dd_p[2][3] - g_dd_m[2][3]) / (2. * delta_num);
        dg[i][3][3] = (g_dd_p[3][3] - g_dd_m[3][3]) / (2. * delta_num);
    }
    // Solve the Christoffel connection equation
    int alpha;

    for (alpha = 0; alpha < 4; alpha++) {
        for (k = 0; k < 4; k++) {
            lconn[alpha][0][0] +=
                g_uu[alpha][k] * (0.5 * (2. * dg[0][k][0] - dg[k][0][0]));
            lconn[alpha][0][1] +=
                g_uu[alpha][k] *
                (0.5 * (dg[1][k][0] + dg[0][k][1] - dg[k][0][1]));
            lconn[alpha][0][2] +=
                g_uu[alpha][k] *
                (0.5 * (dg[2][k][0] + dg[0][k][2] - dg[k][0][2]));
            lconn[alpha][0][3] +=
                g_uu[alpha][k] *
                (0.5 * (dg[3][k][0] + dg[0][k][3] - dg[k][0][3]));

            lconn[alpha][1][1] +=
                g_uu[alpha][k] * (0.5 * (2. * dg[1][k][1] - dg[k][1][1]));
            lconn[alpha][1][2] +=
                g_uu[alpha][k] *
                (0.5 * (dg[2][k][1] + dg[1][k][2] - dg[k][1][2]));
            lconn[alpha][1][3] +=
                g_uu[alpha][k] *
                (0.5 * (dg[3][k][1] + dg[1][k][3] - dg[k][1][3]));

            lconn[alpha][2][2] +=
                g_uu[alpha][k] * (0.5 * (2. * dg[2][k][2] - dg[k][2][2]));
            lconn[alpha][2][3] +=
                g_uu[alpha][k] *
                (0.5 * (dg[3][k][2] + dg[2][k][3] - dg[k][2][3]));

            lconn[alpha][3][3] +=
                g_uu[alpha][k] * (0.5 * (2. * dg[3][k][3] - dg[k][3][3]));
        }
        lconn[alpha][1][0] = lconn[alpha][0][1];

        lconn[alpha][2][0] = lconn[alpha][0][2];
        lconn[alpha][2][1] = lconn[alpha][1][2];

        lconn[alpha][3][0] = lconn[alpha][0][3];
        lconn[alpha][3][1] = lconn[alpha][1][3];
        lconn[alpha][3][2] = lconn[alpha][2][3];
    }
}

/* stopping criterion for geodesic integrator */
/* K not referenced intentionally */

#define RMAX 1e4
#define ROULETTE 1.e3
int stop_criterion(struct of_photon *ph) {
    double wmin, X1min, X1max;

    double r = get_r(ph->X);

    wmin = WEIGHT_MIN; /* stop if weight is below minimum weight */

    X1min = (Rh);         /* this is coordinate-specific; stop
                                at event horizon */
    X1max = (RMAX * 1.1); /* this is coordinate and simulation
                                specific: stop at large distance */

    if (r < X1min) {
        return 1;
    }

    if (r > X1max) {
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
        return 1;
    }

    return (0);
}

/* criterion for recording photon */
int record_criterion(struct of_photon *ph) {
    //  const double X1max = log(RMAX);
    const double X1max = RMAX * 1.1;

    double r = get_r(ph->X);

    if (r > X1max) {
        return (1);
    } else
        return (0);
}

/* EPS doublely ought to be related to the number of
   zones in the simulation. */
#define EPS 0.01
// #define EPS   0.01

double stepsize(double X_u[4], double U_u[4], double dx_local) {
    double r = get_r(X_u);
    double rfac;
#if (CKS)
    rfac = r * r;
#elif (MKS)
    rfac = 1;
#endif

    double dlx1 = EPS / (fabs(U_u[1]) + SMALL * SMALL);
    double dlx2 = EPS / (fabs(U_u[2]) + SMALL * SMALL);
    double dlx3 = EPS / (fabs(U_u[3]) + SMALL * SMALL);

    double idlx1 = 1. / (fabs(dlx1) + SMALL * SMALL);
    double idlx2 = 1. / (fabs(dlx2) + SMALL * SMALL);
    double idlx3 = 1. / (fabs(dlx3) + SMALL * SMALL);

    double maxU = fmax(fmax(fabs(U_u[1]), fabs(U_u[2])), fabs(U_u[3]));

    double step;

    if (r < 40) {
        double step_grid = dx_local / 2. / maxU;

        step = fmin((rfac) / (idlx1 + idlx2 + idlx3), step_grid); //
    }

    else
        step = (rfac) / (idlx1 + idlx2 + idlx3);

    return step;
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
    double r = get_r(ph->X);
    double theta = (acos(ph->X[3] / r));
    double phi = fmod(atan2(ph->X[2], ph->X[1]), 2 * M_PI);
    if (phi < 0)
        phi += 2 * M_PI;
    if (theta < 0 || theta > M_PI)
        fprintf(stderr, "issue with theta %e\n", theta);

    /***************************************************/
    /*  get theta bin, without folding around equator   */
    dx2 = (M_PI) / ((double)N_THBINS);
    ix2 = (int)(theta / dx2);
    dx3 = (2 * M_PI) / ((double)N_PHIBINS);
    ix3 = (int)(phi / dx3);
    /***************************************************/

    /* check limits */
    if (ix2 < 0 || ix2 >= N_THBINS) {
        fprintf(stderr, "outside theta bin %e %d\n", theta, ix2);
        return;
    }
    if (ix3 < 0 || ix3 >= N_PHIBINS) {
        fprintf(stderr, "outside phi bin %e %d\n", phi, ix3);
        return;
    }

    double Xi[4], del[4], gcov[NDIM][NDIM], g;
    Xi[0] = 0;
    Xi[1] = ph->X1i;
    Xi[2] = ph->X2i;
    Xi[3] = ph->X3i;

    gcov_func(Xi, gcov);
    g = gdet_func(gcov, Xi);

    double dix1, dix2, dix3;

    /* get energy bin */
    lE = log(ph->E);
    iE =
        (int)((lE - lE0) / dlE + 2.5) - 2; /* bin is centered on iE*dlE + lE0 */
    /* check limits */
    if (iE < 0 || iE >= N_EBINS) {
        return;
    }

#pragma omp atomic
    N_superph_recorded++;
#pragma omp atomic
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
                MPI_Reduce(&(spect[i][j][k]), &(shared_spect[i][j][k]), 1,
                           tstype, sumstruct, 0, MPI_COMM_WORLD);
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

    sprintf(filename, "spec_kappa_%d.dat", index);
    fp = fopen(filename, "w");
    if (fp == NULL) {
        fprintf(stderr, "trouble opening spectrum file\n");
        exit(0);
    }
#if MPI
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif

    /* output */
    max_tau_scatt = 0.;
    L = 0.;
    for (i = 0; i < N_EBINS; i++) {

        /* output log_10(photon energy/(me c^2)) */
        fprintf(fp, "%10.5g ", exp(i * dlE + lE0) * ME * CL * CL / HPL);

        for (k = 0; k < N_PHIBINS; k++) {
            for (j = 0; j < N_THBINS; j++) {
                /* convert accumulated photon number in each bin
                   to \nu L_\nu, in units of Lsun */
                // to fold the SED
                // dx2 = (stopx[2] - startx[2]) / (2. * N_THBINS);

                dx2 = (1.) / (1. * N_THBINS);
                dx3 = (2 * M_PI) / (1. * N_PHIBINS);

                /* factor of 2 accounts for folding around equator */
                // to fold the sed
                // dOmega = 2. * dOmega_func(j * dx2, (j + 1) * dx2);

                dOmega =
                    dOmega_func(j * dx2, (j + 1) * dx2, k * dx3, (k + 1) * dx3);
                d = 2.6228263e22;
                // 5.061e25; M87
                //                                 d= 2.6228263e22;
                //                                 //SgrA*//RMAX*L_unit;

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
                    (shared_spect[j][k][i].dNdlE / ((double)world_size) +
                     SMALL);
#endif

#if OPENMP
                tau_scatt = shared_spect[j][k][i].tau_scatt /
                            (shared_spect[j][k][i].dNdlE + SMALL);
#endif
                Inu *= nuLnu / (exp(i * dlE + lE0) * ME * CL * CL / HPL);

                fprintf(fp, "%10.5g %10.5g ", nuLnu, Inu);
                if (tau_scatt > max_tau_scatt)
                    max_tau_scatt = tau_scatt;
                L += nuLnu * dOmega * dlE;
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(
        stderr,
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
                // back to the cureent set up, th_len=pi, adn hslope=0 in our
                // case
                th = M_PI * X[2]; // th_len*X[2] + th_beg + hslope*sin(2. * M_PI
                                  // * X[2]);
                phi = fmod(X[3], stopx[3]);
                if (phi < 0.)
                    phi = stopx[3] + phi;

                if (k == 0 || k == N3 - 1) {
                    fprintf(fp, "%g %g %g\n", r * sin(th) * cos(phi),
                            r * cos(th), 0.0);
                    Xi_spec[i][j][k] =
                        (Xi_spec[i][j][N3 - 1] + Xi_spec[i][j][0]) / 2.;
                } else
                    fprintf(fp, "%g %g %g\n", r * sin(th) * cos(phi),
                            r * cos(th), r * sin(th) * sin(phi));
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
