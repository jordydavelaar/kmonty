
#include "decs.h"
#include "BHAC_model.h"

extern double ****p;
extern struct of_geom **geom;

extern struct of_spectrum ***shared_spect;

void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]);
void coord(int i, int j, int k, double *X);
void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B,
                    double *sigma, double *beta, double Ucon[NDIM],
                    double Bcon[NDIM], int *ACCZONE, double *dx_local,
                    int *igrid_c);

/********************************************************************

    Interpolation routines

********************************************************************/

double lnu_min, lnu_max, dlnu;

static void init_linear_interp_weight() {

    lnu_min = log(NUMIN);
    lnu_max = log(NUMAX);
    dlnu = (lnu_max - lnu_min) / (N_ESAMP);
}

static double linear_interp_weight(double nu) {
    int i;
    double di, lnu;

    lnu = log(nu);

    di = (lnu - lnu_min) / dlnu;
    i = (int)di;
    di = di - i;

    return exp((1. - di) * wgt[i] + di * wgt[i + 1]); //* exp(-nu/nu_cutoff);
}

/***********************************************************************************

     End interpolation routines

***********************************************************************************/

#define JCST (M_SQRT2 * EE * EE * EE / (27 * ME * CL * CL))
void init_weight_table(void) {

    int i, j, k, l, lstart, lend, myid, nthreads;
    double Ne, Thetae, B, K2;
    double sum[N_ESAMP + 1], nu[N_ESAMP + 1];
    double fac_th, fac_nth, sfac;
    double Ucon[NDIM], Bcon[NDIM];
    double sigma, beta, dx_local;
    int igrid;
    /*      Set up interpolation */
    init_linear_interp_weight();

    for (i = 0; i <= N_ESAMP; i++) {
        sum[i] = 0.;
        sum[i] = 0.;
        nu[i] = exp(i * dlnu + lnu_min);
    }

    {
        nthreads = 1;
        myid = 0;
        lstart = myid * (N_ESAMP / nthreads);
        lend = (myid + 1) * (N_ESAMP / nthreads);
        if (myid == nthreads - 1)
            lend = N_ESAMP + 1;
        int ACCZONE;
        for (i = 0; i < N1; i++) {
            for (j = 0; j < N2; j++)
                for (k = 0; k < N3; k++) {
                    igrid = i;
                   get_fluid_zone(i, j, k, &Ne, &Thetae, &B, &sigma, &beta,
                                   Ucon, Bcon, &ACCZONE, &dx_local, &igrid);
                    sfac = dx_local * dx_local * dx_local * L_unit * L_unit *
                           L_unit;

                    if (Ne == 0. || Thetae < THETAE_MIN) {
                        continue;
                    }
                    // thermal
#if (THERMAL)
                    K2 = K2_eval(Thetae);
                    fac_th = (JCST * Ne * B * Thetae * Thetae / K2) * sfac *
                             geom[i][j].g;
#elif (KAPPA || POWERLAW)
                    fac_nth = (EE * EE * EE / (2. * M_PI * ME * CL * CL)) * Ne *
                              B * sfac * geom[i][j].g;
#elif (MIXED)
                    K2 = K2_eval(Thetae);
                    fac_th = (JCST * Ne * B * Thetae * Thetae / K2) * sfac *
                             geom[i][j].g;
                    fac_nth = (EE * EE * EE / (2. * M_PI * ME * CL * CL)) * Ne *
                              B * sfac * geom[i][j].g;
#endif
                    for (l = lstart; l < lend; l++) {
#if (THERMAL)
                        sum[l] += fac_th * F_eval_th(Thetae, B, nu[l]);
#elif (KAPPA || POWERLAW)
                        sum[l] += fac_nth * F_eval_nth(Thetae, B, nu[l]);
#elif (MIXED)
                        if (ACCZONE) {
                            sum[l] += perct_thermal * fac_th *
                                          F_eval_th(Thetae, B, nu[l]) +
                                      (1 - perct_thermal) * fac_nth *
                                          F_eval_nth(Thetae, B, nu[l]);
                        } else
                            sum[l] += fac_th * F_eval_th(Thetae, B, nu[l]);
#endif
                        if (isnan(sum[l]))
                            fprintf(stderr,
                                    "fac_th %e K2 %e Thetae %e B %e Ne %e "
                                    "nu[l] %e F_eval %e\n",
                                    fac_th, K2, Thetae, B, Ne, nu[l],
                                    F_eval_th(Thetae, B, nu[l]));
                    }
                }
        }
    }
    // #pragma omp parallel for schedule(static) private(i)
    for (i = 0; i <= N_ESAMP; i++) {
        wgt[i] = log(sum[i] / (HPL * Ns) + WEIGHT_MIN);
    }
    fflush(stderr);
    return;
}

#undef JCST
#if (POWERLAW)
#define BTHSQMIN (1.e-5)
#define BTHSQMAX (1.e-2)
#define NINT (20000)
#else
#define BTHSQMIN (1.e-5)
#define BTHSQMAX (1.e-2)
#define NINT (20000)
#endif

double lb_min, dlb;
double nint[NINT + 1];
double dndlnu_max[NINT + 1];
void init_nint_table(void) {
    int i, j;
    double Bmag, dn = 0;
    static int firstc = 1;

    if (firstc) {
        lb_min = log(BTHSQMIN);
        dlb = log(BTHSQMAX / BTHSQMIN) / NINT;
        firstc = 0;
    }

    for (i = 0; i <= NINT; i++) {
        nint[i] = 0.;
        Bmag = exp(i * dlb + lb_min);
        dndlnu_max[i] = 0.;
        for (j = 0; j < N_ESAMP; j++) {
#if (THERMAL)
            dn = F_eval_th(1., Bmag, exp(j * dlnu + lnu_min)) /
                 (exp(wgt[j]) + 1.e-100);
#elif (KAPPA || POWERLAW)
            dn = F_eval_nth(1., Bmag, exp(j * dlnu + lnu_min)) /
                 (exp(wgt[j]) + 1.e-100);

#endif
            if (dn > dndlnu_max[i]) {
                dndlnu_max[i] = dn;
            }

            nint[i] += dlnu * dn;
        }
#if (THERMAL)
        double fac_th = M_SQRT2 * EE * EE * EE / (27. * ME * CL * CL);
        nint[i] *= L_unit * L_unit * L_unit * fac_th * 1. / HPL;
#elif (KAPPA || POWERLAW)
        double fac_nth = (EE * EE * EE / (2. * M_PI * ME * CL * CL));
        nint[i] *= L_unit * L_unit * L_unit * fac_nth * 1. / HPL;
#endif
        nint[i] = log(nint[i]);
        dndlnu_max[i] = log(dndlnu_max[i]);
    }

    return;
}

static void init_zone(int i, int j, int k, double *nz, double *dnmax) {

    int l;
    double Ne, Thetae, Bmag, lbth, fac_th, fac_nth;
    double sigma, beta, dx_local;
    int igrid;
    double dl, dn, ninterp, K2;
    double Ucon[NDIM], Bcon[NDIM];
    int ACCZONE = 0;
    get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, &sigma, &beta, Ucon, Bcon,
                   &ACCZONE, &dx_local, &igrid);

    if (Ne == 0. || Thetae < THETAE_MIN) {
        *nz = 0.;
        *dnmax = 0.;
        return;
    }
#if (POWERLAW)
    lbth = log(Bmag);
#else
    lbth = log(Bmag * Thetae * Thetae);
#endif
    lb_min = log(BTHSQMIN);
    dlb = log(BTHSQMAX / BTHSQMIN) / NINT;

    dl = (lbth - lb_min) / dlb;
    l = (int)dl;
    dl = dl - l;

    if (l < 0) {
        *dnmax = 0.;
        *nz = 0.;
        fprintf(stderr,
                "l smaller than zero l %d lbth %e dlb %e Bmag %e Thetae %e\n",
                l, lbth, dlb, Bmag, Thetae);
        return;
    } else if (l >= NINT || MIXED || POWERLAW) {
        ninterp = 0.;
        *dnmax = 0.;
        for (l = 0; l <= N_ESAMP; l++) {
#if (KAPPA || POWERLAW)
            dn =
                F_eval_nth(Thetae, Bmag, exp(l * dlnu + lnu_min)) / exp(wgt[l]);
#elif (THERMAL)
            dn = F_eval_th(Thetae, Bmag, exp(l * dlnu + lnu_min)) / exp(wgt[l]);
#elif (MIXED)
            K2 = K2_eval(Thetae);
            if (K2 == 0.) {
                *nz = 0.;
                *dnmax = 0.;
                return;
            }

            if (ACCZONE) {
                dn = perct_thermal *
                         F_eval_th(Thetae, Bmag, exp(l * dlnu + lnu_min)) /
                         (exp(wgt[l] + 1e-100)) +
                     (27. / (M_SQRT2 * 2. * M_PI)) * (K2 / (Thetae * Thetae)) *
                         (1 - perct_thermal) *
                         F_eval_nth(Thetae, Bmag, exp(l * dlnu + lnu_min)) /
                         (exp(wgt[l] + 1e-100));
            } else
                dn = F_eval_th(Thetae, Bmag, exp(l * dlnu + lnu_min)) /
                     (exp(wgt[l]));
#endif

            if (dn > *dnmax) {
                *dnmax = dn;
            }
            ninterp += dlnu * dn;
        }
#if (THERMAL)
        fac_th = M_SQRT2 * EE * EE * EE / (27. * ME * CL * CL);
        ninterp *= dx_local * dx_local * dx_local * L_unit * L_unit * L_unit *
                   fac_th * 1. / HPL;
#elif (KAPPA || POWERLAW)
        fac_nth = (EE * EE * EE / (2. * M_PI * ME * CL * CL));
        ninterp *= dx_local * dx_local * dx_local * L_unit * L_unit * L_unit *
                   fac_nth * 1. / HPL;

#elif (MIXED)
        double fac_mixed = M_SQRT2 * EE * EE * EE / (27. * ME * CL * CL);

        ninterp *= dx_local * dx_local * dx_local * L_unit * L_unit * L_unit *
                   fac_mixed * 1. / HPL;

#endif
    } else {
        if (isinf(nint[l]) || isinf(nint[l + 1])) {
            ninterp = 0.;
            *dnmax = 0.;
        } else {
            ninterp = exp((1. - dl) * nint[l] + dl * nint[l + 1]) * dx_local *
                      dx_local * dx_local;
            *dnmax = exp((1. - dl) * dndlnu_max[l] + dl * dndlnu_max[l + 1]);
        }
    }

#if (THERMAL)
    K2 = K2_eval(Thetae);
    if (K2 == 0.) {
        *nz = 0.;
        *dnmax = 0.;
        fprintf(stderr, "K2 is zero\n");
        return;
    }

    fac_th = Ne * Bmag * Thetae * Thetae / K2;
    *nz = geom[i][j].g * fac_th * ninterp;
#elif (KAPPA || POWERLAW)
    fac_nth = Ne * Bmag;
    *nz = geom[i][j].g * fac_nth * ninterp;

#elif (MIXED)
    double fac_mixed = Ne * Bmag * Thetae * Thetae / K2;
    *nz = geom[i][j].g * fac_mixed * ninterp;
#endif

    if (*nz > Ns * log(NUMAX / NUMIN)) {
        fprintf(stderr,
                "Something very wrong in zone %d %d: \nB=%g  Thetae=%g  "
                "K2=%g  ninterp=%g\n\n",
                i, j, Bmag, Thetae, K2, ninterp);
        *nz = 0.;
        *dnmax = 0.;
    }

    return;
}

int zone_flag;
int get_zone(int *i, int *j, int *k, double *dnmax) {
    /* Return the next zone and the number of superphotons that need to be *
     * generated in it.                                */

    int in2gen;
    double n2gen;
    static int zk = 0;

    static int zi = 0;
    static int zj = -1;

    zone_flag = 1;
    zj++;
    if (zj >= N2) {
        zj = 0;
        zi++;
        if (zi >= N1) {
            in2gen = 1;
            *i = N1;
            return 1;
        }
    }
    init_zone(zi, zj, zk, &n2gen, dnmax);

    if (n2gen == 0. || *dnmax == 0.) {
        in2gen = 0;
    } else if (fmod(n2gen, 1.) > monty_rand()) {
        in2gen = (int)n2gen + 1;
    } else {
        in2gen = (int)n2gen;
    }

    *i = zi;
    *j = zj;
    *k = zk;

    // we only need a fraction of these for every processor

    return in2gen;
}

void sample_zone_photon(int i, int j, int k, double dnmax,
                        struct of_photon *ph) {
    /* Set all initial superphoton attributes */

    int l;
    double K_tetrad[NDIM], tmpK[NDIM], E, Nln;
    double nu, th, cth, sth, phi, sphi, cphi, jmax, weight;
    double Ne, Thetae, Bmag, Ucon[NDIM], Bcon[NDIM], bhat[NDIM];
    static double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];
    double sigma, beta, dx_local;
    int igrid = i;

    calc_coord(j, nx, ndimini, block_info[i].lb, block_info[i].dxc_block,
               ph->X);
    Nln = lnu_max - lnu_min;
    int ACCZONE = 0;
    get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, &sigma, &beta, Ucon, Bcon,
                   &ACCZONE, &dx_local, &igrid);
    double dV = pow(dx_local, 3.);

    /* Sample from superphoton distribution in current sfimulation zone */
    do {
        nu = exp(monty_rand() * Nln + lnu_min);
        weight = linear_interp_weight(nu);
    } while (monty_rand() >
             (F_eval(Thetae, Bmag, nu, ACCZONE) / weight) / dnmax);

    ph->w = weight;
    jmax = jnu_synch(nu, Ne, Thetae, Bmag, M_PI / 2., ACCZONE);
    do {
        cth = 2. * monty_rand() - 1.;
        th = acos(cth);
    } while (monty_rand() >
             jnu_synch(nu, Ne, Thetae, Bmag, th, ACCZONE) / jmax);

    sth = sqrt(1. - cth * cth);
    phi = 2. * M_PI * monty_rand();
    cphi = cos(phi);
    sphi = sin(phi);

    E = nu * HPL / (ME * CL * CL);
    K_tetrad[0] = E;

    K_tetrad[0] = E;
    K_tetrad[1] = E * cth;
    K_tetrad[2] = E * cphi * sth;
    K_tetrad[3] = E * sphi * sth;

    if (zone_flag) {
        /* first photon created in this zone, so make the tetrad */
        if (Bmag > 0.) {
            for (l = 0; l < NDIM; l++)
                bhat[l] = Bcon[l] * B_unit / Bmag;
        } else {
            for (l = 1; l < NDIM; l++)
                bhat[l] = 0.;
            bhat[1] = 0.;
        }
        make_tetrad(ph->X, Ucon, bhat, geom[i][j].gcov, Econ, Ecov);
        zone_flag = 0;
    }

    tetrad_to_coordinate(Econ, K_tetrad, ph->K);

    K_tetrad[0] *= -1.;
    tetrad_to_coordinate(Ecov, K_tetrad, tmpK);

    ph->E = ph->E0 = ph->E0s = -tmpK[0];
    ph->L = tmpK[3];
    ph->tau_scatt = 0.;
    ph->tau_abs = 0.;
    ph->X1i = ph->X[1];
    ph->X2i = ph->X[2];
    ph->nscatt = 0;
    ph->ne0 = Ne;
    ph->b0 = Bmag;
    ph->thetae0 = Thetae;

    return;
}

void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]) {
    double phi;
    /* Map X[3] into sim range, assume startx[3] = 0 */
    phi = fmod(X[3], stopx[3]);
    if (phi < 0.)
        phi = stopx[3] + phi;

    *i = (int)((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
    *j = (int)((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
    *k = (int)(phi / dx[3] - 0.5 + 1000) - 1000;

    if (*i < 0) {
        *i = 0;
        del[1] = 0.;
    } else if (*i > N1 - 2) {
        *i = N1 - 2;
        del[1] = 1.;
    } else {
        del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
    }

    if (*j < 0) {
        *j = 0;
        del[2] = 0.;
    } else if (*j > N2 - 2) {
        *j = N2 - 2;
        del[2] = 1.;
    } else {
        del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
    }

    if (*k < 0) {
        *k = 0;
        del[3] = 0.;
    } else if (*k > N3 - 2) {
        *k = N3 - 2;
        del[3] = 1.;
    } else {
        del[3] = (phi - ((*k + 0.5) * dx[3])) / dx[3];
    }

    return;
}

/* return boyer-lindquist coordinate of point */
void bl_coord(double *X, double *r, double *th) {

    // for cfg data files
    // *r = Rin * exp(X[1]);
    // for hotaka data files
    *r = exp(X[1]) + R0;
    *th = X[2] * M_PI + (hslope / 2.) * sin(2. * M_PI * X[2]);

    return;
}

void coord(int i, int j, int k, double *X) {

    /* returns zone-centered values for coordinates */
    X[0] = startx[0];
    X[1] = Xgrid[i][j][0];
    X[2] = Xgrid[i][j][1];
    X[3] = Xgrid[i][j][2];

    return;
}

void set_units(char *munitstr, char *time) {
    double MBH;
    double t_step;
    //	MBH = 7e7;
    // MBH = 1e8;
    // MBH = 6.2e9;
    // MBH=1.54e8; //ngc1052
    // //MBH=11.;
    MBH = 8.152e39; // 4.1e6;
    //								MBH=4.1e6;
    sscanf(munitstr, "%lf", &M_unit);
    sscanf(time, "%lf", &t_step);
    M_unit = M_unit; // *exp(t_step/4470.46);
    /** input parameters appropriate to Sgr A* **/
    // MBH *=MSUN; //1.251e+43; // MSUN;

    /** from this, calculate units of length, time, mass,
        and derivative units **/
    L_unit = GNEWT * MBH / (CL * CL);
    T_unit = L_unit / CL;

    RHO_unit = M_unit / pow(L_unit, 3);
    U_unit = RHO_unit * CL * CL;
    B_unit = CL * sqrt(4. * M_PI * RHO_unit);
    Ne_unit = RHO_unit / (MP + ME);

    max_tau_scatt = 1e-5 * (6. * L_unit) * RHO_unit * 0.4;
    // Get the rank of the process
    int world_rank;
#if MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0) {
#endif
        fprintf(stderr, "\nUNITS\n");
        fprintf(stderr, "Mbh,L,T,M: %g %g %g %g\n", MBH, L_unit, T_unit,
                M_unit);
        fprintf(stderr, "RHO_unit,U_unit,B_unit: %g %g %g\n", RHO_unit, U_unit,
                B_unit);
        fprintf(stderr, "Ne: %g \n", Ne_unit);
        fprintf(stderr, "max_tau_scatt: %g\n", max_tau_scatt);
#if MPI
    }
#endif
}

/* set up all grid functions */
void init_geometry() {
    int i, j;
    double X[NDIM];

    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {

            /* zone-centered, independent of X[3] */
            // coord(i, j, 0, X);
            calc_coord(j, nx, ndimini, block_info[i].lb,
                       block_info[i].dxc_block, X);
            gcov_func(X, geom[i][j].gcov);

            geom[i][j].g = gdet_func(geom[i][j].gcov, X);

            gcon_func(X, geom[i][j].gcon);
        }
    }
    /* done! */
}

/*

   return solid angle between points x2i, x2f
   and over all x3.

 */

double dOmega_func(double x2i, double x2f, double x3i, double x3f) {
    double dO;
#if(CKS)
    dO = (x3f - x3i) * (-cos(M_PI * x2f) + cos(M_PI * x2i) );
#elif(MKS)             
	dO = (x3f - x3i) *
	    (-cos(M_PI * x2f + 0.5 * (hslope) * sin(2 * M_PI * x2f))
	     + cos(M_PI * x2i + 0.5 * (hslope) * sin(2 * M_PI * x2i))
	    );
#endif                    

    return (dO);
}

void *malloc_rank1(int n1, int size) {
    void *A;
    if ((A = malloc(n1 * size)) == NULL) {
        fprintf(stderr, "malloc failure in malloc_rank1\n");
        exit(123);
    }

    return A;
}

static void **malloc_rank2(int n1, int n2, int size) {

    void **A;
    int i;

    if ((A = (void **)malloc(n1 * sizeof(void *))) == NULL) {
        fprintf(stderr, "malloc failure in malloc_rank2\n");
        exit(124);
    }

    for (i = 0; i < n1; i++) {
        A[i] = malloc_rank1(n2, size);
    }

    return A;
}

double ***malloc_rank3(int n1, int n2, int n3) {
    double ***A;
    double *space;
    int i, j;

    space = malloc_rank1(n1 * n2 * n3, sizeof(double));

    A = malloc_rank1(n1, sizeof(double *));

    for (i = 0; i < n1; i++) {
        A[i] = malloc_rank1(n2, sizeof(double *));
        for (j = 0; j < n2; j++) {
            A[i][j] = &(space[n3 * (j + n2 * i)]);
        }
    }

    return A;
}

double ****malloc_rank4(int n1, int n2, int n3, int n4) {
    double ****A;
    double *space;
    int i, j, k;

    space = malloc_rank1(n1 * n2 * n3 * n4, sizeof(double));

    A = malloc_rank1(n1, sizeof(double *));

    for (i = 0; i < n1; i++) {
        A[i] = malloc_rank1(n2, sizeof(double *));
        for (j = 0; j < n2; j++) {
            A[i][j] = malloc_rank1(n3, sizeof(double *));
            for (k = 0; k < n3; k++) {
                A[i][j][k] = &(space[n4 * (k + n3 * (j + n2 * i))]);
            }
        }
    }

    return A;
}

double *****malloc_rank5(int n1, int n2, int n3, int n4, int n5) {

    double *****A;
    double *space;
    int i, j, k, l;

    space = malloc_rank1(n1 * n2 * n3 * n4 * n5, sizeof(double));

    A = malloc_rank1(n1, sizeof(double *));

    for (i = 0; i < n1; i++) {
        A[i] = malloc_rank1(n2, sizeof(double *));
        for (j = 0; j < n2; j++) {
            A[i][j] = malloc_rank1(n3, sizeof(double *));
            for (k = 0; k < n3; k++) {
                A[i][j][k] = malloc_rank1(n4, sizeof(double *));
                for (l = 0; l < n4; l++) {
                    A[i][j][k][l] =
                        &(space[n5 * (l + n4 * (k + n3 * (j + n2 * i)))]);
                }
            }
        }
    }

    return A;
}

void init_storage(void) {
    int i;

    p = (double ****)malloc_rank1(NPRIM, sizeof(double *));
    for (i = 0; i < NPRIM; i++)
        p[i] = malloc_rank3(N1, N2, N3);

    geom = (struct of_geom **)malloc_rank2(N1, N2, sizeof(struct of_geom));

    int j, k;
    shared_spect = (struct of_spectrum ***)malloc(
        N_THBINS * sizeof(struct of_spectrum **));
    for (i = 0; i < N_THBINS; i++) {
        shared_spect[i] = (struct of_spectrum **)malloc(
            N_PHIBINS * sizeof(struct of_spectrum *));
        for (j = 0; j < N_PHIBINS; j++) {
            shared_spect[i][j] = (struct of_spectrum *)malloc(
                N_EBINS * sizeof(struct of_spectrum));
        }
    }

    return;
}
