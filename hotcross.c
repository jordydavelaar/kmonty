

#include "decs.h"

/*

   given energy of photon in fluid rest frame w, in units of electron rest mass
   energy, and temperature of plasma, again in electron rest-mass units, return
   hot
   cross section in cgs.

   This has been checked against Wienke's Table 1, with some disagreement at
   the one part in 10^{-3} level, see wienke_table_1 in the subdirectory
   hotcross.
   It is not clear what this is due to, but Table 1 does appear to have been
   evaluated
   using Monte Carlo integration (!).

   A better way to do this would be to make a table in w*thetae and w/thetae;
   most
   of the variation is accounted for by w*thetae.

 */

#define MINW 1.e-12
#define MAXW 1.e15
#define MINT 0.0001
#define MAXT 1.e4
#define NW 220
#define NT 80

#define HOTCROSS "hotcross.dat"

double table[NW + 1][NT + 1];
double dlw, dlT, lminw, lmint;

void init_hotcross(void) {
    int i, j, idum, jdum, nread;
    double lw, lT;
    double total_compton_cross_num(double w, double thetae);
    FILE *fp;

    dlw = log10(MAXW / MINW) / NW;
    dlT = log10(MAXT / MINT) / NT;
    lminw = log10(MINW);
    lmint = log10(MINT);

    fp = fopen(HOTCROSS, "r");
    if (fp == NULL || KAPPA || MIXED || THERMAL) {
        for (i = 0; i <= NW; i++)
            for (j = 0; j <= NT; j++) {

                lw = lminw + i * dlw;
                lT = lmint + j * dlT;
                table[i][j] =
                    log10(total_compton_cross_num(pow(10., lw), pow(10., lT)));
                if (isnan(table[i][j])) {
                    fprintf(stderr, "%d %d %g %g\n", i, j, lw, lT);
                    exit(0);
                }
                //                                        fprintf(stderr,"i %d
                //                                        imax %d j %d jmax %d
                //                                        tablevalue
                //                                        %e\n",i,NW,j,NT,table[i][j]);
            }
        fp = fopen(HOTCROSS, "w");
        if (fp == NULL) {
            fprintf(stderr, "couldn't write to file\n");
            exit(0);
        }
        for (i = 0; i <= NW; i++)
            for (j = 0; j <= NT; j++) {
                lw = lminw + i * dlw;
                lT = lmint + j * dlT;
                fprintf(fp, "%d %d %g %g %15.10g\n", i, j, lw, lT, table[i][j]);
            }
    } else {

        for (i = 0; i <= NW; i++)
            for (j = 0; j <= NT; j++) {
                nread = fscanf(fp, "%d %d %lf %lf %lf\n", &idum, &jdum, &lw,
                               &lT, &table[i][j]);
                if (isnan(table[i][j]) || nread != 5) {
                    fprintf(stderr, "error on table read: %d %d\n", i, j);
                    exit(0);
                }
            }
    }

    fclose(fp);

    return;
}

double total_compton_cross_lkup(double w, double thetae) {
    int i, j;
    double lw, lT, di, dj, lcross;
    double total_compton_cross_num(double w, double thetae);
    double hc_klein_nishina(double we);

    /* cold/low-energy: just use thomson cross section */
    if (w * thetae < 1.e-6)
        return (SIGMA_THOMSON);

    /* cold, but possible high energy photon: use klein-nishina */
    if (thetae < MINT)
        return (hc_klein_nishina(w) * SIGMA_THOMSON);

    /* in-bounds for table */
    if ((w > MINW && w < MAXW) && (thetae > MINT && thetae < MAXT)) {

        lw = log10(w);
        lT = log10(thetae);
        i = (int)((lw - lminw) / dlw);
        j = (int)((lT - lmint) / dlT);
        di = (lw - lminw) / dlw - i;
        dj = (lT - lmint) / dlT - j;

        lcross = (1. - di) * (1. - dj) * table[i][j] +
                 di * (1. - dj) * table[i + 1][j] +
                 (1. - di) * dj * table[i][j + 1] +
                 di * dj * table[i + 1][j + 1];

        if (isnan(lcross)) {
            fprintf(stderr, "lcross isnan %g %g %d %d %g %g\n", lw, lT, i, j,
                    di, dj);
        }

        return (pow(10., lcross));
    }

    fprintf(stderr, "out of bounds: %g %g\n", w, thetae);
    return (total_compton_cross_num(w, thetae));
}

#define MAXGAMMA 12.
#define DMUE 0.05
#define DGAMMAE 0.05

double calc_norm(double thetae) {
    // we integrate over beta since this is nicely bound betwoon 0 and 1,
    // indsted of 1 and infinity.
    double kappa_function_int(double beta, void *params);

    double norm;
    double result, error;
    gsl_function integrand;
    integrand.function = &kappa_function_int;
    integrand.params = &thetae;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    /* Integrator parameters */
    double absolute_error = 0.;
    double relative_error = 1.e-6;
    size_t limit = 1000;

    gsl_integration_qag(&integrand, 0, log(1. + 1000. * thetae), absolute_error,
                        relative_error, limit, GSL_INTEG_GAUSS61, w, &result,
                        &error);

    gsl_integration_workspace_free(w);
    norm = 1. / result;

    return norm;
}

double total_compton_cross_num(double w, double thetae) {
    double dmue, dgammae, mue, gammae, f, cross, norm;
    double dNdgammae(double thetae, double gammae);
    double boostcross(double w, double mue, double gammae);
    double hc_klein_nishina(double we);

    if (isnan(w)) {
        fprintf(stderr, "compton cross isnan: %g %g\n", w, thetae);
        return (0.);
    }

    /* check for easy-to-do limits */
    if (thetae < MINT && w < MINW)
        return (SIGMA_THOMSON);
    if (thetae < MINT)
        return (hc_klein_nishina(w) * SIGMA_THOMSON);

    dmue = DMUE;
    dgammae = thetae * DGAMMAE;

/* integrate over mu_e, gamma_e, where mu_e is the cosine of the
   angle between k and u_e, and the angle k is assumed to lie,
   wlog, along the z axis */
#if (KAPPA)
    norm = calc_norm(thetae);
#elif (THERMAL)
    norm = 1;
#endif
    cross = 0.;
    for (mue = -1. + 0.5 * dmue; mue < 1.; mue += dmue)
        for (gammae = 1. + 0.5 * dgammae; gammae < 1. + MAXGAMMA * thetae;
             gammae += dgammae) {

            f = 0.5 * norm * dNdgammae(thetae, gammae);

            cross += dmue * dgammae * boostcross(w, mue, gammae) * f;

            if (isnan(cross)) {
                fprintf(stderr, "cross is nan %g %g %g %g %g %g\n", w, thetae,
                        mue, gammae, dNdgammae(thetae, gammae),
                        boostcross(w, mue, gammae));
            }
        }

    return (cross * SIGMA_THOMSON);
}

double dNdgammae(double thetae, double gammae) {
    double dNdgammae_nth(double thetae, double gammae);
    double dNdgammae_th(double thetae, double gammae);

#if (KAPPA)
    return dNdgammae_nth(thetae, gammae);
#elif (THERMAL)
    return dNdgammae_th(thetae, gammae);
#endif
    return 0;
}

double dNdgammae_th(double thetae, double gammae) {
    double K2f;

    if (thetae > 1.e-2) {
        K2f = gsl_sf_bessel_Kn(2, 1. / thetae) * exp(1. / thetae);
    } else {
        K2f = sqrt(M_PI * thetae / 2.);
    }

    return ((gammae * sqrt(gammae * gammae - 1.) / (thetae * K2f)) *
            exp(-(gammae - 1.) / thetae));
}

// use beta instead of gamma, because its nicely bound between 0 and 1.
double kappa_function_int(double beta, void *params) {
    double gamma = exp(beta);

    double thetae = *(double *)params;
    double kappa = kappa_synch;
    double w = (kappa - 3.) / kappa * thetae;

    // double gamma = 1/sqrt(1-beta*beta);
    // return beta*beta* pow(gamma,5.)*pow(( 1 + ( gamma -
    // 1)/(kappa*thetae)),-(kappa+1))*exp(-gamma/gamma_max);

    return gamma * gamma * sqrt(gamma * gamma - 1.) *
           pow((1 + (gamma - 1) / (kappa * w)), -(kappa + 1)) *
           exp(-gamma / gamma_max);
}

double dNdgammae_nth(double thetae, double gammae) {
    double kappa = kappa_synch;
    double w = (kappa - 3.) / kappa * thetae;
    // This is not yet normalized!
    double dNdgam = gammae * sqrt(gammae * gammae - 1.) *
                    pow((1 + (gammae - 1) / (kappa * w)), -(kappa + 1)) *
                    exp(-gammae / gamma_max);
    return dNdgam;
}

double boostcross(double w, double mue, double gammae) {
    double we, boostcross, v;
    double hc_klein_nishina(double we);

    /* energy in electron rest frame */
    v = sqrt(gammae * gammae - 1.) / gammae;
    we = w * gammae * (1. - mue * v);

    boostcross = hc_klein_nishina(we) * (1. - mue * v);

    if (boostcross > 2) {
        fprintf(stderr, "w,mue,gammae: %g %g %g\n", w, mue, gammae);
        fprintf(stderr, "v,we, boostcross: %g %g %g\n", v, we, boostcross);
        fprintf(stderr, "kn: %g %g %g\n", v, we, boostcross);
    }

    if (isnan(boostcross)) {
        fprintf(stderr, "isnan: %g %g %g\n", w, mue, gammae);
        exit(0);
    }

    return (boostcross);
}

double hc_klein_nishina(double we) {
    double sigma;

    if (we < 1.e-3)
        return (1. - 2. * we);

    sigma = (3. / 4.) *
            (2. / (we * we) +
             (1. / (2. * we) - (1. + we) / (we * we * we)) * log(1. + 2. * we) +
             (1. + we) / ((1. + 2. * we) * (1. + 2. * we)));

    return (sigma);
}
