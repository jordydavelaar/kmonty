

/*

   model-independent radiation-related utilities.

 */

#include "decs.h"
#include <gsl/gsl_sf_gamma.h>

double anu_synch_powerlaw(double nu, double Ne, double Thetae, double B,
                          double theta) {
    double nuc, sth, X, factor;
    double As;
    double p = 3;
    double gmin = 25.;
    double gmax = 1.e7;

    sth = sin(theta);
    nuc = EE * B / (2. * M_PI * ME * CL);
    factor = (Ne * pow(EE, 2.)) / (nu * ME * CL);

    X = nu / (nuc * sth);

    As = pow(3., (p + 1) / 2.) * (p - 1) /
         (4 * (pow(gmin, 1 - p) - pow(gmax, 1 - p)));
    As *= gsl_sf_gamma((3 * p + 2) / 12.) * gsl_sf_gamma((3 * p + 22) / 12.) *
          pow(X, -(p + 2) / 2.);

    return As * factor;
}

double anu_synch_kappa(double nu, double Ne, double Thetae, double B,
                       double theta) {
    // absortion for the kappa distribution function, see Pandya et al. 2016
    double nuc, sth, nus, x, w, X_kappa, factor;
    double A_low, A_high, A_s;
    double kappa = kappa_synch;
    // w = Thetae; //sqrt(2./9./kappa *Thetae * Thetae);
    w = (kappa - 3.) / kappa * Thetae;
    nuc = EE * B / (2. * M_PI * ME * CL);
    sth = sin(theta);

    factor = Ne * EE / (B * sth); // * exp(-nu/5e14);

    nus = nuc * sth * pow(w * kappa, 2);

    //    if (nu > 1.e12 * nus || Thetae < THETAE_MIN){
    //		printf("nu %e nus %e nuc %e B %e Te %e th %e sth
    //%e\n",nu,nus,nuc,B,w, theta,sth);
    X_kappa = nu / nus;

    if (sth < 1e-150 || Thetae < THETAE_MIN || X_kappa > 1e10) {
        return (0.);
    }
    // if(X_kappa>1.e7 || X_kappa<0.002){
    //     return 0;
    // }
    double a = kappa - 1. / 3.;
    double b = kappa + 1.;
    double c = kappa + 2. / 3.;
    double z = -kappa * w;
    double hyp2F1;
    if (fabs(z) == 1.)
        return 0;
    if (fabs(z) < 1)
        hyp2F1 = gsl_sf_hyperg_2F1(a, b, c, z);
    else {
        hyp2F1 = pow(1. - z, -a) * gsl_sf_gamma(c) * tgamma(b - a) /
                     (tgamma(b) * tgamma(c - a)) *
                     gsl_sf_hyperg_2F1(a, c - b, a - b + 1., 1. / (1. - z)) +
                 pow(1. - z, -b) * tgamma(c) * tgamma(a - b) /
                     (tgamma(a) * tgamma(c - b)) *
                     gsl_sf_hyperg_2F1(b, c - a, b - a + 1., 1. / (1. - z));
    }

    A_low = pow(X_kappa, -5. / 3.) * pow(3, 1. / 6.) * (10. / 41.) *
            pow(2 * M_PI, 2) / pow(w * kappa, 16. / 3. - kappa) * (kappa - 2.) *
            (kappa - 1.) * kappa / (3. * kappa - 1.) * tgamma(5. / 3.) * hyp2F1;
    A_high = pow(X_kappa, -(3. + kappa) / 2.) * (2. * pow(M_PI, 5. / 2.) / 3.) *
             ((kappa - 2.) * (kappa - 1.) * kappa / pow(w * kappa, 5.)) *
             (2 * tgamma(2. + kappa / 2.) / (2. + kappa) - 1.) *
             (pow(3. / kappa, 19. / 4.) + 3. / 5.);

    x = pow(-7. / 4. + 8. * kappa / 5., -43. / 50.);

    A_s = pow((pow(A_low, -x) + pow(A_high, -x)), -1. / x);

    if (A_s != A_s)
        fprintf(stderr, "poblems in anu!!!\n");

    return (factor * A_s);
}

double Bnu_inv(double nu, double Thetae) {

    double x;

    x = HPL * nu / (ME * CL * CL * Thetae);

    if (x < 1e-3) /* Taylor expand */
        return ((2. * HPL / (CL * CL)) /
                (x / 24. * (24. + x * (12. + x * (4. + x)))));
    else
        return ((2. * HPL / (CL * CL)) / (exp(x) - 1.));
}

double jnu_inv(double nu, double Thetae, double Ne, double B, double theta) {
    double j;
    int ACCZONE = 0;
    j = jnu_synch(nu, Ne, Thetae, B, theta, ACCZONE);
    if (isnan(j))
        printf("j is nan\n");

    return (j / (nu * nu));
}

/* return Lorentz invariant scattering opacity */
double alpha_inv_scatt(double nu, double Thetae, double Ne) {
    double kappa;
#if COMPTON

    kappa = kappa_es(nu, Thetae);

    return (nu * kappa * Ne * MP);
#else
    return 0;
#endif
}

/* return Lorentz invariant absorption opacity */
double alpha_inv_abs(double nu, double Thetae, double Ne, double B,
                     double theta, int ACCZONE) {
    double a;
#if (THERMAL)
    a = alpha_inv_abs_th(nu, Thetae, Ne, B, theta);
#elif (KAPPA || POWERLAW)
    a = alpha_inv_abs_nth(nu, Thetae, Ne, B, theta) * exp(-nu / nu_cutoff);
#elif (MIXED)
    if (ACCZONE) {
        double a_th = alpha_inv_abs_th(nu, Thetae, Ne, B, theta);
        double a_nth = alpha_inv_abs_nth(nu, Thetae, Ne, B, theta);
        a = perct_thermal * a_th + (1 - perct_thermal) * a_nth;
    } else
        a = alpha_inv_abs_th(nu, Thetae, Ne, B, theta);
#endif
    //	printf("nth %e th %e ratio %e nu %e\n",a, a2, a/a2,nu);
    return a;
}

double alpha_inv_abs_th(double nu, double Thetae, double Ne, double B,
                        double theta) {
    double j, bnu;

    j = jnu_synch_th(nu, Ne, Thetae, B, theta) / (nu * nu);
    bnu = Bnu_inv(nu, Thetae);
    if (j > 0)
        return (j / (bnu + 1.e-100));
    return 0;
}

double alpha_inv_abs_nth(double nu, double Thetae, double Ne, double B,
                         double theta) {
#if (KAPPA)
    return (nu * anu_synch_kappa(nu, Ne, Thetae, B, theta));
#else
    return (nu * anu_synch_powerlaw(nu, Ne, Thetae, B, theta));
#endif
}

/* return electron scattering opacity, in cgs */
double kappa_es(double nu, double Thetae) {
    double Eg;

    /* assume pure hydrogen gas to
       convert cross section to opacity */
    Eg = HPL * nu / (ME * CL * CL);
    if (Eg > 1e75)
        fprintf(stderr, "out of bounds: %g %g %g\n", Eg, Thetae, nu);
    return (total_compton_cross_lkup(Eg, Thetae) / MP);
}

/* get frequency in fluid frame, in Hz */
double get_fluid_nu(double X[4], double K[4], double Ucov[NDIM]) {
    double ener, nu;

    /* this is the energy in electron rest-mass units */
    ener = -(K[0] * Ucov[0] + K[1] * Ucov[1] + K[2] * Ucov[2] + K[3] * Ucov[3]);

    nu = ener * ME * CL * CL / HPL;
    if (nu > 1e75) {
        fprintf(stderr, "problem in fluid nu; %e", nu);
        fprintf(stderr, "problem get_fluid_nu, K: %g %g %g %g\n", K[0], K[1],
                K[2], K[3]);
        fprintf(stderr, "problem get_fluid_nu, X: %g %g %g %g\n", X[0], X[1],
                X[2], X[3]);
        fprintf(stderr, "problem get_fluid_nu, U: %g %g %g %g\n", Ucov[0],
                Ucov[1], Ucov[2], Ucov[3]);
    }
    if (isnan(ener)) {
        fprintf(stderr, "isnan get_fluid_nu, K: %g %g %g %g\n", K[0], K[1],
                K[2], K[3]);
        fprintf(stderr, "isnan get_fluid_nu, X: %g %g %g %g\n", X[0], X[1],
                X[2], X[3]);
        fprintf(stderr, "isnan get_fluid_nu, U: %g %g %g %g\n", Ucov[0],
                Ucov[1], Ucov[2], Ucov[3]);
    }

    return nu;
}

/* return angle between magnetic field and wavevector */
double get_bk_angle(double X[NDIM], double K[NDIM], double Ucov[NDIM],
                    double Bcov[NDIM], double B) {

    double k, mu;
    //    return (M_PI / 2.);

    if (B == 0.) {
        return (M_PI / 2.);
        printf("B is zero\n");
    }
    k = fabs(K[0] * Ucov[0] + K[1] * Ucov[1] + K[2] * Ucov[2] + K[3] * Ucov[3]);

    /* B is in cgs but Bcov is in code units */
    mu = (K[0] * Bcov[0] + K[1] * Bcov[1] + K[2] * Bcov[2] + K[3] * Bcov[3]) /
         (k * B / B_unit);

    if (fabs(mu) > 1.)
        mu /= fabs(mu);

    return (acos(mu));
}
