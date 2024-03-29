

#include "decs.h"
#include "sphere_model.h"
/*

   this is the main photon orbit integrator

 */

#define FAST_CPY(in, out)                                                      \
    {                                                                          \
        out[0] = in[0];                                                        \
        out[1] = in[1];                                                        \
        out[2] = in[2];                                                        \
        out[3] = in[3];                                                        \
    }

#define ETOL 1.e-3
#define MAX_ITER 2
void push_photon(double X[NDIM], double Kcon[NDIM], double dKcon[NDIM],
                 double dl, double *E0, int n) {
    //  double lconn[NDIM][NDIM][NDIM];
    double lconn[NDIM][NDIM][NDIM];
    double Kcont[NDIM], K[NDIM], dK;
    double Xcpy[NDIM], Kcpy[NDIM], dKcpy[NDIM];
    double Gcov[NDIM][NDIM], E1;
    double dl_2, err, errE;
    int i, k, iter;
#if MKS
    if (X[1] < startx[1])
        return;
#endif
    FAST_CPY(X, Xcpy);
    FAST_CPY(Kcon, Kcpy);
    FAST_CPY(dKcon, dKcpy);

    dl_2 = 0.5 * dl;
    /* Step the position and estimate new wave vector */
    for (i = 0; i < NDIM; i++) {
        dK = dKcon[i] * dl_2;
        Kcon[i] += dK;
        K[i] = Kcon[i] + dK;
        X[i] += Kcon[i] * dl;
    }

    get_connection(X, lconn);

    /* We're in a coordinate basis so take advantage of symmetry in the
     * connection
     */
    iter = 0;
    do {
        iter++;
        FAST_CPY(K, Kcont);

        err = 0.;
        for (k = 0; k < 4; k++) {
            dKcon[k] = -2. * (Kcont[0] * (lconn[k][0][1] * Kcont[1] +
                                          lconn[k][0][2] * Kcont[2] +
                                          lconn[k][0][3] * Kcont[3]) +
                              Kcont[1] * (lconn[k][1][2] * Kcont[2] +
                                          lconn[k][1][3] * Kcont[3]) +
                              lconn[k][2][3] * Kcont[2] * Kcont[3]);

            dKcon[k] -= (lconn[k][0][0] * Kcont[0] * Kcont[0] +
                         lconn[k][1][1] * Kcont[1] * Kcont[1] +
                         lconn[k][2][2] * Kcont[2] * Kcont[2] +
                         lconn[k][3][3] * Kcont[3] * Kcont[3]);

            K[k] = Kcon[k] + dl_2 * dKcon[k];
            err += fabs((Kcont[k] - K[k]) / (K[k] + SMALL));
        }
    } while (err > ETOL && iter < MAX_ITER);

    FAST_CPY(K, Kcon);

    gcov_func(X, Gcov);
    E1 = -(Kcon[0] * Gcov[0][0] + Kcon[1] * Gcov[0][1] + Kcon[2] * Gcov[0][2] +
           Kcon[3] * Gcov[0][3]);
    errE = fabs((E1 - (*E0)) / (*E0));

    if (n < 7 && (errE > 1.e-4 || err > ETOL || isnan(err) || isinf(err))) {
        FAST_CPY(Xcpy, X);
        FAST_CPY(Kcpy, Kcon);
        FAST_CPY(dKcpy, dKcon);
        push_photon(X, Kcon, dKcon, 0.5 * dl, E0, n + 1);
        push_photon(X, Kcon, dKcon, 0.5 * dl, E0, n + 1);
        E1 = *E0;
    }

    *E0 = E1;

    /* done! */
}

/* spare photon integrator: 4th order Runge-Kutta */
void push_photon4(double X[], double K[], double dK[], double dl) {

    int k;
    double lconn[NDIM][NDIM][NDIM];
    double Kt[NDIM], Xt[NDIM];
    double f1x[NDIM], f2x[NDIM], f3x[NDIM], f4x[NDIM];
    double f1k[NDIM], f2k[NDIM], f3k[NDIM], f4k[NDIM];
    double dl_2 = 0.5 * dl;

    for (k = 0; k < NDIM; k++)
        f1x[k] = K[k];

    get_connection(X, lconn);
    for (k = 0; k < NDIM; k++) {
        f1k[k] = -2. * (K[0] * (lconn[k][0][1] * K[1] + lconn[k][0][2] * K[2] +
                                lconn[k][0][3] * K[3]) +
                        K[1] * (lconn[k][1][2] * K[2] + lconn[k][1][3] * K[3]) +
                        lconn[k][2][3] * K[2] * K[3]);

        f1k[k] -= (lconn[k][0][0] * K[0] * K[0] + lconn[k][1][1] * K[1] * K[1] +
                   lconn[k][2][2] * K[2] * K[2] + lconn[k][3][3] * K[3] * K[3]);
    }

    for (k = 0; k < NDIM; k++) {
        Kt[k] = K[k] + dl_2 * f1k[k];
        f2x[k] = Kt[k];
        Xt[k] = X[k] + dl_2 * f1x[k];
    }

    get_connection(Xt, lconn);
    for (k = 0; k < NDIM; k++) {
        f2k[k] =
            -2. * (Kt[0] * (lconn[k][0][1] * Kt[1] + lconn[k][0][2] * Kt[2] +
                            lconn[k][0][3] * Kt[3]) +
                   Kt[1] * (lconn[k][1][2] * Kt[2] + lconn[k][1][3] * Kt[3]) +
                   lconn[k][2][3] * Kt[2] * Kt[3]);

        f2k[k] -=
            (lconn[k][0][0] * Kt[0] * Kt[0] + lconn[k][1][1] * Kt[1] * Kt[1] +
             lconn[k][2][2] * Kt[2] * Kt[2] + lconn[k][3][3] * Kt[3] * Kt[3]);
    }

    for (k = 0; k < NDIM; k++) {
        Kt[k] = K[k] + dl_2 * f2k[k];
        f3x[k] = Kt[k];
        Xt[k] = X[k] + dl_2 * f2x[k];
    }

    get_connection(Xt, lconn);
    for (k = 0; k < NDIM; k++) {
        f3k[k] =
            -2. * (Kt[0] * (lconn[k][0][1] * Kt[1] + lconn[k][0][2] * Kt[2] +
                            lconn[k][0][3] * Kt[3]) +
                   Kt[1] * (lconn[k][1][2] * Kt[2] + lconn[k][1][3] * Kt[3]) +
                   lconn[k][2][3] * Kt[2] * Kt[3]);

        f3k[k] -=
            (lconn[k][0][0] * Kt[0] * Kt[0] + lconn[k][1][1] * Kt[1] * Kt[1] +
             lconn[k][2][2] * Kt[2] * Kt[2] + lconn[k][3][3] * Kt[3] * Kt[3]);
    }

    for (k = 0; k < NDIM; k++) {
        Kt[k] = K[k] + dl * f3k[k];
        f4x[k] = Kt[k];
        Xt[k] = X[k] + dl * f3x[k];
    }

    get_connection(Xt, lconn);
    for (k = 0; k < NDIM; k++) {
        f4k[k] =
            -2. * (Kt[0] * (lconn[k][0][1] * Kt[1] + lconn[k][0][2] * Kt[2] +
                            lconn[k][0][3] * Kt[3]) +
                   Kt[1] * (lconn[k][1][2] * Kt[2] + lconn[k][1][3] * Kt[3]) +
                   lconn[k][2][3] * Kt[2] * Kt[3]);

        f4k[k] -=
            (lconn[k][0][0] * Kt[0] * Kt[0] + lconn[k][1][1] * Kt[1] * Kt[1] +
             lconn[k][2][2] * Kt[2] * Kt[2] + lconn[k][3][3] * Kt[3] * Kt[3]);
    }

    for (k = 0; k < NDIM; k++) {
        X[k] +=
            0.166666666666667 * dl * (f1x[k] + 2. * (f2x[k] + f3x[k]) + f4x[k]);
        K[k] +=
            0.166666666666667 * dl * (f1k[k] + 2. * (f2k[k] + f3k[k]) + f4k[k]);
    }

    init_dKdlam(X, K, dK);

    /* done */
}

void init_dKdlam(double X[], double Kcon[], double dK[]) {
    int k;
    double lconn[NDIM][NDIM][NDIM];

    get_connection(X, lconn);

    for (k = 0; k < 4; k++) {

        dK[k] =
            -2. *
            (Kcon[0] * (lconn[k][0][1] * Kcon[1] + lconn[k][0][2] * Kcon[2] +
                        lconn[k][0][3] * Kcon[3]) +
             Kcon[1] * (lconn[k][1][2] * Kcon[2] + lconn[k][1][3] * Kcon[3]) +
             lconn[k][2][3] * Kcon[2] * Kcon[3]);

        dK[k] -= (lconn[k][0][0] * Kcon[0] * Kcon[0] +
                  lconn[k][1][1] * Kcon[1] * Kcon[1] +
                  lconn[k][2][2] * Kcon[2] * Kcon[2] +
                  lconn[k][3][3] * Kcon[3] * Kcon[3]);
    }

    return;
}
