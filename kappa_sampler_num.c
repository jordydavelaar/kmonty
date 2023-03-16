
//
//  kappasampler
//
//

#include "decs.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

// Function that, when zero, gives gamma for which dN/d log gam is maximized
double dfdgam(double ge, void *params)
{
  double Thetae = *(double *)params;

  #if KAPPA
  double kap = kappa_synch;
  double w = (kap-3.)/kap * Thetae;
  return 2. + ge * ( ge / ( ge*ge - 1 ) - 1. / gamma_max - (kap+1)/kap/w / (1. + (ge-1.)/kap/w) );
  #else
  return ge - pow(ge,3.) - 2.*Thetae + 3.*pow(ge,2.)*Thetae;
  #endif
}

// electron distribution function (Maxwell-Juettner or kappa below)
// dN / d log gamma
double fdist(double ge, double Thetae)
{
  #if KAPPA
  double kap = kappa_synch;
  double w = (kap-3.)/kap * Thetae;
  return ge*ge*sqrt(ge*ge - 1.)*pow(1. + (ge - 1.)/(kap * w), - kap - 1.)*exp(-ge/gamma_max);
  #else
  return ge*ge*sqrt(ge*ge-1.)*exp(-ge/Thetae);
  #endif
}

void sample_beta_distr_num(double Thetae, double *gamma_e, double *beta_e)
{

double fdist(double ge, double Thetae);
double dfdgam(double ge, void *params);

  double kap = kappa_synch;
  double w = (kap-3.)/kap * Thetae;
  // Relativistic kappa distribution does not like very small Thetae. Ugly kludge.
  if (w < 0.01) {
    *gamma_e = 1.000001;
          *beta_e = sqrt(1. - 1. / (*gamma_e * *gamma_e));
    return;
  }

  // Get maximum for window
  int status, iter = 0, max_iter = 1000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double ge_max = 1. + Thetae;
  double ge_lo = GSL_MAX(1.000001, 0.01*w); // dfdgam -> +inf as ge -> 1+
  double ge_hi = GSL_MAX(100., 5000.*w);

  //printf("Thetae = %e ge_lo = %e ge_hi = %e\n", Thetae, ge_lo, ge_hi);
  gsl_function F;

  //printf("%e %e\n", dfdgam(ge_lo, &params), dfdgam(ge_hi, &params));
  F.function = &dfdgam;
  F.params = &Thetae;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, ge_lo, ge_hi);
  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    ge_max = gsl_root_fsolver_root(s);
    ge_lo = gsl_root_fsolver_x_lower(s);
    ge_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(ge_lo, ge_hi, 0, 1e-7);
  } while (status == GSL_CONTINUE && iter < max_iter);

  double f_max = fdist(ge_max, Thetae);
  gsl_root_fsolver_free(s);
  //fprintf(stderr, "max is %g at %g for %g\n", f_max, ge_max, Thetae);

  // Sample electron gamma
  double ge_samp;
  do {
    double lge_min = log(GSL_MAX(1., 0.01*w));
    double lge_max = log(GSL_MAX(100., 5000.*w));
    ge_samp = exp(lge_min + (lge_max - lge_min)*monty_rand());
  } while (fdist(ge_samp, Thetae)/f_max < monty_rand());

  *gamma_e = ge_samp;
  *beta_e = sqrt(1. - 1. / (*gamma_e * *gamma_e));
}

