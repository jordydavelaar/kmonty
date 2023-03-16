#include "decs.h"
#include <gsl/gsl_sf_gamma.h>

#pragma omp threadprivate(r)

/*

   "mixed" emissivity formula

   interpolates between Petrosian limit and
   classical thermal synchrotron limit

   good for Thetae > 1

 */

double jnu_synch(double nu, double Ne, double Thetae, double B, double theta,
                 int ACCZONE) {
  double j;
#if (THERMAL)
  j = jnu_synch_th(nu, Ne, Thetae, B, theta);
#elif (KAPPA || POWERLAW)
  j = jnu_synch_nth(nu, Ne, Thetae, B, theta); // *exp(-nu/5e14);
#elif (MIXED)
  if (ACCZONE)
    j = (1 - perct_thermal) * jnu_synch_nth(nu, Ne, Thetae, B, theta) +
        perct_thermal * jnu_synch_th(nu, Ne, Thetae, B, theta);
  else
    j = jnu_synch_th(nu, Ne, Thetae, B, theta);

// *exp(-nu/5e14);

#endif

  return (j);
}

double int_jnu(double Ne, double Thetae, double Bmag, double nu) {
  /* Returns energy per unit time at *
   * frequency nu in cgs */

  double j_nu = 0;
#if (THERMAL)
  j_nu = int_jnu_th(Ne, Thetae, Bmag, nu);
#elif (KAPPA || POWERLAW)
  j_nu = int_jnu_nth(Ne, Thetae, Bmag, nu);
#endif
  return j_nu;
}

double jnu_integrand(double th, void *params) {
  double j_nu = 0;
#if (THERMAL)
  j_nu = jnu_integrand_th(th, params);
#elif (KAPPA || POWERLAW)
  j_nu = jnu_integrand_nth(th, params);
#endif
  return j_nu;
}



double jnu_integrand_kappa(double th, void *params) {

  double jnu;
  double K = *(double *)params;
  double sth = sin(th);
  double X_kappa = K / sth;
  double x, factor;
  double J_low, J_high, J_s;
  double kappa = kappa_synch;

  sth = sin(th);
  if ( X_kappa > 2.e8)
    return 0.;
  //        if (sth < 1.e-40) //|| x > 2.e12)
 // if(th<1e-2)
//	return 0;
  factor = (sth*sth);

 // X_kappa = x;

  J_low = pow(X_kappa, 1. / 3.) * 4. * M_PI * tgamma(kappa - 4. / 3.) /
          (pow(3., 7. / 3.) * tgamma(kappa - 2.));

  J_high = pow(X_kappa, -(kappa - 2.) / 2.)* pow(3., (kappa - 1.) / 2.) *
           (kappa - 2.) * (kappa - 1.) / 4. *
           tgamma(kappa / 4. - 1. / 3.) *
           tgamma(kappa / 4. + 4. / 3.);
  x = 3. * pow(kappa, -3. / 2.);

  J_s = pow((pow(J_low, -x) + pow(J_high, -x)), -1. / x);

  jnu = J_s * factor; // *exp(-X_kappa/1.e7);
  return jnu;
}

double jnu_integrand_powerlaw(double th, void *params) {
 double jnu;
 double K = *(double *)params;
 double sth = sin(th);
 double x = K / sth;

 double factor;
 double Js;
 double p=3.;
 double gmin=25.;
 double gmax=1.e7;

 //if (sth < 1.e-150 || x > 1.e8)
 //   return 0.;

 factor = sth;

 Js = pow(3.,p/2.)*(p-1)*sth/(2*(p+1)*(pow(gmin,1-p)-pow(gmax,1-p)));
 Js *= gsl_sf_gamma((3*p-1)/12.)*gsl_sf_gamma((3*p+19)/12.)*pow(x,-(p-1)/2.);

 return Js*factor;

}

double jnu_integrand_nth(double th, void *params) {
#if(KAPPA)
  return jnu_integrand_kappa(th, params);
#else
  return jnu_integrand_powerlaw(th, params);
#endif
}


#define JCST (EE * EE * EE / (2. * M_PI * ME * CL * CL))
double int_jnu_nth(double Ne, double Thetae, double Bmag, double nu) {
  /* Returns energy per unit time at *
   * frequency nu in cgs
   */
  int ACCZONE = 0;
  double F_eval(double Thetae, double B, double nu, int ACCZONE);

  if (Thetae < THETAE_MIN) {
    return 0.;
  }

  return JCST * Ne * Bmag * F_eval(Thetae, Bmag, nu, ACCZONE);
}

#undef JCST


double jnu_synch_kappa(double nu, double Ne, double Thetae, double B,
                     double theta) {
  // emissivity for the kappa distribution function, see Pandya et al. 2016
  double nuc, sth, nus, x, w, X_kappa, factor;
  double J_low, J_high, J_s;
  double kappa = kappa_synch;
  w = (kappa-3.)/kappa *  Thetae; // Thetae; //sqrt(  2./9./kappa *Thetae * Thetae);
  nuc = EE * B / (2. * M_PI * ME * CL);
  sth = sin(theta);

  factor = (Ne * pow(EE, 2.) * nuc * sth) / CL;

//  nus = nuc * sth * pow(w * kappa, 2.);
  nus = nuc * sth * ( w * kappa) * ( w * kappa);
 // if (Thetae < THETAE_MIN || sth < 1e-150)
 //   return 0.;
  //if (nu > 1.e7 * nus)
  //  return (0.);
  // if (nu > 1.e14 * nus)
  //      return (0.);
 // if(theta<1e-2)
 //       return 0;

  X_kappa = nu / nus;
//  if(X_kappa>1.e10)
  //    return 0;

  J_low = pow(X_kappa, 1. / 3.)* 4. * M_PI * tgamma(kappa - 4. / 3.) /
          (pow(3., 7. / 3.) * tgamma(kappa - 2.));

  J_high = pow(X_kappa, -(kappa - 2.) / 2.) * pow(3., (kappa - 1.) / 2.) *
           (kappa - 2.) * (kappa - 1.) / 4. *
           tgamma(kappa / 4. - 1. / 3.) *
           tgamma(kappa / 4. + 4. / 3.);

  x = 3. * pow(kappa, -3. / 2.);

  J_s = pow((pow(J_low, -x) + pow(J_high, -x)), -1. / x);

  if (J_s != J_s){
    printf("B %e\n", B);
    printf("J_s %e %e %e %e %e %e %e %e\n", X_kappa, B, nuc, nus, J_low, J_high,
           x, kappa);
  }
//fprintf(stderr,"B %e Thetae %e Ne %e\n",B,Thetae,Ne);
//fprintf(stderr,"J_low %e J_high %e x\n",J_low,J_high,x);
//fprintf(stderr,"Js %e Jcgs %e factor %e\n",J_s,J_s*factor,factor);
//exit(1);
  return (J_s * factor)  * exp(-nu / nu_cutoff);
}

double jnu_synch_powerlaw(double nu, double Ne, double Thetae, double B,
                     double theta) {
double nuc,sth,Xs,factor;
double Js;
double p=3.;
double gmin=25.;
double gmax=1.e7;

sth = sin(theta);
nuc = EE * B / (2. * M_PI * ME * CL);
factor = (Ne * pow(EE,2.) * nuc)/CL;

  if (Thetae < THETAE_MIN || sth < 1e-150)
    return 0.;
  if (nu > 1.e8 * nuc)
    return (0.);

Xs = nu/(nuc*sth);

Js = pow(3.,p/2.)*(p-1)*sth/(2*(p+1)*(pow(gmin,1-p)-pow(gmax,1-p)));
Js *= gsl_sf_gamma((3*p-1)/12.)*gsl_sf_gamma((3*p+19)/12.)*pow(Xs,-(p-1)/2.);

return Js*factor;
}

double jnu_synch_nth(double nu, double Ne, double Thetae, double B,
                     double theta) {
#if(KAPPA)
  return jnu_synch_kappa(nu,Ne,Thetae,B,theta);
#elif(POWERLAW)
  return jnu_synch_powerlaw(nu,Ne,Thetae,B,theta);
#endif
  return 0;
}

#define CST 1.88774862536 /* 2^{11/12} */
double jnu_synch_th(double nu, double Ne, double Thetae, double B,
                    double theta) {
  double K2, nuc, nus, x, f, j, sth, xp1, xx;
  double K2_eval(double Thetae);
  // theta = M_PI/2.;
  if (Thetae < THETAE_MIN)
    return 0.;

  K2 = K2_eval(Thetae);

  nuc = EE * B / (2. * M_PI * ME * CL);
  sth = sin(theta);
  nus = (2. / 9.) * nuc * Thetae * Thetae * sth;
 // if (nu > 1.e12 * nus)
 //   return (0.);
  x = nu / nus;
  xp1 = pow(x, 1. / 3.);
  xx = sqrt(x) + CST * sqrt(xp1);
  f = xx * xx;
  j = (M_SQRT2 * M_PI * EE * EE * Ne * nus / (3. * CL * K2)) * f * exp(-xp1);

  return (j);
}

#undef CST

#define JCST (M_SQRT2 * EE * EE * EE / (27 * ME * CL * CL))
double int_jnu_th(double Ne, double Thetae, double Bmag, double nu) {
  /* Returns energy per unit time at
   * *
   * frequency nu in cgs
   */

  double j_fac, K2;
  int ACCZONE = 0;
  double F_eval(double Thetae, double B, double nu, int ACCZONE);
  double K2_eval(double Thetae);

  if (Thetae < THETAE_MIN)
    return 0.;

  K2 = K2_eval(Thetae);
  if (K2 == 0.)
    return 0.;

  j_fac = Ne * Bmag * Thetae * Thetae / K2;

  return JCST * j_fac * F_eval(Thetae, Bmag, nu, ACCZONE);
}

#undef JCST

#define CST 1.88774862536 /* 2^{11/12} */
double jnu_integrand_th(double th, void *params) {

  double K = *(double *)params;
//  th = M_PI/2.;
  double sth = sin(th);
  double x = K / sth;

 // if (sth < 1.e-150 || x > 1.e8)
  //  return 0.;

  return sth * sth * pow(sqrt(x) + CST * pow(x, 1. / 6.), 2.) *
         exp(-pow(x, 1. / 3.));
}

#undef CST

/* Tables */
double F[N_ESAMP + 1], K2[N_ESAMP + 1];
double lK_min, dlK;
double lT_min, dlT;

#define EPSABS 0
#define EPSREL 1.e-6
#if THERMAL
#define KMIN (0.002) //(0.002) //(1e-12) //(0.002) //1e-7 for kappa
#define KMAX (1.e7)  //1e10 for kappa/powerlaw //1e7 for thermal
#else
#define KMIN (1e-7) //(0.002) //(1e-12) //(0.002) //1e-7 for kappa
#define KMAX (1.e10)
#endif
#define TMIN (THETAE_MIN)
#define TMAX (1.e2)
void init_emiss_tables_th(void) {
  int k;
  double result, err, K, T;
  gsl_function func;
  gsl_integration_workspace *w;

  func.function = &jnu_integrand_th;
  func.params = &K;

  lK_min = log(KMIN);
  dlK = log(KMAX / KMIN) / (N_ESAMP);

  lT_min = log(TMIN);
  dlT = log(TMAX / TMIN) / (N_ESAMP);

  /*  build table for F(K) where F(K) is given by
     \int_0^\pi ( (K/\sin\theta)^{1/2} + 2^{11/12}(K/\sin\theta)^{1/6})^2
     \exp[-(K/\sin\theta)^{1/3}]
     so that J_{\nu} = const.*F(K)
   */
  w = gsl_integration_workspace_alloc(5000);
  for (k = 0; k <= N_ESAMP; k++) {
    K = exp(k * dlK + lK_min);
    gsl_integration_qag(&func, 0., M_PI / 2., EPSABS, EPSREL, 5000,
                        GSL_INTEG_GAUSS61, w, &result, &err);

    F_th[k] = log(4 * M_PI * result);
  }
  gsl_integration_workspace_free(w);

  /*  build table for quick evaluation of the bessel function K2 for emissivity
   */
  for (k = 0; k <= N_ESAMP; k++) {
    T = exp(k * dlT + lT_min);
    K2[k] = log(gsl_sf_bessel_Kn(2, 1. / T));
  }

  /* Avoid doing divisions later */
  dlK = 1. / dlK;
  dlT = 1. / dlT;

  return;
}

#define HYPMIN (1e-5)
#define HYPMAX (10000)
#define N_HYP (9000)
#define kappa_min (3.)
#define kappa_max (10.)
#define N_k (70)
#define dkappa (0.1)

double hypergeom[N_k][N_HYP];
void init_emiss_tables_nth(void) {

  int k;
  double result, err, K;
  gsl_function func;
  gsl_integration_workspace *w;
#if(KAPPA)
  func.function = &jnu_integrand_kappa;
  func.params = &K;
#elif(POWERLAW)
  func.function = &jnu_integrand_powerlaw;
  func.params = &K;
#endif
  lK_min = log(KMIN);
  dlK = log(KMAX / KMIN) / (N_ESAMP);

  lT_min = log(TMIN);
  dlT = log(TMAX / TMIN) / (N_ESAMP);

  /*  build table for F(K) where F(K) is given by
     \int_0^\pi ( (K/\sin\theta)^{1/2} + 2^{11/12}(K/\sin\theta)^{1/6})^2
     \exp[-(K/\sin\theta)^{1/3}]
     so that J_{\nu} = const.*F(K)
   */
  w = gsl_integration_workspace_alloc(5000);
  for (k = 0; k <= N_ESAMP; k++) {
    K = exp(k * dlK + lK_min);
    gsl_integration_qag(&func, 0., M_PI / 2., EPSABS, EPSREL, 5000,
                        GSL_INTEG_GAUSS61, w, &result, &err);
 //   gsl_integration_qags(&func, 0.01*M_PI/2., 0.99*M_PI/2. , EPSABS, EPSREL, 10000,
  //                       w, &result, &err);
//	fprintf(stderr,"results %e err %e rel err %e\n",result,err,err/result);
    F_nth[k] = log(4. * M_PI * result);
  }
  gsl_integration_workspace_free(w);

  FILE *input;
  input = fopen("hyper2f1.txt", "r");
  double dummy;
  for (int j = 0; j < N_HYP; j++) {
    for (int i = 0; i < N_k; i++) {
      fscanf(input, "%lf", &dummy);
      hypergeom[i][j] = (dummy);
    }
  }

  /* Avoid doing divisions later */
  dlK = 1. / dlK;
  dlT = 1. / dlT;
  fprintf(stderr, "done reading hypergeom2F1.\n\n");

  return;
}

void init_emiss_tables(void) {
#if (THERMAL)
  init_emiss_tables_th();
#elif (KAPPA || POWERLAW)
  init_emiss_tables_nth();
#elif (MIXED)
  init_emiss_tables_th();
  init_emiss_tables_nth();
#endif
}

/* rapid evaluation of K_2(1/\Thetae) */

double K2_eval(double Thetae) {
  double linear_interp_K2(double);
	double K2;
        if (Thetae < THETAE_MIN)
               return 0.;
        if (Thetae > TMAX)
        	return 2. * Thetae * Thetae;

        return 2. * Thetae * Thetae;;//linear_interp_K2(Thetae);
}

#define KFAC (9 * M_PI * ME * CL / EE)
double F_eval_th(double Thetae, double Bmag, double nu) {

  double K, x;
  double linear_interp_F_th(double);

  K = KFAC * nu / (Bmag * Thetae * Thetae);
  if (K > KMAX) {
    return 0.;
  } else if (K < KMIN) {
    /* use a good approximation */
    x = pow(K, 0.333333333333333333);
    return (x * (37.67503800178 + 2.240274341836 * x));
  } else {
    return linear_interp_F_th(K);
  }
}

double F_eval_kappa(double Thetae, double Bmag, double nu) {

  double K,x;
  double linear_interp_F_nth(double);

  double nuc = EE * Bmag / (2. * M_PI * ME * CL);
  double kappa = kappa_synch;
  double w =(kappa-3.)/kappa *  Thetae;
  double nus = nuc * pow(w * kappa, 2.);

  K = nu / nus;
  if (K > KMAX)
    return 0.;
  if (K < KMIN){
    return (0);
}
  double F_value = linear_interp_F_nth(K) * exp(-nu / nu_cutoff);
  if (isnan(F_value))
    fprintf(stderr, " f_eval %e %e %e %e %e\n", nu, Thetae, nus, Bmag, F_value);
  return F_value;
}

double F_eval_powerlaw(double Thetae, double Bmag, double nu) {

  double K;
  double linear_interp_F_nth(double);
  double nuc = EE * Bmag / (2. * M_PI * ME * CL);

  K = nu / nuc;
  if (K > KMAX)
    return 0.;
  if (K < KMIN)
    return 0.;
  double F_value = linear_interp_F_nth(K) * exp(-nu / nu_cutoff);
  if (isnan(F_value))
    fprintf(stderr, " f_eval %e %e %e %e %e\n", nu, Thetae, nuc, Bmag, F_value);
  return F_value;
}

double F_eval_kappa_exact(double Thetae,double Bmag,double nu){

  double nuc = EE * Bmag / (2. * M_PI * ME * CL);
  double kappa = kappa_synch;
  double wt = (kappa-3.)/kappa * Thetae;
  double nus = nuc * pow(wt * kappa, 2.);
  int k;
  double result, err, K;

  K = nu / nus;
  gsl_function func;
  gsl_integration_workspace *w;

  func.function = &jnu_integrand_kappa;
  func.params = &K;

  w = gsl_integration_workspace_alloc(1500);
  gsl_integration_qag(&func, 0., M_PI / 2., EPSABS, EPSREL, 1500,
                        GSL_INTEG_GAUSS61, w, &result, &err);
  gsl_integration_workspace_free(w);
  return (4. * M_PI * result);

}

double F_eval_nth(double Thetae, double Bmag, double nu) {
 #if(KAPPA)
   return F_eval_kappa(Thetae, Bmag,nu);
 #else
   return F_eval_powerlaw(Thetae, Bmag,nu);
 #endif
}



double hypergeom_eval(double X) {
  double kappa=kappa_synch;
  int i;
  double di,z,hyp2F1;
  double a = 1.;
  double b = -kappa - 1. / 2.;
  double c = 1. / 2.;
/*
     if(fabs(X)<1e-4){
	z=-X;
        hyp2F1=1+a*b*z/c + a*(1+a)*b*(1+b)/(2*c*(1+c))*z*z;
     	return hyp2F1;
     }
     else*/
     if(fabs(X)>1) {
          z = -X;
          hyp2F1 = pow(1.-z, -a) * gsl_sf_gamma(c) * gsl_sf_gamma(b-a)
                                     / (gsl_sf_gamma(b)*gsl_sf_gamma(c-a))
         * gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z))
         + pow(1.-z, -b) * gsl_sf_gamma(c) * gsl_sf_gamma(a-b)
                                     / (gsl_sf_gamma(a) * gsl_sf_gamma(c-b))
         * gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z));
                     return hyp2F1;
             }
     else {
                     //fprintf(stderr,"too small %e %e\n", X,gsl_sf_hyperg_2F1(a,b,c,-X));
                     return gsl_sf_hyperg_2F1(a,b,c,-X);
             }

         /* 
  di = 1000. * log(X / HYPMIN) / log(10);
  i = (int)di;
  double dX = (X - HYPMIN * pow(10, i / 1000.)) /
              (HYPMIN * pow(10, (i + 1) / 1000.) - HYPMIN * pow(10, i / 1000.));
//  if (kappa_synch >= kappa_max) {
//    printf("Kappa=%0.1lf not supported, maximum kappa is %0.1lf\n", kappa_synch,
//           kappa_max - dkappa);
//    exit(1);
//  }
  if (kappa_synch < kappa_min) {
    printf("Kappa=%0.1lf not supported, minimum kappa is %0.1lf\n", kappa_synch,
           kappa_min);
    exit(1);
  }
  int k = (int)((kappa_synch - kappa_min) / dkappa);
  double value = (1 - dX) * hypergeom[k][i] + dX * hypergeom[k][i + 1];

  return value;
*/
	return 0;
}

double F_eval(double Thetae, double Bmag, double nu, int ACCZONE) {

  double F_eval;
#if (THERMAL)
  F_eval = F_eval_th(Thetae, Bmag, nu);
#elif (KAPPA || POWERLAW)
  F_eval = F_eval_nth(Thetae, Bmag, nu);
#elif (MIXED)
  if (ACCZONE) {
    double K2 = K2_eval(Thetae);
    if (K2 == 0.) {
      return 0;
    }

    F_eval = perct_thermal * F_eval_th(Thetae, Bmag, nu) +
             (1 - perct_thermal) * (27. / (M_SQRT2 * 2. * M_PI)) *
                 (K2 / (Thetae * Thetae)) * F_eval_nth(Thetae, Bmag, nu);
  } else
    F_eval = F_eval_th(Thetae, Bmag, nu);

#endif
  //       printf("%e\n",F_eval);
  return F_eval;
}

#undef KFAC
#undef KMIN
#undef KMAX
#undef EPSABS
#undef EPSREL

double linear_interp_K2(double Thetae) {

  int i;
  double di, lT;

  lT = log(Thetae);

  di = (lT - lT_min) * dlT;
  i = (int)di;
  di = di - i;

  return exp((1. - di) * K2[i] + di * K2[i + 1]);
}

double linear_interp_F_th(double K) {

  int i;
  double di, lK;

  lK = log(K);

  di = (lK - lK_min) * dlK;
  i = (int)di;
  di = di - i;

  return exp((1. - di) * F_th[i] + di * F_th[i + 1]);
}

double linear_interp_F_nth(double K) {

  int i;
  double di, lK;

  lK = log(K);

  di = (lK - lK_min) * dlK;
  i = (int)di;
  di = di - i;

  return exp((1. - di) * F_nth[i] + di * F_nth[i + 1]);
}

double linear_interp_F(double K) {
  double F_value = 0;
#if (THERMAL)
  F_value = linear_interp_F_th(K);
#elif (KAPPA || POWERLAW)
  F_value = linear_interp_F_nth(K);

#endif
  return F_value;
}
