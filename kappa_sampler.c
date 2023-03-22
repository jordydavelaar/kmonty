//
//  main.c
//  kappasampler
//
//  Created by Jordy on 13/07/2017.
//  Copyright © 2017 Jordy. All rights reserved.
//

#include "decs.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
// rewrite this

struct f_params {
    double u;
};

double find_y(double u, double (*df)(double), double (*f)(double, void *),
              double w) {
    // Get maximum for window
    int status, steps = 0, max_steps = 1000;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double solution;
    double low = 1e-6; // sqrt(0.001/w);
    double high = 1e6; // sqrt(1e4/w);
    gsl_function F;
    struct f_params params = {u};
    F.function = f;
    F.params = &params;
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, low, high);
    // printf("trying to find %e\n",u);
    do {
        steps++;
        status = gsl_root_fsolver_iterate(s);
        solution = gsl_root_fsolver_root(s);
        low = gsl_root_fsolver_x_lower(s);
        high = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(low, high, 1e-9, 0);
    } while (status == GSL_CONTINUE && steps < max_steps);
    // fprintf(stderr,"steps %d solution %e\n",steps,solution);
    if (steps == max_steps) {
        fprintf(stderr, "max steps with solution %e\n", f(solution, &params));
    }

    return solution;
}

double dF3(double y) {
    double value, denom, num;
    double y2 = y * y;
    double kappa = kappa_synch;

    num = 4. * y2 * pow((kappa + y2) / (kappa), -kappa - 1.) *
          gsl_sf_gamma(kappa);
    denom = sqrt(M_PI) * sqrt(kappa) * gsl_sf_gamma(kappa - 1. / 2.);
    value = num / denom;

    return value;
}

double dF4(double y) {
    double value, denom, num;
    double y2 = y * y;
    double kappa = kappa_synch;

    num = 2. * (kappa - 1) * y2 * y * pow((kappa + y2) / (kappa), -kappa - 1.);
    denom = kappa;
    value = num / denom;

    return value;
}

double dF5(double y) {
    double value, denom, num;
    double y2 = y * y;
    double kappa = kappa_synch;

    num = 8 * y2 * y2 * pow((kappa + y2) / (kappa), -kappa - 1) *
          gsl_sf_gamma(kappa);
    denom =
        3 * sqrt(M_PI) * pow(kappa, 3. / 2.) * gsl_sf_gamma(kappa - 3. / 2.);
    value = num / denom;

    return value;
}

double dF6(double y) {
    double value, denom, num;
    double y2 = y * y;
    double kappa = kappa_synch;

    num = (kappa * kappa - 3. * kappa + 2) * pow(y, 5.) *
          pow((kappa + y2) / (kappa), -kappa - 1);
    denom = kappa * kappa;
    value = num / denom;
    return value;
}

double F3(double y, void *params) {
    struct f_params *p = (struct f_params *)params;
    double u = p->u;
    double value, denom, num, hyp2F1;
    double kappa = kappa_synch;

    double y2 = y * y;
    double z = -y2 / kappa;

    hyp2F1 = hypergeom_eval(-z);

    num = -sqrt(kappa) * pow(((y2 + kappa) / kappa), -kappa) *
          gsl_sf_gamma(kappa) *
          (-kappa * hyp2F1 + y2 * (2 * kappa + 1) + kappa);
    denom = y * sqrt(M_PI) * gsl_sf_gamma(3. / 2. + kappa);

    value = num / denom - u;

    return value;
}

double F4(double y, void *params) {
    struct f_params *p = (struct f_params *)params;
    double u = p->u;
    double value, denom, num;
    double y2 = y * y;
    double kappa = kappa_synch;

    num = 1 - (1 + y2) * pow((kappa + y2) / (kappa), -kappa);
    denom = 1;

    value = num / denom - u;

    return value;
}

double F5(double y, void *params) {
    struct f_params *p = (struct f_params *)params;
    double u = p->u;

    double value, denom, num, hyp2F1;
    double kappa = kappa_synch;

    double y2 = y * y;
    double z = -y2 / kappa;

    hyp2F1 = hypergeom_eval(-z);

    num = pow((y2 + kappa) / kappa, -kappa) * gsl_sf_gamma(kappa) *
          (3 * kappa * kappa * (hyp2F1 - 1) +
           (1. - 4. * kappa * kappa) * y2 * y2 -
           3. * kappa * (2. * kappa + 1.) * y2);
    denom = 3. * pow(kappa, 1. / 2.) * y * sqrt(M_PI) *
            gsl_sf_gamma(3. / 2. + kappa);
    value = num / denom - u;

    return value;
}

double F6(double y, void *params) {
    struct f_params *p = (struct f_params *)params;
    double u = p->u;

    double value, denom, num;
    double y2 = y * y;
    double y4 = y2 * y2;
    double kappa = kappa_synch;

    num = (y4 - (y4 + 2. * y2 + 2.) * kappa) *
          pow((kappa + y2) / (kappa), -kappa);
    denom = 2 * kappa;
    value = num / denom + 1 - u;

    return value;
}

double sample_y_distr_nth(double Thetae) {
    double w = (kappa_synch - 3.) / kappa_synch * Thetae;
    double S_3, pi_3, pi_4, pi_5, pi_6, y = -1, x1, x2, prob;
    double num, den;
    double kappa = kappa_synch;
    yhigh = sqrt((10 * gamma_max - 1) / w);
    // gsl_set_error_handler_off();

    pi_3 = sqrt(kappa) * sqrt(M_PI) * gsl_sf_gamma(-1. / 2. + kappa) /
           (4. * gsl_sf_gamma(kappa));
    pi_4 = kappa / (2. * kappa - 2.) * sqrt(0.5 * w);
    pi_5 = 3. * pow(kappa, 3. / 2.) * sqrt(M_PI) *
           gsl_sf_gamma(-3. / 2. + kappa) / (8. * gsl_sf_gamma(kappa)) * w;
    pi_6 =
        kappa * kappa / (2. - 3. * kappa + kappa * kappa) * w * sqrt(0.5 * w);

    S_3 = pi_3 + pi_4 + pi_5 + pi_6;

    pi_3 /= S_3;
    pi_4 /= S_3;
    pi_5 /= S_3;
    pi_6 /= S_3;
    //    printf("%e\n",pi_3);
    //    printf("%e\n",pi_3+pi_4);
    //    printf("%e\n",pi_3+pi_4+pi_5);
    //    printf("%e\n",pi_3+pi_4+pi_5+pi_6);

    do {
        x1 = monty_rand();
        double u = monty_rand();
        if (x1 < pi_3) {
            y = find_y(u, dF3, F3, w);
        } else if (x1 < pi_3 + pi_4) {
            y = find_y(u, dF4, F4, w);
        } else if (x1 < pi_3 + pi_4 + pi_5) {
            y = find_y(u, dF5, F5, w);
        } else {
            y = find_y(u, dF6, F6, w);
        }

        x2 = monty_rand();
        num = sqrt(1. + 0.5 * w * y * y);
        den = (1. + y * sqrt(0.5 * w));
        prob = (num / den) * exp(-(w * y * y) / gamma_max); //* w*sqrt(2*w)/S_3;
        if (y != y)
            prob = 0;
        //     printf("nth %e %e prob %e
        //     %e\n",(1+w*y*y)/gamma_max,(1+w*y*y),prob,x2 );

    } while (x2 >= prob);

    return (y);
}
