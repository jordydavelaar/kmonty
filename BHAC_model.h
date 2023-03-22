
/***********************************************************************************
    Copyright 2013 Joshua C. Dolence, Charles F. Gammie, Monika Mo\'scibrodzka,
                   and Po Kin Leung

                        GRMONTY  version 1.0   (released February 1, 2013)

    This file is part of GRMONTY.  GRMONTY v1.0 is a program that calculates the
    emergent spectrum from a model using a Monte Carlo technique.

    This version of GRMONTY is configured to use input files from the HARM code
    available on the same site.   It assumes that the source is a plasma near a
    black hole described by Kerr-Schild coordinates that radiates via thermal
    synchrotron and inverse compton scattering.

    You are morally obligated to cite the following paper in any
    scientific literature that results from use of any part of GRMONTY:

    Dolence, J.C., Gammie, C.F., Mo\'scibrodzka, M., \& Leung, P.-K. 2009,
        Astrophysical Journal Supplement, 184, 387

    Further, we strongly encourage you to obtain the latest version of
    GRMONTY directly from our distribution website:
    http://rainman.astro.illinois.edu/codelib/

    GRMONTY is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    GRMONTY is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GRMONTY; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

 

/* Global variables and function prototypes for harm(2d) models */

#include "decs.h"

#ifndef global
#define global extern
#endif

#define NSPIN 3

double ****p;


global double stopx[4], startx[4],dx[4];
global double a;
global double th_len;

global double *neqpar;
global double xprobmin[3],xprobmax[3];
global int ng[3],nxlone[3],nleafs;
global double *eqpar;

global int count_leaf,count_node;
global int *forest; //= NULL;
global struct block *block_info;
global int block_size;
global int forest_size;
int *nx, ndimini;
global int n1, n2, n3;

global double Ladv, dMact;

//global double ** Xcoord;
global double *** Xgrid;
//global double *** Xbar;

/* HARM model internal utilities */
void init_weight_table(void);
void bl_coord(double *X, double *r, double *th);
void make_zone_centered_tetrads(void);
void set_units(char *munitstr, char *time);
void init_geometry(void);
void init_bhac3d_data(char *fname);
void init_bhac3d_grid(char *fname);
void init_nint_table(void);
void init_storage(void);
double dOmega_func(double x2i, double x2f,double x3i, double x3f);

void sample_zone_photon(int i, int j, int k, double dnmax, struct of_photon *ph);
double interp_scalar(double ***var, int i, int j, int k, double coeff[4]);
int get_zone(int *i, int *j, int *k, double *dnamx);
//void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]);
void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B, double *sigma, double *beta, double Ucon[NDIM], double Bcon[NDIM], int *ACCZONE, double *dx_local, int *igrid_c);
void calc_coord(int c, int *nx, int ndimini, double *lb, double *dxc_block, double *X);
void gcov_func(double *X, double gcov[][NDIM]);
void gcon_func(double *X, double gcon[][NDIM]);

double get_r(double X[4]);