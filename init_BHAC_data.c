#include "BHAC_model.h"
#include "decs.h"
#include <stdio.h>
#include <stdlib.h>
/* HDF5 v1.6 API */
// #include <H5LT.h>

/* HDF5 v1.8 API */

extern double ****p;



void new_index(int dims, int i, int *new_i, int *new_j, int *new_k, int ip,
               int jp, int kp) {

    if (dims == 1) {
        *new_i = 2 * (ip) + i;
        *new_j = 0;
        *new_k = 0;
    }

    if (dims == 2) {
        *new_i = 2 * (ip) + i % 2;
        *new_j = 2 * (jp) + i / 2;
        *new_k = 0;
    }

    if (dims == 3) {
        *new_i = 2 * (ip) + i % 2;
        *new_j = 2 * (jp) + ((int)(i / 2.)) % 2;
        *new_k = 2 * (kp) + (int)(i / 4.);
    }
}

void read_node(FILE *file_id, int *igrid, int *refine, int ndimini, int level,
               int ind_i, int ind_j, int ind_k) {
    int buffer_i[1], leaf;
    fread(buffer_i, sizeof(int), 1, file_id);
    leaf = buffer_i[0];

    if (leaf) {
        (*igrid)++;
        int i = (*igrid); //+(*refine);
        if (i > forest_size) {
            forest = (int *)realloc(forest, i * sizeof(int));
            forest_size++;
        }
        if (i > block_size) {
            block_info =
                (struct block *)realloc(block_info, i * sizeof(struct block));
            block_size++;
        }

        forest[forest_size - 1] = 1;

        block_info[block_size - 1].ind[0] = ind_i;
        block_info[block_size - 1].ind[1] = ind_j;
        block_info[block_size - 1].ind[2] = ind_k;
        block_info[block_size - 1].level = level;
    } else {
        (*refine) += 1;
        int i = (*igrid) + (*refine);

        if (i > forest_size) {
            forest = (int *)realloc(forest, i * sizeof(int) * 2);
            forest_size++;
        }

        forest[forest_size - 1] = 0;

        for (int ich = 0; ich < pow(2, ndimini); ich++) {
            int cind_j, cind_i, cind_k;
            new_index(ndimini, ich, &cind_i, &cind_j, &cind_k, ind_i, ind_j,
                      ind_k);
            read_node(file_id, igrid, refine, ndimini, level + 1, cind_i,
                      cind_j, cind_k);
        }
    }
}

double get_detgamma(double x, double y, double z) {

    double g_dd[4][4];
    double X_u[4];

    X_u[0] = 0;
    X_u[1] = x;
    X_u[2] = y;
    X_u[3] = z;

    double rho = X_u[1] * X_u[1] + X_u[2] * X_u[2] + X_u[3] * X_u[3] - a * a;
    double r = sqrt(0.5 * (rho + sqrt(rho * rho + 4. * z * z * a * a)));
    double sqrtgamma;

    if (r < 1)
        sqrtgamma =
            sqrt(1.0 + (1.4142135623730950 *
                        sqrt(rho + sqrt(rho * rho + 4.0 * a * a * z * z))) /
                           sqrt(rho * rho + 4.0 * a * a * z * z + 1.e-3));
    else
        sqrtgamma =
            sqrt(1.0 + (1.4142135623730950 *
                        sqrt(rho + sqrt(rho * rho + 4.0 * a * a * z * z))) /
                           sqrt(rho * rho + 4.0 * a * a * z * z));

    return sqrtgamma;

    /*
       metric_dd(X_u,g_dd);

       double detgamma = g_dd[1][1] * ( g_dd[2][2]*g_dd[3][3] -
       g_dd[3][2]*g_dd[2][3]) - g_dd[1][2] * ( g_dd[2][1] * g_dd[3][3] -
       g_dd[3][1]*g_dd[2][3]) + g_dd[1][3] * (g_dd[2][1]*g_dd[3][2] - g_dd[2][2]
       * g_dd[3][1]); #if(DEBUG) if(isnan(sqrt(detgamma))){ double R2 =
       x*x+y*y+z*z; double a2 = a * a; double r2 = (R2 - a2 + sqrt((R2 - a2)*(R2
       - a2)+4.*a2*z * z))*0.5; double r_current = sqrt(r2);

       fprintf(stderr,"isnan detgam %e %e rc
       %e\n",sqrt(detgamma),detgamma,r_current);
       detgamma=0;
       exit(1);
       }
       #endif
       return sqrt(detgamma);
     */
}

void calc_coord_bar(double X[4], double *dxc_block, double Xbar[4]) {
    double gcov[4][4];
    double coef_1D[3];
    double coef_2D[3][3];
    double coef_3D[3][3][3];

    double is[3], js[3], ks[3];

    coef_1D[0] = 1.;
    coef_1D[1] = 4.;
    coef_1D[2] = 1.;

    is[0] = -1.;
    is[1] = 0.;
    is[2] = 1.;

    js[0] = -1.;
    js[1] = 0.;
    js[2] = 1.;

    ks[0] = -1.;
    ks[1] = 0.;
    ks[2] = 1.;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            coef_2D[i][j] = coef_1D[i] * coef_1D[j];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                coef_3D[i][j][k] = coef_2D[i][j] * coef_1D[k];
            }
        }
    }

    double norm = 0;
    double xbar = 0;
    double ybar = 0;
    double zbar = 0;
    double detgamma;

    for (int k = 0; k < 3; k++) {
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 3; i++) {
                detgamma = get_detgamma(X[0] + dxc_block[0] / 2. * is[i],
                                        X[1] + dxc_block[1] / 2. * js[j],
                                        X[2] + dxc_block[2] / 2. * ks[k]);
                norm += detgamma * coef_3D[i][j][k];
                xbar += detgamma * coef_3D[i][j][k] *
                        (X[0] + dxc_block[0] / 2. * is[i]);
                ybar += detgamma * coef_3D[i][j][k] *
                        (X[1] + dxc_block[1] / 2. * js[j]);
                zbar += detgamma * coef_3D[i][j][k] *
                        (X[2] + dxc_block[2] / 2. * ks[k]);
                //				fprintf(stderr,"x %e %e
                //%e\n",X[0]+dxc_block[0]/2. * is[i],X[1]+dxc_block[1]/2. *
                // js[j],X[2]+dxc_block[2]/2. * ks[k]);
            }
        }
    }

    norm *= (dxc_block[0] / 6.0) * (dxc_block[1] / 6.0) * (dxc_block[2] / 6.0);
    xbar *= (dxc_block[0] / 6.0) * (dxc_block[1] / 6.0) * (dxc_block[2] / 6.0) /
            norm;
    ybar *= (dxc_block[0] / 6.0) * (dxc_block[1] / 6.0) * (dxc_block[2] / 6.0) /
            norm;
    zbar *= (dxc_block[0] / 6.0) * (dxc_block[1] / 6.0) * (dxc_block[2] / 6.0) /
            norm;

    Xbar[1] = xbar;
    Xbar[2] = ybar;
    Xbar[3] = zbar;
#if (DEBUG)
    if (isnan(Xbar[0])) {
        fprintf(stderr, "isnan in calc_coord_bar %e %e %e %e\n", Xbar[1],
                Xbar[2], Xbar[3], norm);
        Xbar[0] = 0;
        Xbar[1] = 0;
        Xbar[2] = 0;
        exit(1);
    }
#endif
}

void calc_coord(int c, int *nx, int ndimini, double *lb, double *dxc_block,
                double *X) {
    int local_ig[3];

    local_ig[0] = (int)((c % nx[0]));
    local_ig[1] = (int)(fmod((((double)c) / ((double)nx[0])), (double)nx[0]));
    if (ndimini == 3) {
        local_ig[2] = (int)(((double)c) / ((double)nx[0] * nx[1]));
    }
    X[0] = 0;
    for (int i = 0; i < ndimini; i++) {
        X[i + 1] = lb[i] + (local_ig[i] + 0.5) * dxc_block[i]; // cell centered
    }

    if (isnan(X[0])) {
        fprintf(stderr, "isnan in calccoord %lf %d %lf %lf\n", lb[0],
                local_ig[0], dxc_block[0], ((double)c) / (nx[0] * nx[1]));
        exit(1);
    }
}

void convert2prim(double prim[8], double **conserved, int c, double X[4],
                  double X_cell[4], double dxc[3]) {

    double r_current = get_r(X);

    double X_u[4];
    double g_dd[4][4], g_uu[4][4];
    // X_u[0] = 0;
    // X_u[1] = Xgrid[0]; // - dxc[0];
    // X_u[2] = Xgrid[1]; // - dxc[0];
    // X_u[3] = Xgrid[2]; // - dxc[0];

    gcov_func(X, g_dd);
    gcon_func(X, g_uu);

    double BS = conserved[S1][c] * conserved[B1][c] +
                conserved[S2][c] * conserved[B2][c] +
                conserved[S3][c] * conserved[B3][c];
    double Bsq = 0;

    double B_d[4];
    double S_u[4];
    S_u[1] = 0;
    S_u[2] = 0;
    S_u[3] = 0;
    B_d[1] = 0;
    B_d[2] = 0;
    B_d[3] = 0;

    double gamma[4][4];
    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
            gamma[i][j] = g_uu[i][j] + g_uu[0][i] * g_uu[0][j] / (-g_uu[0][0]);
        }
    }

    for (int j = 1; j < 4; j++) {
        for (int i = 1; i < 4; i++) {
            S_u[j] += gamma[i][j] * conserved[S1 + i - 1][c];
            B_d[j] += g_dd[i][j] * conserved[B1 + i - 1][c];
        }
    }

    for (int i = 1; i < 4; i++) {
        //		for(int j=1;j<4;j++){
        Bsq += B_d[i] * conserved[B1 + i - 1][c];
    }
    //	}

#if (DEBUG)
    if (isnan(BS) || isnan(Bsq)) {
        fprintf(stderr, "Bsq %e BS %e\n", Bsq, BS);
        fprintf(stderr, "B %e %e %e\n", conserved[B1][c], conserved[B2][c],
                conserved[B3][c]);
        fprintf(stderr, "V %e %e %e\n", conserved[S1][c], conserved[S2][c],
                conserved[S3][c]);
        fprintf(stderr, "rc %e %e X %e %e %e\n", r_current, cutoff_inner, X[0],
                X[1], X[2]);
        LOOP_ij fprintf(stderr, "gij %d %d %e\n", i, j, g_dd[i][j]);
        exit(1);
    }
#endif
    prim[KRHO] = conserved[D][c] / conserved[LFAC][c];
    prim[UU] = (neqpar[0] - 1.) / neqpar[0] *
               (conserved[XI][c] / pow(conserved[LFAC][c], 2.) - prim[KRHO]) /
               (neqpar[0] -
                1.); // need internal energy not pressure so extra 1/(gam-1)

    prim[U1] =
        S_u[1] / (conserved[XI][c] + Bsq) +
        conserved[B1][c] * BS / (conserved[XI][c] * (conserved[XI][c] + Bsq));
    prim[U2] =
        S_u[2] / (conserved[XI][c] + Bsq) +
        conserved[B2][c] * BS / (conserved[XI][c] * (conserved[XI][c] + Bsq));
    prim[U3] =
        S_u[3] / (conserved[XI][c] + Bsq) +
        conserved[B3][c] * BS / (conserved[XI][c] * (conserved[XI][c] + Bsq));

    prim[B1] = conserved[B1][c];
    prim[B2] = conserved[B2][c];
    prim[B3] = conserved[B3][c];
    // #if (DEBUG)
    double VdotV = 0;
    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
            VdotV += g_dd[i][j] * prim[U1 + i - 1] * prim[U1 + j - 1];
        }
    }
    if (VdotV > 1.) {
        fprintf(stderr, "VdotV too large in con2prim %e %e %e %e\n", VdotV,
                X[1], X[2], X[3]);
        exit(1);
    }

    double lor = 1 / sqrt(1 - VdotV);

    if (isnan(lor)) {
        fprintf(stderr, "\n");

        fprintf(stderr, "lor %e vdotv %e lfac %e XI %e Bsq %e BS %e\n", lor,
                VdotV, conserved[LFAC][c], conserved[XI][c], Bsq, BS);
        fprintf(stderr, "xi? %e %e\n",
                conserved[LFAC][c] * conserved[LFAC][c] * prim[KRHO] *
                    (1 + neqpar[0] * prim[UU] / prim[KRHO]),
                conserved[XI][c]);
        fprintf(stderr, "rc %e %e\n", r_current, (1. + sqrt(1. - a * a)));
        fprintf(stderr, "X %e %e %e\n", X[0], X[1], X[2]);
        fprintf(stderr, "dxc %e %e %e\n", dxc[0], dxc[1], dxc[2]);
        //		exit(1);
    }
    if (isnan(lor))
        exit(1);
    // #endif
}

uint64_t mortonEncode(unsigned int ig1, unsigned int ig2, unsigned int ig3) {
    uint64_t answer = 0;
    for (uint64_t i = 0; i < (sizeof(uint64_t) * 64) / ndimini; ++i) {
        answer = answer | ((ig1 & ((uint64_t)1 << i)) << 2 * i) |
                 ((ig2 & ((uint64_t)1 << i)) << (2 * i + 1)) |
                 ((ig3 & ((uint64_t)1 << i)) << (2 * i + 2));
    }
    return answer;
}

int level_one_Morton_ordered(int ***iglevel1_sfc, int **sfc_iglevel1) {
    // first compute how large the block should be assuming its squared
    int ngsq1 = pow(2, ceil(log10(ng[0]) / log10(2.0)));
    int ngsq2 = pow(2, ceil(log10(ng[1]) / log10(2.0)));
    int ngsq3 = pow(2, ceil(log10(ng[2]) / log10(2.0)));

    int ngsqmax = fmax(ngsq1, fmax(ngsq2, ngsq3));
    ngsq1 = ngsqmax;
    ngsq2 = ngsqmax;
    ngsq3 = ngsqmax;

    int gsq_sfc[ngsq1][ngsq2][ngsq3];
    // construct the sfc for the squared grid
    for (uint32_t i = 0; i < ngsq1; i++) {
        for (uint32_t j = 0; j < ngsq2; j++) {
            for (uint32_t k = 0; k < ngsq3; k++) {
                gsq_sfc[i][j][k] = mortonEncode(i, j, k);
            }
        }
    }

    // delete blocks outside of the real grid
    for (int i = 0; i < ngsq1; i++) {
        for (int j = 0; j < ngsq2; j++) {
            for (int k = 0; k < ngsq3; k++) {
                // check which block runs out of the grid
                if (i >= ng[0] || j >= ng[1] || k >= ng[2]) {
                    // if an index is too large, then we need to decrease all
                    // the sfc indices by 1 if they are larger than the sfc
                    // index of that particular block.
                    for (int ic = 0; ic < ngsq1; ic++) {
                        for (int jc = 0; jc < ngsq2; jc++) {
                            for (int kc = 0; kc < ngsq3; kc++) {
                                if (gsq_sfc[ic][jc][kc] > gsq_sfc[i][j][k])
                                    gsq_sfc[ic][jc][kc]--;
                            }
                        }
                    }
                }
            }
        }
    }

    // create maps to go from sfc to normal grid indices
    for (int i = 0; i < ng[0]; i++) {
        for (int j = 0; j < ng[1]; j++) {
            for (int k = 0; k < ng[2]; k++) {
                iglevel1_sfc[i][j][k] = gsq_sfc[i][j][k];
                sfc_iglevel1[gsq_sfc[i][j][k]][0] = i;
                sfc_iglevel1[gsq_sfc[i][j][k]][1] = j;
                sfc_iglevel1[gsq_sfc[i][j][k]][2] = k;
            }
        }
    }
    return 0;
}

void init_bhac_amr_data(char *fname) {

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0)
        fprintf(stderr, "\nReading HEADER...\n");

    ng[0] = 1;
    ng[1] = 1;
    ng[2] = 1;
    nxlone[0] = 96;
    nxlone[1] = 96;
    nxlone[2] = 192;

    xprobmax[0] = 500;
    xprobmax[1] = 500.;
    xprobmax[2] = 1000.;

    xprobmin[0] = -500.;
    xprobmin[1] = -500.;
    xprobmin[2] = -1000.;

    startx[1] = -500;
    stopx[1] = 500;
    startx[2] = -500;
    stopx[2] = 500;
    startx[3] = -1000;
    stopx[3] = 1000;

    a = 0.9375;

    double buffer[1];
    unsigned int buffer_i[1];
    FILE *file_id;

    file_id = fopen(fname, "rb"); // r for read, b for binary

    if (file_id < 0) {
        if (world_rank == 0)
            fprintf(stderr, "file %s does not exist, aborting...\n", fname);
        fflush(stderr);
        exit(1234);
    }
    int i = 1;

    long int offset;
    int j = 0;
    int levmaxini, ndirini, nwini, nws, neqparini, it, t;
    fseek(file_id, 0, SEEK_END);
    offset = -36 - 4; // -4; // 7 int, 1 double = 7 * 4 + 1*8 = 56?
    fseek(file_id, offset, SEEK_CUR);

    fread(buffer_i, sizeof(int), 1, file_id);
    nleafs = buffer_i[0];

    fread(buffer_i, sizeof(int), 1, file_id);
    levmaxini = buffer_i[0];
    fprintf(stderr, "%d\n", levmaxini);
    fread(buffer_i, sizeof(int), 1, file_id);
    ndimini = buffer_i[0];

    fread(buffer_i, sizeof(int), 1, file_id);
    ndirini = buffer_i[0];

    fread(buffer_i, sizeof(int), 1, file_id);
    nwini = buffer_i[0];

    fread(buffer_i, sizeof(int), 1, file_id);
    nws = buffer_i[0];

    fread(buffer_i, sizeof(int), 1, file_id);
    neqparini = buffer_i[0];

    fread(buffer_i, sizeof(int), 1, file_id);
    it = buffer_i[0];

    fread(buffer, sizeof(double), 1, file_id);
    t = buffer[0];

    offset = offset - (ndimini * 4 + neqparini * 8);
    fseek(file_id, offset, SEEK_CUR);

    neqpar = (double *)malloc(neqparini * sizeof(double));
    nx = (int *)malloc(ndimini * sizeof(int));

    for (int k = 0; k < ndimini; k++) {
        fread(buffer_i, sizeof(int), 1, file_id);
        nx[k] = buffer_i[0];
    }
    for (int k = 0; k < neqparini; k++) {
        fread(buffer, sizeof(double), 1, file_id);
        neqpar[k] = buffer[0];
    }
    // a= neqpar[3];

    int cells = 1;
    for (int k = 0; k < ndimini; k++) {
        cells *= nx[k];
    }

    long int size_block = cells * (nwini)*8; // size of a block in bytes

    // Read forest
    offset = nleafs * size_block +
             nleafs * (nx[0] + 1) * (nx[1] + 1) * (nx[2] + 1) * nws * 8;
    // exit(1);
    fseek(file_id, 0, SEEK_SET);
    fseek(file_id, offset, SEEK_CUR);
    int leaf;

    for (int i = 0; i < ndimini; i++) {
        ng[i] = nxlone[i] / nx[i]; // number of blocks in each direction
    }

    int igrid = 0;
    int refine = 0;

    int size_f = 0;
    block_info = (struct block *)malloc(0);

    forest = (int *)malloc(sizeof(int));




    int level = 1;
#if (SFC)
    int ***iglevel1_sfc;
    int **sfc_iglevel1;

    iglevel1_sfc = (int ***)malloc(ng[0] * sizeof(int **));
    for (int i = 0; i < ng[0]; i++) {
        iglevel1_sfc[i] = (int **)malloc(ng[1] * sizeof(int *));
        for (int j = 0; j < ng[1]; j++) {
            iglevel1_sfc[i][j] = (int *)malloc(ng[2] * sizeof(int));
        }
    }

    sfc_iglevel1 = (int **)malloc(ng[0] * ng[1] * ng[2] * sizeof(int *));
    for (int i = 0; i < ng[0] * ng[1] * ng[2]; i++) {
        sfc_iglevel1[i] = (int *)malloc(3 * sizeof(int));
    }

    level_one_Morton_ordered(iglevel1_sfc, sfc_iglevel1);
    int i, j, k;
    for (int sfc_i = 0; sfc_i < ng[2] * ng[1] * ng[0]; sfc_i++) {
        i = sfc_iglevel1[sfc_i][0];
        j = sfc_iglevel1[sfc_i][1];
        k = sfc_iglevel1[sfc_i][2];

        read_node(file_id, &igrid, &refine, ndimini, level, i, j, k);
    }

#else
    for (int k = 0; k < ng[2]; k++) {
        for (int j = 0; j < ng[1]; j++) {
            for (int i = 0; i < ng[0]; i++) {

                read_node(file_id, &igrid, &refine, ndimini, level, i, j, k);
            }
        }
    }
#endif

    double *dx1, *dxc;
    dx1 = (double *)malloc(ndimini * sizeof(double));
    dxc = (double *)malloc(ndimini * sizeof(double));
    // block and cell size on level one.

    for (int i = 0; i < ndimini; i++) {
        dx1[i] = (xprobmax[i] - xprobmin[i]) / ng[i];
        dxc[i] = (xprobmax[i] - xprobmin[i]) / nxlone[i];
    }

    double **values;
    values = (double **)malloc(nwini * sizeof(double *));
    for (int i = 0; i < nwini; i++) {
        values[i] = (double *)malloc(cells * sizeof(double));
    }

    int h_i = 0, h_j = 0, h_k = 0;
    double del[3];

    long int count = 0;
    double Xcent[4];
    double Xbar[4];

    init_storage();

    fseek(file_id, 0, SEEK_SET);
    if (world_rank == 0)
        fprintf(stderr, "\nReading BODY...\n");

    fprintf(stderr, "im here!\n");
    for (int i = 0; i < nleafs; i++) {
        for (int n = 0; n < ndimini; n++) {
            block_info[i].lb[n] =
                (xprobmin[n] + (block_info[i].ind[n]) * dx1[n] /
                                   pow(2., (double)block_info[i].level - 1.));
#if (DEBUG)
            if (isnan(block_info[i].lb[n])) {
                fprintf(stderr, "NaN %d %lf %d", i,
                        pow(2, block_info[i].level - 1), block_info[i].level);
                exit(1);
            }
#endif
            block_info[i].dxc_block[n] =
                (dxc[n] / (pow(2., (double)block_info[i].level - 1.)));
            block_info[i].size[n] = nx[n];
        }

        for (int nw = 0; nw < nwini; nw++) {
            for (int c = 0; c < cells; c++) {
                fread(buffer, sizeof(double), 1, file_id);
                values[nw][c] = buffer[0];
            }
        }

#pragma omp parallel for shared(values, p) schedule(static, 1) private(Xcent)
        for (int c = 0; c < cells; c++) {
            calc_coord(c, nx, ndimini, block_info[i].lb,
                       block_info[i].dxc_block, Xcent);
            // calc_coord_bar(Xcent,block_info[i].dxc_block,Xbar);

#if (DEBUG)
            if (isnan(Xcent[0])) {
                fprintf(stderr, "%d %d", c, i);
                exit(1);
            }
#endif
            double prim[8];
            convert2prim(prim, values, c, Xcent, Xcent,
                         block_info[i].dxc_block);

            p[KRHO][i][c][0] = prim[KRHO];
            p[UU][i][c][0] = prim[UU];

            p[U1][i][c][0] = prim[U1];
            p[U2][i][c][0] = prim[U2];
            p[U3][i][c][0] = prim[U3];

            p[B1][i][c][0] = prim[B1];
            p[B2][i][c][0] = prim[B2];
            p[B3][i][c][0] = prim[B3];
        }

        offset = (nx[0] + 1) * (nx[1] + 1) * (nx[2] + 1) * nws * 8;
        fseek(file_id, offset, SEEK_CUR);
    }

    free(values);
    free(forest);
    // exit(1);
}
