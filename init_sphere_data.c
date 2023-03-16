#include "decs.h"
#include "sphere_model.h"
/*

   get HARM simulation data from fname

   checks for consistency of coordinates in data file with
   values of coordinate parameters

   Uses standard HARM data file format

   CFG 1 Sept 07

 */

void init_sphere_data(char *fname)
{
        double dV, V;
        double RHO0;
        int i,j,k;
        int ACCZONE;

        gam = 13./9.;
        Rin = 0.0;
        Rout = 100.;
        startx[1] = 0.;
        startx[2] = 0;
        startx[3] = 0.;
        dx[1] = (Rout - Rin)/N1;
        dx[2] = (M_PI)/N2;
        dx[3] = 2.*M_PI/N3;
        fprintf(stderr,"%e\n",dx[1]);
        stopx[0] = 1.;
        stopx[1] = startx[1]+N1*dx[1];
        stopx[2] = startx[2]+N2*dx[2];
        stopx[3] = startx[3]+N3*dx[3];

        L_unit = GNEWT*MBH/(CL*CL);
        T_unit = L_unit/CL;
        M_unit = 1.;
        Thetae_unit = MP/ME*(gam-1.)*1./(1. + 3.0);

        // Set remaining units and constants
        RHO_unit = M_unit/pow(L_unit,3);
        U_unit = RHO_unit*CL*CL;
        B_unit = CL*sqrt(4.*M_PI*RHO_unit);
        Ne_unit = RHO_unit/(MP + ME);
        max_tau_scatt = (6.*L_unit)*RHO_unit*0.4;
        //Rh = 1. + sqrt(1. - a * a);

	printf("MBH = %e\n",MBH);
        printf("M_unit = %e\n", M_unit);
        printf("Ne_unit = %e\n", Ne_unit);
        printf("RHO_unit = %e\n", RHO_unit);
        printf("L_unit = %e\n", L_unit);
        printf("T_unit = %e\n", T_unit);
        printf("B_unit = %e\n", B_unit);
        printf("Thetae_unit = %e\n", Thetae_unit);

        init_storage();

//        RHO0 = TAUS/(SIGMA_THOMSON*RSPHERE*L_unit*Ne_unit);
//      printf("RHO0 %e\n",RHO0);

double Ne, Bnet;
#if (1)
        V = dMact = Ladv = 0.;
        dV = dx[1]*dx[2]*dx[3];
        for(int i=0;i<N1;i++){
            for(int j=0;j<N2;j++){
               for(int k=0;k<N3;k++){
                 double X[NDIM], r, th;
                 coord(i, j, X);
                 r = X[1];
                 th = X[2];
               	 Ne = TAUS/(SIGMA_THOMSON*RSPHERE*L_unit);
                 p[KRHO][i][j][k] = Ne/Ne_unit;// *exp(-pow(X[1]/RSPHERE,2));
                 p[UU][i][j][k] = THETAE*p[KRHO][i][j][k]/Thetae_unit;
                 p[U1][i][j][k] = 0.;
                 p[U2][i][j][k] = 0.;
                 p[U3][i][j][k] = 0.;
                 double RHO0 = Ne/Ne_unit;
                 double UU0 = THETAE*RHO0/Thetae_unit;
                 Bnet = sqrt(2.*(gam-1.)*UU0/(BETA));
                 p[B1][i][j][k] = Bnet*cos(th);
                 p[B2][i][j][k] = -Bnet*sin(th)/(r+1e-8);
                 p[B3][i][j][k] = 0.;

	//exit(1);
//        if(j==0){
  //              fprintf(stderr,"r %e theta %e i %d\n",r,th,i);
    //            fprintf(stderr,"rho %e Ne %e Ne_unit %e U %e\n",p[KRHO][i][j][k],Ne,Ne_unit,p[UU][i][j][k]);
      //          fprintf(stderr,"B1 %e B2 %e\n",p[B1][i][j][k],p[B2][i][j][k]);
        //}
              }
            }
          }
        bias_norm /= V;
	fprintf(stderr,"Ne %e Bnet %e Te %e\n",Ne,Bnet*B_unit,THETAE);
#else
        /*      V = dMact = Ladv = 0.;
              dV = dx[1]*dx[2]*dx[3];
              ZLOOP {
                      double X[NDIM];
                      coord(i, j, X);
                      if (X[1] < RSPHERE) {
                              p[KRHO][i][j] = RHO0*exp(-pow(X[1]/RSPHERE,2)); //1.;
                              p[UU][i][j] = THETAE * p[KRHO][i][j]/Thetae_unit;
                              p[U1][i][j] = 0.;
                              p[U2][i][j] = 0.;
                              p[U3][i][j] = 0.;
                              p[B1][i][j] = 0.;
                              p[B2][i][j] = 0.;
                              p[B3][i][j] = 0.;

                      } else {
                              p[KRHO][i][j] = 0.;
                              p[UU][i][j] = 0.;
                              p[U1][i][j] = 0.;
                              p[U2][i][j] = 0.;
                              p[U3][i][j] = 0.;
                              p[B1][i][j]= 0.;
                              p[B2][i][j] = 0.;
                              p[B3][i][j] = 0.;
                      }

                      double Ne = TAUS/(SIGMA_THOMSON*RSPHERE*L_unit);
                      p[KRHO][i][j] = Ne/Ne_unit*exp(-pow(X[1]/RSPHERE,2));
                      p[UU][i][j] = THETAE*p[KRHO][i][j]/Thetae_unit;

                      //                                    printf("%d %d %e %e %e %e\n",i,j,X[1],RSPHERE,p[UU][i][j],Thetae_unit);
                      V += dV*geom[i][j].g;
                      bias_norm += dV*geom[i][j].g*pow(p[UU][i][j]/p[KRHO][i][j]*Thetae_unit,2.);

                      if (i <= 20) {M_PI/2.
                                    double Ne, Thetae, Bmag, Ucon[NDIM], Ucov[NDIM], Bcon[NDIM];
                                    get_fluid_zone(i, j, &Ne, &Thetae, &Bmag, Ucon, Bcon,&ACCZONE);
                                    lower(Ucon, geom[i][j].gcov, Ucov);
                                    dMact += geom[i][j].g*dx[2]*dx[3]*p[KRHO][i][j]*Ucon[1];
                                    Ladv += geom[i][j].g*dx[2]*dx[3]*p[UU][i][j]*Ucon[1]*Ucov[0]; }
              }

              dMact /= 21.;
              Ladv /= 21.;
              bias_norm /= V;
              fprintf(stderr, "dMact: %g, Ladv: %g\n", dMact, Ladv);


              // Add a uniform vertical magnetic field
              ZLOOP {
                      double X[NDIM], r, th;
                      coord(i, j, X);
                      r = X[1];
                      th = X[2];
                      double Ne0 = TAUS/(SIGMA_THOMSON*RSPHERE*L_unit);
                      double RHO0 = Ne0/Ne_unit;
                      double UU0 = THETAE*RHO0/Thetae_unit;
                      //double Bnet = sqrt(2.*(gam-1.)*p[UU][i][j]/(BETA));
                      double Bnet = sqrt(2.*(gam-1.)*UU0/(BETA));
                      p[B1][i][j] = Bnet*cos(th);
                      p[B2][i][j] = -Bnet*sin(th)/r;
                      p[B3][i][j] = 0.;
              }*/
#endif
}
