
/*

   kmonty, based on grmonty

   Using monte carlo method, estimate spectrum of an appropriately
   scaled axisymmetric GRMHD simulation as a function of
   latitudinal viewing angle.

   Input simulation data is assumed to be in dump format provided by
   HARM code.  Location of input file is, at present, hard coded
   (see init_sim_data.c).

   Nph super-photons are generated in total and then allowed
   to propagate.  They are weighted according to the emissivity.
   The photons are pushed by the geodesic equation.
   Their weight decays according to the local absorption coefficient.
   The photons also scatter with probability related to the local
   scattering opacity.

   The electrons are assumed to have a thermal distribution
   function, and to be at the same temperature as the protons.

   CFG 8-17-06

   Implemented synchrotron sampling, 22 Jan 07

   fixed bugs in tetrad/compton scattering routines, 31 Jan 07

   Implemented normalization for output, 6 Feb 07

   Separated out different synchrotron sampling routines
   into separate files, 8 Mar 07

   fixed bug in energy recording; bug used Kcon[0] rather than
   Kcov[0] as energy, 18 Mar 07

   major reorganization to encapsulate problem-dependent parts 5-6 Nov 07

   added accelerated paticles J. Davelaar 2017-2018

 */

#include "decs.h"

/* defining declarations for global variables */
struct of_geom **geom;
struct of_spectrum ***shared_spect;
double ***shared_Xi_spec;
double ***shared_ispec;

double F_th[N_ESAMP + 1], F_nth[N_ESAMP+1], wgt[N_ESAMP + 1];
long int Ns, N_scatt, N_superph_recorded,N_superph_recorded_total;

/* some coordinate parameters */
double a;
double R0, Rin, Rh, Rout, Rms;
double hslope;
double startx[NDIM], stopx[NDIM], dx[NDIM];

double dlE, lE0;
double gam;
double dMsim;
double M_unit, L_unit, T_unit;
double RHO_unit, U_unit, B_unit, Ne_unit, Thetae_unit;
double max_tau_scatt, Ladv, dMact, bias_norm;

double stopx[4], startx[4],dx[4];
double a;
double th_len;

double *neqpar;
double xprobmin[3],xprobmax[3];
int ng[3],nxlone[3],nleafs;
double *eqpar;

int count_leaf,count_node;
int *forest; //= NULL;
struct block *block_info;
int block_size;
int forest_size;

int n1, n2, n3;
int index;

double Ladv, dMact;

double ** Xcoord;
double *** Xgrid;
double *** Xbar;

gsl_rng *r;
gsl_integration_workspace *w;

#include <time.h>


#pragma omp threadprivate(r)

int main(int argc, char *argv[])
{
        double Ntot, N_superph_made_local,N_superph_made;
        int quit_flag, myid;
        struct of_photon ph;
        time_t currtime, starttime;

        if (argc < 4) {
                fprintf(stderr, "usage: grmonty Ns infilename M_unit\n");
                exit(0);
        }
        sscanf(argv[1], "%lf", &Ntot);
        Ns = (long int) Ntot;
        /* initialize random number generator */
       	sscanf(argv[4], "%d", &index);

#if MPI
        MPI_Init(NULL, NULL);

        // Get the number of processes
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        // Get the rank of the process
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        // Get the name of the processor
        char processor_name[MPI_MAX_PROCESSOR_NAME];
        int name_len;
        MPI_Get_processor_name(processor_name, &name_len);

        // Print off a hello world message

        if(world_rank==0) {
                fprintf(stderr,"\nBOOTING UP MPI-spec\n\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);

        fprintf(stderr,"Initialized processor %s, rank %d out of %d processors\n",
                processor_name, world_rank, world_size);
        fflush(stderr);
        MPI_Barrier(MPI_COMM_WORLD);

        if(world_rank==0) {
                fprintf(stderr,"\nMPI-spec running\n\n");
        }

        myid = world_rank;
        init_monty_rand(139 * myid   + time(NULL)*index); /* Arbitrarily picked initial seed */
#endif

#if OPENMP
#pragma omp parallel private(myid)
        {
                myid = omp_get_thread_num();
                init_monty_rand(139 * myid  + time(NULL)*(1337 + index+myid)); /* Arbitrarily picked initial seed */
                fprintf(stderr,"id is %d\n",myid);
       }
#endif 
        long int Ns_total= Ns;

        //   Ns /=(world_size);

        fflush(stderr);

        /* spectral bin parameters */
        dlE = 0.25; /* bin width */
        lE0 = log(1.e-12); /* location of first bin, in electron rest-mass units */

        /* initialize model data, auxiliary variables */
        init_model(argv);
#if MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        /** main loop **/
        N_superph_made_local = 0;
        N_superph_recorded = 0;
        N_scatt = 0;
        starttime = time(NULL);
        quit_flag = 0;
#if MPI
        if(world_rank==0) {
#endif 
               fprintf(stderr, "\nMAIN LOOP\n");
                fflush(stderr);
#if MPI
        }
#endif

#if OPENMP 
#pragma omp parallel private(ph)
        {
#endif
#if MPI
	{
#endif
                while (1) {

                        /* get pseudo-quanta */
#if OPENMP 
#pragma omp critical (MAKE_SPHOT)
			{
#endif
#if MPI
			{
#endif
                                if (!quit_flag)
                                        make_super_photon(&ph, &quit_flag);
                        }
                        if (quit_flag)
                                break;

                        /* push them around */
                        track_super_photon(&ph);

                        /* step */
#if OPENMP
#pragma omp atomic
#endif
                        N_superph_made_local += 1;

                        /* give interim reports on rates */
#if MPI
                        if (((int) (N_superph_made_local)) % 1000 == 0
                            && N_superph_made_local*world_size > 0 && world_rank==0){
                                currtime = time(NULL);
                                fprintf(stderr, "time %g, rate %g ph/s\n",
                                        (double) (currtime - starttime),
                                        N_superph_made_local*world_size / (currtime -
                                                                           starttime));
                                fflush(stderr);

				}
#endif
                        /* give interim reports on rates */
#if OPENMP 
                       if (((int) (N_superph_made_local)) % 10000 == 0
                            && N_superph_made_local > 0) {
                                currtime = time(NULL);
                                fprintf(stderr, "time %g, rate %g ph/s\n",
                                        (double) (currtime - starttime),
                                        N_superph_made_local / (currtime -
                                                          starttime));
                        }
#endif
                }
        }
        currtime = time(NULL);
#if MPI
        MPI_Barrier(MPI_COMM_WORLD);

        if(world_rank==0) {
                fprintf(stderr,"reducing spetra\n");
                fflush(stderr);
        }
#endif

//#ifdef _OPENMP
//#pragma omp parallel
#if MPI
        {
                mpi_reduce_spect();
        }

        if(world_rank==0) {
                fprintf(stderr,"done\n");
                fflush(stderr);
        }


        MPI_Allreduce(&N_superph_made_local, &N_superph_made,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&N_superph_recorded, &N_superph_recorded_total,1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

        if(world_rank==0) {
                fprintf(stderr, "Final time %g, rate %g ph/s\n",
                        (double) (currtime - starttime),
                        N_superph_made / (currtime - starttime));
                report_spectrum(N_superph_made);
        }
        // Finalize the MPI environment.
        MPI_Finalize();
#endif

#if OPENMP
#pragma omp barrier

#pragma omp parallel
        {
                omp_reduce_spect();
        }
        report_spectrum((int) N_superph_made_local);
#endif

        /* done! */
        return (0);

}
