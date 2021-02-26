#
# requires an openmp-enabled version of gcc
#

CC = mpicc #/usr/bin/h5pcc
CCFLAGS = -Wall -std=c99  -fopenmp -g -O0
LDFLAGS = -lm -lgsl -lgslcblas

#CC = gcc
#CFLAGS = -Wall -Ofast -fopenmp
#LDFLAGS = -lm -lgsl -lgslcblas -fopenmp

CC_COMPILE = $(CC) $(CCFLAGS) -c
CC_LOAD = $(CC) $(CCFLAGS)

.c.o:
	$(CC_COMPILE) $*.c

EXE = kmonty
all: $(EXE)

SRCS2 = kmonty_mpi.c compton.c init_geometry.c tetrads.c geodesics.c \
radiation.c jnu_mixed.c hotcross.c track_super_photon.c kappa_sampler.c\
scatter_super_photon.c sphere_model.c init_sphere_data.c sphere_utils.c

OBJS2 = kmonty_mpi.o compton.o init_geometry.o tetrads.o geodesics.o \
radiation.o jnu_mixed.o hotcross.o track_super_photon.o kappa_sampler.o\
scatter_super_photon.o sphere_model.o init_sphere_data.o sphere_utils.o

SRCS3 = kmonty_mpi.c compton.c init_geometry.c tetrads.c geodesics.c \
radiation.c jnu_mixed.c hotcross.c track_super_photon.c kappa_sampler.c\
scatter_super_photon.c BHAC3D_model.c BHAC3D_utils.c init_BHAC3D_data.c

OBJS3 = kmonty_mpi.o compton.o init_geometry.o tetrads.o geodesics.o \
radiation.o jnu_mixed.o hotcross.o track_super_photon.o kappa_sampler.o\
scatter_super_photon.o BHAC3D_model.o BHAC3D_utils.o init_BHAC3D_data.o

INCS = decs.h constants.h sphere_model.h

sphere: $(OBJS2) $(INCS) makefile
	$(CC_LOAD) $(CFLAGS) $(OBJS2) $(LDFLAGS) -o $(EXE)

bhac: $(OBJS3) $(INCS) makefile
	$(CC_LOAD) $(CFLAGS) $(OBJS3) $(LDFLAGS) -o $(EXE)

$(EXE) : $(OBJS) $(INCS) makefile
	$(CC_LOAD) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $(EXE)


clean:
	/bin/rm *.o $(EXE)
