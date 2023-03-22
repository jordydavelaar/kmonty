#
# requires an openmp-enabled version of gcc
#

CC = mpiicc #/usr/bin/h5pcc
CCFLAGS = -Wall -std=c99 -O3 -fopenmp -inline-level=2 -ipo
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


SRCS1 = main.c compton.c init_geometry.c tetrads.c geodesics.c \
radiation.c jnu_mixed.c hotcross.c track_super_photon.c kappa_sampler.c kappa_sampler_num.c\
scatter_super_photon.c sphere_model.c init_sphere_data.c sphere_utils.c

OBJS1 = main.o compton.o init_geometry.o tetrads.o geodesics.o \
radiation.o jnu_mixed.o hotcross.o track_super_photon.o kappa_sampler.o kappa_sampler_num.o\
scatter_super_photon.o sphere_model.o init_sphere_data.o sphere_utils.o

SRCS2 = main.c compton.c init_geometry.c tetrads.c geodesics.c \
radiation.c jnu_mixed.c hotcross.c track_super_photon.c kappa_sampler.c kappa_sampler_num.c\
scatter_super_photon.c BHAC_model.c BHAC_utils.c init_BHAC_data.c

OBJS2 = main.o compton.o init_geometry.o tetrads.o geodesics.o \
radiation.o jnu_mixed.o hotcross.o track_super_photon.o kappa_sampler.o kappa_sampler_num.o\
scatter_super_photon.o BHAC_model.o BHAC_utils.o init_BHAC_data.o

INCS = decs.h constants.h sphere_model.h

sphere: $(OBJS1) $(INCS) makefile
	$(CC_LOAD) $(CFLAGS) $(OBJS1) $(LDFLAGS) -o $(EXE)

bhac: $(OBJS2) $(INCS) makefile
	$(CC_LOAD) $(CFLAGS) $(OBJS2) $(LDFLAGS) -o $(EXE)


$(EXE) : $(OBJS) $(INCS) makefile
	$(CC_LOAD) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $(EXE)


clean:
	/bin/rm *.o $(EXE)
