# kmonty

Public version of kmonty. A General Relativistic Monte Carlo Radation Transport code. The code is currently coupled to the non-uniform data format from the GRMHD code BHAC. kmonty can be used to compute spectra based based on synchrtron emission and Compton upscattering.

When using this code please cite

Dolence et al. 2009 https://arxiv.org/abs/0909.0708

Davelaar et al. submitted 

# Instalation

kmonty uses MPI and gsl libraries.

kmonty is compliled with a makefile.

To compile for BHAC use 

```
make bhac
```

#Usage

To run kmonty use the following command

```
mpirun -n NCORES ./kmonty NPHOTONS FILENAME MSCALING INDEX
```

Where ``` NCORES ``` should be set to the amount of MPI processes. The other run time parameters are

``` NPHOTONS ``` amount of initial super photons

``` FILENAME ``` GRHMHD filename and path

``` MSCALING ``` Mass accretion rate scaling in grams

``` INDEX ``` Additional index for spectrum output files

Additional compilations variables can be found in decs.h, these are, but not limited to;

Space filling curve, needed for some of the BHAC simulations

```#define SFC 1```

Metric type

```
#define MKS 1
#define CKS 0
```

Type of distrubtion function

```
#define THERMAL 1
#define KAPPA 0
#define POWERLAW 0 
```

Frequency range of initial seed photons

```
#define NUMIN 1.e8
#define NUMAX 1.e24
```

Detector properies, amount of energy of angle bins.

```
#define N_ESAMP 200
#define N_EBINS 200
#define N_THBINS 3
#define N_PHIBINS 1
```


# Tutorial

In this short tutorial we will generate an image from a BHAC grmhd model.

Compile the code with 

```
make bhac
```

Download a GRMHD BHAC snapshot

``` wget https://astro.ru.nl/~jordyd/data2000.dat ```

Run the code with

```
mpirun -n 1 ./kmonty 1e5 data2000.dat 1e18 1
```

This should generate the output file ```spec_1.dat``` that contains colums with frequencies, flux and luminosity for 3 different inclination bins in the following way

``` frequencies flux_1 luminosity_1 flux_2 luminosity_2 flux_3 luminosity_3 ```

# Upcoming featerus

1. Restructuring from code in the style of RAPTOR https://github.com/jordydavelaar/raptor
2. Plotting library
3. Polarization
4. GPU acceleration
