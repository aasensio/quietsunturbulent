# Makefile for building the MCMC F90 code BayesME

# Select compiler (gnu/intel)

C_COMPILER = gnu
F_COMPILER = intel

include makefiles/makefile.c$(C_COMPILER)
include makefiles/makefile.f$(F_COMPILER)

HMC_LIB = hmc

LIBS=-L$(HMC_LIB) -lhmc 
# -L$(GEODESIC_LIB) -lgeodesic -L$(BLAS_LIB) -L$(LAPACK_LIB) -llapack_LINUX -lblas_LINUX

.SUFFIXES: .f .o

.f.o:
	$(COMPILER_F90) $(OPTIONS_F90) $< -o $@
	
.SUFFIXES: .f90 .o

.f90.o:
	$(COMPILER_F90) $(OPTIONS_F90) $< -o $@

.SUFFIXES: .c .o

.c.o:
	$(COMPILER_C) $(CFLAGS) -c $< -o $@

OBJECTS_HMC = hmc/guided_hmc.o hmc/mt19937.o hmc/hanson.o 
OBJECTS_GEODESIC = geodesic/Accept.o  geodesic/converge.o  geodesic/destsv.o  geodesic/dgqt.o  geodesic/dpmpar.o  geodesic/fdavv.o \
        geodesic/fdjac.o  geodesic/lambda.o  geodesic/leastsq.o  geodesic/UpdateJac.o
FINAL_EXECUTABLE = turbulentHIERSingle
NVARIABLES = 200

HMC_LIB = hmc
GEODESIC_LIB = geodesic
BLAS_LIB = /scratch/LIBS/BLAS
LAPACK_LIB = /scratch/LIBS/lapack-3.4.2
LIBS=-L$(HMC_LIB) -lhmc -L$(GEODESIC_LIB) -lgeodesic -L$(BLAS_LIB) -L$(LAPACK_LIB) -llapack -lblas

OBJECTS = vars.o maths.o io.o turbulent.o ad.o sampling.o likelihood.o admaths.o geodesiclm.o

turb: $(OBJECTS)
	$(COMPILER_F90) $(OBJECTS) -o $(FINAL_EXECUTABLE) $(LIBS)
	cp $(FINAL_EXECUTABLE) ../
	
all:
	make hmc
	make geodesic
	make turb
	
hmc: $(OBJECTS_HMC)
	ar cr hmc/libhmc.a $(OBJECTS_HMC)
	
geodesic: $(OBJECTS_GEODESIC)
	ar cr geodesic/libgeodesic.a $(OBJECTS_GEODESIC)
	
clean: 
	find . -maxdepth 2 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete ;
	find . -maxdepth 1 -name "*.f90~" -delete ; find . -maxdepth 2 -name "*.a" -delete ;
	find . -maxdepth 1 -name "$(FINAL_EXECUTABLE)" -delete
	find ../ -maxdepth 2 -name "$(FINAL_EXECUTABLE)" -delete
	
vars.o: vars.f90 ad.o
	$(COMPILER_F90) $(OPTIONS_F90) vars.f90

maths.o: maths.f90 vars.o ad.o
	$(COMPILER_F90) $(OPTIONS_F90) maths.f90
	
io.o: io.f90 vars.o maths.o ad.o
	$(COMPILER_F90) $(OPTIONS_F90) io.f90
		
admaths.o: admaths.f90
	$(COMPILER_F90) $(OPTIONS_F90) admaths.f90
	
ad.o: ad.f90 admaths.o
	$(COMPILER_F90) $(OPTIONS_F90) $(PREPROC_F90) -DdualSize=$(NVARIABLES) ad.f90
	
turbulent.o: turbulent.f90 io.o vars.o likelihood.o sampling.o
	$(COMPILER_F90) $(OPTIONS_F90) turbulent.f90
	
sampling.o: sampling.f90 vars.o maths.o likelihood.o geodesiclm.o
	$(COMPILER_F90) $(OPTIONS_F90) sampling.f90

likelihood.o: likelihood.f90 vars.o maths.o
	$(COMPILER_F90) $(OPTIONS_F90) likelihood.f90
	
geodesiclm.o: geodesiclm.f90 vars.o maths.o
	$(COMPILER_F90) $(OPTIONS_F90) geodesiclm.f90
