COMPILER = ifort
FLAGS = 
OPTIONS = -c -O3 -p -g
# -g -traceback
 
FINAL_EXECUTABLE = turbulent
NVARIABLES = 3

OBJECTS = vars.o maths.o io.o partition.o background.o turbulent.o synth.o hydros.o ad.o

turb: $(OBJECTS)
	$(COMPILER) $(OBJECTS) -p -g -o $(FINAL_EXECUTABLE)
	cp $(FINAL_EXECUTABLE) ../
	
clean: 
	find . -maxdepth 2 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete ;
	find . -maxdepth 1 -name "*.f90~" -delete ; find . -maxdepth 2 -name "*.a" -delete ;
	find . -maxdepth 1 -name "hazel" -delete
	find ../ -maxdepth 2 -name "hazel" -delete
	
vars.o: vars.f90
	$(COMPILER) $(OPTIONS)  vars.f90

maths.o: maths.f90 vars.o partition.o ad.o
	$(COMPILER) $(OPTIONS)  maths.f90
	
io.o: io.f90 vars.o maths.o hydros.o partition.o ad.o
	$(COMPILER) $(OPTIONS)  io.f90
	
partition.o: partition.f90
	$(COMPILER) $(OPTIONS)  partition.f90
		
background.o: background.f90 maths.o
	$(COMPILER) $(OPTIONS)  background.f90
	
ad.o: ad.f90 
	$(COMPILER) $(OPTIONS) -fpp -DdualSize=$(NVARIABLES) ad.f90
	
turbulent.o: turbulent.f90 io.o vars.o synth.o
	$(COMPILER) $(OPTIONS)  turbulent.f90
	
hydros.o: hydros.f90 vars.o partition.o maths.o background.o
	$(COMPILER) $(OPTIONS)  hydros.f90
	
synth.o: synth.f90 partition.o vars.o background.o hydros.o
	$(COMPILER) $(OPTIONS)  synth.f90
	