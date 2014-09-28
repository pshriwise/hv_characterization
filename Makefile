
CC = g++

INCLUDE_DIRS = -I/home/shriwise/dagmc_blds/lasso/include \
	       -I/home/shriwise/dagmc_blds/cgm_dev/include/ \
               -I/home/shriwise/dagmc_blds/moabs/include \
	       -I/home/shriwise/meshkit/include/

LIBS = -L/home/shriwise/dagmc_blds/moabs/lib/ -lMOAB -ldagmc \
       -L/home/shriwise/meshkit/lib/ -lMeshKit


all: build scripts

build: sweep hv_mesh


sweep:
	$(CC) sweep.cpp -o sweep

hexmaker.o:
	$(CC) hexmaker.cpp $(INCLUDE_DIRS) $(LIBS) -c -o hexmaker.o

hv_mesh: hexmaker.o

	$(CC) create_hv_mesh.cpp hexmaker.o -o create_hv_mesh $(INCLUDE_DIRS) $(LIBS) 

scripts: sweep

	./sweep
	chmod u+x sweep_param_space

clean: 
	mv cube.h5m cube.saf
	rm -f *.dat sweep create_hv_mesh sweep_param_space *.h5m *.vtk *.stl *~ *.p
	mv cube.saf cube.h5m
