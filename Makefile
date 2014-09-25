
CC = g++

INCLUDE_DIRS = -I/home/shriwise/dagmc_blds/lasso/include \
	       -I/home/shriwise/dagmc_blds/cgm_dev/include/ \
               -I/home/shriwise/dagmc_blds/moabs/include \
	       -I/home/shriwise/meshkit/include/

LIBS = -L/home/shriwise/dagmc_blds/moabs/lib/ -lMOAB \
       -L/home/shriwise/meshkit/lib/ -lMeshKit


all: build scripts

build: sweep hv_mesh


sweep:
	$(CC) sweep.cpp -o sweep

hv_mesh:

	$(CC) create_hv_mesh.cpp -o create_hv_mesh $(INCLUDE_DIRS) $(LIBS) 

scripts: sweep

	./sweep
	chmod u+x sweep_param_space

clean: 
	rm -f *.dat sweep create_hv_mesh sweep_param_space cube_mod.h5m *~ *.p