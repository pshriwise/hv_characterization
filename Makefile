
CC = g++

INCLUDE_DIRS = -I/home/shriwise/dagmc_blds/lasso/include \
	       -I/home/shriwise/dagmc_blds/cgm_dev/include/ \
               -I/home/shriwise/dagmc_blds/moabs/include \
	       -I/home/shriwise/meshkit/include/

LIBS = -L/home/shriwise/dagmc_blds/moabs/lib/ -lMOAB -ldagmc \
       -L/home/shriwise/meshkit/lib/ -lMeshKit


all: build

build: sweep write_obbs hv_cube

sweep: hexmaker.o hv_mesh_gen.o
	$(CC) sweep.cpp hv_mesh_gen.o hexmaker.o -o sweep $(INCLUDE_DIRS) $(LIBS)

write_obbs: hexmaker.o hv_mesh_gen.o
	$(CC) write_obbs.cpp hv_mesh_gen.o hexmaker.o -o write_obbs $(INCLUDE_DIRS) $(LIBS)

hv_cube: hv_mesh_gen.o hexmaker.o
	$(CC) hv_cube.cpp hv_mesh_gen.o hexmaker.o -o hv_cube $(INCLUDE_DIRS) $(LIBS)

hexmaker.o:
	$(CC) hexmaker.cpp $(INCLUDE_DIRS) $(LIBS) -c -o hexmaker.o

hv_mesh_gen.o: hexmaker.o
	$(CC) hv_mesh_gen.cpp $(INCLUDE_DIRS) $(LIBS) -c

plot_datafile:
	paste params.dat data.dat > all.dat

clean: 
	mv cube.h5m cube.saf
	rm -f sweep write_obbs create_hv_mesh hv_cube .#* *# *.h5m *.vtk *.stl *~ *.p *.o
	mv cube.saf cube.h5m
