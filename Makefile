
CC = g++

INCLUDE_DIRS = -I/home/shriwise/dagmc_blds/lasso/include \
	       -I/home/shriwise/dagmc_blds/cgm_dev/include/ \
               -I/home/shriwise/dagmc_blds/moabs/include \
	       -I/home/shriwise/meshkit/include/

LIBS = -L/home/shriwise/dagmc_blds/moabs/lib/ -lMOAB -ldagmc \
       -L/home/shriwise/meshkit/lib/ -lMeshKit


all: build

build: sweep write_obbs hv_cube

test_model: obbhexwriter.o hv_mesh_gen.o ProgOptions.o
	$(CC) test_model.cpp hv_mesh_gen.o obbhexwriter.o ProgOptions.o -o test_model $(INCLUDE_DIRS) $(LIBS)

sweep: obbhexwriter.o hv_mesh_gen.o ProgOptions.o
	$(CC) sweep.cpp hv_mesh_gen.o obbhexwriter.o -o sweep $(INCLUDE_DIRS) $(LIBS)

write_obbs: obbhexwriter.o hv_mesh_gen.o ProgOptions.o
	$(CC) write_obbs.cpp hv_mesh_gen.o obbhexwriter.o -o write_obbs $(INCLUDE_DIRS) $(LIBS)

hv_cube: hv_mesh_gen.o obbhexwriter.o
	$(CC) hv_cube.cpp hv_mesh_gen.o obbhexwriter.o -o hv_cube $(INCLUDE_DIRS) $(LIBS)

obbhexwriter.o:
	$(CC) obbhexwriter.cpp $(INCLUDE_DIRS) $(LIBS) -c -o obbhexwriter.o

hv_mesh_gen.o: obbhexwriter.o
	$(CC) hv_mesh_gen.cpp $(INCLUDE_DIRS) $(LIBS) -c

ProgOptions.o:
	$(CC) ProgOptions.cpp -c -o ProgOptions.o

plot_datafile:
	paste params.dat data.dat > all.dat

clean: 
	rm -f sweep write_obbs create_hv_mesh hv_cube .#* *# *.h5m *.vtk *.stl *~ *.p *.o
