
CC = g++

INCLUDE_DIRS = -I/home/shriwise/dagmc_blds/lasso/include \
	       -I/home/shriwise/dagmc_blds/cgm_dev/include/ \
               -I/home/shriwise/dagmc_blds/moabs/include \
	       -I/home/shriwise/meshkit/include/

LIBS = -L/home/shriwise/dagmc_blds/moabs/lib/ -lMOAB -ldagmc \
       -L/home/shriwise/meshkit/lib/ -lMeshKit


all: build

build: sweep write_obbs hv_cube

test_model: hv_mesh_gen.o ProgOptions.o gen_mb_funcs.o
	$(CC) test_model.cpp hv_mesh_gen.o gen_mb_funcs.o ProgOptions.o -o test_model $(INCLUDE_DIRS) $(LIBS)

sweep: gen_mb_funcs.o hv_mesh_gen.o ProgOptions.o
	$(CC) sweep.cpp hv_mesh_gen.o gen_mb_funcs.o -o sweep $(INCLUDE_DIRS) $(LIBS)

write_obbs: obbhexwriter.o ProgOptions.o gen_mb_funcs.o
	$(CC) write_obbs.cpp gen_mb_funcs.o obbhexwriter.o -o write_obbs $(INCLUDE_DIRS) $(LIBS)

hv_cube: hv_mesh_gen.o gen_mb_funcs.o
	$(CC) hv_cube.cpp hv_mesh_gen.o gen_mb_funcs.o -o hv_cube $(INCLUDE_DIRS) $(LIBS)

obbhexwriter.o:
	$(CC) obbhexwriter.cpp $(INCLUDE_DIRS) $(LIBS) -c -o obbhexwriter.o

hv_mesh_gen.o:
	$(CC) hv_mesh_gen.cpp $(INCLUDE_DIRS) $(LIBS) -c -o hv_mesh_gen.o

ProgOptions.o:
	$(CC) ProgOptions.cpp -c -o ProgOptions.o

gen_mb_funcs.o:
	$(CC) gen_mb_funcs.cpp -c -o gen_mb_funcs.o $(INCLUDE_DIRS) $(LIBS)

clean: 
	rm -f test_model sweep write_obbs create_hv_mesh hv_cube .#* *# *.h5m *.vtk *.stl *~ *.p *.o
