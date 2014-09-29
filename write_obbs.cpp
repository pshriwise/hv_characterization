
#include <iostream>
#include "hv_mesh_gen.hpp"
#include "DagMC.hpp"

int main( int argc, char** argv) 
{

  //hangle arguments
  if( 3 != argc )
    {
      std::cout << "Please include all necessary arguments! Exiting..." << std::endl;
      std::cout << "To generate a high valence region with area ";
      std::cout << "fraction A_f (double) and verts of valency n:" << std::endl;
      std::cout << "$ ./create_hv_mesh <A_f> <n>" << std::endl;
      return 0;
    } 


  double A_f  = atof(argv[1]);

  int valence = atoi(argv[2]);

  //make sure A_f is valid
  if ( (1.0/6.0) < A_f || A_f < 0 ) { std::cout << "Area fraction A_f must be between 0 and 1/6 (for now.)" << std::endl; return 0; }

  //now create the hv mesh and write to file 
  prep_mesh( A_f, valence );

  //start a dagmc instance and load the modified cube file
  moab::DagMC *dag_inst = moab::DagMC::instance();

  moab::ErrorCode result;
  result = dag_inst->load_file( "cube_mod.h5m" );
  if( MB_SUCCESS != result) return MB_FAILURE;
  
  //generate the obb tree
  result = dag_inst->init_OBBTree();
  if( MB_SUCCESS != result) return MB_FAILURE;
  
  moab::Range vols;
  result = get_volumes( dag_inst->moab_instance(), vols ); 

  
  //now write the OBBs to a set of files based on depth in tree
  std::string base_name = "OBBS";
  result = write_obb_mesh( dag_inst, vols[0], base_name );
  if( MB_SUCCESS != result) return MB_FAILURE;

  return 0;
  
}
