
#include <iostream>
#include "hv_mesh_gen.hpp"
#include "DagMC.hpp"

int main( int argc, char** argv) 
{
  //hangle arguments
  if( 2 != argc )
    {
      std::cout << "To write the oriented bounding boxes of a dagmc .h5m file";
      std::cout << "in a .vtk format:"<< std::endl;
      std::cout << "$./write_obbs <filename.h5m>" << std::endl;
      return 0;
    } 


  //assign the filename
  char *filename = argv[1];

  //start a dagmc instance and load the modified cube file
  moab::DagMC *dag_inst = moab::DagMC::instance();

  moab::ErrorCode result;
  result = dag_inst->load_file( filename );
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
