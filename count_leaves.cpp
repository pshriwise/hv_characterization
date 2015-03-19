
#include <iostream>
#include "leafcounter.hpp"
#include "DagMC.hpp"
#include "moab/ProgOptions.hpp"
#include "gen_mb_funcs.hpp"

int main( int argc, char** argv) 
{

  //assign the filename
  std::string filename;
  bool write_tris = false; 

  ProgOptions po("Write_Obbs: a program for writing the OBBs of a mesh as hexes in a .vtk format");

  po.addRequiredArg<std::string>("filename", "Mesh file", &filename);
  
  po.addOpt<void>( "with-tris", "If included, the program will write the triangles contained by each hex to file.", &write_tris); 

  po.parseCommandLine( argc, argv ); 


  //start a dagmc instance and load the modified cube file
  moab::DagMC *dag_inst = moab::DagMC::instance();

  moab::ErrorCode result;
  result = dag_inst->load_file( filename.c_str() );
  if( MB_SUCCESS != result) return MB_FAILURE;
  
  //generate the obb tree
  result = dag_inst->init_OBBTree();
  if( MB_SUCCESS != result) return MB_FAILURE;
  
  moab::Range vols;
  result = get_volumes( dag_inst->moab_instance(), vols );
  
  return 0;
  
}
