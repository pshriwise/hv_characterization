
#include <iostream>
#include "obbhexwriter.hpp"
#include "DagMC.hpp"
#include "ProgOptions.hpp"
#include "gen_mb_funcs.hpp"

moab::ErrorCode write_obb_mesh( moab::DagMC *dag, moab::EntityHandle vol, std::string& base_filename, bool write_tris = false);

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

  
  //now write the OBBs to a set of files based on depth in tree
  std::string base_name = "OBBS";
  result = write_obb_mesh( dag_inst, vols[0], base_name, write_tris );
  if( MB_SUCCESS != result) return MB_FAILURE;

  return 0;
  
}


moab::ErrorCode write_obb_mesh( moab::DagMC *dag, moab::EntityHandle vol, std::string& base_filename, bool write_tris) 
{

  moab::ErrorCode rval; 

  moab::EntityHandle root;

  rval = dag->get_root( vol, root );
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 

  moab::OrientedBoxTreeTool *obbtool = dag->obb_tree();

  rval = obbtool->stats( root, std::cout );
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 

  //make a new moab core for the box hexes
  moab::Core mbi2;

  OBBHexWriter hw( obbtool, &mbi2, write_tris );

  moab::OrientedBoxTreeTool::TrvStats tree_stats;

  rval = obbtool->preorder_traverse( root, hw, &tree_stats );
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 

  rval = hw.write_to_files( base_filename );
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 


}
