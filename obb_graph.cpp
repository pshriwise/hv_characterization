
#include <iostream>
#include "obbgraphwriter.hpp"
#include "DagMC.hpp"
#include "moab/ProgOptions.hpp"
#include "gen_mb_funcs.hpp"

moab::ErrorCode write_obb_mesh( moab::DagMC *dag, moab::EntityHandle vol, std::string& base_filename, bool write_tris = false);

int main( int argc, char** argv) 
{

  //assign the filename
  std::string filename;
  bool write_tris = false; 

  ProgOptions po("Graph_Obbs: a program for writing the OBBs of a mesh in a dot graph format");

  po.addRequiredArg<std::string>("filename", "Mesh file", &filename);
  

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

  OBBGraphWriter gw( obbtool );

  moab::OrientedBoxTreeTool::TrvStats tree_stats;

  gw.open_graph();

  rval = obbtool->preorder_traverse( root, gw, &tree_stats );
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval;
  
  gw.close_graph();

}
