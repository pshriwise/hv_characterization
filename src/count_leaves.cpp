
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
  
  po.parseCommandLine( argc, argv ); 


  //start a dagmc instance and load the modified cube file
  moab::DagMC *dag_inst = new moab::DagMC();

  moab::ErrorCode result;
  result = dag_inst->load_file( filename.c_str() );
  if( MB_SUCCESS != result) return MB_FAILURE;
  
  //generate the obb tree
  result = dag_inst->init_OBBTree();
  if( MB_SUCCESS != result) return MB_FAILURE;
  
  moab::Range vols;
  result = get_volumes( dag_inst->moab_instance(), vols );


  //get the box tree tool 
  moab::OrientedBoxTreeTool *obbtool = dag_inst->obb_tree(); 

  LeafCounter lc( obbtool ); 

  //do a leaf count for each volume
  moab::Range::iterator vol; 
  for( vol = vols.begin(); vol != vols.end(); vol++)
    {
      moab::EntityHandle root; 
      //get the root of the obb for this volume
      result = dag_inst->get_root( *vol, root ); 
      if( MB_SUCCESS != result ) assert(!result); 
      
      result = obbtool->preorder_traverse( root, lc ); 
      if( MB_SUCCESS != result ) assert(!result); 
      

    }

  std::cout << "Minimum number of triangles in a leaf node: " << std::endl
	    << lc.min_leaf_count << std::endl;
  std::cout << "Maximum number of triangles in a leaf node: " << std::endl
	    << lc.max_leaf_count << std::endl;
  std::cout << "Average number of triangles in a leaf node: " << std::endl
	  << lc.get_avg() << std::endl;

  return 0;
  
}
