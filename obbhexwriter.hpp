
#include <assert.h>
#include "moab/Core.hpp"
#include "moab/Types.hpp"
#include "moab/Interface.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include "../src/src/OrientedBox.hpp"


using namespace moab;

class OBBHexWriter : public moab::OrientedBoxTreeTool::Op
{


private:
  OrientedBoxTreeTool *tool;
  EntityHandle new_hex;
  ErrorCode result; 
  //note: this should probably be a different interface than the original instance
  //unless it is desired that the box hexes are added to that instance
  Interface *mbi2;
  std::vector<EntityHandle> hexes;
  std::map< EntityHandle, std::vector<EntityHandle> > tri_map;
  Tag depth_tag;

  ErrorCode transfer_tri_inst( Interface* orig_inst, Interface* new_inst, std::vector<EntityHandle> &tris );

public:

  //constructor declaration
  OBBHexWriter( OrientedBoxTreeTool *tool_ptr, Interface* interface_ptr );

  //destructor declaration
  ~OBBHexWriter();

  // nothing special about the leaves for this op, do nothing
  ErrorCode leaf( EntityHandle node );

  // a visit to a node, will create an OrientedBox object for that node, get a hex from that box and tag it with 
  // an integer tag representing it's depth in the tree, and add it to the list of *hexes* for this tree
  ErrorCode visit( EntityHandle node,
		   int depth, 
		   bool& descend);


  //this function will go through the hexes found for the tree tool of this class
  //and write them to file based on their depth in the tree
  // filename formats: base_filename + "_" + depth_value + ".h5m"
  ErrorCode write_to_files( std::string base_filename );



};
