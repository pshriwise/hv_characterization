


class HexMaker; 


class HexMaker
{


private:
  moab::OrientedBoxTreeTool *tool;
  moab::EntityHandle new_hex;
  moab::ErrorCode result; 
  //note: this should probably be a different interface than the original instance
  //unless it is desired that the box hexes are added to that instance
  moab::Interface *mbi2;
  std::vector<moab::EntityHandle> hexes;
  moab::Tag depth_tag;


public:

  //constructor declaration
  HexMaker( moab::OrientedBoxTreeTool *tool_ptr, moab::Interface* interface_ptr );

  //destructor declaration
  ~HexMaker() {};

  // nothing special about the leaves for this op, do nothing
  moab::ErrorCode leaf( moab::EntityHandle node );

  // a visit to a node, will create an OrientedBox object for that node, get a hex from that box and tag it with 
  // an integer tag representing it's depth in the tree, and add it to the list of *hexes* for this tree
  moab::ErrorCode visit( moab::EntityHandle node,
				   int depth, 
			 bool& descend);


  //this function will go through the hexes found for the tree tool of this class
  //and write them to file based on their depth in the tree
  // filename formats: base_filename + "_" + depth_value + ".h5m"
  moab::ErrorCode write_to_files( std::string base_filename );



}
