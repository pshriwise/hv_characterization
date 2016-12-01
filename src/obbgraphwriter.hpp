
#include <assert.h>
#include "moab/Core.hpp"
#include "moab/Types.hpp"
#include "moab/Interface.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include "../src/src/OrientedBox.hpp"
#include <iostream>
#include <fstream>

using namespace moab; 

class OBBGraphWriter : public moab::OrientedBoxTreeTool::Op
{
  
private:
  //pointer to the tree we're traversing
  OrientedBoxTreeTool* tool;
  
  //output file stream
  std::ofstream output;
  
public:
  
  //constructor declaration
  OBBGraphWriter( OrientedBoxTreeTool *tool_ptr );
  
  //destructor declaration
  ~OBBGraphWriter();

  // nothing special about the leaves for this op, do nothing
  ErrorCode leaf( EntityHandle node );
  //gathers the number of triangles underneath the node, its child relationships and adds them to the graph  
  ErrorCode visit( EntityHandle node,
		   int depth, 
		   bool& descend);

  //opens the graph in the designated output file
  void open_graph();

  //opens the graph in the designated output file
  void close_graph();

};


  
