
#include <assert.h>
#include "moab/Core.hpp"
#include "moab/Types.hpp"
#include "moab/Interface.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include "moab/OrientedBox.hpp"


using namespace moab;

class LeafCounter : public moab::OrientedBoxTreeTool::Op
{

public:

  //constructor declaration
  LeafCounter( OrientedBoxTreeTool *tool_ptr);

  //destructor declaration
  ~LeafCounter();

private:
  OrientedBoxTreeTool *tool;
  EntityHandle new_hex;
 
public:
  int num_leaves; 
  int min_leaf_count; 
  int max_leaf_count; 
  double avg_leaf_count; 

  // nothing special about the leaves for this op, do nothing
  ErrorCode leaf( EntityHandle node );

  // a visit to a node, will create an OrientedBox object for that node, get a hex from that box and tag it with 
  // an integer tag representing it's depth in the tree, and add it to the list of *hexes* for this tree
  ErrorCode visit( EntityHandle node,
		   int depth, 
		   bool& descend);


  double get_avg();


};
