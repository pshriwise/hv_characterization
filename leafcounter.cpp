
#include "leafcounter.hpp"
#include "moab/Types.hpp"
#include "meshkit/MKCore.hpp"

using namespace moab;

LeafCounter::LeafCounter( OrientedBoxTreeTool *tool_ptr )
  : tool(tool_ptr), num_leaves(0), min_leaf_count(0), max_leaf_count(0), avg_leaf_count(0) {}

LeafCounter::~LeafCounter() {};

// nothing special about the leaves for this op, do nothing
ErrorCode LeafCounter::leaf( EntityHandle node ) { return MB_SUCCESS; };

// a visit to a node, will create an OrientedBox object for that node, get a hex from that box and tag it with 
// an integer tag representing it's depth in the tree, and add it to the list of *hexes* for this tree
ErrorCode LeafCounter::visit( EntityHandle node,
			   int depth, 
			   bool& descend)
{
  descend = true; //always descend to get to the leaves

  return MB_SUCCESS;

}
