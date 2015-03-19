
#include "leafcounter.hpp"
#include "moab/Types.hpp"

using namespace moab;

LeafCounter::LeafCounter( OrientedBoxTreeTool *tool_ptr )
  : tool(tool_ptr), num_leaves(0), min_leaf_count(1e8), max_leaf_count(0), avg_leaf_count(0) {}

LeafCounter::~LeafCounter() {};

// nothing special about the leaves for this op, do nothing
ErrorCode LeafCounter::leaf( EntityHandle node ) 
{
  
  ErrorCode rval; 
  //each time we find a leaf do these operations
  num_leaves++; 

  //get the number of triangles in the leaf
  Range tris;
  rval = tool->get_moab_instance()->get_entities_by_type( node, MBTRI, tris);
  assert(MB_SUCCESS == rval );
  if ( MB_SUCCESS != rval ) return rval; 
  
  int tris_in_leaf = int(tris.size());
  if ( tris_in_leaf < min_leaf_count ) min_leaf_count = tris_in_leaf;
  if ( tris_in_leaf > max_leaf_count ) max_leaf_count = tris_in_leaf;
  avg_leaf_count += double(tris_in_leaf);


  
  return MB_SUCCESS; 

};


double LeafCounter::get_avg() { return avg_leaf_count/double(num_leaves); }


// a visit to a node, will create an OrientedBox object for that node, get a hex from that box and tag it with 
// an integer tag representing it's depth in the tree, and add it to the list of *hexes* for this tree
ErrorCode LeafCounter::visit( EntityHandle node,
			   int depth, 
			   bool& descend)
{
  descend = true; //always descend to get to the leaves

  return MB_SUCCESS;

}
