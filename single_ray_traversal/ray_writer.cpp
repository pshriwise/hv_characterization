

#include "ray_writer.hpp"

using namespace moab;

ErrorCode RayTraversalWriter::visit( EntityHandle node,
				     int depth,
				     bool& descend)
{
  //retrieve this node's bounding box
  OrientedBox box;
  ErrorCode rval = tool->box( node, box );
  assert(MB_SUCCESS == rval);
  if (MB_SUCCESS != rval)
    return rval;

  //check if the ray intersects the box
  descend = box.intersect_ray( ray_origin, ray_direction, tol, nonneg_ray_len, 
			       neg_ray_len );

  //if the ray does intersect the box we want to transfer that box and tag it with its depth in the tree
};

//NO DEFINITION FOR LEAF METHOD, WE HAVE NOTHING TO DO THERE
