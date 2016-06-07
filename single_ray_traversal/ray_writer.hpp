





#include <assert.h>
#include "moab/Core.hpp"
#include "moab/Types.hpp"
#include "moab/Interface.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include "moab/CartVect.hpp"



using namespace moab;

class RayTraversalWriter : public moab::OrientedBoxTreeTool::Op
{

private:
  //addition interface used to track and write these elements
  EntityHandle *writeSet;
  Interface *MBI;
  
  OrientedBoxTreeTool* tool;
  const CartVect       ray_origin;
  const CartVect       ray_direction;
  const double*        nonneg_ray_len;  /* length to search ahead of ray origin */
  const double*        neg_ray_len;     /* length to search behind ray origin */
  const double         tol;             /* used for box.intersect_ray, radius of
					   neighborhood for adjacent triangles,
					   and old mode of add_intersection */
  const int            minTolInt;       /* used for old mode of add_intersection */
  const EntityHandle*  rootSet;

public:
  //CONSTRUCTOR
  RayTraversalWriter( OrientedBoxTreeTool*       tool_ptr,
                      const double*              ray_point,
                      const double*              unit_ray_dir,
                      const double*              nonneg_ray_length,
		      const double*              neg_ray_length,
                      double                     tolerance,
                      int                        min_tol_intersections,
		      EntityHandle*              root_set)
    :tool(tool_ptr), ray_origin(ray_point), ray_direction(unit_ray_dir),
     nonneg_ray_len(nonneg_ray_length), neg_ray_len(neg_ray_length),
     tol(tolerance), minTolInt(min_tol_intersections), rootSet(root_set)
  {
    ErrorCode rval;
    // check the limits  
    if(nonneg_ray_len) {
      assert(0 <= *nonneg_ray_len);
    } 
    if(neg_ray_len) {
      assert(0 > *neg_ray_len);
    }
    
    if ( MBI == NULL ) {
      MBI = tool->get_moab_instance();
    }
    
    if ( writeSet == NULL) {
      rval = MBI->create_meshset(MESHSET_SET, *writeSet);
      MB_CHK_ERR_CONT(rval);//, "Could not create the meshset to hold hexes.");
    }
	  
  };

  virtual ErrorCode visit( EntityHandle node,
			   int depth,
			   bool& descend );  
    
};