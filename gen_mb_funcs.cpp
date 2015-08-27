
#include <assert.h>
//#include "MBCore.hpp"
#include "MBTagConventions.hpp"
#include "gen_mb_funcs.hpp"

moab::ErrorCode get_volumes( moab::Interface* mb, moab::Range &volumes)
 {

   moab::ErrorCode rval; 

   //get the GEOM_DIM tag
   moab::Tag geom_dim;
   rval = mb->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1,
			      moab::MB_TYPE_INTEGER, geom_dim);
   assert( moab::MB_SUCCESS == rval);
   if( moab::MB_SUCCESS != rval) return rval; 

   // geom dimension of 3 should indicate a volume
   int dim = 3;
   void* ptr = &dim;
   
   //get all volume meshsets (should only be one)
   rval = mb->get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &geom_dim, &ptr, 1, volumes);
   assert( moab::MB_SUCCESS == rval);
   if( moab::MB_SUCCESS != rval) return rval; 

   return moab::MB_SUCCESS;
 }
