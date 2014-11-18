

#include <assert.h>

#include "moab/ProgOptions.hpp"
#include "MBCore.hpp"
#include "hv_mesh_gen.hpp"

using namespace moab; 

inline void ERR_CHECK( moab::ErrorCode rval )
{

  if (MB_SUCCESS != rval)
    {
      assert(false);
    }
}

int main( int argc, char** argv) 
{

  moab::ErrorCode rval; 
  std::string filename;
 
  Interface *mb = new Core(); 

  ProgOptions po( "Program for running ray_fire tests on a model.");

  po.addRequiredArg<std::string>( "input_file", "Filename of model to test.", &filename); 

  po.parseCommandLine( argc, argv); 

  rval = mb->load_file( filename.c_str() ); 
  ERR_CHECK(rval);

  //get the number of triangles in the model 
  std::vector<EntityHandle> tris;

  rval = mb->get_entities_by_type( 0, MBTRI, tris);
  ERR_CHECK(rval); 
  
  std::cout << "There are " << tris.size() << " triangles in the model." << std::endl; 

  return 0;

}
