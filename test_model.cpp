

#include <assert.h>
#include <fstream>

//sigma includes
#include "moab/ProgOptions.hpp"
#include "MBCartVect.hpp"
#include "MBCore.hpp"
#include "DagMC.hpp"

//local includes
#include "ray_fire.hpp"
#include "gen_mb_funcs.hpp"

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

  //hand model instance to dagmc
  DagMC *dag = DagMC::instance(mb); 

  rval = dag->load_existing_contents(); 
  ERR_CHECK(rval); 

  rval = dag->init_OBBTree();
  ERR_CHECK(rval); 

  //get the first volume in the instance (for now)
  Range vols; 
  rval = get_volumes( dag->moab_instance(), vols); 
  ERR_CHECK(rval); 

  //prep for analysis
  double avg_fire_time; 
  CartVect source_pos; 
  source_pos[0] = 0.0; source_pos[1] = 0.0; source_pos[2] = 0.0; 
  
  fire_rand_rays( dag, vols[0], 100000, avg_fire_time, source_pos); 
  
  std::cout << "The average fire time for this mesh was: " << avg_fire_time << std::endl; 

  std::ofstream datafile; 
  datafile.open("model_data.dat");

  datafile << tris.size() << "\t" << avg_fire_time << std::endl; 

  return 0;

}
