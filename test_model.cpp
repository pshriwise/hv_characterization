

#include <assert.h>
#include <fstream>
#include <sstream>
#include <iomanip>
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
  
  double min_facet_tol = 1e-01;
  double max_facet_tol = 1e-02;
  int facet_tol_intervals = 2;
  
  std::vector<int> tri_numbers;
  std::vector<double> timing;

  for(unsigned int i = 1; i <=facet_tol_intervals; i++)
    {
      //Load the file and facet
      double facet_tol = min_facet_tol + ((max_facet_tol-min_facet_tol)*i)/double(facet_tol_intervals);
      std::stringstream options; options << "FACET_DISTANCE_TOLERANCE=" << std::setprecision(12) << facet_tol;
      std::string opts = options.str();
      rval = mb->load_file( filename.c_str(), 0, opts.c_str() ); 
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

      //clear everything in this moab instance and destroy the DagMC instance 
      rval = mb->delete_mesh();
      ERR_CHECK(rval);
      
      tri_numbers.push_back( int(tris.size()) ); 
      timing.push_back( avg_fire_time ); 
    }

  //Write data to file
  std::ofstream datafile; 
  datafile.open("model_data.dat");

  for( unsigned int i = 0; i < timing.size(); i++)
  datafile << tri_numbers[i] << "\t" << timing[i] << std::endl; 

  return 0;

}
