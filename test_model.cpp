

#include <assert.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>

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
  
  //Argument Vars
  std::string filename;
  double min_facet_tol = 1e-03;
  double max_facet_tol = 0.9;
  int facet_tol_intervals = 4;
 
  Interface *mb = new Core(); 

  ProgOptions po( "Program for running ray_fire tests on a model.");

  po.addRequiredArg<std::string>( "input_file", "Filename of model to test.", &filename); 

  po.addOpt<double>( "min_facet_tol", "Minimum faceting tolerance to test.", &min_facet_tol);

  po.addOpt<double>( "max_facet_tol", "Maximum faceting tolerance allowed. The program will halt if this is reached.", &max_facet_tol);

  po.addOpt<int>( "num_ints", "Number of intervals to add to the minimum faceting tolerance. These are added in multiples of 5 and 2 depending on the iteration.", &facet_tol_intervals); 

  po.parseCommandLine( argc, argv); 
  
  //data storage containers
  std::vector<int> tri_numbers;
  std::vector<double> timing;
  std::vector<double> facet_tols; 
  
  //initalize facet_tol param
  double facet_tol = min_facet_tol;
  for(unsigned int i = 0; i <= facet_tol_intervals; i++)
    {

      if ( 0 != i) facet_tol = (0 == i%2) ? min_facet_tol*pow(5,double(i)/2)*pow(2,double(i)/2) : min_facet_tol*pow(5,((double(i)-1)/2)+1)*pow(2,double((i)-1)/2);

      // stop the program if we've gone over the maximum faceting tolerance
      if (facet_tol > max_facet_tol)
	{
	  std::cout << "Maximum faceting tolerance reached, exiting loop and writing data..." << std::endl;
	  break; 
	}

      //save faceting tolerance for write to datafile
      facet_tols.push_back(facet_tol);
      
      //setup fileoptions string
      std::stringstream options; options << "FACET_DISTANCE_TOLERANCE=" << std::setprecision(12) << facet_tol;
      std::string opts = options.str();

      //Load the file and facet
      rval = mb->load_file( filename.c_str(), 0, opts.c_str() ); 
      ERR_CHECK(rval);
      
      //get the number of triangles in the model 
      std::vector<EntityHandle> tris;
      rval = mb->get_entities_by_type( 0, MBTRI, tris);
      ERR_CHECK(rval); 
      std::cout << "There are " << tris.size() << " triangles in the model." << std::endl; 
      tri_numbers.push_back( int(tris.size()) ); 

      //hand model's moab instance to dagmc
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
      
      timing.push_back( avg_fire_time ); 
    }

  //Write data to file
  std::ofstream datafile; 
  datafile.open("model_data.dat");

  for( unsigned int i = 0; i < timing.size(); i++)
    datafile << tri_numbers[i] << "\t" << timing[i] << "\t" << facet_tols[i] << std::endl; 

  return 0;

}
