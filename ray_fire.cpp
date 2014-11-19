

#include <iostream>

#include "ray_fire.hpp"
#include "DagMC.hpp"
#include "MBCore.hpp"
#include "MBCartVect.hpp"

void fire_rand_rays( moab::DagMC *dagi, moab::EntityHandle vol, int num_rand_rays, double &avg_fire_time, moab::CartVect ray_source)
{

  //always start with the same seed
  srand(randseed);

  moab::CartVect xyz, uvw;

  int random_rays_missed = 0; 

  moab::EntityHandle dum;
  double dist;
  moab::OrientedBoxTreeTool::TrvStats trv_stats;

  //start timer
  std::clock_t start_time, end_time;

  start_time = std::clock();
  for (int j = 0; j < num_random_rays; j++) {
    RNDVEC(uvw, location_az);
    
    xyz = uvw * source_rad + ray_source;
    if (source_rad >= 0.0) {
      RNDVEC(uvw, direction_az);
    }
    
#ifdef DEBUG
    std::cout << "x,y,z,u,v,w,u^2 + v^2 + w^2 = " << xyz 
	      << " " << uvw << " " << uvw%uvw << std::endl;
    uavg += uvw[0]; vavg += uvw[1]; wavg += uvw[2];
#endif
    // added ray orientation
    dagi->ray_fire(vol, xyz.array(), uvw.array(), dum, dist, NULL, 0, 1, &trv_stats );
    
    if( dum == 0){ random_rays_missed++; }
    
  }
  //end timer 
  end_time = std::clock();
  
  if( random_rays_missed ){
    std::cout << "Warning: " << random_rays_missed << " random rays did not hit the target volume" << std::endl;
  }
  
  avg_fire_time = (end_time - start_time) / (double)(CLOCKS_PER_SEC * 100000); 
  std::cout << "Average ray fire time: " << avg_fire_time << "s" << std::endl;
}
