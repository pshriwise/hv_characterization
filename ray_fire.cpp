

#include <iostream>

#include "ray_fire.hpp"
#include "DagMC.hpp"
#include "MBCore.hpp"
#include "MBCartVect.hpp"
#include <sys/resource.h>

void fire_rand_rays( moab::DagMC *dagi, moab::EntityHandle vol, int num_rand_rays, double &avg_fire_time, moab::CartVect ray_source)
{

  std::cout.precision(6); 
  //always start with the same seed
  srand(randseed);

  moab::CartVect xyz, uvw;

  int random_rays_missed = 0; 

  moab::EntityHandle dum;
  double dist;
  moab::OrientedBoxTreeTool::TrvStats trv_stats;

  //set start times
  double ttime_s, utime_s, stime_s, ttime_e, utime_e, stime_e;
  get_time(ttime_s, utime_s, stime_s);

  for (int j = 0; j < num_rand_rays; j++) {

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
  //auto end_epoch = std::chrono::high_resolution_clock::now();
  //auto end_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_epoch-start_time).count();
  
  if( random_rays_missed ){
    std::cout << "Warning: " << random_rays_missed << " random rays did not hit the target volume" << std::endl;
  }

  get_time(ttime_e, utime_e, stime_e);

  double timewith  = ttime_e-ttime_s;
  std::cout << "Total time for ray fires: " << timewith << std::endl; 
  avg_fire_time = timewith/(double)num_rand_rays;
  std::cout << "Average time per ray fire: " << avg_fire_time << std::endl;

}


void get_time( double &tot_time, double &user_time, double& sys_time) 
{

  struct rusage r_usage; 

  getrusage(RUSAGE_SELF, &r_usage); 

  user_time = (double)r_usage.ru_utime.tv_sec +
    ((double)r_usage.ru_utime.tv_usec/1.e6);
  sys_time = (double)r_usage.ru_stime.tv_sec +
    ((double)r_usage.ru_stime.tv_usec/1.e6);
  tot_time = user_time + sys_time;

}
