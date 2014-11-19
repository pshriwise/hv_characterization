
#include <ctime>
#include <limits>
#include <math.h>
#include <time.h>

#include "DagMC.hpp"
#include "MBCore.hpp"
#include "MBCartVect.hpp"

////// Functions for analyzing the hv region \\\\\\\\\\\\\

static const double PI = acos(-1.0);
static const double denom = 1.0 / ((double) RAND_MAX);
static const double denomPI = PI * denom;


static double facet_tol = 1e-4;
static double source_rad = 0;
static int vol_index = 1;
static int num_random_rays = 1000;
static int randseed = 12345;
static bool do_stat_report = false;
static bool do_trv_stats   = false;
static double location_az = 2.0 * PI;
static double direction_az = location_az;
static const char* pyfile = NULL;

inline void RNDVEC(moab::CartVect& uvw, double &az) 
{
  
  double theta = (PI/2)*denom*az*rand()+(PI/4);
  double u = (sqrt(2)/2) - sqrt(2)*denom*rand();
  uvw[0] = sqrt(1-u*u)*cos(theta);
  uvw[1] = sqrt(1-u*u)*sin(theta);
  uvw[2] = u;

}


void fire_rand_rays( moab::DagMC *dagi, moab::EntityHandle vol, int num_rand_rays, double &avg_fire_time, moab::CartVect ray_source);
