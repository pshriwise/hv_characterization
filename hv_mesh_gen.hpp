
#include <iostream>
#include <ctime>
#include <sstream>
#include <iomanip> 
#include <limits> 
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "gen.hpp"
#include "DagMC.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include "hexwriter.hpp"


//timing includes
#include <fstream>
#include <fcntl.h>
#include <unistd.h>
#include <fcntl.h>
#include <cstdlib>

using namespace MeshKit;


////// Functions for generating the hv region \\\\\\\\\\

void prep_mesh( double A_f, int valence );
// creates new facets for the square surface (in const Z-plane) with a high-valence region of size A_f*(surface_area) and a valency of n
void refacet_surface( moab::EntityHandle surf, double A_f, int valence );
// returns the verts of an empty square in the center of the surface and surrounds this square w/ triangles in a watertight fashion
void generate_box_space( moab::EntityHandle surf, double A_f, std::vector<moab::EntityHandle> &box_verts );
// creates the triangles for the high-valence area
void make_hv_region( moab::EntityHandle surf, std::vector<moab::EntityHandle> box_verts, int n );
// returns a surface that is constant in Z

void get_hv_surf( MEntVector surfs, moab::EntityHandle &hv_surf );
// removes and deletes all triangles in the given surface
void tear_down_surface( moab::EntityHandle surf );
// returns the area of a polygon given the ordered verts
double polygon_area( std::vector<moab::EntityHandle> verts );

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

inline void RNDVEC(CartVect& uvw, double &az) 
{
  
  double theta = denom * az * rand();
  double u = 2*denom*rand() - 1;
  uvw[0] = sqrt(1-u*u)*cos(theta);
  uvw[1] = sqrt(1-u*u)*sin(theta);
  uvw[2] = u;

}


moab::ErrorCode write_obb_mesh( moab::DagMC *dag, moab::EntityHandle vol, std::string& base_filename);

moab::ErrorCode get_volumes( moab::Interface* mb, moab::Range &volumes);

void fire_rand_rays( moab::DagMC *dagi, moab::EntityHandle vol, int num_rand_rays, double &avg_fire_time, moab::CartVect ray_source);


moab::ErrorCode test_hv_cube_mesh( double A_f, double valence, double &ray_fire_time );
