

#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "gen.hpp"
#include "DagMC.hpp"
#include "moab/OrientedBoxTreeTool.hpp"

//local includes
#include "obbhexwriter.hpp"
#include "gen_mb_funcs.hpp"

using namespace MeshKit;


////// Functions for generating the hv region \\\\\\\\\\
//creates a cube and meshes it using the SolidSurfaceMesher
void create_cube(); 
// calls functions for setting up the inital mesh, getting the surface we want to make high-valence and then refacets the surface
void prep_mesh( double A_f, int valence );
// creates new facets for the square surface (in const Z-plane) with a high-valence region of size A_f*(surface_area) and a valency of n
void refacet_surface( moab::EntityHandle surf, double A_f, int valence );
// returns the verts of an empty square in the center of the surface and surrounds this square w/ triangles in a watertight fashion
void generate_box_space( moab::EntityHandle surf, double A_f, std::vector<moab::EntityHandle> &box_verts, int axis );
// creates the triangles for the high-valence area
void make_hv_region( moab::EntityHandle surf, std::vector<moab::EntityHandle> box_verts, int n );
// returns a surface that is constant in Z

void get_hv_surf( MEntVector surfs, moab::EntityHandle &hv_surf );
// removes and deletes all triangles in the given surface
void tear_down_surface( moab::EntityHandle surf );
// returns the area of a polygon given the ordered verts
double polygon_area( std::vector<moab::EntityHandle> verts );


moab::ErrorCode test_hv_cube_mesh( double A_f, double valence, double &ray_fire_time, double worst_split_ratio = 0.95 );

//call into save mesh
void save_mesh( std::string filename ); 

void order_corner_verts(std::vector<moab::EntityHandle> corners, std::vector<moab::EntityHandle>& ordered_corners);
