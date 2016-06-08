

#include "moab/Types.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/ProgOptions.hpp"
#include "ray_writer.hpp"


int main (int argc, char** argv)
{

  moab::Interface *MBI = new moab::Core();

  ProgOptions po("Traverse Ray: A program for writing out all OBBs encoutnered by a single ray traversal.");

  int vol_id;
  double ray_origin[3], ray_dir[3];
  std::string fn;
  po.addRequiredArg<double>("x", "Ray origin's x position.", &(ray_origin[0]));
  po.addRequiredArg<double>("y", "Ray origin's y position.", &(ray_origin[1]));
  po.addRequiredArg<double>("z", "Ray origin's z position.", &(ray_origin[2]));
  po.addRequiredArg<double>("u", "Ray direction u value.", &(ray_dir[0]));
  po.addRequiredArg<double>("v", "Ray direction v value.", &(ray_dir[1]));
  po.addRequiredArg<double>("w", "Ray direction w value.", &(ray_dir[2]));
  po.addRequiredArg<int>("vol_id", "GLOBAL ID of the volume to fire the ray in.", &vol_id);
  po.addRequiredArg<std::string>("filename", "DAGMC model.", &fn);

  po.parseCommandLine( argc, argv );

  
  

  
  return 0;

}
