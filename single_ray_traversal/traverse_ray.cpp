

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
  double ray_origin[3];
  po.addRequiredArg<int>("vol_id", "GLOBAL ID of the volume to fire the ray in.", &vol_id);
  //  po.addRequiredArg<double[3]>("Ray Origin", "Origin of the ray to fire.", &ray_origin);
  

  
  return 0;

}
