
#include <iostream>
#include <sstream>
#include <iomanip> 
#include <limits> 
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>


#include "gen.hpp"
#include "meshkit/MKCore.hpp"

using namespace MeshKit;

MKCore *mk;

int main(int argc, char **argv)
{

  mk = new MKCore();
  
  //Load the mesh file
  mk->load_mesh("cube.h5m");

  //Get all of the surface ModelEnts
  
  
  
  return 0;

}
