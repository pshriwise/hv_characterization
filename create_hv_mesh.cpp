
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

MeshKit::MKCore *mk;

int main(int argc, char **argv)
{

  mk = new MeshKit::MKCore();
  
  //Load the mesh file
  mk->load_mesh("cube.h5m");

  
  
  return 0;

}
