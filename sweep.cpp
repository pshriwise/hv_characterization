// The purpose of this file is to creat a bash script which will sweep
// the available parameter spaces of A_f and n

#include <iostream>
#include <fstream>


void gnuplot_script( int A_f, int valence);

int main(int argc, char** argv)
{

  int area_intervals = 10;
  int valence_intervals = 3*10;
  std::ofstream bash_file;
  std::ofstream param_file;

  bash_file.open("sweep_param_space");
  param_file.open("params.dat");

  bash_file << "#!/bin/bash" << std::endl << std::endl;

  bash_file << "# Sweep parameter space of both A_f and n with " 
            << area_intervals << " area fraction intervals " 
	    << "and " << valence_intervals << " valence interfvals." << std::endl;

  double max_A_f = 1;
  
  int max_n = 1e6;

  std::string rayfire_program = "~/dagmc_blds/moabs/bld/tools/dagmc/ray_fire_test";
  for(unsigned int i = 1; i < area_intervals; i++)
    {

      for(unsigned int j = 1; j < valence_intervals; j++)
	{
	  //set values for this run
	  double A_f = (double) i * ( max_A_f / area_intervals);
	  double n = (double) j * ( max_n / valence_intervals );
	  
	  //write these params to our data file
	  param_file << A_f << "\t" << n << std::endl;
	  
	  //create the new cube_mod mesh
	  bash_file << "./create_hv_mesh " << A_f << " " << n << std::endl;
	  
	  //test the mesh
	  bash_file << rayfire_program << " cube_mod.h5m -n 10000 -c 0 0 0" << std::endl;
	  bash_file << std::endl;
	}
    }

  gnuplot_script( area_intervals, valence_intervals );
  return 0;
}


void gnuplot_script( int area_intervals, int valence_intervals)
{

  std::ofstream gp_file; 
  gp_file.open("timing_plot.p");
    
  gp_file << "set dgrid3d " << area_intervals << ", " << valence_intervals << std::endl;
  gp_file << "set style data lines" << std::endl;

  //gp_file << "set hidden3d" << std::endl;

  gp_file << "set xlabel 'Area Fraction (A_f)' " << std::endl;
  gp_file << "set ylabel 'Valence (n)'" << std::endl;
  gp_file << "set zlabel 'Avg Ray Fire Time (micro-s)' " << std::endl;

  gp_file << "splot 'data.dat' u 1:2:($3 * 1000000) notitle" << std::endl;

  gp_file << "pause -1" << std::endl;

}





