// The purpose of this file is to creat a bash script which will sweep
// the available parameter spaces of A_f and n

#include <iostream>
#include <fstream>
#include <string>
//local includes
#include "hv_mesh_gen.hpp"
#include "moab/ProgOptions.hpp"

#include "mpi.h"

void gnuplot_script();

int main(int argc, char** argv)
{


  // start MPI
  MPI_Init(NULL,NULL);

  //
  int cpu_id,size; // cpu id and world size
  MPI_Comm_rank(MPI_COMM_WORLD,&cpu_id);
  MPI_Comm_size(MPI_COMM_WORLD,&size);


  int num_tasks;
  std::vector<double> area_fractions; // list of area fractions
  std::vector<int> valences; // list of valences
  std::vector<double> timing; //vector for results
  std::vector<double> worst_split_ratios; //vector of worst-case split ratios

  // if master task, do the reading and set up. Otherwise wait
  if(cpu_id == 0 )
    {
      int area_intervals = 10;
      int valence_intervals = 4*area_intervals;
      double max_A_f;
      int max_n;

      //Split Ratio Values
      double min_wsr, max_wsr;
      int wsr_intervals;

      
      ProgOptions po("sweep: a program for sweeping two parameter spaces for a manually modified high-valence mesh");
      
      po.addRequiredArg<double>( "af", "Fraction of the cube surface area that should be made high-valence.", &max_A_f);
      
      po.addRequiredArg<int>( "v", "Valence of the high-valence region", &max_n);
      
      po.addOpt<int>( "af_ints", "Number of intervals for the Area fraction param", &area_intervals );
      
      po.addOpt<int>( "v_ints", "Number of intervals for the valence param", &valence_intervals);

      po.addOpt<double>("min_wsr", "Minimum value of the worst-case split ratio. (Used in creation of OBB Trees.)");

      po.addOpt<double>("max_wsr", "Maximum value of the worst-case split ratio. (Used in creation of OBB Trees.)");

      po.addOpt<int>("wsr_ints", "Number of intervals for the worst-case split ratio."); 

      po.parseCommandLine( argc, argv );

      //set defaults for the worst-case split ratio
      if( !po.getOpt( "min_wsr", &min_wsr ) ) min_wsr = 0.95;
      if( !po.getOpt( "max_wsr", &max_wsr ) ) max_wsr = 0.95;
      if( !po.getOpt( "wsr_ints", &wsr_intervals ) ) wsr_intervals = 1;
      
      assert( min_wsr <= max_wsr );
      if (max_wsr == min_wsr)  assert( 1 == wsr_intervals ); 
      
      
      
      double A_f = 0;
      int valence = 0;
      double wsr = 0;
      
      // num_jobs
      num_tasks = (wsr_intervals)*(area_intervals)*valence_intervals;

      // precompute input values      
      for(unsigned int i=1; i <=area_intervals; i++)
	{
	  
	  A_f = (double)i * ( max_A_f / area_intervals);
	  //write this value to params file everytime
	  
	  for(unsigned int j=0; j < valence_intervals; j++)
	    {
	      valence = (double)j * ( max_n / valence_intervals );

	      for(unsigned int k=0; k <wsr_intervals; k++)
		{
		  wsr = min_wsr + k*((max_wsr-min_wsr)/(wsr_intervals));
	      //the first time we go through the inner loop, 
	      //write all of the valence values
	      area_fractions.push_back(A_f);
	      valences.push_back(valence);     
	      worst_split_ratios.push_back(wsr);
	      timing.push_back(0.0);      
		}
	    }
	}
    }//end master-specific commands

  // synchronize
  MPI_Barrier(MPI_COMM_WORLD);
  
  // share num_tasks
  MPI_Bcast(&num_tasks,1,MPI_INT,0,MPI_COMM_WORLD);

  //resize vectors on processors other than slaves
  if (cpu_id != 0 )
    {
      area_fractions.resize(num_tasks);
      valences.resize(num_tasks); 
      worst_split_ratios.resize(num_tasks);
      timing.resize(num_tasks);
    }

  MPI_Barrier(MPI_COMM_WORLD);

  // share area_fractions
  MPI_Bcast(&area_fractions.front(),num_tasks,MPI_DOUBLE,0,MPI_COMM_WORLD);
  // share valences
  MPI_Bcast(&valences.front(),num_tasks,MPI_INT,0,MPI_COMM_WORLD);
  // share ratios
  MPI_Bcast(&worst_split_ratios.front(),num_tasks,MPI_DOUBLE,0,MPI_COMM_WORLD);
  // share timing
  MPI_Bcast(&timing.front(),num_tasks,MPI_DOUBLE,0,MPI_COMM_WORLD);
  //
  MPI_Barrier(MPI_COMM_WORLD);

  // broadcast these vectors
  // broadcast area fractions, valences and timiing
      
  //chop up job
  int num_job_on_proc = num_tasks/size; 
  int num_task_master = (num_tasks/size) + (num_tasks%size);

  int job_id_start = 0;
  int job_id_end = 0;

  if ( cpu_id == 0 ) 
    {
      job_id_start = 0;
      job_id_end = num_task_master - 1;
    }
  else
    {
      job_id_start = num_task_master + ((cpu_id-1)*num_job_on_proc);
      job_id_end = job_id_start + num_job_on_proc - 1;
    }
  
  // run the jobs
  for ( int i = job_id_start ; i <= job_id_end ; i++ )
    {
      test_hv_cube_mesh(area_fractions[i] ,valences[i], timing[i], worst_split_ratios[i]);
    }

  MPI_Barrier(MPI_COMM_WORLD);

  //create a new vector for the results
  std::vector<double>result(num_tasks);

  MPI_Reduce(&timing.front(),&result.front(),num_tasks,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  
  
  // now do output
  if(cpu_id == 0 ) 
    {   
      std::ofstream datafile; 
      datafile.open("data.dat"); 
      
      for( unsigned int i = 0; i < area_fractions.size(); i++)
	datafile << valences[i] << "\t" << area_fractions[i] << "\t" << worst_split_ratios[i] << "\t" << result[i] << std::endl;
      datafile.close();
    }
  
  MPI_Barrier(MPI_COMM_WORLD);

  // end mpi
  MPI_Finalize();


  return 0;

}

