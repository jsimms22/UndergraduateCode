#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"
#include "omp.h"

using std::vector;

#define _cutoff 0.01
#define _density 0.0005

double binSize, gridSize;
int binNum;

void buildBins(vector<bin_t>& bins, particle_t* particles, int n)
{
    gridSize = sqrt(n*_density);
    binSize = _cutoff * 2;
    binNum = int(gridSize / binSize) + 1;
  
    bins.resize(binNum * binNum);

    for (int i = 0; i < n; i++) 
    {
      int x = int(particles[i].x / binSize);
      int y = int(particles[i].y / binSize);
      bins[x*binNum+y].push_back(particles[i]);
    }
}
//
//  benchmarking program
//
int main( int argc, char **argv )
{   
    int navg,nabsavg=0,numthreads; 
    double dmin, absmin=1.0,davg,absavg=0.0;
	
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" ); 
        printf( "-no turns off all correctness checks and particle output\n");   
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;      

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    vector<bin_t> particle_bins;
    vector<bin_t> temp; 
    
    set_size( n );
    init_particles( n, particles );

    buildBins(particle_bins, particles, n);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    #pragma omp parallel private(dmin) 
    {
        #pragma omp master
        {
            numthreads = omp_get_num_threads();
	    temp.resize(numthreads);
        }
        for (int step = 0; step < NSTEPS; step++ )
        {
            navg = 0;
            davg = 0.0;
	    dmin = 1.0;
        //
        //  compute all forces
        //
            #pragma omp for reduction (+:navg) reduction(+:davg)
            for (int i = 0; i < binNum; i++ )
            {
	        for (int j = 0; j < binNum; j++)
	        {
	            bin_t& vec = particle_bins[i*binNum+j];
		    for (int k = 0; k < vec.size(); k++)
		    {
		        vec[k].ax = vec[k].ay = 0.0;
		    }   
		    for (int dx = -1; dx <= 1; dx++)
		    {
		        for (int dy = -1; dy <= 1; dy++)
		        {
			    if (i + dx >= 0 && i + dx < binNum && j + dy >= 0 && j + dy < binNum )
			    {
			        bin_t& vec2 = particle_bins[(i+dx) * binNum + j + dy];
			        for (int k = 0; k < vec.size(); k++)
			        {
	  		            for (int l = 0; l < vec2.size(); l++)
			            {
				        apply_force(vec[k],vec2[l],&dmin,&davg,&navg);
		                    }
			        }
			    }
		        }
		    }
	        } 
            } 	
        //
        //  move particles
        //
            int id = omp_get_thread_num();
	    bin_t& tmp = temp[id];
	    tmp.clear();	
            #pragma omp for
            for (int i = 0; i < binNum; i++ ) 
            {
	        for (int j = 0; j < binNum; j++)
	        {
		    bin_t& vec = particle_bins[i*binNum+j];
		    int tail = vec.size(), k = 0;
		    for (;k < tail;)
		    {
		        move(vec[k]);
		        int x = int(vec[k].x / binSize);
		        int y = int(vec[k].y / binSize);
		        if (x == i && y == j)
		        {
			    k++;
		        } else {
			    tmp.push_back(vec[k]);
			    vec[k] = vec[--tail];
		        }
		    }
	            vec.resize(k);
	        }
	    }
	    #pragma omp master
	    {
	        for (int j = 0; j < numthreads; j++)
	        {
		    bin_t& tmp = temp[j];
		    for (int i = 0; i < tmp.size(); i++)
		    {	
	    		int x = int(tmp[i].x / binSize);
	    		int y = int(tmp[i].y / binSize);
	    		particle_bins[x*binNum+y].push_back(tmp[i]);
		    }
	        }
  	    }
	
    if (find_option( argc, argv, "-no" ) == -1 ) 
    {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) { 
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
	  if (dmin < absmin) absmin = dmin; 
		
          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
	#pragma omp barrier
    }
}
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");
    
    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
