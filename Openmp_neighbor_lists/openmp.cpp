#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"

/*
omp_lock_t lock;
omp_init_lock(&lock);
*/

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
    set_size( n );
    init_particles( n, particles );
    
//stuff for cells structure
	double cutoff = get_cutoff_value();
    	const double size = get_size();
    	const int cells_per_side = (int)((size/cutoff));
	cell *cells = (cell*) malloc(cells_per_side*cells_per_side*sizeof(cell));
	set_cells_per_side(cells_per_side);
	init_cells(cells, n, particles);

//	printf("simulation box size: %f, cutoff radius: %f, number of bins: %d, %d\n",size,cutoff,cells_per_side,cells_per_side*cells_per_side);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    #pragma omp parallel private(dmin)
    {				
    numthreads = omp_get_num_threads();

	for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
	dmin = 1.0;
        //
        //  compute all forces
        //
      

#pragma omp for
for(int i =0;i<n;i++){
    particles[i].ax = particles[i].ay = 0;  	//reset accelerations.
}


#pragma omp for reduction (+:navg) reduction(+:davg) schedule(auto)
for(int i =0;i<cells_per_side*cells_per_side;i++){   
    if(cells[i].particles[0]!=-1){			
	compute_cell(i, cells, particles, dmin, davg, navg);  
    }//END OF IF cell i is empty
}//END OF i LOOP (THROUGH ALL CELLS)
// */
//============================================================================================
		
        //
        //  move particles
        //
       #pragma omp for
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );
        


	if( find_option( argc, argv, "-no" ) == -1 ) 
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
	
//	#pragma omp master
//	update_bins(cells,n,particles); 

	#pragma omp for schedule(auto)
	for(int i=0;i<cells_per_side*cells_per_side;i++){
		for (int k =0; k<MAX_PARTICLES_PER_CELL;k++){
		   cells[i].particles[k] = -1;
		}	
	}	
	

	#pragma omp single
	{
	for(int i =0;i<n;i++){
	int xcell = (particles[i].x/size) *(cells_per_side);
	int ycell = (particles[i].y/size) *(cells_per_side);
		for(int k=0;k<MAX_PARTICLES_PER_CELL;k++){
		   if(cells[xcell+ycell*cells_per_side].particles[k]==-1){
			cells[xcell+ycell*cells_per_side].particles[k]=i;
			break;
		   }	
		}
	}
	}


 	  //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
}
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
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
    free( cells );
    if( fsave )
        fclose( fsave );
    
    return 0;
}








