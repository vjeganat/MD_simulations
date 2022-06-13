#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

//
//  benchmarking program
//



int main( int argc, char **argv )
{    

    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;   



    
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
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

//these three necessary for initializing cell struct.
//should probably rewrite to keep these values in common.cpp
    double cutoff = get_cutoff_value();
    const double size = get_size();
    const int cells_per_side = (int)((size/cutoff));
    set_cells_per_side(cells_per_side);

// not needed for cell struct
 //   const double cell_width = size/cells_per_side;
 //   const int max_particles_per_bin = 3;
    
	cell *cells = (cell*) malloc(cells_per_side*cells_per_side*sizeof(cell));
	make_cell_map();
	init_cells(cells, n, particles);
	sort_particles(cells,particles);

//	printf("simulation box size: %f, cutoff radius: %f, number of bins: %d, %d\n",size,cutoff,cells_per_side,cells_per_side*cells_per_side);

    //
    //  simulate a number of time steps
    //

    double simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ )
    {
	navg = 0;
        davg = 0.0;    
	dmin = 1.0;   

        //
        //  compute forces
        //


//if((step+1)%300==0) sort_particles(cells,particles);


for(int i =0;i<n;i++){
    particles[i].ax = particles[i].ay = 0;  	//reset accelerations. why?
}

for(int i =0;i<cells_per_side*cells_per_side;i++){   
    if(cells[i].particles[0]!=-1){			
	compute_cell(i, cells, particles, dmin, davg, navg);  
    }//END OF IF cell i is empty
}//END OF i LOOP (THROUGH ALL CELLS)

        for( int i = 0; i < n; i++ ) 
            move( particles[i] );		


	update_bins(cells,n,particles);

	//if I want to do correctness checks and write particle positions...
        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          // 
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin; //update minimum radius
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }

    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

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
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particles are not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
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
