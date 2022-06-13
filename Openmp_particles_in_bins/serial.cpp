#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
//#include <vector>
//using namespace std;

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

    const int cells_per_side = get_cells_per_side();

    cell_t *cells = (cell_t*) malloc(cells_per_side*cells_per_side*sizeof(cell_t));
   
    init_cells(cells, n, particles);
	
    free(particles);

//	printf("simulation box size: %f, cutoff radius: %f, number of bins: %d, %d\n",get_size(),get_cutoff(),cells_per_side,cells_per_side*cells_per_side);


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

for(int i =0;i<cells_per_side*cells_per_side;i++){
	cells[i].pq_within_move = cells[i].particle_quantity;
	if(cells[i].particle_quantity!=0){	
	compute_cell(i,cells,dmin,davg,navg);
	}
}

int SHIFT_SIZE=3;

for(int shifty =0;shifty <SHIFT_SIZE;shifty++){			//shiftx

   for(int shiftx=0;shiftx<SHIFT_SIZE;shiftx++){			//shifty

	for(int i =0;i<cells_per_side;i+=SHIFT_SIZE){		//row i --y

	     for(int j=0;j<cells_per_side;j+=SHIFT_SIZE){	//column j --x

		if(((i+shifty)<cells_per_side)&&((j+shiftx)<cells_per_side)){  		 //if shifted i,j are cells in domain

			int current_location= (i+shifty)+(j+shiftx)*cells_per_side;
	
		    for(int k =0; k<cells[current_location].pq_within_move;k++){
// move happens at next line.
		
			int new_particle_location = move(cells[current_location].particles[k]);		
//	cells[current_location].particles[k].ax = cells[current_location].particles[k].ay =  0;		
			if (new_particle_location != current_location){
				//First, move the particle to the correct neighbor location.
				cells[new_particle_location].particles[cells[new_particle_location].particle_quantity] = cells[current_location].particles[k];
				cells[new_particle_location].particle_quantity++;
				// No checks are performed here. Assumes moved to a neighbor cell.
				cells[current_location].particle_quantity--;
				cells[current_location].pq_within_move--;
				if(cells[current_location].pq_within_move>0)  cells[current_location].particles[k]=cells[current_location].particles[cells[current_location].pq_within_move];
				
				if(cells[current_location].particle_quantity>cells[current_location].pq_within_move)	cells[current_location].particles[cells[current_location].pq_within_move] = cells[current_location].particles[cells[current_location].particle_quantity];				
				k--;
			}
		    }// for k
		}// if within domain
             }//for j
	}//for i
    }// for shiftx
}// for shifty
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
              save( fsave, n, cells );
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
    free( cells );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
