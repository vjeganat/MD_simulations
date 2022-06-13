#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"

double size;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
}



double get_size(){
return size;
}

double get_cutoff(){
	return cutoff;
}


int get_cells_per_side(){
return (int) (size/cutoff);
}


//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );
	
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    

    int *shuffle = (int*)malloc( n * sizeof(int) ); 
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;	
    for( int i = 0; i < n; i++ ) 
    {
     
	int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //	
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );  // annihilate shuffle. we're done with it.
}


void compute_cell(int home, cell_t *c, double &dmin, double &davg, int &navg){  

 
  particle_t temp_array[MAX_NEIGHBORS];
   int qp=0;	

   int cps = get_cells_per_side();
   int homex= home/cps;
   int homey= home%cps;
   for(int k = 0;k<c[home].particle_quantity;k++){
	c[home].particles[k].ax = c[home].particles[k].ay = 0;
	for(int shiftx =-1;shiftx<=1;shiftx++){
	    if((homex+shiftx>=0&&homex+shiftx<cps)){
		for(int shifty =-1;shifty<=1;shifty++){
		    if(homey+shifty>=0&&homey+shifty<cps){
			if(shiftx==0&&shifty==0){	
			    for(int l=0; l<c[home].particle_quantity;l++){
   if(k!=l){ temp_array[qp]= c[home].particles[l]; qp++; }
		   	    }
			}else{
			    for(int l=0; l<c[(homey+shifty)+(homex+shiftx)*cps].particle_quantity;l++){
   temp_array[qp] = c[(homey+shifty)+(homex+shiftx)*cps].particles[l]; qp++;
			    }
			}//if-else,(working within cell means avoiding particle forcing itself.
		    }//if - boundary conditions
		}// for shifty
	    } //if - boundary conditions
	}// for shiftx
for(int i=0; i<qp;i++){
apply_force(c[home].particles[k],temp_array[i],&dmin,&davg,&navg);
}
qp=0;

    }// for-k, all particles
	
}//end of method

void init_cells(cell_t *c, int nParticles , particle_t *p){
	int ncells_per_side = get_cells_per_side();
	for(int i=0;i<ncells_per_side*ncells_per_side;i++){
		c[i].particle_quantity = 0;	
		c[i].pq_within_move = 0;
	}	
	for(int i =0;i<nParticles;i++){
	int xpos= (p[i].x/size)*ncells_per_side;
	int ypos= (p[i].y/size)*ncells_per_side;
	int pos = xpos+ypos*ncells_per_side;
	c[pos].particles[c[pos].particle_quantity]=p[i];	
	c[pos].particle_quantity++;
	}
}

//
//  interact two particles
//
//
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{
    
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy; 
    if( r2 > cutoff*cutoff )
        return;
    if (r2 != 0)
    {
	if (r2/(cutoff*cutoff) < *dmin * (*dmin))
	   *dmin = sqrt(r2)/cutoff;
     	(*davg) += sqrt(r2)/cutoff;
     	(*navg) ++;
    }
		
    r2 = fmax( r2, min_r*min_r ); 
    double r = sqrt( r2 ); 
	
    //
    //  very simple short-range repulsive force
    double coef = ( 1 - cutoff / r ) / r2 / mass;  
    particle.ax += coef * dx;	//change acceleration in x
    particle.ay += coef * dy;	//change acceleration in y
}

//==============================================================================

//  integrate the ODE
//
int move(particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //  
    //  velocities changed by calculated accelerations
    //  positions immediately changed by those changed velocities.

    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //  bounce from walls
    //  reflect object position and velocities.
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }

    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }

	int ncells_per_side = get_cells_per_side();

	int xpos= (p.x/size)*ncells_per_side;
	int ypos= (p.y/size)*ncells_per_side;
	int pos = xpos+ypos*ncells_per_side;
	return pos;

}


//===============================================================================

//
//  I/O routines
//  prints to output file, the number of particles and the size,
//  followed by the positions of each particle, at the time this method is called.
//
void save( FILE *f, int n, cell_t *c )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    int nc = get_cells_per_side()*get_cells_per_side();
    for( int i = 0; i < nc; i++ )
	for(int j =0; j<c[i].particle_quantity; j++)
        fprintf( f, "%g %g\n", c[i].particles[j].x, c[i].particles[j].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
