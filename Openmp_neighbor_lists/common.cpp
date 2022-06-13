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
int cells_per_side;
int *cell_map;

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

void set_cells_per_side( int cps){
     cells_per_side = cps;
}

double get_size(){
return size;
}

double get_cutoff_value(){
	return cutoff;
}

int zorder(int x, int y, int sizex, int sizey, int index){
if(sizex==1||sizey==1){
        if(sizex==1&&sizey==1){
        cell_map[index]=y+x*cells_per_side;
        index++;
        }else if(sizey==1){
                if(sizex%2==1){
                	int ssx = sizex/2;
                	int sbx = (sizex+1)/2;
                	index = zorder(x,y,sbx,1,index);
                	index = zorder(x+sbx,y,ssx,1,index);
                }else{
                	int ssx = sizex/2;
                	index = zorder(x,y,ssx,1,index);
                	index = zorder(x+ssx,y,ssx,1,index);
                }
        }else if(sizex==1){
                if(sizey%2==1){
                	int ssy = sizey/2;
                	int sby = (sizey+1)/2;
                	index = zorder(x,y,1,sby,index);
                	index = zorder(x,y+sby,1,ssy,index);
                }else{
                	int ssy = sizey/2;
                	index = zorder(x,y,1,ssy,index);
                	index = zorder(x,y+ssy,1,ssy,index);
                }
        }
}else{
	if(sizex%2==1&&sizey%2==1){
		int ssx = sizex/2;
		int ssy = sizey/2;
		int sbx = (sizex+1)/2;
		int sby = (sizey+1)/2;
 		index = zorder(x,y,sbx,sby,index);
 		index = zorder(x,y+sby,sbx,ssy,index);
 		index = zorder(x+sbx,y,ssx,sby,index);
 		index = zorder(x+sbx,y+sby,ssx,ssy,index);
	}else if(sizex%2==1){
		int ssx = sizex/2;
		int sbx = (sizex+1)/2;
		int ssy = sizey/2;
 		index = zorder(x,y,sbx,ssy,index);
 		index = zorder(x,y+ssy,sbx,ssy,index);
 		index = zorder(x+sbx,y,ssx,ssy,index);
 		index = zorder(x+sbx,y+ssy,ssx,ssy,index);
	}else if(sizey%2==1){
		int ssx = sizex/2;
		int ssy = sizey/2;
		int sby = (sizey+1)/2;
 		index = zorder(x,y,ssx,sby,index);
 		index = zorder(x,y+sby,ssx,ssy,index);
 		index = zorder(x+ssx,y,ssx,sby,index);
 		index = zorder(x+ssx,y+sby,ssx,ssy,index);
	}else{
		int sx = sizex/2;
		int sy = sizey/2;
		index = zorder(x,y,sx,sy,index);
 		index = zorder(x,y+sy,sx,sy,index);
 		index = zorder(x+sx,y,sx,sy,index);
 		index = zorder(x+sx,y+sy,sx,sy,index);
	}
}
return index;

}//zorder end


void make_cell_map(){
cell_map = (int*)malloc(cells_per_side*sizeof(int));
zorder(0,0,cells_per_side,cells_per_side,0);
}






//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );
	
// sx and sy are the division lengths of the sides of the box.
// The two lines below divide the box into enough points, such that 
// the n particles can fit, while taking care to make sure that the box is 
// as square as possible. 
// I believe that sx = sy, or if it cannot, sx = sy + 1
// eg for n = 37,38,...42, sx=6,sy=6 isn't enough, so sx=7,sy=6.  
    
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    

    int *shuffle = (int*)malloc( n * sizeof(int) ); //this creates the space for the shuffle array.
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;			//shuffle array values initialized as their own location in array
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        //Next three lines:
        // random value chosen from array, from first n-i elements.
        // That value is held for assignment later in this loop.
        // value at end of array is copied to the location 
        // of the value just picked.
        // This ensures that the position of each particle is
        // not related to where it is stored. 
        // (The spacial locality in memory isn't indicative of spacial locality in simulation)
        //
     
	int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        //Next two lines:
        //
        //Each particle's position variables x, y are assigned, based on the location picked above.
        //
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //	
        //  assign random velocities within a bound
        //	vx,vy are assigned velocities between -1.0 and 1.0
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );  // annihilate shuffle. we're done with it.
}


void init_cells(cell *c, int nParticles, particle_t *p ){
//column major format

//for every row i
for(int i =0; i< cells_per_side;i++){
	//for every column j
	for(int j =0;j<cells_per_side;j++){
	if(j+1==cells_per_side&&i+1==cells_per_side){
		c[i+j*cells_per_side].E = NULL;
		c[i+j*cells_per_side].NE = NULL;
		c[i+j*cells_per_side].SE = NULL;
		c[i+j*cells_per_side].S = NULL;
		c[i+j*cells_per_side].SW = NULL;
		c[i+j*cells_per_side].W = i+(j-1)*cells_per_side;
		c[i+j*cells_per_side].NW = i-1+(j-1)*cells_per_side;
		c[i+j*cells_per_side].N = i-1+j*cells_per_side;
	}else if(j+1==cells_per_side&&i==0){
		c[i+j*cells_per_side].E = NULL;
		c[i+j*cells_per_side].NE = NULL;
		c[i+j*cells_per_side].SE = NULL;
		c[i+j*cells_per_side].S = i+1+j*cells_per_side;
		c[i+j*cells_per_side].SW = i+1+(j-1)*cells_per_side;
		c[i+j*cells_per_side].W = i+(j-1)*cells_per_side;
		c[i+j*cells_per_side].NW = NULL;
		c[i+j*cells_per_side].N = NULL;
	}else if(j==0&&i+1==cells_per_side){
		c[i+j*cells_per_side].E = i+(j+1)*cells_per_side;
		c[i+j*cells_per_side].NE = i-1+(j+1)*cells_per_side ;
		c[i+j*cells_per_side].SE = NULL;
		c[i+j*cells_per_side].S = NULL;
		c[i+j*cells_per_side].SW = NULL;
		c[i+j*cells_per_side].W = NULL;
		c[i+j*cells_per_side].NW = NULL;
		c[i+j*cells_per_side].N = i-1+j*cells_per_side;
	}else if(j==0&&i==0){
		c[i+j*cells_per_side].E = i+(j+1)*cells_per_side;
		c[i+j*cells_per_side].NE = NULL;
		c[i+j*cells_per_side].SE = i+1+(j+1)*cells_per_side;
		c[i+j*cells_per_side].S = i+1+j*cells_per_side;
		c[i+j*cells_per_side].SW = NULL;
		c[i+j*cells_per_side].W = NULL;
		c[i+j*cells_per_side].NW = NULL;
		c[i+j*cells_per_side].N = NULL;
//just the sides remaining
	}else if(j+1==cells_per_side){
		c[i+j*cells_per_side].E = NULL;
		c[i+j*cells_per_side].NE = NULL;
		c[i+j*cells_per_side].SE = NULL;
		c[i+j*cells_per_side].S = i+1+j*cells_per_side;
		c[i+j*cells_per_side].SW = i+1+(j-1)*cells_per_side;
		c[i+j*cells_per_side].W = i+(j-1)*cells_per_side;
		c[i+j*cells_per_side].NW = i-1+(j-1)*cells_per_side;
		c[i+j*cells_per_side].N = i-1+j*cells_per_side;
	}else if(j==0){
		c[i+j*cells_per_side].E = i+(j+1)*cells_per_side;
		c[i+j*cells_per_side].NE = i-1+(j+1)*cells_per_side;
		c[i+j*cells_per_side].SE = i+1+(j+1)*cells_per_side;
		c[i+j*cells_per_side].S = i+1+j*cells_per_side;
		c[i+j*cells_per_side].SW = NULL;
		c[i+j*cells_per_side].W = NULL;
		c[i+j*cells_per_side].NW = NULL;
		c[i+j*cells_per_side].N = i-1+j*cells_per_side;
	}else if(i==0){
		c[i+j*cells_per_side].E = i+(j+1)*cells_per_side;
		c[i+j*cells_per_side].NE = NULL;
		c[i+j*cells_per_side].SE = i+1+(j+1)*cells_per_side;
		c[i+j*cells_per_side].S = i+1+j*cells_per_side;
		c[i+j*cells_per_side].SW = i+1+(j-1)*cells_per_side;
		c[i+j*cells_per_side].W = i+(j-1)*cells_per_side;
		c[i+j*cells_per_side].NW = NULL;
		c[i+j*cells_per_side].N = NULL;
	}else if(i+1==cells_per_side){
		c[i+j*cells_per_side].E = i+(j+1)*cells_per_side;
		c[i+j*cells_per_side].NE = i-1+(j+1)*cells_per_side;
		c[i+j*cells_per_side].SE = NULL;
		c[i+j*cells_per_side].S = NULL;
		c[i+j*cells_per_side].SW = NULL;
		c[i+j*cells_per_side].W = i+(j-1)*cells_per_side;
		c[i+j*cells_per_side].NW = i-1+(j-1)*cells_per_side;
		c[i+j*cells_per_side].N = i-1+j*cells_per_side;
	}else{	//Except this one. This is for all of the bins not touching a wall.
		//rewrite to place this at beginning of if-else chain?
		c[i+j*cells_per_side].E = i+(j+1)*cells_per_side;
		c[i+j*cells_per_side].NE = i-1+(j+1)*cells_per_side;
		c[i+j*cells_per_side].SE = i+1+(j+1)*cells_per_side;
		c[i+j*cells_per_side].S =  i+1+j*cells_per_side;
		c[i+j*cells_per_side].SW = i+1+(j-1)*cells_per_side;
		c[i+j*cells_per_side].W = i+(j-1)*cells_per_side;
		c[i+j*cells_per_side].NW = i-1+(j-1)*cells_per_side;
		c[i+j*cells_per_side].N = i-1+j*cells_per_side;
	}
	}
}

// Add the particle numbers to the corresponding bins.
	
	for(int i=0;i<cells_per_side*cells_per_side;i++){
		for (int k =0; k<MAX_PARTICLES_PER_CELL;k++){
		   c[i].particles[k] = -1;
		}	
	}	
	for(int i =0;i<nParticles;i++){
	int xcell = (p[i].x/size) *(cells_per_side);
	int ycell = (p[i].y/size) *(cells_per_side);
		for(int k=0;k<MAX_PARTICLES_PER_CELL;k++){
		   if(c[xcell+ycell*cells_per_side].particles[k]==-1){
			c[xcell+ycell*cells_per_side].particles[k]=i;
			break;
		   }	
		}	
	}
}

void update_bins(cell *c, int nParticles, particle_t *p){
	for(int i=0;i<cells_per_side*cells_per_side;i++){
		for (int k =0; k<MAX_PARTICLES_PER_CELL;k++){
		   c[i].particles[k] = -1;
		}	
	}	
	for(int i =0;i<nParticles;i++){
	int xcell = (p[i].x/size) *(cells_per_side);
	int ycell = (p[i].y/size) *(cells_per_side);
		for(int k=0;k<MAX_PARTICLES_PER_CELL;k++){
		   if(c[xcell+ycell*cells_per_side].particles[k]==-1){
			c[xcell+ycell*cells_per_side].particles[k]=i;
			break;
		   }	
		}
	}
}


void compute_cell(int BCL, cell *cells, particle_t *particles, double &dmin, double &davg, int &navg){  
	int temp_particle_array[MAX_NEIGHBORS];
	for(int i =0;i<MAX_NEIGHBORS;i++) temp_particle_array[i]=-1;
	int qp = 0;  //quantity of particles in temp array
	for(int k = 0;k<MAX_PARTICLES_PER_CELL&&cells[BCL].particles[k]!=-1;k++){
		for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[BCL].particles[l]!=-1;l++){
	      	  if(k!=l){ temp_particle_array[qp] = cells[BCL].particles[l]; qp++;}
	   	}
		if(cells[BCL].E!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].E].particles[l]!=-1;l++){
		     temp_particle_array[qp] = cells[cells[BCL].E].particles[l]; qp++;
		  }
		}
		if(cells[BCL].SE!=NULL){ 
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].SE].particles[l]!=-1;l++){
		     temp_particle_array[qp] = cells[cells[BCL].SE].particles[l]; qp++;
		  }
		}
		if(cells[BCL].S!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].S].particles[l]!=-1;l++){
		     temp_particle_array[qp] = cells[cells[BCL].S].particles[l]; qp++;
		  }
		}
		if(cells[BCL].SW!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].SW].particles[l]!=-1;l++){
		     temp_particle_array[qp] = cells[cells[BCL].SW].particles[l]; qp++;
		  }
		}
		if(cells[BCL].W!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].W].particles[l]!=-1;l++){
		     temp_particle_array[qp] = cells[cells[BCL].W].particles[l]; qp++;
		  }
		}
		if(cells[BCL].NW!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].NW].particles[l]!=-1;l++){
		     temp_particle_array[qp] = cells[cells[BCL].NW].particles[l]; qp++;
		  }
		}
		if(cells[BCL].N!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].N].particles[l]!=-1;l++){
		     temp_particle_array[qp] = cells[cells[BCL].N].particles[l]; qp++;
		  }
		}
		if(cells[BCL].NE!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].NE].particles[l]!=-1;l++){
		     temp_particle_array[qp] = cells[cells[BCL].NE].particles[l]; qp++;
		  }
		}
	int kloc = cells[BCL].particles[k];
	for(int i =0; i<qp;i++){
	apply_force(particles[kloc],particles[temp_particle_array[i]],&dmin,&davg,&navg);
	}
	qp = 0;	//<<,^^: reset array to -1's
	}	
}//end of method


// This is much slower than above. collecting list before calculating seems faster than calculating as particles are encountered.
/*
void compute_cell(int BCL, cell *cells, particle_t *particles, double &dmin, double &davg, int &navg){  
	for(int k = 0;k<MAX_PARTICLES_PER_CELL&&cells[BCL].particles[k]!=-1;k++){
	int kloc = cells[BCL].particles[k];
		for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[BCL].particles[l]!=-1;l++){
		apply_force(particles[kloc],particles[cells[BCL].particles[l]],&dmin,&davg,&navg);
	   	}
		if(cells[BCL].E!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].E].particles[l]!=-1;l++){
		apply_force(particles[kloc],particles[cells[cells[BCL].E].particles[l]],&dmin,&davg,&navg);
		  }
		}
		if(cells[BCL].SE!=NULL){ 
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].SE].particles[l]!=-1;l++){
		apply_force(particles[kloc],particles[cells[cells[BCL].SE].particles[l]],&dmin,&davg,&navg);
		  }
		}
		if(cells[BCL].S!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].S].particles[l]!=-1;l++){
		apply_force(particles[kloc],particles[cells[cells[BCL].S].particles[l]],&dmin,&davg,&navg);
		  }
		}
		if(cells[BCL].SW!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].SW].particles[l]!=-1;l++){
		apply_force(particles[kloc],particles[cells[cells[BCL].SW].particles[l]],&dmin,&davg,&navg);
		  }
		}
		if(cells[BCL].W!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].W].particles[l]!=-1;l++){
		apply_force(particles[kloc],particles[cells[cells[BCL].W].particles[l]],&dmin,&davg,&navg);
		  }
		}
		if(cells[BCL].NW!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].NW].particles[l]!=-1;l++){
		apply_force(particles[kloc],particles[cells[cells[BCL].NW].particles[l]],&dmin,&davg,&navg);
		  }
		}
		if(cells[BCL].N!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].N].particles[l]!=-1;l++){
		apply_force(particles[kloc],particles[cells[cells[BCL].N].particles[l]],&dmin,&davg,&navg);
		  }
		}
		if(cells[BCL].NE!=NULL){
		  for(int l=0;l<MAX_PARTICLES_PER_CELL&&cells[cells[BCL].NE].particles[l]!=-1;l++){
		apply_force(particles[kloc],particles[cells[cells[BCL].NE].particles[l]],&dmin,&davg,&navg);
		  }
		}
	}	
}//end of method
// */

//
//  interact two particles
//
//	***Explicitly state what dmin, davg, navg mean***
//
//
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{
	// next three lines calculate distance from neighbor.
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
		
    r2 = fmax( r2, min_r*min_r ); // makes sure radius isn't too large.
    double r = sqrt( r2 ); //division is done after checks, for efficiency
	
    //
    //  very simple short-range repulsive force
    double coef = ( 1 - cutoff / r ) / r2 / mass;  
    particle.ax += coef * dx;	//change acceleration in x
    particle.ay += coef * dy;	//change acceleration in y
 
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //  
    //  velocities changed by calculated accelerations
    //  positions immediately changed by those changed velocities.
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
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

}

//
//  I/O routines
//  prints to output file, the number of particles and the size,
//  followed by the positions of each particle, at the time this method is called.
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
//  initialization, options for particle number, etc.
//
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
