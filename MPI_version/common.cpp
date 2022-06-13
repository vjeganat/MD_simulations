#include <stdlib.h>
#include<iostream>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include<vector>
#include<algorithm>
#include <bits/stdc++.h>
using namespace std;
double size;
int nbin; //number of bins
double  binlen;
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
int set_size( int n )
{
    size = sqrt( density * n );
    nbin = (int) (size/cutoff);
    binlen= size/nbin;
    return nbin;
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p)
{
    srand48( time( NULL ) );
   
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
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
    free( shuffle );
}
void sort_particles( int n, particle_t *p){
particle_t temp;
for (int i=0; i<n; i++)
	for(int j=0; j<n-1; j++){
	if (p[j+1].x<p[j].x){
	temp= p[j];
	p[j]=p[j+1];
	p[j+1]= temp;
	}
	}
}

double get_binlen(){
return binlen;
}
int get_nbin(){
return nbin;
}
double  get_size()
{
return size;
}

void init_bins( int n,  particle_t *p, vector<vector<vector<particle_t> > > &bin){
 for( int i = 0; i < n; i++ ) 
    {
	//binning done
	bin[((int)(p[i].x/binlen))][((int)(p[i].y/binlen))].push_back(p[i]);

    }


}



void rearrange_bins(int n, particle_t *p, vector<vector<vector<particle_t> > > &bin){
for( int i=0; i<nbin; i++)
	for(int j=0; j<nbin; j++)
    		bin[i][j].clear();
   for (int i=0; i<n; i++) 
    bin[((int)(p[i].x/binlen))][((int)(p[i].y/binlen))].push_back(p[i]);

}
void check_bins(int n, int nbin, particle_t *p,  vector<vector<vector<particle_t> > > &bin){
 vector<vector<vector<particle_t> > > cbin(nbin,vector<vector<particle_t> >(nbin,vector <particle_t>(0)));
cbin= bin;

//for( int i=0; i<nbin; i++)
//	for(int j=0; j<nbin; j++)
//   		bin[i][j].clear();

// for (int i=0; i<n; i++){
//     	bin[((int)(p[i].x/binlen))][((int)(p[i].y/binlen))].push_back(p[i]);
//	}
   for (int i=0; i<nbin; i++)
	for (int j=0; j<nbin; j++)
	if (cbin[i][j].size()!=bin[i][j].size())
	cout<<"binnig needs attention"<<endl;	
}
//
//  interact two particles
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
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

void apply_force_reduced( particle_t &particle, vector<vector<vector<particle_t> > > &bin, double *dmin, double *davg, int *navg)
{
 int xbin= (int)(particle.x/binlen);
 int ybin= (int)(particle.y/binlen); 
 for(int x=0; x<=2; x++){
 if ((xbin==0 && x==0) || (xbin==nbin-1 && x==2))
	continue;
	for( int y=0; y<=2; y++){
 	if ((ybin==0 &&y==0)|| (ybin==nbin-1 && y==2))
	continue;
	for(int par=0; par<bin[xbin-1+x][ybin-1+y].size(); par++){	
		apply_force( particle, bin[xbin-1+x][ybin-1+y][par],dmin,davg,navg);
		}
}
 }
}


//
//  integrate the ODE
//
void move( particle_t &p,vector<vector<vector<particle_t> > > &bin )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    particle_t prev= p;
    int xprevbin= (int)(prev.x/binlen);
    int yprevbin= (int)(prev.y/binlen);
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    
    //
    //  bounce from walls
    //
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
    int xcurbin= (int)(p.x/binlen);
    int ycurbin= (int)(p.y/binlen);
for (int find=0; find< bin[xprevbin][yprevbin].size(); find++){
if ((bin[xprevbin][yprevbin][find].x==prev.x) && (bin[xprevbin][yprevbin][find].y==prev.y)){
	bin[xcurbin][ycurbin].push_back(p);
	bin[xprevbin][yprevbin].erase(bin[xprevbin][yprevbin].begin()+find);
	    }
}

}

//
//  I/O routines
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
//MPI EXCLUSIVE ROUTINES
int set_size_mpi( int n )
{
    size = sqrt( density * n );
    nbin = (int) (size/cutoff);
    binlen= size/nbin;
    return nbin;
}

void sort_particles_mpi( int n, particle_t *p){
particle_t temp;
for (int i=0; i<n; i++)
	for(int j=0; j<n-1; j++){
	if (p[j+1].x<p[j].x){
	temp= p[j];
	p[j]=p[j+1];
	p[j+1]= temp;
	}
	}
}

double get_binlen_mpi(){
return binlen;
}
int get_nbin_mpi(){
return nbin;
}
double  get_size_mpi()
{
return size;
}

void init_bins_mpi( int n,  particle_t *p, vector<vector<vector<particle_t> > > &bin){
 for( int i = 0; i < n; i++ ) 
    {
	//binning done
	bin[((int)(p[i].x/binlen))][((int)(p[i].y/binlen))].push_back(p[i]);

    }


}
void init_particles_mpi( int n, particle_t *p)
{
    srand48( time( NULL ) );
   
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
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
    free( shuffle );
}




void init_bins_mpi( int n,  vector <particle_t*> p, vector<vector<vector<particle_t*> > > &bin){
 for( int i = 0; i < n; i++ ) 
    {
	//binning done
	bin[((int)((*p[i]).x/binlen))][((int)((*p[i]).y/binlen))].push_back(p[i]);
    }


}
void clear_bins_mpi(vector<vector<vector<particle_t*> > > &bin){
for( int i=0; i<nbin; i++)
	for(int j=0; j<nbin; j++)
    		bin[i][j].clear();

}

void move_mpi( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    
    //
    //  bounce from walls
    //
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
void apply_force_mpi( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
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
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

void apply_force_reduced_mpi( particle_t &particle, vector<vector<vector<particle_t*> > > &bin, double *dmin, double *davg, int *navg)
{
 int xbin= (int)(particle.x/binlen);
 int ybin= (int)(particle.y/binlen); 
 for(int x=0; x<=2; x++){
 if ((xbin==0 && x==0) || (xbin==nbin-1 && x==2))
	continue;
	for( int y=0; y<=2; y++){
 	if ((ybin==0 &&y==0)|| (ybin==nbin-1 && y==2))
	continue;
	for(int par=0; par<bin[xbin-1+x][ybin-1+y].size(); par++){	
		apply_force_mpi( particle, *bin[xbin-1+x][ybin-1+y][par],dmin,davg,navg);
		}
}
 }
}

void save_mpi( FILE *f, int n, vector<particle_t> &p )
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

