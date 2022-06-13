#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#include<vector>
using namespace std;
inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }
//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
// particle data structure
//
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
int  set_size( int );
double get_binlen();
double  get_size();
int get_nbin();
void init_particles( int n, particle_t *p);
void sort_particles( int n, particle_t *p);
void init_bins(int n, particle_t *p, vector<vector<vector<particle_t> > > &bin);
void rearrange_bins( int n, particle_t *p, vector<vector<vector<particle_t> > > &bin );
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void apply_force_reduced( particle_t &, vector<vector<vector<particle_t> > > &, double *, double *, int *);
void move( particle_t &p,vector<vector<vector<particle_t> > > &bin);
void check_bins( int n, int nbin,  particle_t *p, vector<vector<vector<particle_t> > > &bin );
//MPI exclusive functions
typedef struct
{
int up;
int down;
}exlim;

int  set_size_mpi( int );
double get_binlen_mpi();
double  get_size_mpi();
int get_nbin_mpi();
void clear_bins_mpi(vector<vector<vector<particle_t*> > > &bin);
void init_bins_mpi(int n, vector<particle_t*> p, vector<vector<vector<particle_t*> > > &bin);
void apply_force_reduced_mpi( particle_t &, vector<vector<vector<particle_t*> > > &, double *, double *, int *);
void apply_force_mpi( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move_mpi( particle_t &p);
void init_particles_mpi( int n, particle_t *p);
void sort_particles_mpi( int n, particle_t *p);
void save_mpi( FILE *f, int n, vector<particle_t> &p );
//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
