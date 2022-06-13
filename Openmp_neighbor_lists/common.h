#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;
const int MAX_PARTICLES_PER_CELL = 4;
const int MAX_NEIGHBORS=5;
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


typedef struct
{
  int particles[MAX_PARTICLES_PER_CELL];
  int E;
  int SE;
  int S;
  int SW;
  int W;
  int NW;
  int N;
  int NE;
} cell;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void set_cells_per_side( int cps );
void init_particles( int n, particle_t *p );
void make_cell_map();
void init_cells(cell *c, int nParticles, particle_t *p);
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );
void update_bins(cell *c, int nParticles, particle_t *p);
void sort_particles(cell *c, int nParticles, particle_t *p);
void compute_cell(int BCL, cell *cells, particle_t *particles, double &dmin, double &davg, int &navg);
double get_cutoff_value();
double get_size();
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
