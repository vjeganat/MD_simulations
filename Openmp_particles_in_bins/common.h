#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;
const int MAX_PARTICLES_PER_CELL = 5;
const int MAX_NEIGHBORS=7;
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
  particle_t particles[MAX_PARTICLES_PER_CELL];
  int particle_quantity;
  int pq_within_move;
} cell_t;



//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );
void init_cells(cell_t *c, int nParticles, particle_t *p);
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
int move(particle_t &p );
void update_bins(cell_t *c, int nParticles, particle_t *p);



void compute_cell(int home, cell_t *c, double &dmin, double &davg, int &navg);


double get_cutoff();
double get_size();
int get_cells_per_side();
int particle_quantity_check(int nc, cell_t *c);
//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, cell_t *c );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );



#endif
