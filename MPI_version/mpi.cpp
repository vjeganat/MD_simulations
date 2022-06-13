#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include <vector>
#include<iostream>
using namespace std;
//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg;
        //
    //  process command line parameters
    //
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
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Status status;
    MPI_Status status1;
    MPI_Status status2;
    MPI_Status status3;
    MPI_Request request;

    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;


   particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    int nbin = set_size_mpi(n);
    double  size = get_size_mpi();
    double  binlen= get_binlen_mpi();
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    //  set up the data partitioning across processors
     
    int *particle_per_proc = (int*) malloc (n_proc* sizeof(int));
    int *bin_per_proc = (int*) malloc (n_proc* sizeof(int));
    int *bin_partition_x = (int*) malloc ((n_proc+1)* sizeof(int));
    exlim *ex_bin_partition_x = (exlim*) malloc ((n_proc)* sizeof(exlim));

	if( rank == 0 ){
        init_particles_mpi( n, particles );
	sort_particles_mpi(n, particles);
	}

//Calcualting number of bins per process
	int temp= (int)(nbin/n_proc);
	if (n_proc>1){
	while(1)
	{
	if ((nbin-temp)%(n_proc-1)==0)
		break;
	else
		temp++;
	}
	bin_per_proc[n_proc-1]=temp;
	for (int i=0; i< n_proc-1; i++)
		bin_per_proc[i]= (nbin-temp)/(n_proc-1);
	}
	else
	bin_per_proc[0]=temp;


//Calculating number the partition of bins
	for (int i= 0; i< n_proc +1; i++){
	bin_partition_x[i]=0;
		for (int j=0; j<i; j++){
				if (i==0)
					continue;
				else
					bin_partition_x[i]+= bin_per_proc[j];
	}
	}
//`Calculating extended number of partitions in bins
//
	for (int i= 0; i< n_proc; i++){
	if (i==0)
		ex_bin_partition_x[0].up=0;
	else
	ex_bin_partition_x[i].up=bin_partition_x[i]-1;
	if (i==n_proc-1)
	ex_bin_partition_x[i].down=bin_partition_x[i+1];
	else
	ex_bin_partition_x[i].down=bin_partition_x[i+1]+1;

	}

	if (rank ==0){
	for (int i=0; i<n_proc; i++)
      {   
          
          int count=0;
          for (int j=0; j<n; j++){
          if ((particles[j].x <=bin_partition_x[i+1]*binlen ) && (particles[j].x> bin_partition_x[i]*binlen ))
                  count ++;
          }
          particle_per_proc[i]=count;

      }
     }
    if(n_proc>1)
    MPI_Bcast(particle_per_proc, n_proc, MPI_INT, 0, MPI_COMM_WORLD);
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
        for( int i = 0; i < n_proc+1; i++ ){
	partition_offsets[i]=0;
	if (i==0){
		continue;
		}
	
	for(int j=0; j<i; j++){
        partition_offsets[i] +=  particle_per_proc[j];
	}
    }
//    if (rank==0)
//		cout<<"no of bins"<<nbin<<endl;
//    if (rank==0)
//	for (int i=0; i<n_proc; i++)
//		cout<<"particle_per_proc "<<particle_per_proc[i]<<endl;
//    if (rank==0)
//	for (int i=0; i<n_proc+1; i++)
//		cout<<"particles_offset "<<partition_offsets[i]<<endl; 
//     if (rank==0)
//	for (int i=0; i<n_proc+1; i++)
//		cout<<"bin_ partitions "<<bin_partition_x[i]<<endl;
//     if (rank==0)
//         for (int i=0; i<n_proc; i++)
//                 cout<<"extended bin partition"<< ex_bin_partition_x[i].up<< ex_bin_partition_x[i].down<<endl;
//     if (rank==0)
//         for (int i=0; i<n_proc; i++)
//                 cout<<"bin_per_proc"<<bin_per_proc[i]<<endl;
//     if (rank==0)
//	cout<<"Total particles:"<< n<<endl;
//     if(rank==0)
//	cout<<"Length of bin is"<<binlen<<endl;
//     if (rank==0)
//	cout<<"Size of domain"<<size<<endl;
	
    //
    //  allocate storage for local partition
    //
    int nlocal = particle_per_proc[rank];
    int recvdownamt=0;
    int  recvupamt=0;
   particle_t *local = (particle_t*) malloc( nlocal * sizeof(particle_t) );
   particle_t* send= (particle_t*) malloc( n * sizeof(particle_t) );
   particle_t* recvup=(particle_t*) malloc( n * sizeof(particle_t) );
   particle_t* recvdown=(particle_t*) malloc( n * sizeof(particle_t) );
   vector < particle_t*> exlocal(0);
   vector < particle_t> vsend(0);
   vector < particle_t> vrecvup(0);
   vector < particle_t> vrecvdown(0);
   vector < particle_t> vlocal(0);
   vector<vector<vector<particle_t*> > > localbin(nbin,vector<vector<particle_t*> >(nbin,vector <particle_t*>(0)));
   MPI_Scatterv( particles, particle_per_proc, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
   for (int i=0; i<nlocal; i++)
	vlocal.push_back(local[i]);
   //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    	    
    //
    //  simulate a number of time steps
    //

   double simulation_time = read_timer( );
        
   for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;
    // Sending rank to rank-1. 0 won't send anything 
   if (n_proc>1){
	if (rank>0)	
	for (int i=0; i< vlocal.size(); i++){
			if ((vlocal[i].x> bin_partition_x[rank]*binlen)&& (vlocal[i].x< ex_bin_partition_x[rank-1].down*binlen))
			vsend.push_back(vlocal[i]);}
	send= &vsend[0];
		if(rank%2==1){
	MPI_Send(send, vsend.size(), PARTICLE, (rank+n_proc-1)%n_proc, 4*step, MPI_COMM_WORLD);
	MPI_Probe((rank+1)%n_proc, 4*step, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, PARTICLE, &recvdownamt);
	MPI_Recv(recvdown, recvdownamt, PARTICLE, (rank+1)%n_proc, 4*step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for(int i=0; i<recvdownamt; i++)
		vrecvdown.push_back(recvdown[i]);
	}

	else{
	MPI_Probe((rank+1)%n_proc, 4*step, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, PARTICLE, &recvdownamt);
	MPI_Recv(recvdown, recvdownamt, PARTICLE, (rank+1)%n_proc, 4*step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for(int i=0; i<recvdownamt; i++)
		vrecvdown.push_back(recvdown[i]);
	MPI_Send(send, vsend.size(), PARTICLE, (rank+n_proc-1)%n_proc, 4*step, MPI_COMM_WORLD);
	}
	

	vsend.resize(0);
//Sending particles from rank to rank+1; rank+1 wont' send anything
	
	if (rank< n_proc -1)
	for (int i=0; i< vlocal.size(); i++)
			if ((vlocal[i].x<bin_partition_x[rank+1]*binlen) && (vlocal[i].x> ex_bin_partition_x[rank+1].up*binlen))
			vsend.push_back(vlocal[i]);

	send = &vsend[0];

	if(rank%2==1){
	MPI_Send(send, vsend.size(), PARTICLE, (rank+1)%n_proc, 4*step+ 1, MPI_COMM_WORLD);
	MPI_Probe((rank+n_proc-1)%n_proc, 4*step+ 1, MPI_COMM_WORLD,  &status1);
	MPI_Get_count(&status1, PARTICLE, &recvupamt);
	MPI_Recv(recvup, recvupamt, PARTICLE, (rank+n_proc-1)%n_proc, 4*step+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for(int i=0; i<recvupamt; i++)
		vrecvup.push_back(recvup[i]);
	}

	else{
	MPI_Probe((rank+n_proc-1)%n_proc, 4*step+1, MPI_COMM_WORLD, &status1);
	MPI_Get_count(&status1, PARTICLE, &recvupamt);	
	MPI_Recv(recvup, recvupamt, PARTICLE, (rank+n_proc-1)%n_proc, 4*step+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for(int i=0; i<recvupamt; i++)
		vrecvup.push_back(recvup[i]);
	MPI_Send(send, vsend.size(), PARTICLE, (rank+1)%n_proc, 4*step+1, MPI_COMM_WORLD);
	}
	}

	for (int i=0; i<vlocal.size(); i++)
		exlocal.push_back(&vlocal[i]);
	for(int i=0; i<recvdownamt; i++)
		exlocal.push_back(&vrecvdown[i]);
	for(int i=0; i<recvupamt; i++)
		exlocal.push_back(&vrecvup[i]);

	vsend.resize(0);

	
	clear_bins_mpi(localbin);			
	init_bins_mpi(exlocal.size(), exlocal, localbin);

       if( find_option( argc, argv, "-no" ) == -1 )
          if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
        
        //
        //  compute all forces
        //
        for( int i = 0; i < vlocal.size(); i++ )
        {
            vlocal[i].ax = vlocal[i].ay = 0;
	    apply_force_reduced_mpi( vlocal[i], localbin, &dmin, &davg, &navg );
        }
	        

	if( find_option( argc, argv, "-no" ) == -1 )
        {
          if (n_proc>1){
          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
          }
	else{
	rdavg=davg;
	rnavg=navg;
	rdmin=dmin;
	}
	if (rank == 0){
            //
            // Computing statistical data
            //
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin) absmin = rdmin;
          }
        }

        //
        //  move particles
        //
        for( int i = 0; i < vlocal.size(); i++ )
            move_mpi( vlocal[i] );

	vrecvup.resize(0);
	vrecvdown.resize(0);
   if(n_proc>1){
	if (rank>0) 
	for (int i=0; i<vlocal.size(); i++)
		if (vlocal[i].x< bin_partition_x[rank]*binlen){
		vsend.push_back(vlocal[i]);
		vlocal.erase(vlocal.begin()+i);
		i--;		
		}
	send= &vsend[0];

	if(rank%2==1){
	MPI_Send(send, vsend.size(), PARTICLE, (rank+n_proc-1)%n_proc,4*step + 2, MPI_COMM_WORLD);
	MPI_Probe((rank+1)%n_proc, 4*step+ 2, MPI_COMM_WORLD, &status2);
	MPI_Get_count(&status2, PARTICLE, &recvdownamt);
	MPI_Recv(recvdown, recvdownamt, PARTICLE, (rank+1)%n_proc, 4*step+ 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for(int i=0; i<recvdownamt; i++){
		vrecvdown.push_back(recvdown[i]);}
	}

	else{
	MPI_Probe((rank+1)%n_proc,4*step+ 2, MPI_COMM_WORLD, &status2);
	MPI_Get_count(&status2, PARTICLE, &recvdownamt);
	MPI_Recv(recvdown, recvdownamt, PARTICLE, (rank+1)%n_proc, 4*step+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for(int i=0; i<recvdownamt; i++)
		vrecvdown.push_back(recvdown[i]);
	MPI_Send(send, vsend.size(), PARTICLE, (rank+n_proc-1)%n_proc, 4*step+2, MPI_COMM_WORLD);
	}

	vsend.resize(0);

//Sending particles from rank to rank+1; rank+1 wont' send anything	

	if (rank< n_proc -1)
		 for (int i=0; i<vlocal.size(); i++)
                 	if (vlocal[i].x> bin_partition_x[rank+1]*binlen){
                 	vsend.push_back(vlocal[i]);
                 	vlocal.erase(vlocal.begin()+i);
			i--;
                 	}

	send = &vsend[0];

	if(rank%2==1){
	MPI_Send(send, vsend.size(), PARTICLE, (rank+1)%n_proc,4*step+ 3, MPI_COMM_WORLD);
	MPI_Probe((rank+n_proc-1)%n_proc,4*step+ 3, MPI_COMM_WORLD,  &status3);
	MPI_Get_count(&status3, PARTICLE, &recvupamt);
	MPI_Recv(recvup, recvupamt, PARTICLE, (rank+n_proc-1)%n_proc,4*step+ 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for(int i=0; i<recvupamt; i++)
		vrecvup.push_back(recvup[i]);
	}

	else{
	MPI_Probe((rank+n_proc-1)%n_proc,4*step+ 3, MPI_COMM_WORLD, &status3);
	MPI_Get_count(&status3, PARTICLE, &recvupamt);
	MPI_Recv(recvup, recvupamt, PARTICLE, (rank+n_proc-1)%n_proc, 4*step +3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for(int i=0; i<recvupamt; i++)
		vrecvup.push_back(recvup[i]);
	MPI_Send(send, vsend.size(), PARTICLE, (rank+1)%n_proc,4*step+ 3, MPI_COMM_WORLD);
	}

	for(int i=0; i<recvdownamt; i++)
		vlocal.push_back(vrecvdown[i]);
	for(int i=0; i<recvupamt; i++)
		vlocal.push_back(vrecvup[i]);
}	
	vsend.resize(0);
	vrecvdown.resize(0);
	vrecvup.resize(0);
	exlocal.resize(0);
}
    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
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
      if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
      if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
      }
      printf("\n");     
        
      //  
      // Printing summary data
      //  
      if( fsum)
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
  
    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    free( partition_offsets );
    free( local );
    free( particles );

    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
