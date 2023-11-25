  #include <iostream>
  #include "Timer.h"
  #include <stdlib.h>   // atoi
  #include <math.h>
  #include <mpi.h>
  #include <stdio.h>
  #include <sstream>
  #include <omp.h>

  int default_size = 100;  // the default system size
  int defaultCellWidth = 8;
  double c = 1.0;      // wave speed
  double dt = 0.1;     // time quantum
  double dd = 2.0;     // change in system

  using namespace std;

  int main( int argc, char *argv[] ) {
  
    int my_rank = 0;            // used by MPI
    int mpi_size;           // used by MPI
    // verify arguments
    if ( argc != 5 ) {
      cerr << "usage: Wave2D size max_time interval" << endl;
      return -1;
    }
    int size = atoi( argv[1] );
    int max_time = atoi( argv[2] );
    int interval  = atoi( argv[3] );
    int nThreads = atoi( argv[4] );
    

    if ( size < 100 || max_time < 3 || interval < 0 || nThreads <=0) {
      cerr << "usage: Wave2D size max_time interval" << endl;
      cerr << "       where size >= 100 && time >= 3 && interval >= 0" << endl;
      return -1;
    }


  //MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Bcast(&size,1,MPI_INT,0,MPI_COMM_WORLD);
  // change # of threads
    omp_set_num_threads(nThreads);


  // Calculate stripe boundaries
    int stripe = size / mpi_size;
    int remainder = size % mpi_size;
    int stripe_begins[mpi_size];      // stripe_begins[i] is rank i’s
    int stripe_ends[mpi_size];        // stripe_begins[i] is rank i’s
    
  for (int my_rank = 0; my_rank < mpi_size; my_rank++) {   
      if (my_rank < remainder) {       
          // If the rank is less than the remainder, it gets extra remainder size
          stripe_begins[my_rank] = stripe * my_rank + my_rank;
          stripe_ends[my_rank] =  stripe_begins[my_rank]+ stripe;
      } else {
          // For the last rank, it accounts for the remainder work
          stripe_begins[my_rank] = stripe * my_rank + remainder;
          stripe_ends[my_rank] = stripe_begins[my_rank] + stripe - 1;
      }   
  }


    // create a simulation space
    double z[3][size][size];
    for ( int p = 0; p < 3; p++ ) 
      for ( int i = 0; i < size; i++ )
        for ( int j = 0; j < size; j++ )
    z[p][i][j] = 0.0; // no wave

    // start a timer
    Timer time;
    time.start( );

    // time = 0;
    // initialize the simulation space: calculate z[0][][]
    int weight = size / default_size;
    for( int i = 0; i < size; i++ ) {
      for( int j = 0; j < size; j++ ) {
        if( i > 40 * weight && i < 60 * weight  &&
      j > 40 * weight && j < 60 * weight ) {
    z[0][i][j] = 20.0;
        } else {
    z[0][i][j] = 0.0;
        }
      }
    }
  
    // time = 1
    // calculate z[1][][] 
    // cells not on edge 
    // OMP parallelization 
    #pragma omp parallel for 
      for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == 0 || j == 0 || i == size - 1 || j == size - 1) {
                z[1][i][j] = 0;
            } else {
                z[1][i][j] = z[0][i][j] + ((c * c)/2 ) * pow(dt / dd, 2) * (z[0][i + 1][j] + z[0][i - 1][j] + z[0][i][j + 1] + z[0][i][j - 1] - 4.0 * z[0][i][j]);
            }
        }
    }
    

  for (int t = 2; t < max_time; t++) {
      int p = t % 3;
      int q = (t + 2) % 3;
      int r = (t + 1) % 3;

      //Even Odd communication for deadlock prevention
      //Even rank send , Odd rank receive
      
          if (my_rank % 2 == 0) {
            
              // Even-rank processes send data to their left and right neighbors .Left neighbor means my_rank-1 . Right neighbor means my_rank+1
        
              if (my_rank == 0 && my_rank !=mpi_size-1) {            
                  MPI_Send(&z[q][stripe_ends[my_rank]], size, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD);
                
              } else if (my_rank == mpi_size - 1) {
                  MPI_Send(&z[q][stripe_begins[my_rank]], size, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD);
                  
              } else {
                  MPI_Send(&z[q][stripe_begins[my_rank]], size, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD);
                  MPI_Send(&z[q][stripe_ends[my_rank]], size, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD);
                  
              }

          } else {
            
              // Odd-rank processes receive data from their right and left neighbors

              if (my_rank == 0 && my_rank !=mpi_size-1) {
                  MPI_Status status; 
                  MPI_Recv(&z[q][stripe_begins[my_rank+1]], size, MPI_DOUBLE, my_rank +1, 0, MPI_COMM_WORLD, &status);
                    
              } else if (my_rank == mpi_size - 1) {
                  MPI_Status status; 
                  MPI_Recv(&z[q][stripe_ends[my_rank-1]], size, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, &status);
                    
              } else {
                  MPI_Status status; 
                  MPI_Recv(&z[q][stripe_begins[my_rank+1]], size, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD, &status);
                  MPI_Recv(&z[q][stripe_ends[my_rank-1]], size, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, &status);
                
              
              }
          }
      
            //Even rank receive , Odd rank send
    
      if (my_rank % 2 == 0) {            
            // Even-rank processes receive data from their left and right neighbors

              if (my_rank == 0 && my_rank !=mpi_size-1) {
                  MPI_Status status; 
                  MPI_Recv(&z[q][stripe_begins[my_rank+1]], size, MPI_DOUBLE, my_rank +1, 0, MPI_COMM_WORLD, &status);
                    
              } else if (my_rank == mpi_size - 1) {
                  MPI_Status status; 
                  MPI_Recv(&z[q][stripe_ends[my_rank-1]], size, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, &status);
                    
              } else {
                  MPI_Status status; 
                  MPI_Recv(&z[q][stripe_ends[my_rank-1]], size, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, &status);
                  MPI_Recv(&z[q][stripe_begins[my_rank+1]], size, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD, &status);
              
              }
            
          } else {
            // Odd-rank processes send data to their right and left neighbors

              if (my_rank == 0 && my_rank !=mpi_size-1) {
                  MPI_Send(&z[q][stripe_ends[my_rank]], size, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD);
                
              } else if (my_rank == mpi_size - 1) {
                  MPI_Send(&z[q][stripe_begins[my_rank]], size, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD);
              
              } else {
                
                  MPI_Send(&z[q][stripe_ends[my_rank]], size, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD);
                  MPI_Send(&z[q][stripe_begins[my_rank]], size, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD);
                  
              }
          
          }
      

  if (interval != 0 && ( t % interval == 0)) {
        //Aggregate all results from all ranks
              if (my_rank == 0) {
                  for (int rank = 1; rank < mpi_size; rank++) {
                      MPI_Status status;
                      MPI_Recv(&z[p][stripe_begins[rank]], (stripe_ends[rank] - stripe_begins[rank] + 1) * size, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD,&status);
                      
              }
            cout << t << endl;
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                  cout << z[p][j][i] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }

        else {
          
            MPI_Send(&z[p][stripe_begins[my_rank]], (stripe_ends[my_rank] - stripe_begins[my_rank] + 1) * size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
          
            }
          
    }
        //Parallelization for the Schroedinger's formula
      
      int i;
      #pragma omp parallel for private (i)
      for (int i = stripe_begins[my_rank]; i <= stripe_ends[my_rank]; i++) {
            for (int j = 0; j < size; j++) {
                if (i == 0 || i == size - 1 || j == 0 || j == size - 1) {
                    z[p][i][j] = 0;
                } else {
                    z[p][i][j] = 2.0 * z[q][i][j] - z[r][i][j] + c * c * pow( dt/dd, 2) * (z[q][i + 1][j] + z[q][i - 1][j] + z[q][i][j + 1] + z[q][i][j - 1] - 4.0 * z[q][i][j]);
                }
            }
        }

  }
  

  
  // end of simulation   
  MPI_Finalize(); // shut down MPI
  // Print rank range
  std::stringstream range;

      if (my_rank == 0) {
          for (int rank = 0; rank < mpi_size; rank++) {
              range.str("");  // Clear the stringstream
              if (rank < mpi_size) {
                  range << stripe_begins[rank] << " ~ " << stripe_ends[rank];
              }

              if (rank < mpi_size) {
                  std::cerr << "rank[" << rank << "]'s range = " << range.str() << std::endl;
              }
          }
      }
  
    // finish the timer
    if (my_rank == 0) {
    cerr << "Elapsed time = " << time.lap( ) << endl;
    }
    return 0;
  }
