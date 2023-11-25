# Wave2D
Parallelizing Wave Diffusion with MPI and OpenMP


1. Purpose
In this programming assignment, we will parallelize a sequential version of a two-dimensional wave diffusion program, using a hybrid form of MPI and OpenMPI.
2. Schroedinger's Wave Dissemination
Assume the water surface in a two-dimensional square bucket. To simulate wave dissemination over this water surface, let’s partition this square in mesh and thus into N-by-N cells. Each cell(i, j) where 0 < i, j < N-1 maintains the height of its water surface. A wave is disseminated north, east, south, and west of each cell, and therefore cell(i, j) computes its new surface height from the previous height of itself and its four neighboring cells: cell(i+1, j), cell(i-1, j), cell(i, j+1) and cell(i, j-1). Let Zt_i,j, Zt-1_i,j, and Zt-2_i,j be the surface height of cell(i, j) at time t, time t-1, and time t-2 respectively. No wave at cell(i, j) at time t means Zt_i,j = 0.0. The water surface can go up and down between 20.0 and -20.0 through the wave dissemination.
Schroedinger’s wave formula computes Zt_i,j (where t >= 2 ) as follows:
Zt_i,j = 2.0 * Zt-1_i,j – Zt-2_i,j + c2 * (dt/dd)2 * (Zt-1_i+1,j + Zt-1_i-1,j + Zt-1_i,j+1 + Zt-1_i,j-1 – 4.0 * Zt-1_i,j)
where
c is the wave speed and should be set to 1.0,
dt is a time quantum for simulation, and should be set to 0.1, and
dd is a change of the surface, and should be set to 2.0.
Note that, if a cell is on an edge, (i.e., i = 0, i = N -1, j = 0, or j = N – 1), Zt_i,j = 0.0
The above formula does not work when t = 1. Zt_i,j (at t == 1) should be computed as: Zt_i,j = Zt-1_i,j + c2 / 2 * (dt/dd)2 * (Zt-1_i+1,j + Zt-1_i-1,j + Zt-1_i,j+1 + Zt-1_i,j-1 – 4.0 * Zt-1_i,j)
Note that, if a cell is on an edge, (i.e., i = 0, i = N -1, j = 0, or j = N – 1), Zt_i,j = 0.0
 i–1, j–1
 i, j–1
 i+1, j–1
 i–1, j
 i, j
 i+1, j
 i–1, j+1
 i, j+1
 i+1, j+1
How about t == 0? This is an initialization of the water surface. Let’s create a huge tidal wave in the middle of this square bucket. Set all cells(i, j) to 20.0 where 0.4 * N < i < 0.6 * N, 0.4 * N < j < 0.6 * N.
0.4N 0.6N
N-1
0 ..... 0.4N ~ 0.6N ..... N-1
  0
   Your simulation now starts with t == 0 (initialization), increments t by one, and computes the surface height of all cells(i, j), (i.e., Zt_i,j) at each time t, based on the above formulae (See examples of simulation outputs below).
   
   Look at ~css534/prog2/Wave2D_template.cpp. The main( ) function first reads three parameters: (1) size: the edge length of a 2D simulated square; (2) max_time: # steps to simulate wave dissemination where max_time >= 2; and (3) interval: # simulation steps needed each time the current simulation space is printed out, (e.g., interval == 1 means simulation status to be displayed every single step, interval == 2 prints out the space at time 2, 4, ..., whereas interval == 0 means no simulation output that is necessary to
remove any I/O overheads when measuring execution performance.)
int main( int argc, char *argv[] ) { // verify arguments
if ( argc != 4 ) {
cerr << "usage: Wave2D size max_time interval" << endl;
return -1; }
int size = atoi( argv[1] );
int max_time = atoi( argv[2] ); int interval = atoi( argv[3] );
if ( size < 100 || max_time < 3 || interval < 0 ) {
cerr << "usage: Wave2D size max_time interval" << endl;
cerr << " where size >= 100 && time >= 3 && interval >= 0" << endl; return -1;
}
Thereafter, the main( ) function creates a simulation space and starts a timer.
// create a simulation space
double z[3][size][size];
for ( int p = 0; p < 3; p++ )
for ( int i = 0; i < size; i++ ) for ( int j = 0; j < size; j++ )
z[p][i][j] = 0.0; // no wave
  // start a timer
  Timer time;
  time.start( );
After that, main( ) initializes the water surface, z[0][][] at time = 0.
// time = 0;
// initialize the simulation space: calculate z[0][][] int weight = size / default_size;
for( int i = 0; i < size; i++ ) {
for( int j = 0; j < size; j++ ) { if(i>40*weight&&i<60*weight &&
j > 40 * weight && j < 60 * weight ) { z[0][i][j] = 20.0;
      } else {
        z[0][i][j] = 0.0;
} }

}
We now have to simulate the wave diffusion at time = 1, 2, and all the way to max_time - 1. You must implement the rest of main( ) by yourself. Don't forget to insert the code to print out the simulation space every interval steps. The printing format should be:
t
z[t][0][0] z[t][1][0] ... z[t][size-1][0]
z[t][0][1] z[t][1][1] ... z[t][size-1][1]
...
z[t][0][size-1] z[t][1][size-1] ... z[t][size-1][size-1]
For example, given z[3][3][3], when t == 2, we can print as follows:
2
z[2][0][0] z[2][1][0] z[2][2][0] z[2][0][1] z[2][1][1] z[2][2][1] z[2][0][2] z[t][1][2] z[2][2][2]
Note that, when t == 3, we cannot print z[3][i][j]. This in turn means that you have to shift down values from z[2][][] to z[1][][] and from z[1][][] to z[0][][]. Don't copy values, which slows down the execution. Instead, rotate z[2][][], z[1][][], and z[0][][].
Running this program is not really impressive unless its standard outputs are graphically displayed. For
this purpose, use Wout.java. To see Wave2's outputs graphically, redirect the standard outputs to Wout.
For example, to see a 100 x 100 simulation space every 10 steps from t = 0 to 499, type:
Wave2D 100 500 10 | java Wout 100
3. Parallelization
Follow the parallelization strategies described below:
(1) CopyWave2D_template.cppintoWave2D.cpp,andcompleteitssequentialimplementation.Checkif
it works correctly, using Wout.java.
(2) Copy Wave2D.cpp into Wave2D_mpi.cpp. Start parallelization with MPI first. Divide each square
into small stripes along the i-axis. For instance, if you use four processes, z[3][100][100] should be
divided and allocated to different processors as follows:
rank 0: z[0][0][j] ~ z[0][24][j], z[1][0][j] ~ z[1][24][j], z[2][0][j] ~ z[2][24][j] rank 1: z[0][25][j] ~ z[0][49][j], z[1][25][j] ~ z[1][49][j], z[2][25][j] ~ z[2][49][j] rank 2: z[0][50][j] ~ z[0][74][j], z[1][50][j] ~ z[1][74][j], z[2][50][j] ~ z[2][74][j] rank 3: z[0][75][j] ~ z[0][99][j], z[1][75][j] ~ z[1][99][j], z[2][75][j] ~ z[2][99][j]
Note 0 <= j < 99. For simplicity, each rank may allocate an entire z[3][size][size] array but just use
only the above stripe.
(3) Ineachiterationofsimulationloopt=2throughtomax_time-1,yourorderofoperationsshouldbe
(a) printing out an intermediate simulation if necessary, (b) exchanging data among ranks, and (c)
computing Schroedinger's formula in parallel.
(4) Rank0isresponsibletoprintoutanintermediatestatustothestandardoutput.Forthispurpose,rank
0 must receive all strips from the other ranks 1 ~ 3 before printing out the status.
(5) Two neighboring ranks must exchange boundary data. For instance, rank 1 must send its z[p][25][j] to rank 0 as well as z[p][49][j] to rank 2. At the same time, rank 1 must receive z[p][24][j] from rank 0 as well as z[p][50][j] from rank 2. At time == 2, p is 1. However, beyond time == 2, p will repeatedly change into 0, 1, and back to 2. Note that rank 0 has no left neighbor and rank N -1 has no
right neighbor.
(6) Schroedinger's formula is the most computation intensive part that should be parallelized. Each rank
computes only its own stripe.
(7) AfterverifyingthecorrectnessofyourMPIparallelization,keepworkingonyourparallelizationwith
OpenMP. Focus on Schroedinger's formula.
