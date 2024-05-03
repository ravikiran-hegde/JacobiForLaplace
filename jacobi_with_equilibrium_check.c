#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>

#define TOLERANCE 0.00000001
/*** function declarations ***/

// save matrix to file
void save_gnuplot( double *M, int dim );

// evolve Jacobi
void evolve( double * matrix, double *matrix_new, int n_loc, int dimension);

// return the elapsed time
double seconds( void );

/*** end function declaration ***/

int global_i(int i, int me, int n_loc, int offset){
        int i_g = i + me*n_loc + offset;
        return i_g;
}

void print_matrix( double * mat, int n_loc, int dimension ){

  int i, j;

  for( i = 0; i < n_loc + 2; i++ ){
    for( j = 0; j < dimension + 2; j++ ){
      printf("%.3g ", mat[ i *( dimension + 2) + j ] );
    }
    printf("\n");
  }
}

int matricesEqual(double *matrix1, double *matrix2, int n_loc, int dimension) {
    for (int i = 1; i < (n_loc + 1); i++) {
      for (int j = 1; j < (dimension + 1); j++ ){
        if (fabs( matrix1[i*(dimension + 2) + j ] - matrix2[i*(dimension + 2) + j ]) > TOLERANCE) {
            return 0; // Matrices are not equal
        }
      }
    }
    return 1; // Matrices are equal
}


int main(int argc, char* argv[]){

  // timing variables
  double t_start, t_end, increment;

  // indexes for loops
  int i, j, it, count, i_g;

  // initialize matrix
  double *matrix, *matrix_new, *tmp_matrix, *tmpry_matrix;

  int me, npes, n_loc, offset = 0, rest = 0, n_loc2;

  int dimension = 0, iterations = 0, row_peek = 0, col_peek = 0, freq = 0;
  int byte_dimension = 0;
  bool check_eq;

  MPI_Init( & argc, & argv );
  MPI_Comm_size( MPI_COMM_WORLD, &npes );
  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  ///////////////////////////
  // initial input processing

  // check on input parameters
  dimension = atoi(argv[1]);
  iterations = atoi(argv[2]);
  row_peek = atoi(argv[3]);
  col_peek = atoi(argv[4]);
  freq = atoi(argv[5]);

  if (iterations < 0){
    check_eq = true;
    iterations = abs(iterations);
  } else check_eq = false;

  n_loc = dimension / npes;
  rest = dimension % npes;
  if( me < rest ) n_loc += 1;
  else offset = rest;

  byte_dimension = sizeof(double) * ( n_loc + 2 ) * ( dimension + 2 );
  matrix = ( double* )malloc( byte_dimension );
  matrix_new = ( double* )malloc( byte_dimension );
  tmp_matrix = ( double* )malloc( byte_dimension );


  memset( matrix, 0, byte_dimension );
  memset( matrix_new, 0, byte_dimension );

  int row_start = 0, row_end = n_loc + 2;

  //fill initial values
  if (!me) row_start = 1;
  if (me==npes -1) row_end = n_loc + 1;

  for( i = row_start; i < row_end; i++ ){
    for( j = 1; j < dimension + 1; j++ ){
      matrix[ ( i * ( dimension + 2 ) ) + j ] = 0.5;
      matrix_new[ ( i * ( dimension + 2 ) ) + j ] = 0.5;
    }
  }

  // set up borders
  increment = 100.0 / ( dimension + 1.0 );
  //bottom row
  if( me == npes - 1 ){
    for( j = 1; j < dimension + 1; j++){
      matrix[ (n_loc + 1) * ( dimension + 2 ) + j] = 100 - j * increment;
      matrix_new[ (n_loc + 1) * ( dimension + 2 ) + j] = 100 - j * increment;
      }
  }

  //left column
  for(i = row_start; i < n_loc + 2; i++){
      i_g = global_i(i, me, n_loc, offset);
      matrix[ i*( dimension + 2 ) ] = i_g*increment;
      matrix_new[ i*( dimension + 2 ) ] = i_g*increment;
    }


  //visualising the initialisation

    if (dimension <= 32){

    if( !me ){
    memcpy(tmp_matrix, matrix, byte_dimension);
    printf("\n");
    print_matrix( tmp_matrix, n_loc, dimension );
    printf("\n");
    for( count = 1; count < npes; count += 1 ){
      MPI_Recv(tmp_matrix, (dimension + 2) * (n_loc + 2), MPI_DOUBLE, count,count, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      
      if( count == rest ){
         n_loc2 = n_loc - 1;
      } else n_loc2 = n_loc;
      print_matrix( tmp_matrix, n_loc2, dimension );
      printf("\n");
    }
  printf("\n\n");
  } else MPI_Send(matrix, (dimension + 2) * (n_loc + 2), MPI_DOUBLE, 0, me, MPI_COMM_WORLD);

  }

  // start algorithm
  if (!me)  t_start = seconds();

  for( it = 0; it < iterations; it++ ){

    evolve( matrix, matrix_new, n_loc, dimension );

    // swap the pointers
    tmpry_matrix = matrix;
    matrix = matrix_new;
    matrix_new = tmpry_matrix;

    //saving iteration output
    if(it%freq == 0 ){

      if( !me ){

      const double h = 1;
      FILE *file;
      // double *tmpr_matrix;

      file = fopen( "solution.dat", "a" );

        memcpy(tmp_matrix, matrix, byte_dimension);
        for( i = 0; i < n_loc + 1; i++ ){
          i_g = global_i(i, me, n_loc, offset);
          for( j = 0; j < dimension + 2; j++ ){
            fprintf(file, "%d\t%f\t%f\t%f\n", it, h * j, -h * i_g, tmp_matrix[ ( i * ( dimension + 2 ) ) + j ] );
          }
        }
        n_loc2 = n_loc;
        for( count = 1; count < npes; count++ ){
          if( count == rest ){
            n_loc2 = n_loc2 - 1;
          }
          MPI_Recv(tmp_matrix, (dimension + 2) * (n_loc + 2), MPI_DOUBLE, count,count, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
          //printf("%d/n", n_loc2);
          if (count == npes-1){
             row_end = n_loc2 + 2 ;
          } else row_end = n_loc2 + 1;
          for( i = 1; i < row_end; i++ ){
            i_g = global_i(i, count, n_loc2, offset);
            for( j = 0; j < dimension + 2; j++ ){
              fprintf(file, "%d\t%f\t%f\t%f\n", it, h * j, -h * i_g, tmp_matrix[ ( i * ( dimension + 2 ) ) + j ] );
            }
          }
        }
        fclose( file );
       } else MPI_Send(matrix, (dimension + 2) * (n_loc + 2), MPI_DOUBLE, 0, me, MPI_COMM_WORLD);
  }

  // check equilibrium

  if(it%freq == 0 ){

    if(check_eq){
    int equal = matricesEqual(matrix, matrix_new, n_loc, dimension);

    MPI_Allreduce(MPI_IN_PLACE, &equal, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  
    if (equal) {
            if (!me) {
                printf("Reached equilibrium after %d iterations.\n", it+1);
            }
            break;
        }
    }
  }
  //sharing top row and rcv btm row
  if (me != 0){
  MPI_Sendrecv( &matrix[(dimension + 2) + 1 ], dimension , MPI_DOUBLE, me - 1, 100*me, &matrix[1], dimension, MPI_DOUBLE, me -1, me - 1 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //sharing bottom row and rcv top row
  if(me != npes -1){
  MPI_Sendrecv( &matrix[ (n_loc) * (dimension + 2) + 1], dimension , MPI_DOUBLE, me + 1, me, &matrix[(n_loc + 1)*(dimension + 2) + 1 ], dimension, MPI_DOUBLE, me + 1, 100*(me+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

}



  if(!me){
  t_end = seconds();

  printf( "%d %f\n", npes, t_end - t_start );
  // printf( "\nelapsed time = %f seconds\n", t_end - t_start );
  // printf( "\nmatrix[%d,%d] = %f\n", row_peek, col_peek, matrix[ ( row_peek + 1 ) * ( dimension + 2 ) + ( col_peek + 1 ) ] );
  
  // save_gnuplot( matrix, dimension );
  }

  //Visualising result
    if (dimension <= 32){

    if( !me ){
    memcpy(tmp_matrix, matrix, byte_dimension);
    printf("\nIteration : %d, Me : %d \n\n", it+1 , me);
    printf("\n");
    print_matrix( tmp_matrix, n_loc, dimension );
    printf("\n");
    for( count = 1; count < npes; count += 1 ){
      MPI_Recv(tmp_matrix, (dimension + 2) * (n_loc + 2), MPI_DOUBLE, count,count, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      if( count == rest ){
         n_loc2 = n_loc - 1;
      } else n_loc2 = n_loc;
      print_matrix( tmp_matrix, n_loc2, dimension );
      printf("\n");

    }
  printf("\n\n");
  } else MPI_Send(matrix, (dimension + 2) * (n_loc + 2), MPI_DOUBLE, 0, me, MPI_COMM_WORLD);

  }

  free( matrix );
  free( matrix_new );
//   free( tmp_matrix );
//   free( tmpry_matrix );


  MPI_Finalize();

  return 0;
}

void evolve( double * matrix, double *matrix_new, int n_loc, int dimension){

  int i , j;

  //This will be a row dominant program.
  for( i = 1 ; i <= n_loc; i++ )
    for( j = 1; j <= dimension; j++ )
      matrix_new[ ( i * ( dimension + 2 ) ) + j ] = ( 0.25 ) *
            ( matrix[ ( ( i - 1 ) * ( dimension + 2 ) ) + j ] +
              matrix[ ( i * ( dimension + 2 ) ) + ( j + 1 ) ] +

              matrix[ ( ( i + 1 ) * ( dimension + 2 ) ) + j ] +
              matrix[ ( i * ( dimension + 2 ) ) + ( j - 1 ) ] );
}

void save_gnuplot( double *M, int dimension ){

  int i , j;
  const double h = 0.1;
  FILE *file;

  file = fopen( "solution.dat", "w" );

  for( i = 0; i < dimension + 2; ++i )
    for( j = 0; j < dimension + 2; ++j )
      fprintf(file, "%f\t%f\t%f\n", h * j, -h * i, M[ ( i * ( dimension + 2 ) ) + j ] );

  fclose( file );

}

// A Simple timer for measuring the walltime
double seconds(){

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}