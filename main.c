#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

double dot_prod(double *restrict a, double *restrict b, int count)
{
  double sum=0.0;

  for (int i = 0; i < count; i++)
  {
    sum+=a[i] * b[i];
  }
  return sum;
}
int main(int argc, char const *argv[])
{
  int rank, size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  
  MPI_Finalize();
  return 0;
}
