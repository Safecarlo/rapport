#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

double reduc_somme(double in)
{
  int rang,taille,root;
  double sum,*all_in;

  MPI_Comm_rank(MPI_COMM_WORLD, &rang);
  MPI_Comm_size(MPI_COMM_WORLD, &taille);

  root=0;
  if(rang == root)
  {
    all_in= (double*)malloc(taille*sizeof(double));
  }else
  {
    all_sum=NULL;
  }

  sum=0.0;

  for (int i = 0; i < taille; i++)
  {
    sum+=all_in[i];
  }
  free(all_in);
}

int main(int argc, char const *argv[])
{

  return 0;
}
