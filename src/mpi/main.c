#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <mpi.h>

/*lecture du fichier de donnee*/
void read_file(const char *file, double *vec, int n)
{
  FILE *fd;
  double val;

  fd=fopen(file, "r");
  vec=malloc(n*sizeof(double));
  if(fd==NULL)
  {
    puts("Erreur d'ouverture du fichier");
    exit(0);
  }else
  {
    /*recuperation des valeur du tableau*/
    for (int i = 0; i < n; i++)
    {
      fscanf(fd, "%lf", &val);
      vec[i]=val;
    }
  }

  for (int i = 0; i < n; i++)
  {
    printf("%lf\n",vec[i] );
  }

  free(vec);
  fclose(fd);
}
/*Creation de la fonction de reduction*/
double reduc_somme(double *restrict all_in, int n)
{
  double sum;

    all_in = NULL;
    all_in = (double*)malloc(n*sizeof(double));

      if(all_in==NULL)
      {
        puts("Allocation impossible");
        exit(0);
      }

  sum=0;

  for (int i = 0; i < n; i++)
  {
    sum+=all_in[i];
  }

  free(all_in);
  return sum;
}

int main(int argc, char **argv)
{
  int nelmt,rang;
  char *name;
  double *x,sum;

  name = argv[1];
  nelmt = atoi(argv[2]);

  read_file(name, x, nelmt);

  sum=reduc_somme(x,nelmt);
  printf("Somme = %lf\n",sum);

  return 0;
}
