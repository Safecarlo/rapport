//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a parrallel mpi C version of the NPB CG code. This C        //
//  version is developed by the Center for Manycore Programming at Seoul   //
//  National University and derived from the serial Fortran versions in    //
//  "NPB3.3-SER" developed by NAS.                                         //
//                                                                         //
//  Permission to use, copy, distribute and modify this software for any   //
//  purpose with or without fee is hereby granted. This software is        //
//  provided "as is" without express or implied warranty.                  //
//                                                                         //
//  Information on NPB 3.3, including the technical report, the original   //
//  specifications, source code, results and information on how to submit  //
//  new results, is available at:                                          //
//                                                                         //
//           http://www.nas.nasa.gov/Software/NPB/                         //
//                                                                         //
//  Send comments or suggestions for this C version to cmp@aces.snu.ac.kr  //
//                                                                         //
//          Center for Manycore Programming                                //
//          School of Computer Science and Engineering                     //
//          Seoul National University                                      //
//          Seoul 151-744, Korea                                           //
//                                                                         //
//          E-mail:  cmp@aces.snu.ac.kr                                    //
//                                                                         //
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//                                              //
//-------------------------------------------------------------------------//

//---------------------------------------------------------------------
// NPB CG serial version
//---------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "globals.h"
#include "randdp.h"
#include "timers.h"
#include "print_results.h"

//---------------------------------------------------------------------
/* common / main_int_mem / */
static int colidx[NZ];
static int rowstr[NA+1];
static int iv[NA];
static int arow[NA];
static int acol[NAZ];

/* common / main_flt_mem / */
static double v[NA+1];
static double aelt[NAZ];
static double a[NZ];
static double x[NA+2];
static double z[NA+2];
static double p[NA+2];
static double q[NA+2];
static double r[NA+2];
static double w[NA+2];
/* common / partit_size / */
static int naa;
static int nzz;
static int firstrow;
static int lastrow;
static int firstcol;
static int lastcol;
static int npcols;
static int nprows;
static int exch_proc;
static int exch_recv_legth;
static int send_start;
static int send_len;

/* common /urando/ */
static double amult;
static double tran;

/* common /timers/ */
static logical timeron;

/*work arrays for reduction*/
int l2npcols;

/*process info*/
int num_procs;
int num_proc_rows;
int num_proc_col;
int rank;
int nprocs;
int proc_row;
int proc_col;
int *reduce_exch_proc;
int *reduce_recv_starts;
int *reduce_send_starts;
int *reduce_recv_lengths;
int *reduce_send_legths;
//---------------------------------------------------------------------


//---------------------------------------------------------------------
static void setup_proc_info();

static void setup_submatrix_info();

static void conj_grad(double *rnorm);

static void makea(int n,
                  int nz,
                  double a[],
                  int colidx[],
                  int rowstr[],
                  int nonzer,
                  int firstrow,
                  int lastrow,
                  int firstcol,
                  int lastcol,
                  double rcond,
                  int arow[],
                  int acol[NONZER+1],
                  double aelt[NONZER+1],
                  double v[],
                  int iv[],double shift);
static void sparse(double a[],
                   int colidx[],
                   int rowstr[],
                   int n,
                   int arow[],
                   int acol[NONZER+1],
                   double aelt[NONZER+1],
                   int firstrow,
                   int lastrow,
                   double x[],
                   double mark[],
                   int nzloc[],
                   int nnza);
static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[],int mark[]);
static int icnvrt(double x, int ipwr2);
static void vecset(int n, double v[], int iv[], int *nzv, int i, double val);
//---------------------------------------------------------------------


int main(int argc, char **argv)
{
  int i, j, k, it, root;

  double zeta;
  double rnorm;
  double norm_temp1[2], norm_temp2[2];

  double t, mflops, tmax;
  char Class;
  logical verified;
  double zeta_verify_value, epsilon, err;

  char *t_names[T_last];

  MPI_Init(&argc, &argv);
  MPI_Request request;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  root = 0;

  /*!---------------------------------------------------------------------
!  Set up processor info, such as whether sq num of procs, etc
!---------------------------------------------------------------------*/
  setup_proc_info();

  FILE *fp;
  if ((fp = fopen("timer.flag", "r")) != NULL)
  {
    timeron = true;
    t_names[T_init] = "init";
    t_names[T_bench] = "benchmk";
    t_names[T_conj_grad] = "conjgd";
    fclose(fp);
  } else {
    timeron = false;
  }
   timer_start(T_init);

  if (NA == 1400 && NONZER == 7 && NITER == 15 && SHIFT == 10) {
    Class = 'S';
    zeta_verify_value = 8.5971775078648;
  } else if (NA == 7000 && NONZER == 8 && NITER == 15 && SHIFT == 12) {
    Class = 'W';
    zeta_verify_value = 10.362595087124;
  } else if (NA == 14000 && NONZER == 11 && NITER == 15 && SHIFT == 20) {
    Class = 'A';
    zeta_verify_value = 17.130235054029;
  } else if (NA == 75000 && NONZER == 13 && NITER == 75 && SHIFT == 60) {
    Class = 'B';
    zeta_verify_value = 22.712745482631;
  } else if (NA == 150000 && NONZER == 15 && NITER == 75 && SHIFT == 110) {
    Class = 'C';
    zeta_verify_value = 28.973605592845;
  } else if (NA == 1500000 && NONZER == 21 && NITER == 100 && SHIFT == 500) {
    Class = 'D';
    zeta_verify_value = 52.514532105794;
  } else if (NA == 9000000 && NONZER == 26 && NITER == 100 && SHIFT == 1500) {
    Class = 'E';
    zeta_verify_value = 77.522164599383;
  } else {
    Class = 'U';
  }
  if (rank == root)
  {
    printf("\n\n NAS Parallel Benchmarks (NPB3.3-PRALLEL-MPI-C) - CG Benchmark\n\n");
    printf(" Size: %11d , CLASS :'%c'\n", NA, Class);
    printf(" Iterations: %5d\n", NITER);
    printf("\n");
  }
  /*!---------------------------------------------------------------------
!  Set up partition's submatrix info: firstcol, lastcol, firstrow, lastrow
!---------------------------------------------------------------------*/
  setup_submatrix_info();

  for (i = 0; i < T_last; i++)
  {
    timer_clear(i);
  }

  naa = NA;
  nzz = NZ;
  //---------------------------------------------------------------------
  // Inialize random number generator
  //---------------------------------------------------------------------
  tran    = 314159265.0;
  amult   = 1220703125.0;
  zeta    = randlc(&tran, amult);

  //---------------------------------------------------------------------
  // Set up partition's sparse random matrix for given class size
  //---------------------------------------------------------------------
  makea(naa, nzz, a, colidx, rowstr,NONZER, firstrow, lastrow, firstcol, lastcol, RCOND, arow, acol, aelt, v,iv, SHIFT);

  //---------------------------------------------------------------------
  // Note: as a result of the above call to makea:
  //      values of j used in indexing rowstr go from 0 --> lastrow-firstrow
  //      values of colidx which are col indexes go from firstcol --> lastcol
  //      So:
  //      Shift the col index vals from actual (firstcol --> lastcol )
  //      to local, i.e., (0 --> lastcol-firstcol)
  //---------------------------------------------------------------------
  for (j = 0; j < lastrow - firstrow + 1; j++)
  {
    for (k = rowstr[j]; k < rowstr[j+1]; k++)
    {
      colidx[k] = colidx[k] - firstcol;
    }
  }

  //---------------------------------------------------------------------
  // set starting vector to (1, 1, .... 1)
  //---------------------------------------------------------------------
  for (i = 0; i < NA+1; i++)
  {
    x[i] = 1.0;
  }
  zeta = 0.0;

  //---------------------------------------------------------------------
  //---->
  // Do one iteration untimed to init all code and data page tables
  //---->                    (then reinit, start timing, to niter its)
  //---------------------------------------------------------------------
  for (it = 1; it <= 1; it++)
  {
    //---------------------------------------------------------------------
    // The call to the conjugate gradient routine:
    //---------------------------------------------------------------------
    conj_grad( &rnorm);

    //---------------------------------------------------------------------
    // zeta = shift + 1/(x.z)
    // So, first: (x.z)
    // Also, find norm of z
    // So, first: (z.z)
    //---------------------------------------------------------------------
    norm_temp1[1] = 0.0;
    norm_temp1[2]= 0.0;
    for (j = 0; j < lastcol - firstcol + 1; j++)
    {
      norm_temp1[1] = norm_temp1[1] + x[j] * z[j];
      norm_temp1[2] = norm_temp1[2] + z[j] * z[j];
    }

    for ( i = 0; i < l2npcols; i++)
    {
      MPI_Irecv(norm_temp2, 2, MPI_DOUBLE, reduce_exch_proc[i], i,
                  MPI_COMM_WORLD, &request);
      MPI_Send(norm_temp1, 2, MPI_DOUBLE, reduce_exch_proc[i], i,
                MPI_COMM_WORLD);
      MPI_Wait(&request, MPI_STATUS_IGNORE);

      norm_temp1[1] =norm_temp1[1] + norm_temp2[2];
      norm_temp1[1] =norm_temp1[1] + norm_temp2[2];

    }

    norm_temp1[1] = 1.0 / sqrt(norm_temp1[2]);

    if (rank == root)
    {
      zeta = SHIFT + 1.0/sqrt(norm_temp1[1]);
    }

    //---------------------------------------------------------------------
    // Normalize z to obtain x
    //---------------------------------------------------------------------
    for (j = 0; j < lastcol - firstcol + 1; j++)
    {
      x[j] = norm_temp1[2] * z[j];
    }
  } // end of do one iteration untimed


  //---------------------------------------------------------------------
  // set starting vector to (1, 1, .... 1)
  //---------------------------------------------------------------------
  for (i = 0; i < NA+1; i++) {
    x[i] = 1.0;
  }

  zeta = 0.0;

  timer_stop(T_init);

  printf(" Initialization time = %15.3f seconds\n", timer_read(T_init));

  timer_start(T_bench);

  //---------------------------------------------------------------------
  //---->
  // Main Iteration for inverse power method
  //---->
  //---------------------------------------------------------------------
  for (it = 1; it <= NITER; it++) {
    //---------------------------------------------------------------------
    // The call to the conjugate gradient routine:
    //---------------------------------------------------------------------
    if (timeron) timer_start(T_conj_grad);
    conj_grad( &rnorm);
    if (timeron) timer_stop(T_conj_grad);

    //---------------------------------------------------------------------
    // zeta = shift + 1/(x.z)
    // So, first: (x.z)
    // Also, find norm of z
    // So, first: (z.z)
    //---------------------------------------------------------------------
    norm_temp1[1] = 0.0;
    norm_temp1[2] = 0.0;
    for (j = 0; j < lastcol - firstcol + 1; j++)
    {
      norm_temp1[1] = norm_temp1[1] + x[j]*z[j];
      norm_temp1[2] = norm_temp1[2] + z[j]*z[j];
    }

    for ( i = 0; i < l2npcols; i++)
    {
      MPI_Irecv(norm_temp2, 2, MPI_DOUBLE, reduce_exch_proc[i], i,
                  MPI_COMM_WORLD, &request);
      MPI_Send(norm_temp1, 2, MPI_DOUBLE, reduce_exch_proc[i], i,
                MPI_COMM_WORLD);
      MPI_Wait(&request, MPI_STATUS_IGNORE);

      norm_temp1[1] =norm_temp1[1] + norm_temp2[2];
      norm_temp1[1] =norm_temp1[1] + norm_temp2[2];

    }

    norm_temp1[1] = 1.0 / sqrt(norm_temp1[2]);

    if (rank == root)
    {
      zeta = SHIFT + 1.0 / norm_temp1[1];
      if (it == 1)
      {
        printf("\n   iteration           ||r||                 zeta\n");
        printf("    %5d       %20.14E%20.13f\n", it, rnorm, zeta);
      }
    }


    //---------------------------------------------------------------------
    // Normalize z to obtain x
    //---------------------------------------------------------------------
    for (j = 0; j < lastcol - firstcol + 1; j++)
    {
      x[j] = norm_temp2[2] * z[j];
    }
  } // end of main iter inv pow meth

  timer_stop(T_bench);

  //---------------------------------------------------------------------
  // End of timed section
  //---------------------------------------------------------------------

  t = timer_read(T_bench);
if (rank == root)
{
  printf(" Benchmark completed\n");

  epsilon = 1.0e-10;
  if (Class != 'U') {
    err = fabs(zeta - zeta_verify_value) / zeta_verify_value;
    if (err <= epsilon) {
      verified = true;
      printf(" VERIFICATION SUCCESSFUL\n");
      printf(" Zeta is    %20.13E\n", zeta);
      printf(" Error is   %20.13E\n", err);
    } else {
      verified = false;
      printf(" VERIFICATION FAILED\n");
      printf(" Zeta                %20.13E\n", zeta);
      printf(" The correct zeta is %20.13E\n", zeta_verify_value);
    }
  } else {
    verified = false;
    printf(" Problem size unknown\n");
    printf(" NO VERIFICATION PERFORMED\n");
  }
}


  if (t != 0.0) {
    mflops = (double)(2*NITER*NA)
                   * (3.0+(double)(NONZER*(NONZER+1))
                     + 25.0*(5.0+(double)(NONZER*(NONZER+1)))
                     + 3.0) / t / 1000000.0;
  } else {
    mflops = 0.0;
  }

  print_results("CG", Class, NA, 0, 0,
                NITER, t,
                mflops, "          floating point",
                verified, NPBVERSION, COMPILETIME,
                CS1, CS2, CS3, CS4, CS5, CS6, CS7);

  //---------------------------------------------------------------------
  // More timers
  //---------------------------------------------------------------------
  if (timeron) {
    tmax = timer_read(T_bench);
    if (tmax == 0.0) tmax = 1.0;
    printf("  SECTION   Time (secs)\n");
    for (i = 0; i < T_last; i++) {
      t = timer_read(i);
      if (i == T_init) {
        printf("  %8s:%9.3f\n", t_names[i], t);
      } else {
        printf("  %8s:%9.3f  (%6.2f%%)\n", t_names[i], t, t*100.0/tmax);
        if (i == T_conj_grad) {
          t = tmax - t;
          printf("    --> %8s:%9.3f  (%6.2f%%)\n", "rest", t, t*100.0/tmax);
        }
      }
    }
  }
  MPI_Finalize();
  return 0;
}

/*Setup proc info*/
static void setup_proc_info()
{
  int root;
/*  !---------------------------------------------------------------------
!  set up dimension parameters after partition
!  num_proc_rows & num_proc_cols are set by get_active_nprocs
!---------------------------------------------------------------------*/

/*!---------------------------------------------------------------------
!  num_procs must be a power of 2, and num_procs=num_proc_cols*num_proc_rows.
!  num_proc_cols and num_proc_cols are to be found in npbparams.h.
!  When num_procs is not square, then num_proc_cols must be = 2*num_proc_rows.
!---------------------------------------------------------------------*/
      num_procs = num_proc_cols * num_proc_rows;

/*!---------------------------------------------------------------------
!  num_procs must be a power of 2, and num_procs=num_proc_cols*num_proc_rows
!  When num_procs is not square, then num_proc_cols = 2*num_proc_rows
!---------------------------------------------------------------------
!  First, number of procs must be power of two.
!--------------------------------------------------------------------->*/
  root=0;
  if (nprocs != num_procs)
  {
    if (rank == root)
    {
      printf("ERROR:Number of process\n");
    }
  }

  npcols = num_proc_cols;
  nprows = num_proc_rows;
}

/*submatrix*/
static void setup_submatrix_info()
{
  int col_size, row_size;
  int i, j;
  int div_factor;

  proc_row = rank/npcols;
  proc_col = rank - proc_row*npcols;

  /*!---------------------------------------------------------------------
!  If na evenly divisible by npcols, then it is evenly divisible
!  by nprows
!---------------------------------------------------------------------
*/
  if (NA/npcols*npcols == NA)
  {
    col_size = NA /npcols;
    firstcol = proc_col*col_size + 1;
    lastcol = firstcol -1 +col_size;

    row_size = NA/nprows;
    firstrow = proc_row * row_size +1;
    lastrow = firstrow -1 + row_size;
  }/*!---------------------------------------------------------------------
!  If na not evenly divisible by npcols, then first subdivide for nprows
!  and then, if npcols not equal to nprows (i.e., not a sq number of procs),
!  get col subdivisions by dividing by 2 each row subdivision.
!---------------------------------------------------------------------0*/
  else
  {
    if (proc_row < NA-NA/nprows*nprows)
    {
      row_size = NA/nprows+ 1;
      firstrow = proc_row*row_size + 1;
      lastrow  = firstrow - 1 + row_size;
    }
    else
    {
      row_size = NA/nprows;
      firstrow = (NA - NA/nprows*nprows)*(row_size+1)+
                  (proc_row - (NA/NA/nprows*nprows))*
                  row_size +1;
      lastrow =firstrow -1 + row_size;
    }
    if(npcols == nprows)
    {
      if(proc_col < NA -NA/npcols*npcols)
      {
        col_size = NA/npcols+1;
        firstcol = proc_col * col_size +1;
        lastcol = firstcol - 1 +col_size;
      }
      else
      {
        col_size = NA/npcols;
        firstcol = (NA-NA/npcols*npcols)*(col_size+1)+
                    (proc_col -(NA-NA/npcols*npcols))*col_size+1;
        lastcol = firstcol -1 +col_size;
      }
    }else
    {
      if ((proc_col/2) < NA-NA/(npcols/2)*(npcols/2))
      {
        col_size = NA/(npcols/2) + 1;
        firstcol = (npcols/2)*col_size+1;
        lastcol =firstcol -1 +col_size;
      }
      else
      {
        col_size = NA/(npcols/2) + 1;
        firstcol = (NA -NA/(npcols/2)*(npcols/2))*(col_size + 1)
                    + ((proc_col/2)-(NA-NA/(npcols/2)*(npcols/2)))*col_size+1;
        lastcol = firstcol - 1 +col_size;
      }
      if(rank%2 == 0)
      {
        lastcol = firstcol + (col_size-1)/2+1;

      }
      else
      {
        firstcol =firstcol +(col_size-1)/2 +1;
        lastcol =firstcol -1 + col_size/2;
      }
    }
  }

  if (npcols == nprows)
  {
    send_start = 1;
    send_len = lastrow-firstrow+1;
  }
  else
  {
    if(rank %2 ==0)
    {
      send_start = 1;
      send_len = (1 + lastrow-firstrow+1)/2;
    }
    else
    {
      send_start = 1;
      send_len = (lastrow-firstrow+1)/2;
    }
  }

  /*!---------------------------------------------------------------------
!  Transpose exchange processor
!---------------------------------------------------------------------*/
  if (npcols == nprows)
  {
    exch_proc = (rank%nprows)*nprows + rank/nprows;
  }
  else
  {
    exch_proc = 2*(((rank/2)%nprows)*nprows + rank/nprows)+(rank%2);
  }

  i = npcols/2;
  l2npcols =0;
  do{
    l2npcols = l2npcols +1;
    i =i/2;
  } while(i>0);

/*  !---------------------------------------------------------------------
!  Set up the reduce phase schedules...
!---------------------------------------------------------------------*/
  div_factor = npcols;

  for ( i = 0; i < l2npcols; i++)
  {
    j = ((proc_col+div_factor/2)%div_factor)+ proc_col/div_factor*div_factor;

    reduce_exch_proc[i] = proc_row*npcols +j;

    div_factor = div_factor/2;
  }

  for ( i = l2npcols; i > 0; i--)
  {
    if (nprows == npcols)
    {
      reduce_send_starts[i] = send_start;
      reduce_send_legths[i] = send_len;
      reduce_recv_lengths[i] = lastrow - firstrow +1;
    }
    else
    {
      reduce_recv_lengths[i] = send_len;
      if (i == l2npcols)
      {
        reduce_send_legths[i] = lastrow - firstrow+1 -send_len;
        if (rank/2*2 == rank)
        {
          reduce_send_starts[i] = send_start + send_len;
        }
        else
        {
          reduce_send_starts[i] =1;
        }
      }
      else
      {
        reduce_send_legths[i] = send_len;
        reduce_send_starts[i] = send_start;
      }
    }
  }
  exch_recv_legth = lastcol -firstcol+1;
}


//---------------------------------------------------------------------
// Floaging point arrays here are named as in NPB1 spec discussion of
// CG algorithm
//---------------------------------------------------------------------
static void conj_grad(double *rnorm)
{
  int j, k, root;
  int cgit, cgitmax = 25;
  double d, sum, rho, rho0, alpha, beta;

  /*Variable for MPI*/
  MPI_Request request;

  rho = 0.0;
  root=0;

  /*allocation for rnorm*/
  if (rank==root)
  {
    rnorm = malloc(sizeof(double));
  }


  //---------------------------------------------------------------------
  // Initialize the CG algorithm:
  //---------------------------------------------------------------------
  for (j = 0; j < naa+1; j++) {
    q[j] = 0.0;
    z[j] = 0.0;
    r[j] = x[j];
    p[j] = r[j];
    w[j] =0.0;
  }

  //---------------------------------------------------------------------
  // rho = r.r
  // Now, obtain the norm of r: First, sum squares of r elements locally...
  //---------------------------------------------------------------------
  sum = 0.0;
  for (j = 0; j < lastcol - firstcol + 1; j++)
  {
    sum+= r[j]*r[j];
  }

  /*---------------------------------------------------------------------
!  Exchange and sum with procs identified in reduce_exch_proc
!  (This is equivalent to mpi_allreduce.)
!  Sum the partial sums of rho, leaving rho on all processors
!---------------------------------------------------------------------*/
  for (int i = 0; i < l2npcols; i++)
  {
    MPI_Irecv(&rho, 1, MPI_DOUBLE, reduce_exch_proc[i], i, MPI_COMM_WORLD, &request);
    MPI_Send(&sum, 1, MPI_DOUBLE, reduce_exch_proc[i], i, MPI_COMM_WORLD);
    MPI_Wait(&request, MPI_STATUS_IGNORE);

    sum= sum+ rho;
  }

  rho =sum;
  //---------------------------------------------------------------------
  //---->
  // The conj grad iteration loop
  //---->
  //---------------------------------------------------------------------
  for (cgit = 1; cgit <= cgitmax; cgit++)
  {
    //---------------------------------------------------------------------
    // q = A.p
    // The partition submatrix-vector multiply: use workspace w
    //---------------------------------------------------------------------

    for (j = 0; j < lastrow - firstrow + 1; j++)
    {
      sum = 0.0;
      for (k = rowstr[j]; k < rowstr[j+1]; k++)
      {
        sum = sum + a[k]*p[colidx[k]];
      }
      w[j] = sum;
    }
    /*---------------------------------------------------------------------
!  Sum the partition submatrix-vec A.p's across rows
!  Exchange and sum piece of w with procs identified in reduce_exch_proc
!---------------------------------------------------------------------*/
    for (int i = l2npcols; i > 0; i--)
    {
      MPI_Irecv(&(q[reduce_recv_starts[i]]), reduce_recv_lengths[i],
                MPI_DOUBLE, reduce_exch_proc[i], i, MPI_COMM_WORLD,&request);
      MPI_Send(&(w[reduce_recv_starts[i]]),reduce_send_legths[i],
                MPI_DOUBLE, reduce_exch_proc[i], i, MPI_COMM_WORLD);
      MPI_Wait(&request, MPI_STATUS_IGNORE);

      for ( j = send_start; j < send_start + reduce_send_legths[i]; j++)
      {
        w[j] = w[j] + q[j];
      }
    }
    /*!---------------------------------------------------------------------
!  Exchange piece of q with transpose processor:
!---------------------------------------------------------------------*/
    if (l2npcols!=0)
    {
      MPI_Irecv(q, exch_recv_legth, MPI_DOUBLE, exch_proc, 1,
                MPI_COMM_WORLD, &request);
      MPI_Send(&(w[send_start]), send_len, MPI_DOUBLE, exch_proc
                ,1 ,MPI_COMM_WORLD);
      MPI_Wait(&request, MPI_STATUS_IGNORE);
    }
    else
    {
      for (j = 0; j < exch_recv_legth; j++)
      {
        w[j] = q[j];
      }
    }

    /*!---------------------------------------------------------------------
!  Clear w for reuse...
!---------------------------------------------------------------------*/
    for ( j = 0; j < max(lastrow - firstrow, lastcol - firstcol); j++)
    {
      w[j] = 0.0;
    }
    //---------------------------------------------------------------------
    // Obtain p.q
    //---------------------------------------------------------------------
    sum = 0.0;
    for (j = 0; j < lastcol - firstcol + 1; j++)
    {
      sum = sum + p[j]*q[j];
    }
    /*!---------------------------------------------------------------------
    !  Obtain d with a sum-reduce
    !---------------------------------------------------------------------*/
    for (int i = 0; i < l2npcols; i++)
    {
      MPI_Irecv(&d, 1, MPI_DOUBLE, reduce_exch_proc[i], i, MPI_COMM_WORLD,&request);
      MPI_Send(&sum, 1, MPI_DOUBLE, reduce_exch_proc[i], i, MPI_COMM_WORLD);
      MPI_Wait(&request, MPI_STATUS_IGNORE);

      sum=sum+d;
    }

    d = sum;
    //---------------------------------------------------------------------
    // Obtain alpha = rho / (p.q)
    //---------------------------------------------------------------------
    alpha = rho / d;

    //---------------------------------------------------------------------
    // Save a temporary of rho
    //---------------------------------------------------------------------
    rho0 = rho;

    //---------------------------------------------------------------------
    // Obtain z = z + alpha*p
    // and    r = r - alpha*q
    //---------------------------------------------------------------------
    rho = 0.0;
    for (j = 0; j < lastcol - firstcol + 1; j++)
    {
      z[j] = z[j] + alpha*p[j];
      r[j] = r[j] - alpha*q[j];
    }

    //---------------------------------------------------------------------
    // rho = r.r
    // Now, obtain the norm of r: First, sum squares of r elements locally...
    //---------------------------------------------------------------------
    sum=0.0;
    for (j = 0; j < lastcol - firstcol + 1; j++)
    {
      sum = sum + r[j]*r[j];
    }
    /*!---------------------------------------------------------------------
!  rho = r.r
!  Now, obtain the norm of r: First, sum squares of r elements locally...
!---------------------------------------------------------------------*/
    for (int i = 0; i < l2npcols; i++)
    {
      MPI_Irecv(&rho, 1, MPI_DOUBLE, reduce_exch_proc[i], i, MPI_COMM_WORLD,&request);
      MPI_Send(&sum, 1, MPI_DOUBLE, reduce_exch_proc[i], i,MPI_COMM_WORLD);
      MPI_Wait(&request, MPI_STATUS_IGNORE);
      sum = sum + rho;
    }
    rho = sum;
    //---------------------------------------------------------------------
    // Obtain beta:
    //---------------------------------------------------------------------
    beta = rho / rho0;

    //---------------------------------------------------------------------
    // p = r + beta*p
    //---------------------------------------------------------------------
    for (j = 0; j < lastcol - firstcol + 1; j++)
    {
      p[j] = r[j] + beta*p[j];
    }
  } // end of do cgit=1,cgitmax

  //---------------------------------------------------------------------
  // Compute residual norm explicitly:  ||r|| = ||x - A.z||
  // First, form A.z
  // The partition submatrix-vector multiply
  //---------------------------------------------------------------------
  for (j = 0; j < lastrow - firstrow + 1; j++)
  {
    sum = 0.0;
    for (k = rowstr[j]; k < rowstr[j+1]; k++) {
      sum = sum + a[k]*z[colidx[k]];
    }
    w[j] = sum;
  }

  /*!---------------------------------------------------------------------
!  Sum the partition submatrix-vec A.z's across rows
!---------------------------------------------------------------------*/
  for (int i = l2npcols; i > 0; i--)
  {
    MPI_Irecv(&(r[reduce_recv_starts[i]]), reduce_recv_lengths[i],
              MPI_DOUBLE, reduce_exch_proc[i], i, MPI_COMM_WORLD, &request);
    MPI_Send(&(w[reduce_send_starts[i]]), reduce_send_legths[i],
              MPI_DOUBLE, reduce_exch_proc[i], i, MPI_COMM_WORLD);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    for (j = send_start; j < send_start - reduce_recv_starts[i]; j++)
    {
      w[j] = w[j] + r[j];
    }
  }

  /*!---------------------------------------------------------------------
!  Exchange piece of q with transpose processor:
!---------------------------------------------------------------------*/
  if (l2npcols!=0)
  {
    MPI_Irecv(r, exch_recv_legth, MPI_DOUBLE, exch_proc, 1,
              MPI_COMM_WORLD, &request);
    MPI_Send(&(w[send_start]), send_len, MPI_DOUBLE, exch_proc,
              1, MPI_COMM_WORLD);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
  }
  else
  {
    for ( j = 0; j < exch_recv_legth; j++)
    {
      r[j] = w[j];
    }
  }
  //---------------------------------------------------------------------
  // At this point, r contains A.z
  //---------------------------------------------------------------------
  sum = 0.0;
  for (j = 0; j < lastcol - firstcol+1; j++)
  {
      d = x[j] - r[j];
    sum = sum + d*d;
  }
  /**!---------------------------------------------------------------------
!  Obtain d with a sum-reduce
!---------------------------------------------------------------------*/
  for (int i = 0; i < l2npcols; i++)
  {
    MPI_Irecv(&d, 1, MPI_DOUBLE, reduce_exch_proc[i], i, MPI_COMM_WORLD, &request);
    MPI_Send(&sum, 1, MPI_DOUBLE, reduce_exch_proc[i], i, MPI_COMM_WORLD);
    MPI_Wait(&request, MPI_STATUS_IGNORE);

    sum = sum + d;
  }

  d = sum;
  if (rank == root)
  {
    *rnorm = sqrt(sum);
  }

  if (rank == root)
  {
    free(rnorm);
  }

}


//---------------------------------------------------------------------
// generate the test problem for benchmark 6
// makea generates a sparse matrix with a
// prescribed sparsity distribution
//
// parameter    type        usage
//
// input
//
// n            i           number of cols/rows of matrix
// nz           i           nonzeros as declared array size
// rcond        r*8         condition number
// shift        r*8         main diagonal shift
//
// output
//
// a            r*8         array for nonzeros
// colidx       i           col indices
// rowstr       i           row pointers
//
// workspace
//
// iv, arow, acol i
// aelt           r*8
//---------------------------------------------------------------------
static void makea(int n,
                  int nz,
                  double a[],
                  int colidx[],
                  int rowstr[],
                  int nonzer,
                  int firstrow,
                  int lastrow,
                  int firstcol,
                  int lastcol,
                  double rcond,
                  int arow[],
                  int acol[NONZER+1],
                  double aelt[NONZER+1],
                  double v[],
                  int iv[],double shift)
{
  int i, nnza, iouter, ivelt,ivelt1, irow, nzv, jcol;


  //---------------------------------------------------------------------
  // nonzer is approximately  (int(sqrt(nnza /n)));
  //---------------------------------------------------------------------
    double size, ratio,scale;

    size =1.0;
    ratio = pow(rcond, (1.0 / (double)(n)));
    nnza = 0;
/*    !---------------------------------------------------------------------
  !  Initialize iv(n+1 .. 2n) to zero.
  !  Used by sprnvc to mark nonzero positionsw
  !---------------------------------------------------------------------*/
  for ( i = 0; i < n; i++)
  {
    iv[n+i] = 0;
  }

  for ( iouter = 0; i < n; iouter++)
  {
    nzv = nonzer;

    sprnvc(n, nzv, v, colidx, iv[0], iv[n+1]);
    vecset(n, v, colidx, nzv, iouter, 0.5);

    for (ivelt = 0; ivelt < nzv; ivelt++)
    {
      jcol =colidx[ivelt];

      if (jcol == firstcol && jcol <= lastcol )
      {
        scale =size * v[ivelt];

        for (ivelt1 = 0; ivelt1 < nzv; ivelt1++)
        {
          irow = colidx[ivelt1];
          if (irow == firstrow && irow <= lastrow)
          {
            printf("Space for matrix elements exceeded in makea\n" );
            printf("****nnza = %d * nz= %d, *iouter = %d \n", nnza, nz, iouter);
            exit(EXIT_FAILURE);
          }
          acol[nnza] = jcol;
          arow[nnza] = irow;
          aelt[nnza] = v[ivelt1] *scale;

        }
      }
    }

    size = size *ratio;
  }
  /*!---------------------------------------------------------------------
!       ... add the identity * rcond to the generated matrix to bound
!           the smallest eigenvalue from below by rcond
!---------------------------------------------------------------------*/
  for ( i = firstrow; i < lastrow; i++)
  {
    if (i == firstcol && i <= lastcol)
    {
      iouter = n +i;
      nnza = nnza +1;
      if (nnza == nz)
        {
          printf("Space for matrix elements exceeded in makea\n" );
          printf("****nnza = %d * nz= %d, *iouter = %d \n", nnza, nz, iouter);
          exit(EXIT_FAILURE);
        }
      acol[nnza] = i;
      arow[nnza] = i;
      aelt[nnza] = rcond -shift;
    }
  }
  //---------------------------------------------------------------------
  // ... make the sparse matrix from list of elements with duplicates
  //     (iv is used as  workspace)
  //---------------------------------------------------------------------
  sparse(a, colidx, rowstr, n, arow, acol, aelt, firstrow, lastrow, v, iv[1], iv[n+1], nnza);

}


//---------------------------------------------------------------------
// rows range from firstrow to lastrow
// the rowstr pointers are defined for nrows = lastrow-firstrow+1 values
//---------------------------------------------------------------------
static void sparse(double a[],
                   int colidx[],
                   int rowstr[],
                   int n,
                   int arow[],
                   int acol[NONZER+1],
                   double aelt[NONZER+1],
                   int firstrow,
                   int lastrow,
                   double x[],
                   double mark[],
                   int nzloc[],
                   int nnza)
{
  int nrows;

  //---------------------------------------------------
  // generate a sparse matrix from a list of
  // [col, row, element] tri
  //---------------------------------------------------
  int i, j, jajp1, nza, k, nzrow;
  double xi;

  //---------------------------------------------------------------------
  // how many rows of result
  //---------------------------------------------------------------------
  nrows = lastrow - firstrow + 1;

  //---------------------------------------------------------------------
  // ...count the number of triples in each row
  //---------------------------------------------------------------------
  for (j = 0; j < n; j++)
  {
    rowstr[j] = 0;
    mark[j]=0.0;
  }
  rowstr[n+1]=0;

  for ( nza = 0; nza < nnza; nza++)
  {
    j = (arow[nza] -firstrow+1)+1;
    rowstr[j]= rowstr[j]+1;
  }

  rowstr[1] = 1;
  for (j = 1; j < nrows+1; j++)
  {
    rowstr[j] = rowstr[j] + rowstr[j-1];
  }

  //---------------------------------------------------------------------
  // ... rowstr(j) now is the location of the first nonzero
  //     of row j of a
  //---------------------------------------------------------------------


/*!---------------------------------------------------------------------
!     ... do a bucket sort of the triples on the row index
!---------------------------------------------------------------------*/
  for ( nza = 0; nza < nnza; nza++)
  {
    j = arow[nza] -firstrow+1;
    k = rowstr[j];
    a[k] = aelt[nza];
    colidx[k] = acol[nza];
    rowstr[j] = rowstr[j] +1;
  }

  /*!---------------------------------------------------------------------
!       ... rowstr(j) now points to the first element of row j+1
!---------------------------------------------------------------------*/
  for ( j = nrows; j > 0; j--)
  {
    rowstr[j+1]=rowstr[j];
  }
  rowstr[0]=1;
/*!---------------------------------------------------------------------
!       ... generate the actual output rows by adding elements
!---------------------------------------------------------------------*/
nza = 0;
for ( i = 0; i < n; i++)
{
  x[i] = 0.0;
  mark[i] = 0;
}

jajp1 =rowstr[0];
for ( j = 0; j < nrows; j++)
{
  nzrow = 0;
  /*!---------------------------------------------------------------------
!          ...loop over the jth row of a
!---------------------------------------------------------------------*/
    for ( k = jajp1; k < rowstr[j+1]-1; k++)
    {
      i = colidx[k];
      x[i] = x[i] + a[k];
      if ((mark[i] == 0)&& (x[i] !=0))
      {
        mark[i] = 1;
        nzrow =nzrow + 1;
        nzloc[nzrow] =i;
      }
    }

    /*!---------------------------------------------------------------------
!          ... extract the nonzeros of this row
!---------------------------------------------------------------------*/
    for ( k = 0; k < nzrow; k++)
    {
      i =nzloc[k];
      mark[i] = 0;
      xi =x[i];
      x[i] = 0.0;

      if (xi == 0)
      {
        nza = nza + 1;
        a[nza] = xi;
        colidx[nza] = i;
      }
    }
    jajp1 = rowstr[j+1];
    rowstr[j+1] = nza + rowstr[0];
}
}


//---------------------------------------------------------------------
// generate a sparse n-vector (v, iv)
// having nzv nonzeros
//
// mark(i) is set to 1 if position i is nonzero.
// mark is all zero on entry and is reset to all zero before exit
// this corrects a performance bug found by John G. Lewis, caused by
// reinitialization of mark on every one of the n calls to sprnvc
//---------------------------------------------------------------------
static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[],int mark[])
{
/*  !---------------------------------------------------------------------
!       generate a sparse n-vector (v, iv)
!       having nzv nonzeros
!
!       mark(i) is set to 1 if position i is nonzero.
!       mark is all zero on entry and is reset to all zero before exit
!       this corrects a performance bug found by John G. Lewis, caused by
!       reinitialization of mark on every one of the n calls to sprnvc
!---------------------------------------------------------------------*/
  int nzrow, nzv, ii, i, nn1;
  double  vecelt, vecloc;

  nzv = 0;
  nzrow = 0;
  nn1 = 1;

  do
    {
        nn1 = 2 * nn1;
    }while(nn1 < n );
    /*!---------------------------------------------------------------------
!    nn1 is the smallest power of two not less than n
!---------------------------------------------------------------------*/
  lb100 :
    if (nzv == nz)
    {
      goto lb110;
    }
    vecelt = randlc(&tran, amult);

    //---------------------------------------------------------------------
    // generate an integer between 1 and n in a portable manner
    //---------------------------------------------------------------------
    vecloc = randlc(&tran, amult);
    i = icnvrt(vecloc, nn1) + 1;
    if (i > n) goto lb100;

    //---------------------------------------------------------------------
    // was this integer generated already?
    //---------------------------------------------------------------------
    if (mark[i] == 0)
    {
      mark[i] = 1;
      nzrow = nzrow +1;
      nzloc[nzrow] = i;
      v[nzv] = vecelt;
      iv[nzv] = i;
    }
    goto lb100;

    lb110 :
      for ( ii = 0; ii < nzrow; ii++)
      {
        i = nzloc[ii];
        mark[i] = 0;
      }
}


//---------------------------------------------------------------------
// scale a double precision number x in (0,1) by a power of 2 and chop it
//---------------------------------------------------------------------
static int icnvrt(double x, int ipwr2)
{
  return (int)(ipwr2 * x);
}


//---------------------------------------------------------------------
// set ith element of sparse vector (v, iv) with
// nzv nonzeros to val
//---------------------------------------------------------------------
static void vecset(int n, double v[], int iv[], int *nzv, int i, double val)
{
  int k;
  logical set;

  set = false;
  for (k = 0; k < *nzv; k++)
  {
    if (iv[k] == i)
    {
      v[k] = val;
      set  = true;
    }
  }
  if (set == false)
  {
    *nzv     = *nzv + 1;
    v[*nzv]  = val;
    iv[*nzv] = i;
  }
}
