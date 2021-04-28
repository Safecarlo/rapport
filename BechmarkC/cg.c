/*--------------------------------------------------------------------
	Information on NAS Parallel Benchmarks is available at:
	http://www.nas.nasa.gov/Software/NPB/
	Authors: M. Yarrow
           C. Kuszmaul
	CPP version:
  			Dalvan Griebler <dalvangriebler@gmail.com>
			Júnior Löff <loffjh@gmail.com>
--------------------------------------------------------------------*/

/*
c---------------------------------------------------------------------
c  Note: please observe that in the routine conj_grad three
c  implementations of the sparse matrix-vector multiply have
c  been supplied.  The default matrix-vector multiply is not
c  loop unrolled.  The alternate implementations are unrolled
c  to a depth of 2 and unrolled to a depth of 8.  Please
c  experiment with these to find the fastest for your particular
c  architecture.  If reporting timing results, any of these three may
c  be used without penalty.
c---------------------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>

#include "npbparams.h"

#define NZ  NA*(NONZER+1)*(NONZER+1)+NA*(NONZER+2)
/*global variables*/

/*common /partit_size/*/
static int naa;
static int nzz;
static int firstrow;
static int lastrow;
static int firtcol;
static int lastcol;

/* common /main_int_mem/ */
static int colidx[NZ+1];	/* colidx[1:NZ] */
static int rowstr[NA+1+1];	/* rowstr[1:NA+1] */
static int iv[2*NA+1+1];	/* iv[1:2*NA+1] */
static int arow[NZ+1];		/* arow[1:NZ] */
static int acol[NZ+1];		/* acol[1:NZ] */

/* common /main_flt_mem/ */
static double v[NA+1+1];	/* v[1:NA+1] */
static double aelt[NZ+1];	/* aelt[1:NZ] */
static double a[NZ+1];		/* a[1:NZ] */
static double x[NA+2+1];	/* x[1:NA+2] */
static double z[NA+2+1];	/* z[1:NA+2] */
static double p[NA+2+1];	/* p[1:NA+2] */
static double q[NA+2+1];	/* q[1:NA+2] */
static double r[NA+2+1];	/* r[1:NA+2] */
static double w[NA+2+1];	/* w[1:NA+2] */

/* common /urando/ */
static double amult;
static double tran;

/* function declarations */
static void conj_grad (int colidx[], int rowstr[], double x[],
	double z[], double a[], double p[], double q[], double r[],
	double w[], double *rnorm);
static void makea(int n, int nz, double a[], int colidx[], int rowstr[],
	int nonzer, int firstrow, int lastrow, int firstcol,
	int lastcol, double rcond, int arow[], int acol[],
	double aelt[], double v[], int iv[], double shift );
static void sparse(double a[], int colidx[], int rowstr[], int n,
	int arow[], int acol[], double aelt[],
	int firstrow, int lastrow,
	double x[], boolean mark[], int nzloc[], int nnza);
static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[], int mark[]);
static int icnvrt(double x, int ipwr2);
static void vecset(int n, double v[], int iv[], int *nzv, int i, double val);

/*--------------------------------------------------------------------
      program cg
--------------------------------------------------------------------*/

int main(int argc, char const *argv[])
{
  int	i, j, k, it;
	double zeta;
	double rnorm;
	double norm_temp11;
	double norm_temp12;
	double t, mflops;
	char class_npb;
	boolean verified;
	double zeta_verify_value, epsilon;

	firstrow = 1;
	lastrow  = NA;
	firstcol = 1;
	lastcol  = NA;
  if (NA == 1400 && NONZER == 7 && NITER == 15 && SHIFT == 10.0) {
		class_npb = 'S';
		zeta_verify_value = 8.5971775078648;
	} else if (NA == 7000 && NONZER == 8 && NITER == 15 && SHIFT == 12.0) {
		class_npb = 'W';
		zeta_verify_value = 10.362595087124;
	} else if (NA == 14000 && NONZER == 11 && NITER == 15 && SHIFT == 20.0) {
		class_npb = 'A';
		zeta_verify_value = 17.130235054029;
	} else if (NA == 75000 && NONZER == 13 && NITER == 75 && SHIFT == 60.0) {
		class_npb = 'B';
		zeta_verify_value = 22.712745482631;
	} else if (NA == 150000 && NONZER == 15 && NITER == 75 && SHIFT == 110.0) {
		class_npb = 'C';
		zeta_verify_value = 28.973605592845;
	} else {
		class_npb = 'U';
	}

	printf("NAS Parallel Benchmarks 4.0 PARALLEL version with MPI" " - CG Benchmark\n");
  	printf("Developed by:  Hery ANDRIANANTENAINA\n\n");
  	printf(" Size: %10d\n", NA);
	printf(" Iterations: %5d\n", NITER);

	naa = NA;
	nzz = NZ;
  return 0;
}
