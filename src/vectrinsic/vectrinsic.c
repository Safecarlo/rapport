#include <stdio.h>
#include <stdlib.h>
#include <emmintrin.h>

// Vector type
typedef double double2 __attribute__((ext_vector_type(2)));

// Simple test for vector type and intrinsics
int main()
{
  // Init
  double2 a = {0.56125898, 0.6797};
  double2 b = {0.96769, 0.66456343};

  // Compute
  __m128d res = _mm_add_pd(*(double2 *)(void *)&a, *(double2 *)(void *)&b);

  // Print value
  printf("first  : %lf\n", res[0]);
  printf("second : %lf\n", res[1]);
  
  // Print difference
  printf("diff first  : %lf\n", res[0] - (a[0] + b[0]));
  printf("diff second : %lf\n", res[1] - (a[1] + b[1]));
  
  return 0;
}
