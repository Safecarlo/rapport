#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

// Vector type
typedef double double4 __attribute__((ext_vector_type(4)));

// Simple test for vector type and intrinsics
int main()
{  
  // Init
  double4 a = {0.56125898, 0.1};
  double4 b = {0.96769, 0.66456343};

  // Compute
#ifdef __AVX__
  __m256d res = _mm256_add_pd(*(__m256d *)(void *)&a, *(__m256d *)(void *)&b);
#else
  __m256d res;
  res[0] = a[0] + b[0];
  res[1] = a[1] + b[1];
#endif

  // Print value
  printf("first  : %lf\n", res[0]);
  printf("second : %lf\n", res[1]);
  
  // Print difference
  printf("diff first  : %lf\n", res[0] - (a[0] + b[0]));
  printf("diff second : %lf\n", res[1] - (a[1] + b[1]));
  
  return 0;
}
