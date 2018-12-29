#include "clenshaw.h"

#include <stdio.h>
#include <time.h>

static double cheb_coeffs[128] = { 1.0, 0.5, 7.0, 2.0, 3.0, 4.0, 5.0, 6.0,
    8.80101171e-01, 8.80101171e-01, 8.80101171e-01, 8.80101171e-01, 8.80101171e-01, 8.80101171e-01, 8.80101171e-01, 8.80101171e-01, 
    1.00000000e+00, 5.00000000e-01, 7.00000000e+00, 2.00000000e+00, 3.00000000e+00, 4.00000000e+00, 5.00000000e+00, 6.00000000e+00,              
   -3.91267080e-02, -3.91267080e-02, -3.91267080e-02, -3.91267080e-02, -3.91267080e-02, -3.91267080e-02, -3.91267080e-02, -3.91267080e-02, 
    1.01573120e-16, 4.65373831e-17, 9.04151333e-16, 2.60345101e-16, 5.44075991e-16, 6.76275359e-16, 8.14303449e-16, 1.08964014e-15, 
    4.99515460e-04, 4.99515460e-04, 4.99515460e-04, 4.99515460e-04, 4.99515460e-04, 4.99515460e-04, 4.99515460e-04, 4.99515460e-04, 
    1.66195964e-17, 4.52495641e-17, -2.49029035e-16, -2.89749554e-17, 7.16295185e-17, 1.15597650e-16, 1.30198693e-16, 5.54405952e-17, 
   -3.00465163e-06, -3.00465163e-06, -3.00465163e-06, -3.00465163e-06, -3.00465163e-06, -3.00465163e-06, -3.00465163e-06, -3.00465163e-06, 
   -3.81578178e-17, -1.92583670e-17, -7.15226489e-16, -1.19679883e-16, -4.68070633e-16, -5.01924168e-16, -4.81380344e-16, -7.72880747e-16, 
    1.04985004e-08, 1.04985003e-08, 1.04985003e-08, 1.04985005e-08, 1.04985005e-08, 1.04985004e-08, 1.04985004e-08, 1.04985005e-08, 
    1.48029737e-17, -2.22044605e-17, -3.55271368e-16, 5.92118946e-17, 2.96059473e-17, -3.55271368e-16, -6.51330841e-16, -5.92118946e-17, 
   -2.39601459e-11, -2.39601629e-11, -2.39606408e-11, -2.39601734e-11, -2.39601769e-11, -2.39604196e-11, -2.39602854e-11, -2.39604759e-11, 
   -2.96059473e-17, -3.70074342e-17, -2.66453526e-16, -1.85037171e-17, 7.40148683e-18, -4.81096644e-16, -3.92278802e-16, -2.22044605e-16, 
    3.84877315e-14, 3.84285196e-14, 3.88430029e-14, 3.83693077e-14, 3.86061553e-14, 3.88430029e-14, 3.88430029e-14, 3.83693077e-14, 
   -1.18423789e-16, 0.00000000e+00, -8.28966525e-16, -2.36847579e-16, -3.55271368e-16, -5.92118946e-16, -7.10542736e-16, -4.73695157e-16, 
   -1.77635684e-16, -5.92118946e-17, 0.00000000e+00, -1.18423789e-16, -2.36847579e-16, 0.00000000e+00, 0.00000000e+00, -4.73695157e-16 };

/* Allocate an aligned block of memory using whatever API is known
 * to the compiler. */
static inline double *
alloc_tmp(size_t size)
{
    double * ptr = NULL;
#ifdef HAVE_POSIX_ALLOC
    posix_memalign((void **)&ptr, ALIGNMENT, ALIGN_VECTOR(size));
#elif defined HAVE_C11_ALLOC
    ptr = aligned_alloc(ALIGNMENT, ALIGN_VECTOR(size));
#elif defined HAVE_WINDOWS_ALLOC
    ptr = _aligned_malloc(ALIGN_VECTOR(size), ALIGNMENT);
#else
#error "Your C compiler can't deal with aligned allocation. Sorry!"
#endif
    return ptr;
}

int main(int argc, char ** argv)
{
  const size_t ylen = 8;
  const size_t clen = 16;
  double x = -1;
  double * out = alloc_tmp(ylen);
  double * coeffs = alloc_tmp(ylen * clen);

  printf("aligned ptr: %p stack ptr: %p\n",  coeffs, &cheb_coeffs);
  for (int i = 0; i < (ylen * clen); ++i)
  {
    coeffs[i] = cheb_coeffs[i];
    /* printf("coeff %d: %.10e\n", i, coeffs[i]); */
  }

  clock_t start, end;
  start = clock();
  int call_count = 0;
  while (x < 1)
  {
    clenshaw(out, ylen, &x, clen, coeffs);
    x = x + 0.000001;
    call_count++;
  }
  end = clock();
  double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("time taken: %.3fs with stride: %d\n", cpu_time_used, STRIDE);


  printf("called %d times with x = %.3f\n", call_count, x);
  for (int i = 0; i < ylen; ++i)
  {
    printf("out %d: %.10e\n", i, out[i]);
  }

  free(coeffs);
  free(out);

  return 0;
}
