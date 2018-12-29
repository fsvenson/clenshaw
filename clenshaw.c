/*=================================================================
 * clenshaw.c
 *
 * Evaluate a polynomial in Chebyshev basis over the interval
 * [-1, 1] at a vector of evaluation points. Called from matlab
 * with two inputs X and C, representing the vector of points and
 * a matrix of coefficients, with one polynomial per column.
 * Returns a vector Y with the same number of rows as X and the
 * same number of columns as C.
 *
 * Input:   vector X, matrix C
 * Output:  matrix Y
 *
 *=================================================================*/
#include "clenshaw.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/**
 * @param[out] y: output point
 * @param[in]  x: input point
 * @param[in]  c: coefficients
 *
 * Both px and py should be pointers to SIMD-vector-size aligned arrays of
 * doubles, of length xlen. The array of coefficients pc need to be aligned.
 * Use function alloc_tmp to ensure correct alignment.
 *
 * @return:
 */
void
clenshaw(double *py,
         size_t ylen, double *px,
         size_t clen, double *pc)
{
    int64_t i, k;
    unsigned d;
    vdouble x, x2;
    vdouble b1, b2;
    vdouble ck, ck1;
    double zero = 0.0;
    vdouble b = broadcast_val(zero); 
    double b0 = 0;
    const double two = 2.0;

    /* Degree is one less than the number of coefficients. */
    d = clen - 1;
    bool odd = false;
    if (d & 1) {
        /* For odd degree polynomial, do one assignment outside the loop. */
        d--;
        odd = true;
    }

    x = broadcast_val(*px);
    x2 = mul_pd(x, broadcast_val(two));
    /* Loop over the polynomials, each time taking as many coeffs as
     * fit in a SIMD vector. */
    for (i = 0; i < (int64_t) ylen; i += STRIDE) {
        b2 = broadcast_val(b0);
        if (odd) {
          /* For odd degree polynomial, do one assignment outside the loop. */
          b1 = load_vector(pc[ylen * d + i]);
        }
        else {
          b1 = b;
        }

        for (k = d; k > 0; k -= 2) {
            /* b2 = 2*x * b1 - b2 + c[k]
             * b1 = 2*x * b2 - b1 + c[k-1] */
            ck = load_vector(pc[ylen * k + i]);
            ck1 = load_vector(pc[ylen * (k - 1) + i]);

            b2 = add_pd(fmsub_pd(x2, b1, b2), ck);
            b1 = add_pd(fmsub_pd(x2, b2, b1), ck1);
        }
        b2 = add_pd(fmsub_pd(x, b1, b2),
                    load_vector(pc[i]));
        /* b2 = fmsub_pd(x, b1, b2); */

        store_vector(py[i], b2);
    }
}

/* void */
/* clenshaw_complex(double *pyr, double *pyi, */
/*                  size_t xlen, double *pxr, double *pxi, */
/*                  size_t clen, double *pcr, double *pci) */
/* { */
/*     int i, k; */
/*     unsigned d; */
/*     vdouble xr, xi, x2r, x2i; */
/*     vdouble b1r, b1i, b2r, b2i; */

/*     double br, bi, b0 = 0; */
/*     const double two = 2.0; */

/*     d = clen - 1; */
/*     if (d & 1) { */
/*         br = pcr[d]; */
/*         bi = pci[d]; */
/*         d--; */
/*     } */

/*     for (i = 0; i < (int64_t) xlen; i += STRIDE) { */
/*         xr = load_vector(pxr[i]); */
/*         xi = load_vector(pxi[i]); */
/*         x2r = mul_pd(xr, broadcast_val(two)); */
/*         x2i = mul_pd(xi, broadcast_val(two)); */
/*         b2r = broadcast_val(b0); */
/*         b2i = broadcast_val(b0); */
/*         b1r = broadcast_val(br); */
/*         b1i = broadcast_val(bi); */

/*         for (k = d; k > 0; k -= 2) { */
/*            /1* b2 = 2*x * b1 - b2 + c[k] */
/*             * b1 = 2*x * b2 - b1 + c[k-1] */
/*             * */
/*             * In complex arithmetic, we have: */
/*             * re(b2) = re(2*x)*re(b1) - im(2*x)*im(b1) - re(b2) + re(c[k]) */
/*             * im(b2) = re(2*x)*im(b1) + im(2*x)*re(b1) - im(b2) + im(c[k]) */
/*             * re(b1) = re(2*x)*re(b2) - im(2*x)*im(b2) - re(b1) + re(c[k-1]) */
/*             * im(b1) = re(2*x)*im(b2) + im(2*x)*re(b2) - im(b1) + im(c[k-1]) *1/ */
/*             b2r = sub_pd(fmadd_pd(x2r, b1r, broadcast_val(pcr[k])), */
/*                          fmadd_pd(x2i, b1i, b2r)); */
/*             b2i = add_pd(fmadd_pd(x2r, b1i, broadcast_val(pci[k])), */
/*                          fmsub_pd(x2i, b1r, b2i)); */
/*             b1r = sub_pd(fmadd_pd(x2r, b2r, broadcast_val(pcr[k-1])), */
/*                          fmadd_pd(x2i, b2i, b1r)); */
/*             b1i = add_pd(fmadd_pd(x2r, b2i, broadcast_val(pci[k-1])), */
/*                          fmsub_pd(x2i, b2r, b1i)); */
/*         } */
/*         b2r = sub_pd(fmadd_pd(xr, b1r, broadcast_val(pcr[0])), */
/*                      fmadd_pd(xi, b1i, b2r)); */
/*         b2i = add_pd(fmadd_pd(xr, b1i, broadcast_val(pci[0])), */
/*                      fmsub_pd(xi, b1r, b2i)); */

/*         store_vector(pyr[i], b2r); */
/*         store_vector(pyi[i], b2i); */
/*     } */
/* } */
