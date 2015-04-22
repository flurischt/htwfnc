/**
 *      _________   _____________________  ____  ______
 *     / ____/   | / ___/_  __/ ____/ __ \/ __ \/ ____/
 *    / /_  / /| | \__ \ / / / /   / / / / / / / __/
 *   / __/ / ___ |___/ // / / /___/ /_/ / /_/ / /___
 *  /_/   /_/  |_/____//_/  \____/\____/_____/_____/
 *
 *  http://www.inf.ethz.ch/personal/markusp/teaching/
 *  How to Write Fast Numerical Code 263-2300 - ETH Zurich
 *  Copyright (C) 2015  Alen Stojanov      (astojanov@inf.ethz.ch)
 *                      Daniele Spampinato (daniele.spampinato@inf.ethz.ch)
 *                      Gagandeep Singh    (gsingh@inf.ethz.ch)
 *	                    Markus Pueschel    (pueschel@inf.ethz.ch)
 *	Copyright (C) 2013  Georg Ofenbeck     (ofenbeck@inf.ethz.ch)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see http://www.gnu.org/licenses/.
 */

#include <stdlib.h>
#include <stdio.h>
#include <x86intrin.h>
#include "comp.h"

#define USE_UNARY_MINUS_VERSION 1

void simd_conjugate_transpose (complex_t in[2][2], complex_t out[2][2])
{
    //Load    
    __m128 row1 = _mm_loadu_ps((float*)in[0]);
    __m128 row2 = _mm_loadu_ps((float*)in[1]);
#if USE_UNARY_MINUS_VERSION
    // try to avoid a 5-cycle _mm_mul call. use shuffle and unary-minus
    __m128 re = _mm_shuffle_ps(row1, row2, _MM_SHUFFLE(2, 0, 2, 0)); //[re1 re2 re3 re4]
    __m128 im = _mm_shuffle_ps(row1, row2, _MM_SHUFFLE(3, 1, 3, 1));// [im1 im2 im3 im4]
    im = -im;
    __m128 conj_row1 = _mm_shuffle_ps(re, im, _MM_SHUFFLE(2, 0, 2, 0)); // [re1 re3 im1 im3]
    __m128 conj_row2 = _mm_shuffle_ps(re, im, _MM_SHUFFLE(3, 1, 3, 1)); // [re2, re4, im2, im4]
    conj_row1 = _mm_shuffle_ps(conj_row1, conj_row1, _MM_SHUFFLE(3, 1, 2, 0)); // [re1, im1, re3, im3]
    conj_row2 = _mm_shuffle_ps(conj_row2, conj_row2, _MM_SHUFFLE(3, 1, 2, 0));
#else
    // straightforward implementation using a multiplier
    __m128 multiplier = _mm_set_ps(-1.0, 1.0, -1.0, 1.0);
    //Compute
    row1 = _mm_mul_ps(row1, multiplier);
    row2 = _mm_mul_ps(row2, multiplier);
    //Shuffle/Store
    __m128 conj_row1 = _mm_shuffle_ps(row1, row2, _MM_SHUFFLE(1, 0, 1, 0));
    __m128 conj_row2 = _mm_shuffle_ps(row1, row2, _MM_SHUFFLE(3, 2, 3, 2));
#endif
    _mm_storeu_ps((float*)out[0], conj_row1);
    _mm_storeu_ps((float*)out[1], conj_row2);
}

// a small helper to debug __m128 datatypes
// they are all output as floats
void debug_m128(char* msg, __m128 out)
{
    float *a = malloc((sizeof(float)*4));
    printf(msg);
    printf("[");
    _mm_store_ps(a, out);
    int i;
    for(i=0;i<4;i++)
    {
        printf(" %f", a[i]);
    }
    printf(" ]\n");
}

void simd_pairs_multiplications (float * x, float * y, float * z, size_t n)
{
    __m128 ab, cd, dc, re_part, im_part, result;
    __m128 multiplier = _mm_set_ps(2.0, 1.0, 2.0, 1.0);
    size_t i;
    for(i=0;i<=2*n-4;i+=4)
    {
        //Load two complex numbers from x and two from y
        ab = _mm_loadu_ps(x+i); // LSB[a1, b1, a2, b2]
        cd = _mm_loadu_ps(y+i); // LSB[c1, d1, c2, d2]
        dc = _mm_shuffle_ps(cd, cd, _MM_SHUFFLE(2, 3, 0, 1)); // LSB[d1, c1, d2, c2]
        //Compute
        re_part = _mm_mul_ps(ab, cd); // LSB[a1*c1, b1*d1, a2*c2, b2*d2]
        im_part = _mm_mul_ps(ab, dc); // LSB[a1*d1, b1*c1, a2*d2, b2*c2]
        // TODO there's a dependency between the two re_part = multiplicatoins...
        re_part = _mm_mul_ps(re_part, multiplier); // LSB[a1c1, 2*b1d1, a2c2, 2*b2d2]
        result = _mm_hadd_ps(re_part, im_part);
        result = _mm_shuffle_ps(result, result, _MM_SHUFFLE(3, 1, 2,0));
        //Store
        _mm_storeu_ps(z+i, result);
    }
}

void simd_ceil (float * m, size_t n)
{
    __m128 input, ones, twos, threes, combined;
    __m128 one_constants = _mm_set1_ps(1);
    __m128 two_constants = _mm_set1_ps(2);
    __m128 three_constants = _mm_set1_ps(3);
    size_t i;

    // check for alignment, peel the loop if necessary
    size_t peel = ((unsigned long)m) & 0xf;
    if(peel != 0)
    {
        // if more than 4 elements, try to vectorize, but unaligned
        peel = 16 - peel;
        peel = (peel < n*n) ? peel : n*n;
        for(i=0;i<=peel-4;i+=4)
        {
            input = _mm_loadu_ps(m+i);
            ones = _mm_cmple_ps(input, one_constants);
            twos = _mm_cmple_ps(input, two_constants);
            combined = _mm_blendv_ps(two_constants, one_constants, ones);
            threes = _mm_or_ps(ones, twos);
            combined = _mm_blendv_ps(three_constants, combined, threes);
            _mm_storeu_ps(m+i, combined);
        }
        // finish remaining elements
        for(;i<peel;i++)
            *(m+i) = (*(m+i) < 1) ? 1 : (*(m+i) < 2) ? 2 : 3;
        peel--;
    }
    // now compute aligned
    float * next_elems = m+peel;
    size_t num_elements = n*n;
    size_t num_vectorizable_elements = num_elements -4;
    for(i=peel;i<=num_vectorizable_elements;i+=4,next_elems+=4)
    {
        //printf("i=%i, n=%i\n", i, num_vectorizable_elements);
        input = _mm_load_ps(next_elems);
        //printf("load_success!\n");
        ones = _mm_cmple_ps(input, one_constants);
        twos = _mm_cmple_ps(input, two_constants);
        combined = _mm_blendv_ps(two_constants, one_constants, ones);
        threes = _mm_or_ps(ones, twos);
        combined = _mm_blendv_ps(three_constants, combined, threes);
        _mm_store_ps(next_elems, combined);
        //printf("store_success!\n");
    }
    // and finish remaining 1,2 or 3 elements
    for(;i<num_elements;i++)
        *(m+i) = (*(m+i) < 1) ? 1 : (*(m+i) < 2) ? 2 : 3;
}

