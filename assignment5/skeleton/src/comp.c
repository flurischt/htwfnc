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
#include "xmmintrin.h"
#include "comp.h"

void simd_conjugate_transpose (complex_t in[2][2], complex_t out[2][2])
{
    //Load
    __m128 row1 = _mm_load_ps((float*)in[0]);
    __m128 row2 = _mm_load_ps((float*)in[1]);
    __m128 multiplier = _mm_set_ps(-1.0, 1.0, -1.0, 1.0);
    //Compute
    row1 = _mm_mul_ps(row1, multiplier);
    row2 = _mm_mul_ps(row2, multiplier);
    //Shuffle/Store
    __m128 conj_row = _mm_shuffle_ps(row1, row2, _MM_SHUFFLE(1, 0, 1, 0));
    _mm_store_ps((float*) out[0], conj_row);
    conj_row = _mm_shuffle_ps(row1, row2, _MM_SHUFFLE(3, 2, 3, 2));
    _mm_store_ps((float*) out[1], conj_row);
}

void simd_pairs_multiplications (float * x, float * y, float * z, size_t n)
{
    size_t i;
    for(i=0;i<=n-4;i+=4)
    {
        //Load two complex numbers from x and two from y
        __m128 ab = _mm_load_ps(x+i); // LSB[a1, b1, a2, b2]
        __m128 cd = _mm_load_ps(y+i); // LSB[c1, d1, c2, d2]
    
        //Compute
        __m128 prod1 = _mm_mul_ps(ab, cd); // LSB[a1*c1, b1*d1, a2*c2, b2*d2]

    }
}

void simd_ceil (float * m, size_t n)
{
	// your code goes here ...
}

