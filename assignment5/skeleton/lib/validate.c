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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "complex.h"
#include "comp.h"
#include "validate.h"

int validate (float * a, float * b, size_t size)
{
	int i, result = 1;
	for (i = 0; i < size; i++) {
		if (fabs((a[i] - b[i]) / a[i]) > EPSILON) {
			result = 0;
			break;
		}
	}
	return result;
}


void sisd_conjugate_transpose (complex_t in[2][2], complex_t out[2][2])
{
	out[0][0] = complex_conjugate(in[0][0]);
	out[0][1] = complex_conjugate(in[1][0]);
	out[1][0] = complex_conjugate(in[0][1]);
	out[1][1] = complex_conjugate(in[1][1]);
}


void sisd_pairs_multiplications (float * x, float * y, float * z, size_t n)
{
	int i;
	for (i = 0; i < 2 * n; i += 2) {
		z[i+0] = x[i+0] * y[i+0] + x[i+1] * y[i+1] * 2;
		z[i+1] = x[i+0] * y[i+1] + x[i+1] * y[i+0];
	}
}

void sisd_ceil(float *m, size_t n)
{
	int i, j;
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++) {
			if (m[n*i + j] <= 1)
				m[n*i + j] = 1;
			else if (m[n*i + j] <= 2)
				m[n*i + j] = 2;
			else
				m[n*i + j] = 3;
		}
	}
}
