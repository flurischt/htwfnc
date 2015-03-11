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
 *                      Singh Gagandeep    (gsingh@inf.ethz.ch)
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

#include <math.h>
#include "smat.h"
#include "function.h"

#define EPS (1e-8)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* f(x, i, j) multiplies or divides by a sine expression depending
   on i + j being odd or even
 */
double ver_f (smat_t * a, int i, int j) {
	double x = get_elt(a, i, j);
	int t = inc_smat_counter();
	if ( ((i + j) % 3) & 0x1 ) {
		return x / (t + sin(M_PI/(i+1)));
	} else {
		return x * sin(M_PI/(i+1));
	}
}

/* This is the function you need to optimize. It takes one
   square matrix as input
 */
void verfunc(smat_t *a)
{
	int i, j;
	double x,x2;

	reset_smat_counter();

	// i is the column of a we're computing right now
	for(i = 0; i < a->n; i++)
	{
		// j is the row of a we're computing right now
		for(j = 0; j < a->n; j++)
		{
			// First, compute f(A) for the element of a in question
			x = get_elt(a, i, j);
			x = ver_f(a, i, j);

			// Add this to the value of a we're computing and store it
			x2 = get_elt(a, i, j);
			x = x * x2;
			set_elt(a, i, j, x);

		}
	}
}

/*
 * Checks whether or not the results of two computation functions are the
 * same (within a tolerance defined by EPS, to deal with rounding issues
 */
int check_valid(comp_func checkFunc, int chk_n)
{
	int i;

	int valid = 1;
	smat_t *a1 = random_SMat(chk_n);

	double *aMat1 = a1->mat;

	smat_t *a2 = copy_SMat(a1);

	double *aMat2 = a2->mat;

	verfunc(a1 );
	checkFunc(a2);

	for(i = 0; i < chk_n * chk_n; i++) {
		if (isnan(aMat2[i])) {
			valid = 0;
			break;
		}

		if (! (fabs(aMat1[i] - aMat2[i]) < EPS)) {
			valid = 0;
			break;
		}
	}

	free_SMat(a1);
	free_SMat(a2);
	return valid;
}
