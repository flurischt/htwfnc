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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"

void print_complex (complex_t complex) {
	if ( complex.im >= 0 ) {
		printf("(%.2f + i%.2f)", complex.re, complex.im);
	} else {
		printf("(%.2f - i%.2f)", complex.re, -complex.im);
	}
}

complex_t complex_conjugate (complex_t complex) {
	complex_t result = {complex.re, complex.im * (-1)};
	return result;
}

void print_matrix(complex_t m[2][2]) {
	print_complex(m[0][0]); printf(" ");
	print_complex(m[0][1]); printf("\n");
	print_complex(m[1][0]); printf(" ");
	print_complex(m[1][1]); printf("\n");
}
