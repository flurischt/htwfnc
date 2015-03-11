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

#include "comp.h"
#include "perf.h"
#include "performance.h"
#include "main.h"

static unsigned int test_sizes[] = {300, 600, 800, 1000};
static unsigned int test_iters[] = {100, 80, 80, 80};

#define NUM_TESTS (sizeof(test_sizes) / sizeof(unsigned int))
#define NUM_ITERS (sizeof(test_iters) / sizeof(unsigned int))

/*
 * Main driver routine - calls register_funcs to get student functions, then
 * tests all functions registered, and reports the best performance
 */
int main(int argc, char **argv)
{
	double perf, maxPerf = 0;
	int i, maxInd = 0;

	print_CPU_info ();

	if (argc >= 2) {
		MHz = atoi(argv[1]);
		printf("CPU frequency fixed at %d MHz\n", MHz);
	} else {
		get_CPU_freq ();
		printf("CPU frequency estimated at %d MHz\n", MHz);
		printf("If that is not the correct value for your CPU, ");
		printf("please specify it as an input argument in MHz\n\n");
	}

	register_functions();

	if (numFuncs == 0) {
		printf("No functions registered - nothing for driver to do\n");
		printf("Register functions by calling register_func(f, name)\n");
		printf("in register_funcs()\n");
		return 0;
	}

	if (NUM_TESTS != NUM_ITERS) {
		printf("The number of tests cases must correspond to the number of iterations\n");
		return 0;
	}

	for (i = 0; i < NUM_ITERS; i++) {
		if (test_iters[i] > MAX_NUM_MEASUREMENTS) {
			printf("The number of iterations per test case can not exceed MAX_NUM_MEASUREMENTS\n");
			return 0;
		}
	}

	for (i = 0; i < numFuncs; i++) {
		perf = perf_test(userFuncs[i], funcNames[i], test_sizes, test_iters, NUM_TESTS);
		if (perf > maxPerf) {
			maxInd = i;
			maxPerf = perf;
		}
	}

	if (maxPerf > 0) {
		printf("Best: %s\nPerf: %.4f F/C\n", funcNames[maxInd], maxPerf);
	} else {
		printf("No valid functions registered\n");
	}

	return 0;
}
