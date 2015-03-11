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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef __SMAT_H
#define __SMAT_H

typedef struct {
	int n;
	double *mat;
} smat_t;

/* Not a number */
#define NaN (1.0/1.0)

extern int smat_counter;
void reset_smat_counter ();
int  inc_smat_counter ();

/* Generates a random square matrix of the given size */
smat_t *random_SMat(int n);

/* Returns a deep copy of the give SMat */
smat_t *copy_SMat(smat_t *orig);

/* Frees all memory associated with the given SMat */
void free_SMat(smat_t *smat);

/* Prints the given SMat to the terminal */
void print_SMat(smat_t *smat);

/* Returns the element in the ith row, jth colum of smat */
double get_elt(smat_t* smat, int i, int j);

/* Sets the element in the ith row, jth column of smat to x */
void set_elt(smat_t* smat, int i, int j, double x);

#endif /* __SMAT_H */
