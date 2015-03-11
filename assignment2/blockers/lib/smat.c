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

#include "smat.h"

double d_rand()  {
    return (double)rand() / (double)RAND_MAX;
}

int smat_counter = 0;

void reset_smat_counter () {
	smat_counter = 0;
}

int inc_smat_counter () {
	smat_counter += 2;
	return smat_counter++;
}

/* returns a random SMat with entries in [-1, 1] */
smat_t *random_SMat(int n)
{
    int i, j;
    smat_t *smat = malloc(sizeof(smat_t));
    double *mat = malloc(n * n * sizeof(double));
    for (i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            mat[i * n + j] = (d_rand() * 2) - 1;
        }
    }
    smat->mat = mat;
    smat->n = n;
    return smat;
}

/* Returns a deep copy of the give SMat */
smat_t *copy_SMat(smat_t *orig)
{
    int n = orig->n;
    size_t matSize = n * n * sizeof(double);
    smat_t *new = malloc(sizeof(smat_t));
    new->n = orig->n;
    new->mat = malloc(matSize);
    memcpy(new->mat, orig->mat, matSize);
    return new;
}

/* Frees all memory associated with the given SMat */
void free_SMat(smat_t *smat)
{
    free(smat->mat);
    free(smat);
}

/* Prints the given SMat to the terminal */
void print_SMat(smat_t *smat)
{
    int i, j;
    int n = smat->n;
    double *mat = smat->mat;
    
    printf("%d x %d Matrix:\n\n", n, n);

    for ( i= 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%1.3f\t", mat[i * n + j]);
        }
        printf("\n");
    }
}

/* Returns the element in the ith row, jth colum of smat */
double get_elt(smat_t* smat, int i, int j)
{
    double *mat = smat->mat;
    int n = smat->n;
    if (i >= n || i < 0 || j >= n || j < 0)
        return NaN;
    return mat[i * n + j];
}

/* Sets the element in the ith row, jth column of smat to x */
void set_elt(smat_t* smat, int i, int j, double x)
{
    double *mat = smat->mat;
    int n = smat->n;
    
    if(i >= 0 && i < n && j >= 0 && j < n)
        mat[i * n + j] = x;
    else
        printf("Set failed\n");
}
