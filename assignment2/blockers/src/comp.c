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
#include <stdio.h>
#include "comp.h"
#include "perf.h"
#include "performance.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/**
 * f(x, i, j) multiplies or divides by a sine expression depending
 * on (i + j) % 3 being odd or even. It uses the incremental smat
 * counter to perform the devision.
 */
double f (smat_t * a, int i, int j) {
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
void superslow(smat_t *a)
{
    int i, j;
    double x,x2;
    reset_smat_counter ();
    // i is the column of a we're computing right now
    for(i = 0; i < a->n; i++) {
        // j is the row of a we're computing right now
        for(j = 0; j < a->n; j++) {
            // First, compute f(A) for the element of a in question
            x = f (a, i, j);
            // Add this to the value of a we're computing and store it
            x2 = get_elt(a, i, j);
            x = x * x2;
            set_elt(a, i, j, x);
        }
    }
}

extern int smat_counter;

//
// changes to superslow() 
//  - no call to f()
//  - no calls to get_elt and set_elt 
//        (and therefore unnecessary no boundary checks)
//
void superslow_inlined(smat_t *a)
{
    int i, j;
    double x,x2;
    reset_smat_counter ();
    // i is the column of a we're computing right now
    for(i = 0; i < a->n; i++) {
        // j is the row of a we're computing right now
        for(j = 0; j < a->n; j++) {
            // First, compute f(A) for the element of a in question
            // #####BEGIN previous f() call
            //x = f (a, i, j);
            //x = get_elt(a, i, j);
            x = a->mat[i * a->n + j];
            int t = inc_smat_counter();
            if ( ((i + j) % 3) & 0x1 ) {
                x = x / (t + sin(M_PI/(i+1)));
            } else {
                x = x * sin(M_PI/(i+1));
            }
            // #####END f() call
            // Add this to the value of a we're computing and store it
            //x2 = get_elt(a, i, j);
            x2 = a->mat[i * a->n + j];
            x = x * x2;
            //set_elt(a, i, j, x);
            a->mat[i * a->n + j] = x;
        }
    }
}

extern int smat_counter;

// changes between superslow_inlined and trigo():
// - compute sin() only once for each i
// - remove *_smat_counter() calls
//
void trigo(smat_t *a)
{
    int i, j;
    double x,x2;
    smat_counter = 0;
    double sine_i = 0;
    // i is the column of a we're computing right now
    for(i = 0; i < a->n; i++) {
        sine_i = sin(M_PI/(i+1));
        // j is the row of a we're computing right now
        for(j = 0; j < a->n; j++) {
            // First, compute f(A) for the element of a in question
            x2 = x = a->mat[i * a->n + j];
            smat_counter += 2;
            int t = smat_counter++;
            if ( ((i + j) % 3) & 0x1 ) {
                x = x / (t + sine_i);
            } else {
                x = x * sine_i;
            }
            // Add this to the value of a we're computing and store it
            x = x * x2;
            a->mat[i * a->n + j] = x;
        }
    }
}

void unroll1(smat_t *a)
{
    int i, j;
    double x,x2;
    double x_2, x2_2; // I know, stupid variable naming
    double x_3, x2_3; // even worse... (it was late...)
    smat_counter = 0;
    double sine_i = 0;
    // i is the column of a we're computing right now
    for(i = 0; i < a->n; i++) {
        sine_i = sin(M_PI/(i+1));
        // j is the row of a we're computing right now
        for(j = 0; j < a->n-2; j+=3) {
            // First, compute f(A) for the element of a in question
            x2 = x = a->mat[i * a->n + j];
            x2_2 = x_2 = a->mat[i * a->n + j + 1];
            x2_3 = x_3 = a->mat[i * a->n + j + 2];
            smat_counter += 2;
            int t = smat_counter++;
            // x
            if ( ((i + j) % 3) & 0x1 ) {
                x = x / (t + sine_i);
                x_2 = x_2 * sine_i;
                x_3 = x_3 * sine_i;
                smat_counter += 6;
            } else if ( ((i + j + 1) % 3) & 0x1 ) {
                x = x * sine_i;
                smat_counter += 2;
                t = smat_counter++;
                x_2 = x_2 / (t + sine_i);
                x_3 = x_3 * sine_i;
                smat_counter += 3;
            } else {
                x = x * sine_i;
                x_2 = x_2 * sine_i;
                smat_counter += 5;
                t = smat_counter;
                x_3 = x_3 / (t + sine_i);
                smat_counter++;
            }
            // Add this to the value of a we're computing and store it
            x = x * x2;
            x_2 = x_2 * x2_2;
            x_3 = x_3 * x2_3;
            a->mat[i * a->n + j] = x;
            a->mat[i * a->n + j+1] = x_2;
            a->mat[i * a->n + j+2] = x_3;
        }
        // calculate remaining elements
        for(; j < a->n; j++) {
            // First, compute f(A) for the element of a in question
            x2 = x = a->mat[i * a->n + j];
            smat_counter += 2;
            int t = smat_counter++;
            if ( ((i + j) % 3) & 0x1 ) {
                x = x / (t + sine_i);
            } else {
                x = x * sine_i;
            }
            // Add this to the value of a we're computing and store it
            x = x * x2;
            a->mat[i * a->n + j] = x;
        }
    }
}

/**
 * Called by the driver to register your functions
 * Use add_function(func, description) to add your own functions
 */
void register_functions()
{
    // Registers comp_superslow with the driver
    add_function(&superslow, "superslow: original function");
    // my functions
    add_function(&superslow_inlined, "inlined f() and direct array access");
    add_function(&trigo, "trigo(), compute sin() only n times, and stop using *smat_counter()");
    add_function(&unroll1, "unroll");
}

