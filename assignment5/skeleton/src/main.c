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
#include <string.h>

#include "main.h"
#include "comp.h"
#include "complex.h"
#include "validate.h"

void test_conjugate_transpose ()
{
	printf("Testing SIMD Conjugate Transpose function\n");
	printf("-----------------------------------------\n");
	int i; complex_t m[2][2], m_simd[2][2], m_sisd[2][2];
	for (i = 0; i < 10; i++) {
		m[0][0] = (complex_t) {RANDOM, RANDOM};
		m[0][1] = (complex_t) {RANDOM, RANDOM};
		m[1][0] = (complex_t) {RANDOM, RANDOM};
		m[1][1] = (complex_t) {RANDOM, RANDOM};
		simd_conjugate_transpose(m, m_simd);
		sisd_conjugate_transpose(m, m_sisd);
		if (validate((float *)&m_simd[0][0],  (float *)&m_sisd[0][0], 8)) {
			printf("Test %02d: OK\n", i+1);
		} else {
			printf("Test %02d: FAIL\n", i+1);
			// print_matrix(m_simd);
			// print_matrix(m_sisd);
		}
	}
}


void test_paris_multiplications ()
{
	printf("Testing SIMD Pairs Multiplication function\n");
	printf("-----------------------------------------\n");
	int i, j, n = 64;
	float x[2 * n], y[2 * n], z_simd[2 * n], z_sisd[2 * n];

	for (i = 0; i < 10; i++) {
		for (j = 0; j < 2 * n; ++j)  {
			x[j] = (float)(rand())/RAND_MAX;
			y[j] = (float)(rand())/RAND_MAX;;
		}
		simd_pairs_multiplications(x, y, z_simd, n);
		sisd_pairs_multiplications(x, y, z_sisd, n);
		if (validate(z_simd,  z_sisd, 2 * n)) {
			printf("Test %02d: OK\n", i+1);
		} else {
			printf("Test %02d: FAIL\n", i+1);
		}

	}
}

void test_ceil ()
{
	printf("Testing SIMD Ceil function\n");
	printf("-----------------------------------------\n");
	int i, j, n = 32;
	float m_simd[n * n], m_sisd[n * n], temp, max = RAND_MAX / 3;
	for(i = 0; i < 10; i++){
		for(j = 0; j < n * n; j++){
			temp = (float) (rand()) / max;
			m_simd[j] = temp;
			m_sisd[j] = temp;
		}
		simd_ceil(m_simd, n);
		sisd_ceil(m_sisd, n);
		if (validate(m_simd, m_sisd, n * n)) {
			printf("Test %02d: OK\n", i+1);
		} else {
			printf("Test %02d: FAIL\n", i+1);
		}
	}
}

void test_ceil_unaligned()
{
	printf("Testing SIMD unaligned function\n");
	printf("-----------------------------------------\n");
	int i, j, n = 1000, k;
    char mem[4*n*n+16];
	float *m_simd, m_sisd[n * n], temp, max = RAND_MAX / 3;
    k=1;
    for(k=0;k<=16;k++)
    {
        // simulate any possible alignment
        m_simd = (float*)&mem[k];
        for(i = 0; i < 10; i++){
            for(j = 0; j < n * n; j++){
                temp = (float) (rand()) / max;
                m_simd[j] = temp;
                m_sisd[j] = temp;
            }
            simd_ceil(m_simd, n);
            sisd_ceil(m_sisd, n);
            if (validate(m_simd, m_sisd, n * n)) {
//                printf("Test %02d: OK\n", i+1);
            } else {
                printf("Test %02d: FAIL\n", i+1);
            }
        }
    }
}

int main(int argc, char **argv)
{
	srand (21169);
	test_conjugate_transpose ();
	printf("\n");
	test_paris_multiplications ();
	printf("\n");
	test_ceil();
	printf("\n");
    test_ceil_unaligned();
    printf("\n");
	return 0;
}
