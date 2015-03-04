#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "tsc_x86.h"
#define FREQUENCY 2.6e9
#define CYCLES_REQUIRED 1e8
#define min(a,b) a < b ? a
#define NUM_RUNS 1
#define CALIBRATE

#define NUM_MEASUREMENTS_PER_N 21

int n=0;
double *x = 0;
double *y = 0;
double *z = 0;
double *w = 0;

// allocates and initalizes the arrays
// x,y and w are filled with random values
// z is zeroed
// sets the global n to size
void arr_init(int size)
{
    if(n != size)
    {
        n = size;
        if(x != 0)
            free(x);
        if(y != 0)
            free(y);
        if(z != 0)
            free(z);
        if(w != 0)
            free(w);

        x = (double*)malloc(size*sizeof(double));
        y = (double*)malloc(size*sizeof(double));
        z = (double*)malloc(size*sizeof(double));
        w = (double*)malloc(size*sizeof(double));
    }
    int i;
    for(i=0;i<n;i++)
    {
        x[i] = rand();
        y[i] = rand();
        w[i] = rand() / ((double) RAND_MAX); // get a rand out of the range [0, 1]
        z[i] = 0;
    }
#ifdef DEBUG
    //output a generated table
    printf("x \t y \t w \t z\n");
    for(i=0;i<n;i++)
        printf("%lf \t %lf \t %lf \t %lf\n", x[i], y[i], w[i], z[i]);
#endif
}

void compute()
{
    int i;
    for(i=0;i<n;i++)
    {
        z[i] = w[i] * x[i] + (1-w[i])*y[i];
    }
}

/* 
 * (COPYED FROM THE FLOYD WARSHALL EXERCISE!)
 * Timing function based on the TimeStep Counter of the CPU.
 *
 */
double rdtsc() 
{
    int i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = NUM_RUNS;

    /* 
     *
     * The CPUID instruction serializes the pipeline.
     * Using it, we can create execution barriers around the code we want to time.
     * The calibrate section is used to make the computation large enough so as to 
     * avoid measurements bias due to the timing overhead.
     */
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i) {
            compute();
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
        compute();
    }

    cycles = stop_tsc(start)/num_runs;
    return cycles;
}

int compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

int main(int argc, char **argv)
{
    time_t t;
    srand((unsigned) time(&t));

    printf("n \t cycles \t flops/cycle\n");
    int i,j;
    for(i=2;i<5e6;i+=i)
    {
        double flops = 4*i; // with n=i compute() takes 4*n flops
        double *cycles = malloc(NUM_MEASUREMENTS_PER_N*sizeof(double));
        for(j=0;j<NUM_MEASUREMENTS_PER_N;j++)
        {
            arr_init(i);
            cycles[j] = rdtsc();

        }
        qsort(cycles, NUM_MEASUREMENTS_PER_N, sizeof(double), compare_doubles);
        double median = cycles[((int)NUM_MEASUREMENTS_PER_N/2)+1];
        printf("%i \t %lf \t %lf\n", i, median, flops/median);
        free(cycles);
    }
}

