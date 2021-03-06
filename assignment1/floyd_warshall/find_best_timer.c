//#error Please comment out the next two lines under linux, then comment this error
//#include "stdafx.h" //Visual studio expects this line to be the first one, comment out if different compiler
//#include <windows.h> // Include if under windows

#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "tsc_x86.h"

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define FREQUENCY 2.6e9
#define min(a,b) a < b ? a : b

#define CALIBRATE

int n;
double *A;

/******
	Initialize the input
******/

void init_mat(){
  int i,j;
  for(i =0; i<n;i++) A[n*i+i]=0;
  for(i =0; i<n;i++) {
	for(j=0; j<n;j++){
		if(i!=j){
			A[n*i+j] =  lrand48()%n;
		} 
	}
  }
}


/* 
 * Straightforward implementation of Floyd-Warshall
 * 
 */

void compute() {
  int i,j,k;
  for(k = 0; k < n; k++) {
    for(i = 0; i < n; i++) {
      for(j = 0; j < n; j++) {
        A[n*i+j] = min(A[n*i+j],A[n*i+k] + A[n*k+j]);
      }
    }
  }
}

/* 
 * Timing function based on the TimeStep Counter of the CPU.
 */

double rdtsc() {
  int i, num_runs;
  myInt64 cycles;
  myInt64 start;
  num_runs = NUM_RUNS;

/* 
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

double c_clock() {
  int i, num_runs;
  double cycles;
  clock_t start, end;

  num_runs = NUM_RUNS;
#ifdef CALIBRATE
  while(num_runs < (1 << 14)) {
      start = clock();
      for (i = 0; i < num_runs; ++i) {
          compute();
      }
      end = clock();

      cycles = (double)(end-start);

      // Same as in c_clock: CYCLES_REQUIRED should be expressed accordingly to the order of magnitude of CLOCKS_PER_SEC
      if(cycles >= CYCLES_REQUIRED/(FREQUENCY/CLOCKS_PER_SEC)) break;

      num_runs *= 2;
  }
#endif

  start = clock();
  for(i=0; i<num_runs; ++i)
    { compute(); }
  end = clock();

  return (double)(end-start)/num_runs;
}

#ifndef WIN32
double timeofday() {
  int i, num_runs;
  double cycles;
  struct timeval start, end;

  num_runs = NUM_RUNS;
#ifdef CALIBRATE
  while(num_runs < (1 << 14)) {
      gettimeofday(&start, NULL);
      for (i = 0; i < num_runs; ++i) {
          compute();
      }
      gettimeofday(&end, NULL);

      cycles = (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6)*FREQUENCY;

      if(cycles >= CYCLES_REQUIRED) break;

      num_runs *= 2;
  }
#endif

  gettimeofday(&start, NULL);
  for(i=0; i < num_runs; ++i) {
    compute(); 
  }
  gettimeofday(&end, NULL);
  
  return (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6)/ num_runs;
}

#else

double gettickcount() {
  int i, num_runs;
  double cycles, start, stop;

  num_runs = NUM_RUNS;
#ifdef CALIBRATE
  while(num_runs < (1 << 14)) {
      start = (double)GetTickCount();
      for (i = 0; i < num_runs; ++i) {
          compute();
      }
      end = (double)GetTickCount();

      cycles = (end-start)*FREQUENCY/1e3; // end-start provides a measurement in the order of milliseconds

      if(cycles >= CYCLES_REQUIRED) break;

      num_runs *= 2;
  }
#endif

  start = (double)GetTickCount();
  for(i=0; i < num_runs; ++i) {
    compute(); 
  }
  end = (double)GetTickCount();

  return (start-end)/num_runs;
}

double queryperfcounter(LARGE_INTEGER f) {
  int i, num_runs;
  double cycles;
  LARGE_INTEGER start, end;

  num_runs = NUM_RUNS;
#ifdef CALIBRATE
  while(num_runs < (1 << 14)) {
      QueryPerformanceCounter(&start);
      for (i = 0; i < num_runs; ++i) {
          compute();
      }
      QueryPerformanceCounter(&end);

      cycles = (double)(end.QuadPart - start.QuadPart);
      
      // Same as in c_clock: CYCLES_REQUIRED should be expressed accordingly to the order of magnitude of f
      if(cycles >= CYCLES_REQUIRED/(FREQUENCY/f)) break; 

      num_runs *= 2;
  }
#endif

  QueryPerformanceCounter(&start);
  for(i=0; i < num_runs; ++i) {
    compute(); 
  }
  QueryPerformanceCounter(&end);

  return (double)(end.QuadPart - start.QuadPart)/num_runs;
}

#endif

int
compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

//
// run both timers "run" times and then find the median.
// for clock and rdtsc the sum of differences to their median is calculated
//
int main(int argc, char **argv){

  int i,j;
  double r;
  double c;
  double c_seconds;
  double *rdtsc_arr;
  double *clock_arr;
  int runs = 20;

  rdtsc_arr = (double*)malloc(runs*sizeof(double));
  clock_arr = (double*)malloc(runs*sizeof(double));

  n=700;

  for(i=0;i<runs;i++)
  {
    A = (double*)malloc(n*n*sizeof(double));
    init_mat();
    rdtsc_arr[i] = rdtsc()/FREQUENCY;
    init_mat();
    clock_arr[i] = c_clock()/CLOCKS_PER_SEC;
    free(A);
  }
  qsort (rdtsc_arr, runs, sizeof (double), compare_doubles);
  qsort (clock_arr, runs, sizeof (double), compare_doubles);

  float clock_median = clock_arr[(int)runs/2];
  float rdtsc_median = rdtsc_arr[(int)runs/2];

  float clock_median_difference = 0;
  float rdtsc_median_difference = 0;

  for(i=0;i<runs;i++)
  {
  //    printf("%lf, %lf\n", rdtsc_arr[i], clock_arr[i]);
      clock_median_difference += fabs(clock_arr[i] - clock_median);
      rdtsc_median_difference += fabs(rdtsc_arr[i] - rdtsc_median);
  }

  printf("Clock median difference: %lf, rdtsc median difference: %lf\n", clock_median_difference, rdtsc_median_difference);
  return 0;
}

