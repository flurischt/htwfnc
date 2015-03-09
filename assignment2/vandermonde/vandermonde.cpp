#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <string>

#define FLOAT double
#define CYCLES_REQUIRED 1e7
#define REP 15

#include "tsc_x86.h"

using namespace std;

void dumpList(vector< double > const & l, string const & filename, unsigned x)
{
  ofstream outfile;
  outfile.open(filename.c_str(), ios_base::app);

  vector< double >::const_iterator it = l.begin();
  outfile << x << " " << *it;
  for (++it; it != l.end(); ++it)
	  outfile << " " << *it;
  outfile << "\n";
  outfile.close();
}

double median(vector< double > const & l) {
  double res;
  size_t n = l.size();
  vector< double > ordl = l;
  
  if(n>1) {
    sort(ordl.begin(), ordl.end());
    res = n%2 == 0 ? ordl[n/2] : (ordl[n/2] + ordl[n/2+1])/2.;
  } else {
    res = l[0];
  }
  
  return res;
}

/*
 * Substitute with the opcount of method Vandermonde::det
 */
  
#define OPCOUNT ((n * n) - n)

template<typename F> 
class V {
  unsigned _n;
  F * _r1;

  void _randInit()
  {
    for (unsigned i = 0; i < _n; ++i) _r1[i] = static_cast<F>(rand()+1)/RAND_MAX;
  }

public:
  V(unsigned n): _n(n) {
    _r1 = new double[n];
    _randInit();
  }

  ~V() {
    delete [] _r1;
  }
  
  void print() const {
    for (unsigned i = 0; i < _n; i++) {
      cout << "| ";
      for (unsigned j = 0; j < _n; j++) {
        cout << setprecision(2) << pow(_r1[j], static_cast< double >(i)) << "\t";
      }
      cout << " |" << endl;
    }      
  }
  
  F det() const {
    F d = 1;
    for (unsigned i = 1; i < _n; i++)
      for (unsigned j = 0; j < i; j++)
        d *= (_r1[i] - _r1[j]);
    
    return d;
  }

  F det_opt() const {
    
    // Your optimized code here
    
  }

  bool validate(double threshold) const {
    
    F d_ref = det(), d = det_opt();
    
    double err = fabs(d - d_ref);
    bool succ = true;
    
	  if(err > threshold) {
	    succ = false;
		  cerr << "Error: det_ref = " << d_ref << "\t-- det = " << d << "\t-- Err = " << err << endl;
    }
    
    return succ;
  }
  
};

void test(size_t n)
{

  V<FLOAT> v(n);
  FLOAT * va;
  FLOAT r;
  myInt64 start, end;
  ofstream ns("/dev/null");

  double cycles = 0.;
  unsigned num_runs = 2, multiplier = 1;

  if( !v.validate(1e-5) )
    exit(EXIT_FAILURE);
  else
    cout << "Test passed." << endl;
    
  init_tsc();

  do{
    num_runs = num_runs * multiplier;
    start = start_tsc();
    for(size_t i = 0; i < num_runs; i++) {
      r += v.det_opt();
    }
    end = stop_tsc(start);
        
    cycles = (double) end;
    multiplier = ceil (  (CYCLES_REQUIRED) / (cycles)  + 1.0 );

  }while (multiplier > 2);

  vector< double > cyclesList, perfList;

  for (size_t j = 0; j < REP; j++) {

    start = start_tsc();
    for (size_t i = 0; i < num_runs; ++i) {
      r += v.det_opt();
    }
    end = stop_tsc(start);
     
    cycles = ((double) end) / num_runs;

    cyclesList.push_back(cycles);
    perfList.push_back(OPCOUNT/cycles);

  }
  
  if (r>0) ns << r;
  
  cout << "Computation with N =" << n <<  " - Performance [flops/cycle]: " << median(perfList) << endl;
  cout << "Op count [flops]: " << OPCOUNT << endl;
  cout << "Runtime [cycle]: " << median(cyclesList) << endl;
  
}

int main(int argc, char **argv)
{
  srand(time(NULL));
  
  for(int i = 4; i <= 1024; i*=2)
    test(5+i);
  
  return 0;
}
