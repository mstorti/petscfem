#include <time.h>
#include <sys/time.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Manage chronometers.

*/ 
/// 
class Chrono {
private:
  /// structures for `libc' time function calls.
  clock_t start_time;
public:
#if 0
  /// create a chronometer
  Chrono();
  /// destroy chronometer
  ~Chrono() {};
#endif
  /// return elapsed CPU time from start
  double elapsed() const {return ((double) (clock() - start_time)) / CLOCKS_PER_SEC;};
  /// reset start time to actual time
  void start() {start_time = clock();};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** High precision chronometer
 */ 
/// 
class HPChrono {
  double gettod() const {
    struct timeval tv;
    gettimeofday(&tv,0);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
  }
private:
  /// structures for `libc' time function calls.
  double start_time;
public:
  /// return elapsed CPU time from start
  double elapsed() const {return gettod()-start_time;};
  /// reset start time to actual time
  void start() {start_time = gettod();};
};
