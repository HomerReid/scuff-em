/* Readers/writers lock class for libscuff.  Steven G. Johnson (2012) */

/* For both pthreads and OpenMP, we just use pthread_rwlock_t.  

   Perhaps some future version of OpenMP will add rwlocks, and
   we can eliminat the pthread dependency that way. */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#if defined(USE_OPENMP) || defined(USE_PTHREAD)
#  include <pthread.h>
#endif

namespace scuff {

class rwlock {
public:
  rwlock();
  ~rwlock();

  void read_lock();
  void read_unlock();
  void write_lock();
  void write_unlock();

private:
  // make these private, so that trying to copy a lock gives compile-time error
  rwlock(const rwlock &L) { (void) L; }
  void operator=(const rwlock &L) { (void) L; }

#if defined(USE_OPENMP) || defined(USE_PTHREAD)
  pthread_rwlock_t lock;
#endif
};

}
