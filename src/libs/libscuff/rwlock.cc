/* Readers/writers lock class for libscuff.  Steven G. Johnson (2012) */

/* For both pthreads and OpenMP, we just use pthread_rwlock_t.  

   Perhaps some future version of OpenMP will add rwlocks, and
   we can eliminat the pthread dependency that way. */

#include <stdlib.h>
#include "rwlock.h"

namespace scuff {

rwlock::rwlock() {
#if defined(USE_OPENMP) || defined(USE_PTHREAD)
  pthread_rwlock_init(&lock,0);
#endif
}

rwlock::~rwlock() {
#if defined(USE_OPENMP) || defined(USE_PTHREAD)
  pthread_rwlock_destroy(&lock);
#endif
}

void rwlock::read_lock() {
#if defined(USE_OPENMP) || defined(USE_PTHREAD)
  pthread_rwlock_rdlock(&lock);
#endif
}

void rwlock::read_unlock() {
#if defined(USE_OPENMP) || defined(USE_PTHREAD)
  pthread_rwlock_unlock(&lock);
#endif
}

void rwlock::write_lock() {
#if defined(USE_OPENMP) || defined(USE_PTHREAD)
  pthread_rwlock_wrlock(&lock);
#endif
}

void rwlock::write_unlock() {
#if defined(USE_OPENMP) || defined(USE_PTHREAD)
  pthread_rwlock_unlock(&lock);
#endif
}

} // namespace
