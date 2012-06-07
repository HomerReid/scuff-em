/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

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
