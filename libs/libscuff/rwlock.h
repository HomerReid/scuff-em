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
