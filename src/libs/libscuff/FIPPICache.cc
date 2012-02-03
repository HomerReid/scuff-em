/*
 * FIPPICache.cc -- implementation of the FIPPICache class for libscuff
 * 
 * homer reid    -- 11/2005 -- 1/2011
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tr1/unordered_map>

#include <libhrutil.h>

#include "libscuff.h"
#include "libscuffInternals.h"

#define KEYLEN 15

int Found;
int NotFound;

/*--------------------------------------------------------------*/
/*- note: i found this on wikipedia ... ------------------------*/
/*--------------------------------------------------------------*/
long JenkinsHash(char *key, size_t len)
{
    long hash; 
    int i;
    for(hash = i = 0; i < len; ++i)
    {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

long HashFunction(const double *Key)
{ return JenkinsHash( (char *)Key, KEYLEN*sizeof(double) );
} 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
typedef struct
 { double Key[KEYLEN];
 } KeyStruct;

typedef std::pair<KeyStruct, QIFIPPIData *> KeyValuePair;

struct KeyHash
 {
   long operator() (const KeyStruct &K) const { return HashFunction(K.Key); }
 };

typedef struct 
 { 
   bool operator()(const KeyStruct &K1, const KeyStruct &K2) const
    { int nd;
      for(nd=0; nd<KEYLEN; nd++)
       if ( fabs(K1.Key[nd] - K2.Key[nd]) > 1.0e-8 * fabs(K1.Key[nd]) )
        return false;
      return true;
    };

 } KeyCmp;

typedef std::tr1::unordered_map< KeyStruct,
                                 QIFIPPIData *, 
                                 KeyHash, 
                                 KeyCmp> KeyValueMap;

/*--------------------------------------------------------------*/
/*- class constructor ------------------------------------------*/
/*--------------------------------------------------------------*/
FIPPICache::FIPPICache()
{
  DoNotCompute=0;

  pthread_mutex_init(&FCMutex,0);

  KeyValueMap *KVM=new KeyValueMap;
  opTable = (void *)KVM;
}

/*--------------------------------------------------------------*/
/*- class destructor  ------------------------------------------*/
/*--------------------------------------------------------------*/
FIPPICache::~FIPPICache()
{
  KeyValueMap *KVM=(KeyValueMap *)opTable;
  delete KVM;
} 

/*--------------------------------------------------------------*/
/*- routine for fetching a FIPPI data record from a FIPPIDT:    */
/*- we look through our table to see if we have a record for    */
/*- this panel pair, we return it if we do, and otherwise we    */
/*- compute a new FIPPI data record for this panel pair and     */
/*- add it to the table.                                        */
/*- important note: the vertices are assumed to be canonically  */
/*- ordered on entry.                                           */
/*--------------------------------------------------------------*/
QIFIPPIData *FIPPICache::GetQIFIPPIData(double **OVa, double **OVb, int ncv)
{
  /***************************************************************/
  /* the search key is kinda stupid, just a string of 15 doubles */
  /* as follows:                                                 */
  /* 0--2    VMed  - VMin [0..2]                                 */
  /* 3--5    VMax  - VMin [0..2]                                 */
  /* 6--8    VMinP - VMin [0..2]                                 */
  /* 9--11   VMedP - VMin [0..2]                                 */
  /* 12--14  VMaxP - VMin [0..2]                                 */
  /***************************************************************/
  KeyStruct K;
  VecSub(OVa[1], OVa[0], K.Key+0 );
  VecSub(OVa[2], OVa[0], K.Key+3 );
  VecSub(OVb[0], OVa[0], K.Key+6 );
  VecSub(OVb[1], OVa[0], K.Key+9 );
  VecSub(OVb[2], OVa[0], K.Key+12);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  KeyValueMap *KVM=(KeyValueMap *)opTable;
  KeyValueMap::iterator p=KVM->find(K);
  if ( p != (KVM->end()) )
   { Found++;
     return (QIFIPPIData *)(p->second);
   };
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  NotFound++;
  KeyStruct *K2 = (KeyStruct *)malloc(sizeof(*K2));
  memcpy(K2->Key, K.Key, KEYLEN*sizeof(double));
  QIFIPPIData *QIFD=(QIFIPPIData *)malloc(sizeof *QIFD);
  if (DoNotCompute==0)
   ComputeQIFIPPIData(OVa, OVb, ncv, QIFD);
   
  pthread_mutex_lock(&FCMutex);
  KVM->insert( KeyValuePair(*K2, QIFD) );
  pthread_mutex_unlock(&FCMutex);

  return QIFD;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
FIPPICache GlobalFIPPICache;
