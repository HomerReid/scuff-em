/*
 * FIPPITable.cc -- implementation of the FIPPITable class for libscuff
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

long HashFunction(const double *Key)
{ return 0; 
} 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
typedef std::pair<const double *, QIFIPPIData *> KeyValuePair;

struct KeyHash
 {
   long operator() (const double *Key) const { return HashFunction(Key); }
 };

typedef struct 
 { 
   bool operator()(const double *D1, const double *D2) const
    { int nd;
      for(nd=0; nd<KEYLEN; nd++)
       if ( fabs(D1[nd]-D2[nd]) > 1.0e-8 * fabs(D1[nd]) )
        return false;
      return true;
    };

 } KeyCmp;

typedef std::tr1::unordered_map< const double *, 
                                 QIFIPPIData *, 
                                 KeyHash, 
                                 KeyCmp> KeyValueMap;

KeyValueMap *MyKeyValueMap=0;

/*--------------------------------------------------------------*/
/*- class constructor ------------------------------------------*/
/*--------------------------------------------------------------*/
FIPPITable::FIPPITable()
{
  MyKeyValueMap=new KeyValueMap;
  DoNotCompute=0;
}

/*--------------------------------------------------------------*/
/*- class destructor  ------------------------------------------*/
/*--------------------------------------------------------------*/
FIPPITable::~FIPPITable()
{

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
QIFIPPIData *FIPPITable::GetQIFIPPIData(double **OVa, double **OVb, int ncv)
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
  double Key[15];
  VecSub(OVa[1], OVa[0], Key+0 );
  VecSub(OVa[2], OVa[0], Key+3 );
  VecSub(OVb[0], OVa[0], Key+6 );
  VecSub(OVb[1], OVa[0], Key+9 );
  VecSub(OVb[2], OVa[0], Key+12);

#if 0
KeyCmp MyKeyCmp;
if (MyKeyCmp(Key,Key))
 printf("yes\n");
double Key2[15];
memcpy(Key2,Key,15*sizeof(double));
Key[0]+=1.0;
if (MyKeyCmp(Key,Key2))
 printf("yes\n");
else
 printf("no\n");
#endif

  KeyValueMap::iterator p=MyKeyValueMap->find(Key);
  if ( p == (MyKeyValueMap->end()) )
   { QIFIPPIData *QIFD=(QIFIPPIData *)malloc(sizeof *QIFD);
     if (DoNotCompute==0)
      ComputeQIFIPPIData(OVa, OVb, ncv, QIFD);
     MyKeyValueMap->insert( KeyValuePair(Key, QIFD) );
     return QIFD;
   }
  else
{
printf(" found it!\n");
   return p->second;
}

}
