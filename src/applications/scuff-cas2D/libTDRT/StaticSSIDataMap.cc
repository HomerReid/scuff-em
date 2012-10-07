/*
 * StaticSSIDataTable.cc -- libTDRT routines for handling the data structures
 *                       -- that store collections of segment-segment integrals
 *
 * homer reid -- 6/2009 -- 10/2010
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "libhrutil.h"
#include "libTDRT.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

/*******************************************************************/
/* function that returns 1 if the two edges are nearby each other **/
/*******************************************************************/
int Nearby(double *Xs, double *Xe, double *Xsp, double *Xep)
{ 
   double MinD2, rRel, l, lp;

   l=VecDistance(Xs,Xe);
   lp=VecDistance(Xsp,Xep);
   MinD2=VecD2(Xs,Xsp);
   MinD2=fmin(MinD2,VecD2(Xs,Xep));
   MinD2=fmin(MinD2,VecD2(Xe,Xsp));
   MinD2=fmin(MinD2,VecD2(Xe,Xep));
   rRel=sqrt(MinD2)/fmax(l,lp);
   
   if ( rRel<DESINGULARIZATION_RADIUS )
    return 1;
   else
    return 0;

}

/*******************************************************************/
/* Given a pair of line segments, compute a unique key into the    */
/* StaticSSIDataTable for the given pair of segments.              */
/* 20101017: NOTE: because i am in a hurry to get this working, i  */
/* will temporarily use an inefficient scheme in which i fail to   */
/* take advantage of certain obvious symmetries; as a result the   */
/* resulting data tables are 4 or 8 times larger than they need to */
/* be and take that much more time to compute; FIXME FIXME         */
/*******************************************************************/
unsigned long GetStaticSSIDataTableKey(TDRTObject *Oa, int iXs, int iXe, 
                                       TDRTObject *Ob, int iXsp, int iXep)
{ 
  int NMax;

  NMax=Oa->NumVertices;
  if (Ob->NumVertices > NMax) 
   NMax=Ob->NumVertices;

  return iXep + NMax*(iXsp + NMax*(iXe + NMax*iXs));

} 

/*******************************************************************/
/* Attempt to retrieve static SSI data for the line segments       */
/* with the given endpoints. (The endpoints are given as indices   */
/* within the lists of vertices for the given objects).            */
/* If a corresponding data record is found in the data table, its  */
/* contents are copied into the buffer to which OutputBuffer       */ 
/* points, and OutputBuffer is returned.                           */ 
/* If no corresponding data record is found, NULL is returned.     */
/*******************************************************************/
StaticSSIDataRecord *GetStaticSSIData(StaticSSIDataTable *SSSIDT,
                                      TDRTObject *Oa, int iXs, int iXe, 
                                      TDRTObject *Ob, int iXsp, int iXep, 
                                      StaticSSIDataRecord *OutputBuffer)
{ 
  unsigned long Key;
  StaticSSIDataRecord *SSSIDR;
  StaticSSIDataMap::const_iterator it;
  StaticSSIDataMap *Map;

  if(SSSIDT==0) return 0;
   
  Key=GetStaticSSIDataTableKey(Oa, iXs, iXe, Ob, iXsp, iXep);
  Map=SSSIDT->Map;

  it=Map->find(Key); 

  if ( it==(Map->end()) ) /* no matching data record was found */
   return 0;

  /* a matching data record was found; copy the data into the output */
  /* buffer and then massage them as necessary before returning      */
  SSSIDR=it->second;
  memcpy(OutputBuffer,SSSIDR,sizeof(*SSSIDR));
 
  return SSSIDR;

}

/***************************************************************/
/* This is an alternative version of the previous routine that */
/* returns a pointer into the actual table data corresponding  */
/* to a key in the table (if one is found, or NULL otherwise). */
/* This routine is intended to be used only by the thread      */
/* routine below when the table is being initialized.          */
/***************************************************************/
StaticSSIDataRecord *FindStaticSSIData(StaticSSIDataTable *SSSIDT,
                                       TDRTObject *Oa, int iXs, int iXe, 
                                       TDRTObject *Ob, int iXsp, int iXep)
{ 
  unsigned long Key;
  StaticSSIDataMap::const_iterator it;
  StaticSSIDataMap *Map;

  if(SSSIDT==0) return 0;
   
  Key=GetStaticSSIDataTableKey(Oa, iXs, iXe, Ob, iXsp, iXep);
  Map=SSSIDT->Map;

  it=Map->find(Key);

  if ( it==(Map->end()) ) /* no matching data record was found */
   return 0;
  return it->second;

}

/***************************************************************/
/* data structure passed to thread routine below               */
/***************************************************************/
typedef struct ThreadData
 { 
   TDRTObject *Oa, *Ob;
   StaticSSIDataTable *SSSIDT;

   int nt, NumTasks;

 } ThreadData;

void *CreateStaticSSIDataTable_Thread(void *pTD)
{ 
  ThreadData *TD=(ThreadData *)pTD;

  /* local copies of fields in ThreadData structure */
  TDRTObject *Oa              = TD->Oa;
  TDRTObject *Ob              = TD->Ob;
  StaticSSIDataTable *SSSIDT  = TD->SSSIDT;

  /* other local variables */
  int niva, nivb;
  int NeedDerivatives;
  int iCVa, iEV1a, iEV2a;
  double *CVa, *EV1a, *EV2a;
  int iCVb, iEV1b, iEV2b;
  double *CVb, *EV1b, *EV2b;
  StaticSSIDataRecord *SSSIDR;
  int nt;

  NeedDerivatives = (Oa==Ob ? 0 : 1);
  
  nt=0;
  for(niva=0; niva<Oa->NumIVs; niva++)
   for(nivb=(Oa==Ob ? niva : 0); nivb<Ob->NumIVs; nivb++)
    {
      nt++;
      if (nt==TD->NumTasks) nt=0;
      if (nt!=TD->nt) continue;

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      iCVa=Oa->IVs[niva]            ;  // index of center vertex of BF A 
      iEV1a=Oa->Neighbors[2*niva]   ;  // index of end-vertex 1 of BF A
      iEV2a=Oa->Neighbors[2*niva+1] ;  // index of end-vertex 2 of BF A 

      iCVb=Ob->IVs[nivb]            ;  // index of center vertex of BF B 
      iEV1b=Ob->Neighbors[2*nivb]   ;  // index of end-vertex 1 of BF B
      iEV2b=Ob->Neighbors[2*nivb+1] ;  // index of end-vertex 2 of BF B

      CVa=Oa->Vertices   + 2*iCVa;     // center vertex of BF A
      EV1a=Oa->Vertices  + 2*iEV1a;    // end-vertex 1 of BF A 
      EV2a=Oa->Vertices  + 2*iEV2a;    // end-vertex 2 of BF A

      CVb=Ob->Vertices   + 2*iCVb;     // center vertex of BF B
      EV1b=Ob->Vertices  + 2*iEV1b;    // end-vertex 1 of BF B 
      EV2b=Ob->Vertices  + 2*iEV2b;    // end-vertex 2 of BF B
 
      if ( (SSSIDR=FindStaticSSIData(SSSIDT,Oa,iEV1a,iCVa,Ob,iEV1b,iCVb)) ) 
       ComputeStaticSSIData(EV1a,CVa,EV1b,CVb,NeedDerivatives, SSSIDR);
      if ( (SSSIDR=FindStaticSSIData(SSSIDT,Oa,iEV1a,iCVa,Ob,iEV2b,iCVb)) )
       ComputeStaticSSIData(EV1a,CVa,EV2b,CVb,NeedDerivatives, SSSIDR);
      if ( (SSSIDR=FindStaticSSIData(SSSIDT,Oa,iEV2a,iCVa,Ob,iEV1b,iCVb)) )
       ComputeStaticSSIData(EV2a,CVa,EV1b,CVb,NeedDerivatives, SSSIDR);
      if ( (SSSIDR=FindStaticSSIData(SSSIDT,Oa,iEV2a,iCVa,Ob,iEV2b,iCVb)) )
       ComputeStaticSSIData(EV2a,CVa,EV2b,CVb,NeedDerivatives, SSSIDR);

if (TDRTGeometry::LogLevel>=2)
 { if ( nivb==0 && niva ==   (Oa->NumIVs)/5 ) Log(" 20%%");
   if ( nivb==0 && niva == 2*(Oa->NumIVs)/5 ) Log(" 40%%");
   if ( nivb==0 && niva == 3*(Oa->NumIVs)/5 ) Log(" 60%%");
   if ( nivb==0 && niva == 4*(Oa->NumIVs)/5 ) Log(" 80%%");
 };

    };
  return 0;
}

/*******************************************************************/
/* Create a new SSIDataTable containing static segment-segment     */
/* integrals for all pairs of nearby segments on the given object  */
/* or pair of objects.                                             */
/*******************************************************************/
StaticSSIDataTable *CreateStaticSSIDataTable(TDRTObject *O, int NumThreads)
 { return CreateStaticSSIDataTable(O, O, NumThreads); }

StaticSSIDataTable *CreateStaticSSIDataTable(TDRTObject *Oa, TDRTObject *Ob, int NumThreads)
{  
  int niva, nivb;
  unsigned long Key;
  int NNearby, nNearby;
  int NeedDerivatives;

  int iCVa, iEV1a, iEV2a;    /* indices of center and end vertices for BF A */
  double *CVa, *EV1a, *EV2a; /* center and end vertices for BF A */

  int iCVb, iEV1b, iEV2b;    /* indices of center and end vertices for BF B */ 
  double *CVb, *EV1b, *EV2b; /* center and end vertices for BF B */

  StaticSSIDataTable *SSSIDT;
  StaticSSIDataRecord *SSSIDR;

  NeedDerivatives = (Oa==Ob ? 0 : 1);

  /*******************************************************************/
  /* first pass to count how many nearby segments there are.         */
  /*******************************************************************/
  NNearby=0;
  for(niva=0; niva<Oa->NumIVs; niva++)
   for(nivb=0; nivb<Ob->NumIVs; nivb++)
    { 
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      iCVa=Oa->IVs[niva]            ;  // index of center vertex of BF A 
      iEV1a=Oa->Neighbors[2*niva]   ;  // index of end-vertex 1 of BF A
      iEV2a=Oa->Neighbors[2*niva+1] ;  // index of end-vertex 2 of BF A 

      iCVb=Ob->IVs[nivb]            ;  // index of center vertex of BF B 
      iEV1b=Ob->Neighbors[2*nivb]   ;  // index of end-vertex 1 of BF B
      iEV2b=Ob->Neighbors[2*nivb+1] ;  // index of end-vertex 2 of BF B

      CVa=Oa->Vertices   + 2*iCVa;     // center vertex of BF A
      EV1a=Oa->Vertices  + 2*iEV1a;    // end-vertex 1 of BF A 
      EV2a=Oa->Vertices  + 2*iEV2a;    // end-vertex 2 of BF A

      CVb=Ob->Vertices   + 2*iCVb;     // center vertex of BF B
      EV1b=Ob->Vertices  + 2*iEV1b;    // end-vertex 1 of BF B 
      EV2b=Ob->Vertices  + 2*iEV2b;    // end-vertex 2 of BF B

      if ( Nearby(EV1a, CVa, EV1b, CVb) )
       NNearby++;
      if ( Nearby(EV1a, CVa, EV2b, CVb) )
       NNearby++;
      if ( Nearby(EV2a, CVa, EV1b, CVb) )
       NNearby++;
      if ( Nearby(EV2a, CVa, EV2b, CVb) )
       NNearby++;

    };

  if (TDRTGeometry::LogLevel>=2)
   Log(" %i/%i nearby pairs (%.2f %%)",
           NNearby,4*Oa->NumIVs*Ob->NumIVs, 
           ((double)NNearby)/((double)(4*Oa->NumIVs*Ob->NumIVs)));

  /*******************************************************************/
  /* allocate a new StaticSSIDataTable with enough room to store     */
  /* NNearby data records.                                           */
  /*******************************************************************/
  SSSIDT=(StaticSSIDataTable *)malloc(sizeof(*SSSIDT));
  SSSIDT->Map=new StaticSSIDataMap;
  // the following two commands are needed if we are using the 
  // google dense_hash_map implementation of map, but not for  
  // the std::tr1 implementation
  // SSSIDT->Map->set_empty_key(-1); 
  // SSSIDT->Map->set_deleted_key(-2);
  SSSIDT->Buffer=(StaticSSIDataRecord *)malloc(NNearby*sizeof(StaticSSIDataRecord));
  if (SSSIDT->Buffer==0)
   { if (TDRTGeometry::LogLevel>=2)
      Log("...insufficient memory for static SSI data table");
     delete SSSIDT->Map;
     free(SSSIDT);
     return 0;
   };

  /*******************************************************************/
  /* second pass to assign (key,data record) pairs.  *****************/
  /*******************************************************************/
  nNearby=0;
  for(niva=0; niva<Oa->NumIVs; niva++)
   for(nivb=0; nivb<Ob->NumIVs; nivb++)
    { 
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      iCVa=Oa->IVs[niva]            ;  // index of center vertex of BF A 
      iEV1a=Oa->Neighbors[2*niva]   ;  // index of end-vertex 1 of BF A
      iEV2a=Oa->Neighbors[2*niva+1] ;  // index of end-vertex 2 of BF A 

      iCVb=Ob->IVs[nivb]            ;  // index of center vertex of BF B 
      iEV1b=Ob->Neighbors[2*nivb]   ;  // index of end-vertex 1 of BF B
      iEV2b=Ob->Neighbors[2*nivb+1] ;  // index of end-vertex 2 of BF B

      CVa=Oa->Vertices   + 2*iCVa;     // center vertex of BF A
      EV1a=Oa->Vertices  + 2*iEV1a;    // end-vertex 1 of BF A 
      EV2a=Oa->Vertices  + 2*iEV2a;    // end-vertex 2 of BF A

      CVb=Ob->Vertices   + 2*iCVb;     // center vertex of BF B
      EV1b=Ob->Vertices  + 2*iEV1b;    // end-vertex 1 of BF B 
      EV2b=Ob->Vertices  + 2*iEV2b;    // end-vertex 2 of BF B

      if ( Nearby(EV1a, CVa, EV1b, CVb) )
       { Key=GetStaticSSIDataTableKey(Oa, iEV1a, iCVa, Ob, iEV1b, iCVb);
         SSSIDR=SSSIDT->Buffer + nNearby++;
         SSSIDT->Map->insert( std::pair<unsigned long, StaticSSIDataRecord *>(Key,SSSIDR) );
       };
         
      if ( Nearby(EV1a, CVa, EV2b, CVb) )
       { Key=GetStaticSSIDataTableKey(Oa, iEV1a, iCVa, Ob, iEV2b, iCVb);
         SSSIDR=SSSIDT->Buffer + nNearby++;
         SSSIDT->Map->insert( std::pair<unsigned long, StaticSSIDataRecord *>(Key,SSSIDR) );
       };
         
      if ( Nearby(EV2a, CVa, EV1b, CVb) )
       { Key=GetStaticSSIDataTableKey(Oa, iEV2a, iCVa, Ob, iEV1b, iCVb);
         SSSIDR=SSSIDT->Buffer + nNearby++;
         SSSIDT->Map->insert( std::pair<unsigned long, StaticSSIDataRecord *>(Key,SSSIDR) );
       };
         
      if ( Nearby(EV2a, CVa, EV2b, CVb) )
       { Key=GetStaticSSIDataTableKey(Oa, iEV2a, iCVa, Ob, iEV2b, iCVb);
         SSSIDR=SSSIDT->Buffer + nNearby++;
         SSSIDT->Map->insert( std::pair<unsigned long, StaticSSIDataRecord *>(Key,SSSIDR) );
       };

    };

  /*******************************************************************/
  /* finally, a multithreaded third pass to compute the data records.*/
  /*******************************************************************/
#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[NumThreads], *TD;
  pthread_t *Threads = new pthread_t[NumThreads];
  for(int nt=0; nt<NumThreads; nt++)
   { 
     TD=&(TDS[nt]);
     TD->nt=nt;
     TD->NumTasks=NumThreads;

     TD->Oa=Oa;
     TD->Ob=Ob;
     TD->SSSIDT=SSSIDT;

     if (nt+1 == NumThreads)
      CreateStaticSSIDataTable_Thread((void *)TD);
     else
      pthread_create( &(Threads[nt]), 0, CreateStaticSSIDataTable_Thread, (void *)TD);
   };

  for(int nt=0; nt<NumThreads; nt++)
   pthread_join(Threads[nt],0);

  delete[] Threads;
  delete[] TDs;

#else 
#ifndef USE_OPENMP
  int NumTasks=NumThreads=1;
#else
  int NumTasks=NumThreads*100;
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nt=0; nt<NumTasks; nt++)
   { 
     ThreadData TD1;
     TD1.nt=nt;
     TD1.NumTasks=NumTasks;

     TD1.Oa=Oa;
     TD1.Ob=Ob;
     TD1.SSSIDT=SSSIDT;

     CreateStaticSSIDataTable_Thread((void *)&TD1);
   };
#endif

  if (TDRTGeometry::LogLevel>=2) 
   Log(" done!");

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  return SSSIDT;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void DestroyStaticSSIDataTable(StaticSSIDataTable *SSSIDT)
{
  if (!SSSIDT || !(SSSIDT->Buffer) )
   return;
  free(SSSIDT->Buffer); 
  free(SSSIDT);

}
