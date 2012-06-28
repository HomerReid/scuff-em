/***************************************************************/
/*                                                             */
/***************************************************************/
void GetGCEdgeHalfEdge(OpenSurface *OSA, int ntea,
                       OpenSurface *OSB, int nteb,
                       k, cdouble GC[2])
{
  /*--------------------------------------------------------------*/
  /*- panel-panel interactions -----------------------------------*/
  /*--------------------------------------------------------------*/

  /*--------------------------------------------------------------*/
  /*- edge-panel interactions  -----------------------------------*/
  /*--------------------------------------------------------------*/
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ThreadData
 { 
   AOCMBArgStruct *Args;
   int nt, nThread;
 } ThreadData;

typedef struct ACCMBArgStruct 
{
  RWGGeometry *G;
  RWGComposite *CA, *CB;

} AOCMBArgStruct;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void *AOCMBThread(void *data)
{
  /***************************************************************/
  /* local copies of fields in argument structure                */
  /***************************************************************/
  ThreadData *TD = (ThreadData *)data;
  AOCMBArgStruct *Args  = TD->Args;
  RWGGeometry *G        = Args->G;
  RWGComposite *CA      = Args->CA;
  RWGComposite *CB      = Args->CB;
  int RowOffset         = Args->RowOffset;
  int ColOffset         = Args->ColOffset;
  HMatrix *M            = Args->M;

  /***************************************************************/
  /* other local variables ***************************************/
  /***************************************************************/
  int nosa, NOSA = CA->NumOpenSurfaces;
  int nosb, NOSB = CB->NumOpenSurfaces;
  int ntea, NTEA;
  int nteb, NTEB; 

  int RegionA1, RegionA2, RegionB1, RegionB2;
  int CommonRegion[2], NumCommonRegions;

  cdouble Eps1, Mu1, k1, EEPreFac1, MEPreFac1, MMPreFac1;
  cdouble Eps2, Mu2, k2, EEPreFac2, MEPreFac2, MMPreFac2;

  cdouble GC[2];

  for(nosa=0; nosa<NOSA; nosa++)
   for(nosb=0; nosb<NOSB; nosb++)
    { 
      /*--------------------------------------------------------------*/
      /* figure out if open surfaces #nosa and #nosb share 0, 1, or 2 */
      /* regions in common. if the answer is 0 then basis functions   */
      /* on the two surfaces do not interact with one another.        */
      /*--------------------------------------------------------------*/
      RegionA1 = Regions[2*nosa+0];
      RegionA2 = Regions[2*nosa+1];
      RegionB1 = Regions[2*nosb+0];
      RegionB2 = Regions[2*nosb+1];
      NumCommonRegions=0;
      if ( RegionA1==RegionB1 || RegionA1==RegionB2 )
       CommonRegion[NumCommonRegions++]=RegionA1;
      if ( RegionA2==RegionB1 || RegionA2==RegionB2 )
       CommonRegion[NumCommonRegions++]=RegionA2;
      if (NumCommonRegions==0) 
       continue;

      /*--------------------------------------------------------------*/
      /* if the two RWG composites in question are distinct ... ------*/
      /*--------------------------------------------------------------*/
      if ( (CA!=CB) )
       ErrExit("%s:%i: feature not yet implemented",__FILE__,__LINE__);

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      Eps1=EpsTF[ CommonRegion[0] ];
      Mu1=EpsTF[ CommonRegion[0] ];
      k1=csqrt2(Eps1*Mu1)*Omega;
      EEPreFac1 = Sign*II*Mu1*Omega;
      EMPreFac1 = -Sign*II*k1;
      MMPreFac1 = -1.0*Sign*II*Eps1*Omega;

      if (NumCommonRegions==2)
       { Eps2=EpsTF[ CommonRegion[1] ]; 
         Mu2=EpsTF[ CommonRegion[1] ]; }
         k2=csqrt2(Eps2*Mu2)*Omega;
         EEPreFac2 = Sign*II*Mu2*Omega;
         EMPreFac2 = -Sign*II*k2;
         MMPreFac2 = -1.0*Sign*II*Eps2*Omega;
       };

      /***************************************************************/
      /* preinitialize an argument structure to be passed to         */
      /* GetEdgeEdgeInteractions() below                             */
      /***************************************************************/
      GetEEIArgStruct MyGetEEIArgs, *GetEEIArgs=&MyGetEEIArgs;
      InitGetEEIArgs(GetEEIArgs);

      GetEEIArgs->PanelsA=OSA->Panels;
      GetEEIArgs->VerticesA=CA->Vertices;
      GetEEIArgs->PanelsB=OSB->Panels;
      GetEEIArgs->VerticesB=CB->Vertices;
      GetEEIArgs->NumGradientComponents = GradB ? 3 : 0;
      GetEEIArgs->NumTorqueAxes=NumTorqueAxes;
      GetEEIArgs->GammaMatrix=GammaMatrix;
    
      /* pointers to arrays inside the structure */
      cdouble *GC=GetEEIArgs->GC;
      cdouble *GradGC=GetEEIArgs->GradGC;
      cdouble *dGCdT=GetEEIArgs->dGCdT;
          
      /*--------------------------------------------------------------*/
      /* now loop over all basis functions (both full and half RWG    */
      /* functions) on open surfaces #nosa and #nosb.                 */
      /*--------------------------------------------------------------*/
      NTEA=OpenSurfaces[nosa]->NumTotalEdges;
      OffsetA = RowOffset + BFIndexOffset[nosa];
      NTEB=OpenSurfaces[nosb]->NumTotalEdges;
      OffsetB = ColOffset + BFIndexOffset[nosb];
      for(ntea=0; ntea<NTEA; ntea++)
       for(nteb=0; nteb<NTEB; nteb++)
        { 
          nt++;
          if (nt==TD->nThread) nt=0;
          if (nt!=TD->nt) continue;

          if (G->LogLevel>=SCUFF_VERBOSELOGGING)
           for(int PerCent=0; PerCent<9; PerCent++)
            if ( nteb==Symmetric*ntea &&  (ntea == (PerCent*NTEA)/10) )
             MutexLog("%i0 %% (%i/%i)...",PerCent,ntea,NTEA);
          
          GetGCEdgeHalfEdge(OSA, ntea, OSB, nteb, k1, GC);
          B->SetEntry(OffsetA + 2*ntea+0, OffsetB + 2*nteb+0, EEPreFac1*GC[0]);
          B->SetEntry(OffsetA + 2*ntea+0, OffsetB + 2*nteb+1, EMPreFac1*GC[1]);
          B->SetEntry(OffsetA + 2*ntea+1, OffsetB + 2*nteb+0, EMPreFac1*GC[1]);
          B->SetEntry(OffsetA + 2*ntea+1, OffsetB + 2*nteb+1, MMPreFac1*GC[0]);
          
          if (NumCommonRegions==2)
           { GetGCEdgeHalfEdge(OSA, ntea, OSB, nteb, k2, GC);
             B->AddEntry(OffsetA + 2*ntea+0, OffsetB + 2*nteb+0, EEPreFac2*GC[0]);
             B->AddEntry(OffsetA + 2*ntea+0, OffsetB + 2*nteb+1, EMPreFac2*GC[1]);
             B->AddEntry(OffsetA + 2*ntea+1, OffsetB + 2*nteb+0, EMPreFac2*GC[1]);
             B->AddEntry(OffsetA + 2*ntea+1, OffsetB + 2*nteb+1, MMPreFac2*GC[0]);

           };
  
        };

    }; // for (nosa=0 ...) for nosb=0 ...)
}

/***************************************************************/
/* this is a version of the AssembleBEMMatrixBlock() routine   */
/* for use when one (or both) of the two objects in question   */
/* is an RWGComposite instead of an RWGObject.                 */
/***************************************************************/
/***************************************************************/  
/***************************************************************/  
/***************************************************************/
void AssembleCCMatrixBlock(ACCMBArgStruct *Args)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGGeometry  *G=Args->G;
  RWGComposite *CA=Args->CA;
  RWGComposite *CB=Args->CB;
  cdouble Omega=Args->Omega;

  /***************************************************************/
  /* precompute material properties of all regions on both       */
  /* composites                                                  */
  /***************************************************************/
  int nr;
  for(nr=0; nr<CA->NumRegions; nr++)
   CA->RegionMPs[nr] -> GetEpsMu(Omega, CA->EpsTF + nr, CA->MuTF + nr);
  if (CB!=CA)
   { for(nr=0; nr<CB->NumRegions; nr++)
      CB->RegionMPs[nr] -> GetEpsMu(Omega, CB->EpsTF + nr, CB->MuTF + nr);
   };

  /***************************************************************/
  /* fire off threads ********************************************/
  /***************************************************************/
  GlobalFIPPICache.Hits=GlobalFIPPICache.Misses=0;

  int nt, nThread=Args->nThread;

#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[nThread], *TD;
  pthread_t *Threads = new pthread_t[nThread];
  for(nt=0; nt<nThread; nt++)
   { 
     TD=&(TDs[nt]);
     TD->nt=nt;
     TD->nThread=nThread;
     TD->Args=Args;
     if (nt+1 == nThread)
       ACCMBThread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, ACCMBThread, (void *)TD);
   }
  for(nt=0; nt<nThread-1; nt++)
   pthread_join(Threads[nt],0);
  delete[] Threads;
  delete[] TDs;

#else 
#ifndef USE_OPENMP
  nThread=1;
#else
#pragma omp parallel for schedule(dynamic,1), num_threads(nThread)
#endif
  for(nt=0; nt<nThread*100; nt++)
   { 
     ThreadData TD1;
     TD1.nt=nt;
     TD1.nThread=nThread*100;
     TD1.Args=Args;
     ACCMBThread((void *)&TD1);
   };
#endif

  if (G->LogLevel>=SCUFF_VERBOSELOGGING)
   Log("  %i/%i cache hits/misses",GlobalFIPPICache.Hits,GlobalFIPPICache.Misses);

}
