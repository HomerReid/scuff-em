
/***************************************************************/
/***************************************************************/
/***************************************************************/
struct OpenSurface
{  
  int NumEdges;
  RWGEdge  *Edges;  // internal edges to which we assign a full RWG function

  int NumHEdges;
  HRWGEdge *HEdges; // external edges to which we assign a half-RWG function

  int NumTotalEdges;

};

class RWGComposite
 {
   int NumRegions;
   int NumOpenSurfaces;
   int TotalBFs;

   // MP[0] = material properties of EXTERIOR medium
   // MP[r] = material properties of internal region r for r=1,...,NumRegions
   MatProp **MP;

   // OpenSurface structures for each open subregion of  
   // the object bounding surface 
   OpenSurface *OpenSurfaces;

   // BFIndexOffset[nosa] is the index of the first basis function
   // on open surface #nosa within the vector of all basis functions 
   // on this RWGComposite.
   int *BFIndexOffset;

   // for nos = 0 , 1, ..., NumOpenSurfaces-1, 
   // Regions[2*nos+0] and Regions[2*nos+1] are the indices of 
   // the two regions bounded by open surface #nos
   int *Regions;


   // EpsTF[nr] = Epsilon for region #nr at the present frequency
   //             (nr==0 for exterior medium)                   
   //             (note: TF stands for 'this frequency')
   cdouble *EpsTF;
   cdouble *MuTF;

 };

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

typedef struct AOCMBArgStruct 
{
  RWGGeometry *G;
  RWGComposite *C1, *C2;

} AOCMBArgStruct;


/***************************************************************/
/***************************************************************/
/***************************************************************/
void *AOCMBThread(void *data)
{
  /***************************************************************/
  /* extract local copies of fields in argument structure        */
  /***************************************************************/
  ThreadData *TD = (ThreadData *)data;
  AOCMBArgStruct *Args  = TD->Args;
  RWGGeometry *G        = Args->G;
  RWGComposite *C       = Args->C;

  
  int NOSA, NOSB;
  int nosa, nosb;
  int NTEA, NTEB;
  int ntea, nteb; // 'num total edges' for object a/b 

  int RegionA1, RegionA2, RegionB1, RegionB2; 
  int CommonRegion[2], NumCommonRegions;

  cdouble Eps1, Mu1, k1, EEPreFac1, MEPreFac1, MMPreFac1;
  cdouble Eps2, Mu2, k2, EEPreFac2, MEPreFac2, MMPreFac2;

  cdouble GC[2];
 
  B->Zero();

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
void AssembleBEMMatrixBlock_Composite()
{ 
 
  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 

}
