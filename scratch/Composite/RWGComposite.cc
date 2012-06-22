
/***************************************************************/
/* an **********************************************************/
/***************************************************************/
struct OpenSurface
{  
    int NumEdges;
    RWGEdge  *Edges;  // internal edges to which we assign a full RWG function

    int NumHEdges;
    HRWGEdge *HEdges; // external edges to which we assign a half-RWG function

};

class RWGComposite
 {
   int NumRegions;
   int NumOpenSurfaces;
   int TotalBFs;

   // MP[0] = material properties of EXTERIOR medium
   // MP[r] = material properties of internal region r for r=1,...,NumRegions
   MatProp **MP;

   OpenSurface *OpenSurfaces;

   // for nos = 0 , 1, ..., NumOpenSurfaces-1, 
   // Regions[2*nos+0] and Regions[2*nos+1] are the indices of 
   // the two regions bounded by open surface #nos
   int *Regions;

   // EpsTF[nr] = Epsilon for region #nr at the present frequency
   //             (nr==0 for exterior medium)                   
   cdouble *EpsTF;
   cdouble *MuTF;

 };

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOpenSurfaceInteraction()
{
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void *ABMBThread()
{
  int NOSA, NOSB;
  int nosa, nosb;
  int ntea, nteb; // 'num total edges' for object a/b 

  int CommonRegion1, CommonRegion2;
  for(nosa=0; nosa<NOSA; nosa++)
   for(nosb=0; nosb<NOSB; nosb++)
    { 
      /*--------------------------------------------------------------*/
      /* figure out if open surfaces #nosa and #nosb share 2, 1, or 0 */
      /* regions in common. if the answer is 0 then basis functions   */
      /* on the two surfaces do not interact with one another.        */
      /*--------------------------------------------------------------*/
      if (nosa==nosb)
       {
         CommonRegion1=Regions[2*nosa+0];
         CommonRegion2=Regions[2*nosa+1];
       }
      else 
       { CommonRegion2=-1;
         if (     Regions[2*nosa+0] == Regions[2*nosb+0]
               || Regions[2*nosa+0] == Regions[2*nosb+1] 
            ) CommonRegion1=Regions[2*nosa+0];
         else if (     Regions[2*nosa+1] == Regions[2*nosb+0]
                    || Regions[2*nosa+1] == Regions[2*nosb+1] 
                 ) CommonRegion1=Regions[2*nosa+1];
         else
          CommonRegion1=-1;
       };

      if (CommonRegion1==-1)
       continue;

      /*--------------------------------------------------------------*/
      /* now loop over all basis functions (both full and half RWG    */
      /* functions) on open surfaces nosa and nosb.                   */
      /*--------------------------------------------------------------*/
  for(ntea=0; nbfa<a; nea++)
   for(nbf=Symmetric*nea; neb<NEb; neb++)
    { 
      nt++;
      if (nt==TD->nThread) nt=0;
      if (nt!=TD->nt) continue;

      if (G->LogLevel>=SCUFF_VERBOSELOGGING)
       for(int PerCent=0; PerCent<9; PerCent++)
        if ( neb==Symmetric*nea &&  (nea == (PerCent*NEa)/10) )
         MutexLog("%i0 %% (%i/%i)...",PerCent,nea,NEa);

      if 

    };

    }; // for (nosa=0 ...) for nosb=0 ...)
}

void AssembleMatrixBlock()
{ 
 
  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 

}
