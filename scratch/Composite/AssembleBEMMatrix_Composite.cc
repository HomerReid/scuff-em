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

/*--------------------------------------------------------------*/
/*-  RWGComposite constructor that begins reading an open      -*/
/*-  file immediately after a line like COMPOSITE MyLabel.     -*/
/*-                                                            -*/
/*-  On return, if an object was successfully created, its     -*/
/*-   ErrMsg field is NULL, and the file read pointer points   -*/
/*-   the line immediately following the ENDCOMPOSITE line.    -*/
/*-                                                            -*/
/*-  If there was an error in parsing the COMPOSITE section,   -*/
/*-   ErrMsg field points to an error message string.          -*/
/*-                                                            -*/
/*- syntax:                                                    -*/
/*- COMPOSITE MyCompositeLabel                                 -*/
/*-                                                            -*/
/*-           MESHFILE MyMesh.msh                              -*/
/*-                                                            -*/
/*-           REGION 1 SiO2                                    -*/
/*-           REGION 2 Gold                                    -*/
/*-                                                            -*/
/*-           OPENSURFACE 1 REGIONS 0 1                        -*/
/*-           OPENSURFACE 2 REGIONS 0 2                        -*/
/*-           OPENSURFACE 3 REGIONS 1 2                        -*/
/*-                                                            -*/
/*- ENDCOMPOSITE                                               -*/
/*--------------------------------------------------------------*/
RWGComposite::RWGComposite(FILE *f, const char *pLabel, int *LineNum)
{ 
  /***************************************************************/
  /* initialize class fields *************************************/
  /***************************************************************/
  Label = strdup(pLabel);
  ErrMsg=0;

  OpenSurfaces = 0 ;
  NumOpenSurfaces = 0;
  OpenSurfaceIDs = 0;

  // we start with one region, namely, the exterior medium
  //  (the 0th slot of the RegionMPs array is filled in later)
  NumRegions = 1; 
  RegionMPs = (MatProp *) mallocEC( 1*sizeof(MatProp *) );

  Log("Processing COMPOSITE %s",pLabel);

  /***************************************************************/
  /* read lines from the file one at a time **********************/
  /***************************************************************/
  char Line[MAXSTR], LineCopy[MAXSTR];
  int NumTokens, TokensConsumed;
  char *Tokens[MAXTOK];
  int ReachedTheEnd=0;
  char *pMeshFileName=0;
  GTransformation *OTGT=0; // 'one-time geometrical transformation'
  int RegionID, RegionID1, RegionID2, OpenSurfaceID;
  MatProp *MP;$
  while ( ReachedTheEnd==0 && fgets(Line, MAXSTR, f) )
   { 
     (*LineNum)++;
     strcpy(LineCopy,Line);
     NumTokens=Tokenize(LineCopy, Tokens, MAXTOK);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

     /*--------------------------------------------------------------*/
     /*- switch off based on first token on the line ----------------*/
     /*--------------------------------------------------------------*/
     if ( !strcasecmp(Tokens[0],"MESHFILE") )
      { if (NumTokens!=2)
         { ErrMsg=strdup("MESHFILE keyword requires one argument");
           return;
         };
        pMeshFileName=strdup(Tokens[1]);
      }
     else if ( !strcasecmp(Tokens[0],"REGION") )
      { 
        if (NumTokens!=3)
         { ErrMsg=strdup("REGION keyword requires two arguments ");
           return;
         };

        if ( 1!=sscanf(Tokens[1],"%i", &RegionID) || RegionID != (NumRegions+1) )
         { ErrMsg=vstrdup("REGIONs must be numbered sequentially starting from 1");
           return; 
         };

        MP=new MatProp(Tokens[2]);
        if (MP->ErrMsg)
         { ErrMsg=vstrdup("material %s: %s",Tokens[2],MP->ErrMsg);
           return; 
         };

        // add a new region to the list of regions for this composite
        NumRegions++;
        RegionMPs=(MatProp *)realloc( RegionMPs, NumRegions*sizeof(MatProp *));
        RegionMPs[NumRegions-1]=MP;
        Log(" Adding new region (%i) for material %s",MP->Name);

      }
     else if ( !strcasecmp(Tokens[0],"OPENSURFACE") )
      { 
        if (NumTokens!=5 || strcasecmp(Tokens[2],"REGIONS") )
         ErrMsg=strdup("invalid syntax for OPENSURFACE keyword");
        if ( !ErrMsg && 1!=sscanf(Tokens[1],"%i",&OpenSurfaceID) )
         ErrMsg=vstrdup("syntax error");
        if ( !ErrMsg && 1!=sscanf(Tokens[3],"%i",&RegionID1) )
         ErrMsg=vstrdup("syntax error");
        if ( !ErrMsg && 1!=sscanf(Tokens[4],"%i",&RegionID2) )
         ErrMsg=vstrdup("syntax error");

        if ( RegionID1<0 || RegionID1>NumRegions ) 
         ErrMsg=vstrdup("invalid region ID (%i)",RegionID1);
        if ( RegionID2<0 || RegionID2>NumRegions ) 
         ErrMsg=vstrdup("invalid region ID (%i)",RegionID1);

        if (ErrMsg)
         return;

        NumOpenSurfaces++;
        OpenSurfaceIDs=(int *)realloc(OpenSurfaceIDs, NumOpenSurfaces * sizeof(int) );
        OpenSurfaceID[NumOpenSurfaces-1]=OpenSurfaceID;
        Regions = (int *)realloc(Regions, 2*NumOpenSurfaces*sizeof(int) );
        Regions[ 2*(NumOpenSurfaces-1) + 0 ] = RegionID1;
        Regions[ 2*(NumOpenSurfaces-1) + 1 ] = RegionID2;
        Log(" Adding new open surface (mesh physical surface %i) bounding %s and %s", OpenSurfaceID, 
              RegionMPs[RegionID1] ? RegionMPs[RegionID1]->Name : "exterior",
              RegionMPs[RegionID2] ? RegionMPs[RegionID2]->Name : "exterior");

      }
     else if ( !strcasecmp(Tokens[0],"DISPLACED") || !strcasecmp(Tokens[0],"ROTATED") )
      { 
        // try to parse the line as a geometrical transformation.
        // note that OTGT is used as a running GTransformation that may
        // be augmented by multiple DISPLACED ... and/or ROTATED ...
        // lines within the OBJECT...ENDOBJECT section, and which is 
        // applied to the object at its birth and subsequently discarded.
        // in particular, OTGT is NOT stored as the 'GT' field inside 
        // the Object class, which is intended to be used for 
        // transformations that are applied and later un-applied 
        // during the life of the object. 
	OTGT = new GTransformation(Tokens, NumTokens, &ErrMsg, &TokensConsumed);
        if (ErrMsg)
         return;
        if (TokensConsumed!=NumTokens) 
         { ErrMsg=strdup("junk at end of line");
           return;
         };
      }
     else if ( !strcasecmp(Tokens[0],"ENDCOMPOSITE") )
      { 
        ReachedTheEnd=1;
      }
     else
      { ErrMsg=vstrdup("unknown keyword %s in COMPOSITE section",Tokens[0]);
        return;
      };
   }; 

  if (pMeshFileName==0)
   ErrMsg=vstrdup("COMPOSITE section must include a MESHFILE specification",Tokens[0]);

  InitRWGComposite(pMeshFileName, OTGT);
  
  free(pMeshFileName);

}
