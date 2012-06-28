/*
 * RWGComposite.cc -- class constructor and other misc functions
 *                 -- for the RWG composite class 
 *                 
 * Homer Reid      -- 6/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libscuff.h>

#include "RWGComposite.h"

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

  OpenSurfaces=0 ;
  NumOpenSurfaces=0;
  OpenSurfaceIDs=0;

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
        RegionMPs=(MatProp *)realloc(RegionMPs, (NumRegions+1)*sizeof(MatProp *));
        RegionMPs[NumRegions-1]=MP;
        Log(" Adding new region (%i) for material %s",MP->Name);

      }
     else if ( !strcasecmp(Tokens[0],"OPENSURFACE") )
      { 
        if (NumTokens!=5 || strcasecmp(Tokens[2],"REGIONS") )
         ErrMsg=strdup("invalid syntax for OPENSURFACE keyword");
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

        // add a new open surface to this composite.
        //  (actually for now we just make a note of its existence; the
        //   actual OpenSurface data structure is not created until later)
        NumOpenSurfaces++;
        OpenSurfaceLabels=(char *)realloc(OpenSurfaceLabels, NumOpenSurfaces * sizeof(char *) );
        OpenSurfaceLabels[NumOpenSurfaces-1]=vstrdup(Tokens[1]);
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
   ErrMsg=vstrdup("COMPOSITE section must include a MESHFILE specification");
  if (NumRegions==0)
   ErrMsg=vstrdup("COMPOSITE section must include one or more REGION designations");
  if (NumOpenSurfaces==0)
   ErrMsg=vstrdup("COMPOSITE section must include one or more OPENSURFACE designations");

  InitRWGComposite(pMeshFileName, OTGT);
  
  free(pMeshFileName);

}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void RWGComposite::InitRWGComposite(char *pMeshFileName, GTransformation *OTGT)
{
  ErrMsg=0;
  //kdPanels = NULL;

  /*------------------------------------------------------------*/
  /*- try to open the mesh file.                                */
  /*------------------------------------------------------------*/
  FILE *MeshFile=fopen(pMeshFileName,"r");
  if (!MeshFile)
   ErrExit("could not open file %s",pMeshFileName);
   
  /*------------------------------------------------------------*/
  /*- initialize simple fields ---------------------------------*/
  /*------------------------------------------------------------*/
  NumEdges=NumPanels=NumVertices=NumRefPts=0;
  ContainingObject=0;
  MeshFileName=strdup(pMeshFileName);

  NumPanelsPerOpenSurface=(int *)mallocEC(NumOpenSurfaces*sizeof(int));

  /*------------------------------------------------------------*/
  /*- note: the 'OTGT' parameter to this function is distinct   */
  /*- from the 'GT' field inside the class body. the former is  */
  /*- an optional 'One-Time Geometrical Transformation' to be   */
  /*- applied to the object once at its creation. the latter    */
  /*- is designed to store a subsequent transformation that may */
  /*- be applied to the object, and is initialized to zero.     */
  /*------------------------------------------------------------*/
  GT=0;

  /*------------------------------------------------------------*/
  /*- Switch off based on the file type to read the mesh file:  */
  /*-  1. file extension=.msh    --> ReadGMSHFile              -*/
  /*-  2. file extension=.mphtxt --> ReadComsolFile            -*/
  /*------------------------------------------------------------*/
  char *p=GetFileExtension(MeshFileName);
  if (!p)
   ErrExit("file %s: invalid extension",MeshFileName);
  else if (!strcasecmp(p,"msh"))
   ReadGMSHFile(MeshFile,MeshFileName,OTGT);
  //else if (!strcasecmp(p,"mphtxt"))
  //  ReadComsolFile(MeshFile,MeshFileName,OTGT);
  else
   ErrExit("file %s: unknown extension %s",MeshFileName,p);

  /*------------------------------------------------------------*/
  /*- create OpenSurface structures ----------------------------*/
  /*------------------------------------------------------------*/
  OpenSurfaces=(OpenSurface **)mallocEC(NumOpenSurfaces*sizeof(OpenSurface *));
  BFIndexOffset=(int *)mallocEC(NumOpenSurfaces*sizeof(int));
  OpenSurface *OS;
  int np, npTOS;
  NumBFs=0;
  for(nos=0; nos<NumOpenSurfaces; nos++)
   { 
     // create a new OpenSurface structure for this open surface 
     OS = OpenSurfaces[nos]=(OpenSurface *)mallocEC(sizeof(OpenSurface *));

     // initialize the list of panels that reside on this open
     // surface. Note that the entries of the Panels[] array  
     // within the OpenSurface structure are just pointers to 
     // the same RWGPanel structures that were created by the 
     // meshfile input routine. 
     // note: npTOS = 'num panels, this open surface'
     OS->NumPanels = NumPanelsPerOpenSurface[nos];
     OS->Panels=(RWGPanel *)mallocEC(OS->NumPanels*sizeof(RWGPanel *));
     for(npTOS=np=0; np<NumPanels; np++)
      if ( Panels[np]->SurfaceIndex==nos )
       OpenSurfaces[nos]->Panels[npTOS++]= Panels[np];

     // analyze edge connectivity for this open surface 
     InitEdgeList();

     NumBFs += 2*(OS->NumEdges + OS->NumHEdges);

     BFIndexOffset[nos] = (nos==0) ? 0 : NumBFs + BFIndexOffset[nos-1];

   };

} 
