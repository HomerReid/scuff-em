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

/*
 * RWGGeometry.cc -- implementation of some methods in the RWGGeometry class
 *
 * homer reid      -- 3/2007 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#include <fenv.h>

#include <libhrutil.h>
#include <BZIntegration.h> // needed for GetRLBasis

#include "libscuff.h"

namespace scuff {

#define MAXSTR 1000
#define MAXTOK 50
#define MAXMMJ 10

/***************************************************************/
/* initialization of static class variables                    */
/***************************************************************/
bool RWGGeometry::UseHRWGFunctions=false;
bool RWGGeometry::UseHighKTaylorDuffy=true;
bool RWGGeometry::UseTaylorDuffyV2P0=true;
bool RWGGeometry::DisableCache=false;

/***********************************************************************/
/* subroutine to parse the MEDIUM...ENDMEDIUM section in a .scuffgeo   */
/* file. (currently the only keyword supported for this section is     */
/* MATERIAL xx).                                                       */
/* note: this is legacy syntax for changing the material properties of */
/* the exterior region; the modern way to do it would be to use a      */
/* line like REGION Exterior MATERIAL SiC.                             */
/***********************************************************************/
void RWGGeometry::ProcessMEDIUMSection(FILE *f, char *FileName, int *LineNum)
{
  char Line[MAXSTR];
  while( fgets(Line,MAXSTR,f) )
   { 
     (*LineNum)++;
     char *Tokens[MAXTOK];
     int NumTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

     if ( !StrCaseCmp(Tokens[0],"MATERIAL") )
      {
        if (NumTokens!=2)
         ErrExit("%s:%i: syntax error",FileName,*LineNum);
        if (RegionMPs[0]) delete RegionMPs[0];
        RegionMPs[0] = new MatProp(Tokens[1]); 
        if (RegionMPs[0]->ErrMsg)
         ErrExit("%s:%i: %s",FileName,*LineNum,RegionMPs[0]->ErrMsg);
        Log("Setting material properties of exterior region to %s.",Tokens[1]);
      }
     else if ( !StrCaseCmp(Tokens[0],"ENDMEDIUM") )
      return;
     else
      ErrExit("%s:%i: unknown keyword %s",FileName,*LineNum,Tokens[0]);
   }
  ErrExit("%s: unexpected end of file",FileName);
}

/***************************************************************/
/* this is a hopefully temporary routine that asserts that the */
/* user-specified lattice comports with our current limitations*/
/* regarding the types of lattices that are supported;         */
/* this routine will hopefully vanish in the future as we      */
/* increase support for all types of lattices.                 */
/***************************************************************/
void CheckLattice(HMatrix *LBasis)
{
  int LDim = LBasis->NC;

  if (LDim==1)
   { if (    LBasis->GetEntryD(1,0)!=0.0
          || LBasis->GetEntryD(2,0)!=0.0
        )
      ErrExit("1D lattice vector must be parallel to x axis");
   }
  else if (LDim==2)
   { if (    LBasis->GetEntryD(2,0)!=0.0
          || LBasis->GetEntryD(2,1)!=0.0
        )
      ErrExit("2D lattice vectors must lie in xy plane");

     if (    LBasis->GetEntryD(1,0)!=0.0
          || LBasis->GetEntryD(0,1)!=0.0
       )
      ErrExit("non-rectangular 2D lattices not yet supported");
   }
  else if (LDim==3)
   ErrExit("3D lattices not yet supported");
}

/***********************************************************************/
/* subroutine to parse the LATTICE...ENDLATTICE section in a .scuffgeo */
/* file.                                                               */
/***********************************************************************/
void RWGGeometry::ProcessLATTICESection(FILE *f, char *FileName, int *LineNum)
{
  char Line[MAXSTR];
  double LBV[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  LDim=0;
  while( fgets(Line,MAXSTR,f) )
   { 
     (*LineNum)++;

     char *Tokens[MAXTOK];
     int NumTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

     if ( !StrCaseCmp(Tokens[0],"VECTOR") )
      {
        if (LDim==MAXLDIM)
         ErrExit("%s:%i: too many lattice vectors",FileName,*LineNum);
        
        if (NumTokens<2 || NumTokens>4) 
         ErrExit("%s:%i: lattice vectors must have 1, 2, or 3 components",FileName,*LineNum);
        for(int nc=0; nc<(NumTokens-1); nc++)
         if ( 1 != sscanf(Tokens[nc+1],"%le",&(LBV[LDim][nc])))
          ErrExit("%s:%i: invalid vector component %s",FileName,*LineNum,Tokens[nc+1]);

        Log("Adding lattice basis vector (%g,%g,%g).",LBV[LDim][0],LBV[LDim][1],LBV[LDim][2]);
        LDim++;
      }
     else if ( !StrCaseCmp(Tokens[0],"ENDLATTICE") )
      { 
        LBasis=new HMatrix(3, LDim);
        for(int nd=0; nd<LDim; nd++)
         for(int i=0; i<3; i++)
          LBasis->SetEntry(i,nd,LBV[nd][i]);

        RLBasis=GetRLBasis(LBasis, &LVolume, &RLVolume);
        if (LVolume==0.0)
         ErrExit("%s:%i: lattice has empty unit cell",FileName,*LineNum);

        return; 
      }
     else
      ErrExit("%s:%i: unknown keyword %s",FileName,*LineNum,Tokens[0]);
   }
  ErrExit("%s: unexpected end of file",FileName);
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
void RWGGeometry::AddRegion(char *RegionLabel, char *MaterialName, int LineNum)
{
  int RegionID=0; 

  /*--------------------------------------------------------------*/
  /*- first check to see if the region label is 'Exterior', in   -*/
  /*- which case we aren't adding a new region but just redefining*/
  /*- the properties of the exterior region.                      */
  /*--------------------------------------------------------------*/
  if ( !StrCaseCmp(RegionLabel,"EXTERIOR") )
   { 
     RegionID = 0; 
   }
  else if ( !StrCaseCmp(RegionLabel,"PEC") )
   { 
     ErrExit("%s:%i: REGIONs cannot be PEC (use an OBJECT instead)",GeoFileName,LineNum);
   }
  else
   { 
     RegionLabels = (char **)reallocEC( RegionLabels, (NumRegions+1)*sizeof(char *));
     RegionLabels[NumRegions] = strdupEC(RegionLabel);
     RegionMPs = (MatProp **)reallocEC( RegionMPs, (NumRegions+1)*sizeof(MatProp *));
     NumRegions++;

     RegionID = NumRegions-1;
   }

  RegionMPs[RegionID] = new MatProp(MaterialName);
  if ( RegionMPs[RegionID]->ErrMsg )
   ErrExit("%s:%i: %s\n",GeoFileName,LineNum,RegionMPs[RegionID]->ErrMsg);
}

/***************************************************************/
/* 20150818 convenience hack to support special file names of  */
/* the following forms                                         */
/*  MESH__MyMesh                                               */
/*  MESH__MyMesh__MAT__MatProp                                 */
/*  MESH__MyMesh__MAT__MatProp__MED__MatProp                   */
/***************************************************************/
bool ProcessSpecialSCUFFGeoFileName(const char *GeoFileName)
{
  // check for a special geometry file name, which
  //  begins with MESH__
  int N=strlen(GeoFileName);
  if ( N<7 || strncmp(GeoFileName,"MESH__",6) )
   return false;

  char GFNCopy[MAXSTR];
  strncpy(GFNCopy, GeoFileName+6, MAXSTR);

  char *MeshStr=GFNCopy;
  char *MatStr=strstr(GFNCopy,"__MAT__");
  char *MedStr=strstr(GFNCopy,"__MED__");

  if (MatStr)
   { *MatStr=0;
     MatStr+=7;
   }
  if (MedStr)
   { *MedStr=0;
     MedStr+=7;
   }

  Log("Processing special file name: mesh=%s",MeshStr);
  if (MatStr)
   Log("           special file name: mat=%s",MatStr);
  if (MedStr)
   Log("           special file name: med=%s",MedStr);

  FILE *f=vfopen("%s.scuffgeo","w",GeoFileName);
  fprintf(f,"OBJECT %s\n",GetFileBase(MeshStr));
  fprintf(f," MESHFILE %s.msh\n",GetFileBase(MeshStr));
  if (MatStr) 
   fprintf(f," MATERIAL %s\n",MatStr);
  fprintf(f,"ENDOBJECT\n");
  if (MedStr) 
   fprintf(f,"\nMEDIUM\n MATERIAL %s\nENDMEDIUM\n",MedStr);
  fclose(f);

  return true;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MMJTHRESH 1.0e-4
void RWGGeometry::DetectMultiMaterialJunctions()
{
  NumMMJs=0;
  MultiMaterialJunctions=0;

  // loop over all edges on all surfaces
  for(int ns=0; ns<NumSurfaces; ns++)
   for(int ne=0; ne<Surfaces[ns]->NumEdges; ne++)
    { 
      RWGSurface *S = Surfaces[ns];
      RWGEdge *E    = S->Edges[ne];
      if (E->iQM!=-1) continue;

      // check that this edge not already part of an MMJ
      bool InMMJ=false;
      if (ns>0)
       for(int nMMJ=0; !InMMJ && nMMJ<NumMMJs; nMMJ++)
        { MMJData *MMJ=MultiMaterialJunctions[nMMJ];
          for(int n=1; !InMMJ && n<MMJ->NumEdges; n++)
           if ( ns==MMJ->SurfaceIndices[n] && ne==MMJ->EdgeIndices[n])
            InMMJ=true;
        };
      if (InMMJ) continue;

      // count number of edges on other surfaces that coincide with this edge
      double *V1    = S->Vertices + 3*E->iV1;
      double *V2    = S->Vertices + 3*E->iV2;
      double Length = E->Length;
      int NumEdges=0;
      int SurfaceIndices[MAXMMJ];
      int EdgeIndices[MAXMMJ];
      for(int nsp=ns+1; nsp<NumSurfaces; nsp++)
       for(int nep=0; nep<Surfaces[nsp]->NumEdges; nep++)
        { 
          RWGSurface *SP = Surfaces[nsp];
          RWGEdge *EP    = SP->Edges[nep];
          if (EP->iQM!=-1) continue;

          double *V1P   = SP->Vertices + 3*EP->iV1;
          double *V2P   = SP->Vertices + 3*EP->iV2;
   
          bool Match1  = (     ( VecDistance(V1,V1P)<MMJTHRESH*Length )
                            && ( VecDistance(V2,V2P)<MMJTHRESH*Length )
                         );
          bool Match2  = (     ( VecDistance(V1,V2P)<MMJTHRESH*Length )
                            && ( VecDistance(V2,V1P)<MMJTHRESH*Length )
                         );
          if (Match1==false && Match2==false) continue;

          if (NumEdges==MAXMMJ)
           ErrExit("Surface %i edge %i: too many surfaces meeting at multi-material junction",ns,ne);

          if ( S->IsPEC != SP->IsPEC)
           ErrExit("surfaces in a MMJ must be all PEC or all non-PEC "
                   "{surface %s edge %i <> surface %s edge %i}",
                   S->Label,ne,SP->Label,nep);

          if (NumEdges==0)
           { SurfaceIndices[NumEdges] = ns;
             EdgeIndices[NumEdges++]  =ne;
           };
          SurfaceIndices[NumEdges] = nsp;
          EdgeIndices[NumEdges++]  = nep;

        };
      if (NumEdges==0) continue;

      // we have detected a new multi-material junction
      MMJData *Data = (MMJData *)mallocEC(sizeof(MMJData));
      Data->NumEdges       = NumEdges;
      Data->SurfaceIndices = (int *)memdup(SurfaceIndices, NumEdges*sizeof(int));
      Data->EdgeIndices    = (int *)memdup(EdgeIndices, NumEdges*sizeof(int));

      MultiMaterialJunctions =
       (MMJData **)reallocEC(MultiMaterialJunctions, (NumMMJs+1) * sizeof(MMJData *));
      MultiMaterialJunctions[NumMMJs++] = Data;

      if (LogLevel>=SCUFF_VERBOSE2)
       { Log(" MMJ #03i: ",NumMMJs);
         for(int n=0; n<NumEdges; n++)
          LogC(" {%i,%04i} ",SurfaceIndices[n],EdgeIndices[n]);
       };

    };
  Log("Detected %i multi-material junctions.",NumMMJs);
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
RWGGeometry::RWGGeometry(const char *pGeoFileName, int pLogLevel)
{ 
  /***************************************************************/
  /* NOTE: i am not sure where to put this. put it here for now. */
  /***************************************************************/
  MatProp::SetLengthUnit(1.0e-6);
   
  /***************************************************************/
  /* initialize simple fields ************************************/
  /***************************************************************/
  LogLevel=pLogLevel;
  NumSurfaces=TotalBFs=TotalEdges=TotalPanels=0;
  GeoFileName=strdupEC(pGeoFileName);
  Surfaces=0;
  AllSurfacesClosed=1;
  Substrate=0;

  /***************************************************************/
  /* initially assume compact (non-periodic) geometry            */
  /***************************************************************/
  LDim=0;
  LBasis=RLBasis=0;
  LVolume=RLVolume=0.0;
  for(int nd=0; nd<MAXLDIM; nd++)
   { NumStraddlers[nd]=NULL;
     RegionIsExtended[nd]=NULL;
   }
  tolVecClose=0.0; // to be updated once mesh is read in
  TBlockCacheNameAddendum=0;

  // we always start with a single Region, for the exterior,
  // taken to be vacuum by default
  NumRegions=1;
  RegionLabels=(char **)mallocEC(1*sizeof(char *));
  RegionLabels[0]=strdupEC("EXTERIOR");
  RegionMPs=(MatProp **)mallocEC(1*sizeof(MatProp *));
  RegionMPs[0] = new MatProp("VACUUM");

  /***************************************************************/
  /* check for various environment variables *********************/
  /***************************************************************/
  char *s;
  if ( (s=getenv("SCUFF_LOGLEVEL")) )
   {      if ( !strcasecmp(s, "NONE"     ) ) LogLevel=SCUFF_NOLOGGING;
     else if ( !strcasecmp(s, "TERSE"    ) ) LogLevel=SCUFF_TERSELOGGING;
     else if ( !strcasecmp(s, "VERBOSE"  ) ) LogLevel=SCUFF_VERBOSELOGGING;
     else if ( !strcasecmp(s, "VERBOSE2" ) ) LogLevel=SCUFF_VERBOSE2;
   }

  RWGGeometry::DisableCache=CheckEnv("SCUFF_DISABLE_CACHE");
  RWGGeometry::UseHRWGFunctions=CheckEnv("SCUFF_HALF_RWG");

  if (CheckEnv("SCUFF_ABORT_ON_FPE"))
   {
#ifndef __APPLE__
     feenableexcept(FE_INVALID | FE_OVERFLOW);
     Log("Enabling abort-on-floating-point-exception.");
#else
     Log("Can't enable abort-on-floating-point-exception on OS X because feenableexcept is not available.");
#endif
   }

  /***************************************************************/
  /* try to open input file **************************************/
  /***************************************************************/
  if (ProcessSpecialSCUFFGeoFileName(GeoFileName))
   GeoFileName=vstrappend(GeoFileName,".scuffgeo");
  char *Dir;
  FILE *f=fopenPath(getenv("SCUFF_GEO_PATH"),GeoFileName,"r",&Dir);
  if (!f)
   ErrExit("could not open %s",GeoFileName);
  if (strcmp(Dir,"."))
   Log("Found file %s in directory %s",GeoFileName,Dir);
  char *GeoFileDir=strdup(Dir);
  AppendEnv("SCUFF_MESH_PATH",":","%s/../mshFiles",GeoFileDir);

  /***************************************************************/
  /* read and process lines from input file one at a time        */
  char Line[MAXSTR];
  int LineNum=0; 
  while( fgets(Line,MAXSTR,f) )
   { 
     LineNum++;

     /*--------------------------------------------------------------*/
     /*- break up line into tokens; skip blank lines and comments ---*/
     /*--------------------------------------------------------------*/
     char *Tokens[MAXTOK];
     int nTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( nTokens==0 || Tokens[0][0]=='#' )
      continue; 
    
     /*--------------------------------------------------------------*/
     /*- switch off based on first token on the line ----------------*/
     /*--------------------------------------------------------------*/
     if ( !StrCaseCmp(Tokens[0],"MEDIUM") )
      { 
        ProcessMEDIUMSection(f,GeoFileName,&LineNum);
      }
     else if ( !StrCaseCmp(Tokens[0],"LATTICE") )
      { 
        ProcessLATTICESection(f,GeoFileName,&LineNum);
        //UseHRWGFunctions=false;
      }
     else if (    !StrCaseCmp(Tokens[0],"MESHPATH")
               || !StrCaseCmp(Tokens[0],"MESH_PATH")
             )
      { 
        if ( nTokens!=2 )
         ErrExit("%s:%i: invalid MESHPATH specification",GeoFileName,LineNum);

        // if MESHPATH starts with a '/' then it is an absolute path; otherwise
        // it is taken relative to the path in which the .scuffgeo file was found
        if (Tokens[1][0]=='/')
         AppendEnv("SCUFF_MESH_PATH",":","%s",Tokens[1]);
        else 
         AppendEnv("SCUFF_MESH_PATH",":","%s/%s",GeoFileDir,Tokens[1]);
      }
     else if ( !StrCaseCmp(Tokens[0],"MATERIAL") )
      {
        /*--------------------------------------------------------------*/
        /* hand off to MatProp class constructor to parse this section  */
        /*--------------------------------------------------------------*/
        if ( nTokens==1 )
         ErrExit("%s:%i: no name given for MATERIAL ",GeoFileName,LineNum);
        else if ( nTokens>2 )
         ErrExit("%s:%i: syntax error",GeoFileName,LineNum);
         
        char *ErrMsg=AddMaterialToMatPropDataBase(f, GeoFileName, Tokens[1], &LineNum);
        if (ErrMsg)
         ErrExit("%s:%i: %s",GeoFileName,LineNum,ErrMsg);
      }
     else if ( !StrCaseCmp(Tokens[0],"REGION") )
      {
        /*--------------------------------------------------------------*/
        /*- add a new region to our geometry ---------------------------*/
        /*--------------------------------------------------------------*/
        if ( nTokens!=4 || StrCaseCmp(Tokens[2],"MATERIAL") )
         ErrExit("%s:%i: syntax error",GeoFileName,LineNum);
        AddRegion(Tokens[1], Tokens[3], LineNum);
      }
     else if ( !StrCaseCmp(Tokens[0],"OBJECT") || !StrCaseCmp(Tokens[0],"SURFACE") )
      { 
        RWGSurface *S=0;
        if ( nTokens==2 )
         S=new RWGSurface(f,Tokens[1],&LineNum,Tokens[0]);
        else if (nTokens!=1)
         ErrExit("%s:%i: syntax error",GeoFileName,LineNum);
        else
         { char Label[MAXSTR];
           snprintf(Label,MAXSTR,"Surface_%i",NumSurfaces+1);
           S=new RWGSurface(f,Label,&LineNum,Tokens[0]);
         }

        if (S->ErrMsg)
         ErrExit("%s:%i: %s",GeoFileName,LineNum,S->ErrMsg); 

        if ( S->IsObject )
         { 
            /* for an OBJECT, if it is non-PEC, we need to add a new Region for  */
            /* the region interior to the object. In this case, for now, we set  */
            /* S->RegionLabels[0] = EXTERIOR and S->RegionLabels[1] = the label  */
            /* specified for the object in question. We may subsequently need to */
            /* change S->RegionLabels[0] (and S->RegionIndices[0]) to another    */
            /* region if we detect that the object is embedded in another object.*/
            S->RegionIndices[0] = 0;
            if ( S->IsPEC )
             S->RegionIndices[1] = -1;
            else
             { AddRegion(S->RegionLabels[1], S->MaterialName, S->MaterialRegionsLineNum );
               S->RegionIndices[1] = NumRegions - 1 ;
             }
         }
        else
         { 
           AllSurfacesClosed=0;

           /* On the other hand, for a SURFACE, we need to check that the       */
           /* REGIONS specified in the SURFACE...ENDSURFACE description         */
           /* are regions that have been previously declared.                   */
           S->RegionIndices[0] = GetRegionByLabel(S->RegionLabels[0]);
           if (S->RegionIndices[0]==-1)
            ErrExit("%s:%i: unknown region %s",GeoFileName,S->MaterialRegionsLineNum,S->RegionLabels[0]);

           if( S->IsPEC )
            S->RegionIndices[1] = -1;
           else
            { S->RegionIndices[1] = GetRegionByLabel(S->RegionLabels[1]);
              if (S->RegionIndices[1]==-1)
               ErrExit("%s:%i: unknown region %s",GeoFileName,S->MaterialRegionsLineNum,S->RegionLabels[1]);
            };
         };

        NumSurfaces++;
        Surfaces=(RWGSurface **)realloc(Surfaces, NumSurfaces*sizeof(RWGSurface *) );
        Surfaces[NumSurfaces-1]=S;
        S->Index=NumSurfaces-1;
      }
     else if ( !StrCaseCmp(Tokens[0],"SUBSTRATE") )
      {
        /*--------------------------------------------------------------*/
        /* hand off to LayeredSubstrate class constructor to parse this section  */
        /*--------------------------------------------------------------*/
        if ( nTokens>1 )
         ErrExit("%s:%i: syntax error",GeoFileName,LineNum);
         
        Substrate = new LayeredSubstrate(f, &LineNum);
        if (Substrate->ErrMsg)
         ErrExit("%s:%i: %s",GeoFileName,LineNum,Substrate->ErrMsg);
      }
     else 
      { 
        /*--------------------------------------------------------------*/
        /* unknown keyword                                              */
        /*--------------------------------------------------------------*/
        ErrExit("%s:%i: syntax error",GeoFileName,LineNum);
      }

   } // while( fgets(Line,MAXSTR,f) )

  /*******************************************************************/
  /* Before doing any vector comparisons with VecClose (e.g.
     in InitPBCData), initialize tolVecClose to 1e-3*(min edge len)  */
  /*******************************************************************/
  tolVecClose = 1./0.; // infinity
  for(int ns=0; ns<NumSurfaces; ns++)
    for(int ne=0; ne<Surfaces[ns]->NumEdges; ne++) 
      tolVecClose = fmin(tolVecClose, Surfaces[ns]->Edges[ne]->Length);
  tolVecClose *= 1e-3;
  for(int ns=0; ns<NumSurfaces; ns++)
    Surfaces[ns]->tolVecClose = tolVecClose;
    
  /*******************************************************************/
  /* if a lattice is present, then we need to do some preliminary    */
  /* setup                                                           */
  /*******************************************************************/
  if (LBasis)
   { LDim=LBasis->NC;
     InitPBCData();
   }

  /*******************************************************************/
  /* Autodetect nesting relationships & topologically sort           */
  /* (so that if A contains B, then B comes after A).                */
  /* Note that the Contains() function and the autodetection of      */
  /* containership relations are only applicable to CLOSED           */
  /* surfaces; for OPEN surfaces there is no autodetection, and      */
  /* we rely instead on the user to specify explicitly which         */
  /* surfaces bound which regions.                                   */
  /*******************************************************************/
  for (int ns = 0; ns < NumSurfaces; ++ns)
   Surfaces[ns]->InitkdPanels(false, LogLevel);
  for (int ns = 1; ns < NumSurfaces; ++ns) 
   { 
     RWGSurface *S=Surfaces[ns];
     int NewSlotForS=-1;
     for (int nsp = ns-1; nsp >= 0; --nsp) // innermost to outermost
      { 
        RWGSurface *SP = Surfaces[nsp];
        if ( SP->Contains(S) ) 
         {
           if (SP->IsPEC) 
            ErrExit("%s: PEC object %s cannot contain object %s",GeoFileName,SP->Label,S->Label);

	   // if SP contains S, then we set the exterior region for S to be
	   // the interior region for SP.
	   S->RegionIndices[0] = SP->RegionIndices[1];
	   free(S->RegionLabels[0]);
	   S->RegionLabels[0]=strdupEC(SP->RegionLabels[1]);
	   break;
         } 
        else if ( S->Contains(SP) )
         NewSlotForS = nsp;
      }

     if ( NewSlotForS != -1 )
      {  
        // for each surface between the current occupant of slot #NewSlotForS  
        // and slot #ns, move the surface one slot forward in the array and
        // adjust its region ID as necessary
        for(int nsp = ns; nsp > NewSlotForS; --nsp) 
         { 
           Surfaces[nsp] = Surfaces[nsp-1];

           RWGSurface *SP = Surfaces[nsp];
           if (    SP->RegionIndices[0] == S->RegionIndices[0] 
                && S->Contains(SP) 
              ) 
            { 
              if (S->IsPEC) 
               ErrExit("%s: PEC object %s cannot contain object %s",GeoFileName,S->Label,SP->Label);
              SP->RegionIndices[0] = S->RegionIndices[1];
	      free(SP->RegionLabels[0]);
	      SP->RegionLabels[0]=strdupEC(S->RegionLabels[1]);
            };
         }
        Surfaces[NewSlotForS] = S;
      }
    } // for(ns = 1 ...)

  /*******************************************************************/
  /* make sure that all panel normals point in the correct direction.*/
  /*******************************************************************/
  int NumFlipped=0;
  for(int ns=0; ns<NumSurfaces; ns++)
   { 
     RWGSurface *S=Surfaces[ns];
     int ExteriorRegion=S->RegionIndices[0]; 
     int InteriorRegion=S->RegionIndices[1]; 
     double X[3];
     RWGPanel *P;
     for(int np=0; np<S->NumPanels; np++)
      {
        /* how it works: if we start at the panel centroid */
        /* and travel a short distance in the direction of */
        /* the panel normal, then that should place us in  */
        /* the exterior region (since the panel normal     */
        /* should point into the exterior region). if,     */
        /* instead, we wind up in the interior, then the   */
        /* panel normal is backwards and must be flipped.  */
        /* (if we wind up in neither the interior nor the  */
        /* exterior then something more serious is wrong.) */
        P=S->Panels[np];
        VecScaleAdd(P->Centroid, 0.5*P->Radius, P->ZHat, X);
        if (PointInRegion(ExteriorRegion,X))
         { 
           ; // don't need to do anything
         }
        else if (PointInRegion(InteriorRegion,X))
         { 
           P->ZHatFlipped=true;
           VecScale(P->ZHat, -1.0);
           NumFlipped++;
         }
/*
        else
         ErrExit("%s:%i: internal error (%i,%i)",__FILE__,__LINE__,ns,np);
*/

      }
   }
  Log("Flipped %i panel normals to comport with region definitions.",NumFlipped);
 
  /*******************************************************************/
  /* compute average panel area for statistical bookkeeping purposes */
  /*******************************************************************/
  AveragePanelArea=0.0; 
  for(int ns=0; ns<NumSurfaces; ns++)
   for(int np=0; np<Surfaces[ns]->NumPanels; np++)
    AveragePanelArea+=Surfaces[ns]->Panels[np]->Area;
  AveragePanelArea/=((double) TotalPanels);

  /***************************************************************/
  /* initialize arrays of basis-function, edge, and panel offsets*/
  /***************************************************************/
  BFIndexOffset    = (int *)mallocEC(NumSurfaces*sizeof(int));
  EdgeIndexOffset  = (int *)mallocEC(NumSurfaces*sizeof(int));
  PanelIndexOffset = (int *)mallocEC(NumSurfaces*sizeof(int));
  TotalBFs    = Surfaces[0]->NumBFs;
  TotalEdges  = Surfaces[0]->NumEdges;
  TotalPanels = Surfaces[0]->NumPanels;
  BFIndexOffset[0]=EdgeIndexOffset[0]=PanelIndexOffset[0]=0;
  for(int ns=1; ns<NumSurfaces; ns++)
   { 
     TotalBFs    += Surfaces[ns]->NumBFs;
     TotalEdges  += Surfaces[ns]->NumEdges;
     TotalPanels += Surfaces[ns]->NumPanels;

     BFIndexOffset[ns]    = BFIndexOffset[ns-1]    + Surfaces[ns-1]->NumBFs;
     EdgeIndexOffset[ns]  = EdgeIndexOffset[ns-1]  + Surfaces[ns-1]->NumEdges;
     PanelIndexOffset[ns] = PanelIndexOffset[ns-1] + Surfaces[ns-1]->NumPanels;
   }

  /***************************************************************/
  /* allocate space for cached epsilon and mu values *************/
  /***************************************************************/
  StoredOmega=0.0;
  EpsTF = (cdouble *)mallocEC(NumRegions * sizeof(cdouble));
  MuTF  = (cdouble *)mallocEC(NumRegions * sizeof(cdouble));

  /***************************************************************/
  /* initialize Mate[] array.                                    */
  /*                                                             */
  /* how it works:                                               */
  /*                                                             */
  /* (1) two surfaces are considered identical if                */
  /*     (a) they were read in from the same physical region of  */
  /*         the same mesh file, and                             */
  /*     (b) the two regions they bound have the same material   */
  /*         properties.                                         */
  /*                                                             */
  /* (2) Mate[] array: If surfaces i, j, k, ... are identical and*/
  /*                   i<j<k<..., then we set                    */
  /*                   Mate[i] = -1                              */
  /*                   Mate[j] = i                               */
  /*                   Mate[k] = i                               */
  /***************************************************************/
  Mate=(int *)mallocEC(NumSurfaces*sizeof(int));
  Mate[0]=-1;
  for(int ns=1; ns<NumSurfaces; ns++)
   { RWGSurface *S=Surfaces[ns];
     Mate[ns]=-1;
     int nr1=S->RegionIndices[0];
     int nr2=S->RegionIndices[1];
     for(int nsp=0; nsp<ns && Mate[ns]==-1; nsp++)
      { RWGSurface *SP=Surfaces[nsp];
        int nr1p=SP->RegionIndices[0];
        int nr2p=SP->RegionIndices[1];
        if (    ( !strcmp(S->MeshFileName, SP->MeshFileName) )
             && ( S->MeshTag == SP->MeshTag )
             && ( !strcmp(RegionMPs[nr1]->Name, RegionMPs[nr1p]->Name) )
             && (   (S->IsPEC && SP->IsPEC)
                 || (!S->IsPEC && !SP->IsPEC && !strcmp(RegionMPs[nr2]->Name, RegionMPs[nr2p]->Name) )
                )
           ) { Mate[ns]=nsp;
               Log("Noting that surface %i (%s) is a duplicate of surface %i (%s)...",ns,S->Label,nsp,SP->Label);
             }
      }
   }

  /***************************************************************/
  /* initialize SurfaceMoved[] array.                            */
  /* the values of this array are only defined after             */
  /* a call to the RWGGeometry::Transform() function, when we    */
  /* have SurfaceMoved[i]=1 if the ith surface was modified by   */
  /* the transformation.                                         */
  /***************************************************************/
  SurfaceMoved=(int *)mallocEC(NumSurfaces*sizeof(int));

  /***************************************************************/
  /* initialize all equivalent edge-pair tables to NULL (they    */
  /* are created lazily)                                         */
  /***************************************************************/
  if (!CheckEnv("SCUFF_IGNORE_EEPS"))
   for(int ns=0; ns<NumSurfaces; ns++)
    EEPTables.push_back( vector<EquivalentEdgePairTable *>(NumSurfaces, NULL) );

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  NumMMJs=0;
  MultiMaterialJunctions=0;
  if (UseHRWGFunctions)
   DetectMultiMaterialJunctions();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FIBBICaches = (void **)mallocEC(NumSurfaces*sizeof(void *));
  for(int ns=0; ns<NumSurfaces; ns++)
   if (Mate[ns]!=-1)
    FIBBICaches[ns] = FIBBICaches[ Mate[ns] ];
   else
    FIBBICaches[ns] = CreateFIBBICache(Surfaces[ns]->MeshFileName);
  fclose(f);
}

/***************************************************************/
/* RWGGeometry class destructor *******************************/
/***************************************************************/
RWGGeometry::~RWGGeometry()
{
  for(int ns=0; ns<NumSurfaces; ns++)
   delete Surfaces[ns];

  free(Surfaces);

  for(int nr=0; nr<NumRegions; nr++)
   { delete RegionMPs[nr];
     free(RegionLabels[nr]); 
   }
  free(RegionMPs);
  free(RegionLabels);

  free(BFIndexOffset);
  free(EdgeIndexOffset);
  free(PanelIndexOffset);
  free(EpsTF);
  free(MuTF);
  free(Mate);
  free(SurfaceMoved);
  free(GeoFileName);

  for(int ns=0; ns<NumSurfaces; ns++)
   if (Mate[ns]==-1)
    DestroyFIBBICache(FIBBICaches[ns]);
  free(FIBBICaches);

}

/***************************************************************/
/* return the index of the region whose label is Label. (This  */
/* index is between 0 and NumRegions-1, and is the index into  */
/* the RegionLabels[] and RegionMPs[] arrays for the object    */
/* in question.)                                               */
/*                                                             */
/* if there was no such region, -1 is returned.                */
/***************************************************************/
int RWGGeometry::GetRegionByLabel(const char *Label)
{
  if (Label==0)
   return -1;
  
  for(int nr=0; nr<NumRegions; nr++)
   if ( !StrCaseCmp(Label, RegionLabels[nr]) ) 
    return nr;

  return -1;
}

/***************************************************************/
/* return the RWGSurface whose label is Label. if pns is       */
/* non-NULL on entry, then on return it is set to the index of */
/* the RWGSurface. If no corresponding RWGSurface was found,   */
/* the return value is NULL and *pns is set to -1.             */
/***************************************************************/
RWGSurface *RWGGeometry::GetSurfaceByLabel(const char *Label, int *pns)
{
  if (pns) *pns=-1;
 
  if (Label==0)
   return NULL;
  
  for(int ns=0; ns<NumSurfaces; ns++)
   if ( !StrCaseCmp(Label, Surfaces[ns]->Label) ) 
    { if (pns) *pns = ns;
      return Surfaces[ns];
    }
  return NULL;
}

/***************************************************************/
/* Apply the specified GTComplex to transform the geometry.    */
/* (Note that a 'GTComplex' is a list of GTransformations, each*/
/* of which is applied to one specific surface in the geometry.)*/
/***************************************************************/
void RWGGeometry::Transform(GTComplex *GTC)
{ 
  // assume that no objects will be modified by this operation
  memset(SurfaceMoved, 0, NumSurfaces*sizeof(int));

  // loop over the individual transformations in the complex
  for(unsigned ns=0; ns<GTC->GTs.size(); ns++)
   { 
     // find the surface corresponding to the label for this transformation
     RWGSurface *S = 0;
     char *Label = GTC->SurfaceLabels[ns];
     if (!Label)
      Warn("no OBJECT/SURFACE label specified for GTransformation %s (ignoring)",GTC->Tag);
     else 
     { S=GetSurfaceByLabel(Label);
       if (!S) Warn("GTransformation %s refers to unknown object/surface %s (ignoring)",GTC->Tag,Label);
     }

     // apply the transformation to that surface
     if (S)
      { S->Transform(GTC->GTs[ns]);
        SurfaceMoved[S->Index]=1;
      }
   }
}

/***************************************************************/
/* Undo transformations. ***************************************/
/***************************************************************/
void RWGGeometry::UnTransform()
{ 
  for(int ns=0; ns<NumSurfaces; ns++)
   Surfaces[ns]->UnTransform();
}

/***************************************************************/
/* Quick sanity check to make sure that a given list of        */
/* GTComplex structures actually makes sense for the given     */
/* geometry, which is to say that it doesn't request           */
/* transformations on any surfaces that don't exist in the     */
/* geometry.                                                   */
/* Returns 0 if the check passed, or an error message if not.  */
/***************************************************************/
char *RWGGeometry::CheckGTCList(GTCList GTCs)
{
  for(unsigned nt=0; nt<GTCs.size(); nt++)
   for (unsigned ns=0; ns<GTCs[nt]->GTs.size(); ns++)
    { char *Label = GTCs[nt]->SurfaceLabels[ns];
      if (!Label) 
       return vstrdup("no OBJECT/SURFACE specified for GTransformation %s",GTCs[nt]->Tag);
      if (!GetSurfaceByLabel(Label))
       return vstrdup("unknown OBJECT/SURFACE %s specified for GTransformation %s",Label,GTCs[nt]->Tag);
    }

  return 0;
}

/***************************************************************/
/* Return the dimension of the linear system. ******************/
/***************************************************************/
int RWGGeometry::GetDimension()
{ return TotalBFs; }

/***************************************************************/
/* Given a vector of surface-current coefficients obtained as  */
/* the solution to a scattering problem, extract the           */
/* coefficients of the electric and magnetic RWG functions     */
/* associated with edge #ne of surface #ns.                    */
/***************************************************************/
void RWGGeometry::GetKNCoefficients(HVector *KN, int ns, int ne,
                                    cdouble *pKAlpha, cdouble *pNAlpha)
{
  if (KN==0)
   { *pKAlpha=0.0;
     if (pNAlpha) *pNAlpha=0.0;
     return;
   }

  RWGSurface *S=Surfaces[ns];
  int Offset = BFIndexOffset[ns];
  cdouble KAlpha, NAlpha;

  if (S->IsPEC)
   { KAlpha = KN->GetEntry(Offset + ne);
     NAlpha = 0.0;
   }
  else
   { KAlpha =       KN->GetEntry(Offset + 2*ne + 0);
     NAlpha = -ZVAC*KN->GetEntry(Offset + 2*ne + 1);
   };
  *pKAlpha=KAlpha;
  if (pNAlpha) *pNAlpha=NAlpha;
}

/***************************************************************/
/* given an interior triangle edge specified by an index into  */
/* the overall list of internal edges in a geometry, determine */
/* the surface to which the edge belongs and the index of the  */
/* edge within that surface; also determine the index of the   */
/* electric surface-current coefficient for this edge within   */
/* the overall list of surface-current coefficients for the    */
/* entire geometry                                             */
/***************************************************************/
RWGSurface *RWGGeometry::ResolveEdge(int neFull, int *pns, int *pne, int *pKNIndex)
{
  int ns=0, NSm1=NumSurfaces - 1;
  while( (ns < NSm1) && (neFull >= EdgeIndexOffset[ns+1]) )
   ns++;

  int ne  = neFull - EdgeIndexOffset[ns];
  
  int Mult    = Surfaces[ns]->IsPEC ? 1 : 2;
  int KNIndex = BFIndexOffset[ns] + Mult*ne;

  if (pns) *pns=ns;
  if (pne) *pne=ne;
  if (pKNIndex) *pKNIndex=KNIndex;
  return Surfaces[ns];
}

int RWGGeometry::UnResolveEdge(int ns, int ne)
{ return EdgeIndexOffset[ns] + ne; }

/***************************************************************/
/* like the previous routine, except the input is an index into*/
/* the overall BEM system.                                     */
/***************************************************************/
RWGSurface *RWGGeometry::ResolveBF(int nbfFull, int *pns, int *pne, bool *pIsMagnetic)
{
  int ns=0, NSm1=NumSurfaces - 1;
  while( (ns < NSm1) && (nbfFull >= BFIndexOffset[ns+1]) )
   ns++;

  int nbf  = nbfFull - BFIndexOffset[ns];
  
  bool IsMagnetic=false;
  int ne=0;
  if (Surfaces[ns]->IsPEC)
   { ne = nbf;
     IsMagnetic=false;
   }
  else
   { ne = nbf/2;
     IsMagnetic = (nbf%2)==1;
   };

  if (pns) *pns=ns;
  if (pne) *pne=ne;
  if (pIsMagnetic) *pIsMagnetic=IsMagnetic;
  return Surfaces[ns];
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::SetLogLevel(int NewLogLevel)
 { LogLevel=NewLogLevel; }

void RWGGeometry::SetLogLevel(const char *Level)
{
  if (Level==0) return;

  if ( !strcasecmp(Level, "NONE") ) 
   LogLevel=SCUFF_NOLOGGING;
  else if ( !strcasecmp(Level, "TERSE") ) 
   { LogLevel=SCUFF_TERSELOGGING;
     Log("Setting log level to terse.");
   }
  else if ( !strcasecmp(Level, "VERBOSE") ) 
   { LogLevel=SCUFF_VERBOSELOGGING;
     Log("Setting log level to verbose.");
   }
  else if ( !strcasecmp(Level, "VERBOSE2") ) 
   { LogLevel=SCUFF_VERBOSE2;
     Log("Setting log level to verbose2.");
   }
  else
   Warn("unknown log level %s",Level);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::SetEpsMu(const char *Label, cdouble Eps, cdouble Mu)
{ 
  int nr=GetRegionByLabel(Label);
  if (nr==-1)
   Warn("unknown object %s specified in SetEpsMu() (ignoring)",Label);
  else
   RegionMPs[nr]->SetEpsMu(Eps, Mu);
} 

void RWGGeometry::SetEpsMu(cdouble Eps, cdouble Mu)
{ SetEpsMu(0, Eps, Mu); }

void RWGGeometry::SetEps(const char *Label, cdouble Eps)
{ SetEpsMu(Label, Eps, 1.0); }

void RWGGeometry::SetEps(cdouble Eps)
{ SetEpsMu(0, Eps, 1.0); }

/***************************************************************/
/* update the internally stored cache of epsilon and mu values */
/* for a given frequency                                       */
/***************************************************************/
void RWGGeometry::UpdateCachedEpsMuValues(cdouble Omega)
{
  /*--------------------------------------------------------------*/
  /*- update cached epsilon and mu values at this frequency if    */
  /*- necessary.                                                  */
  /*-                                                             */
  /*- Note (20120924): Previously I looked to see if              */
  /*- Omega==StoredOmega (i.e. the frequency is the same as it    */
  /*- was the last time we updated the EpsMu values), and if so   */
  /*- I bypassed the update. However, this doesn't work, because  */
  /*- one or more of the RegionMPs may have been zeroed or        */
  /*- unzeroed since the last call. So now I just do the updating */
  /*- in all cases no matter what.                                */
  /*--------------------------------------------------------------*/
  //if (Omega != StoredOmega )
  if (1)
   { StoredOmega=Omega;
     for(int nr=0; nr<NumRegions; nr++)
      RegionMPs[nr]->GetEpsMu(Omega, &(EpsTF[nr]), &(MuTF[nr]) );
   };
}

} // namespace scuff
