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
 * libSubstrate.cc -- implicit handling of multilayered material substrates
 *                 -- this file: stuff common to static and full-wave cases
 *
 * homer reid   -- 3/2017
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include "libhrutil.h"
#include "libMDInterp.h"
#include "libMatProp.h"
#include "libSGJC.h"
#include "libSubstrate.h"

/***************************************************************/
/* constructor entry point 1: construct from a .substrate file */
/***************************************************************/
LayeredSubstrate::LayeredSubstrate(const char *FileName)
{
  char *Dir=0;
  FILE *f=fopenPath(getenv("SCUFF_SUBSTRATE_PATH"), FileName,"r",&Dir);
  if (!f)
   { ErrMsg=vstrdup("could not open file %s",FileName);
     return; 
   };
  Log("Reading substrate definition from %s/%s.",Dir ? Dir : ".",FileName);

  Initialize(f,FileName,0);
  fclose(f);
}

/***************************************************************/
/* constructor entry point 2: construct starting from the      */
/* second line of a SUBSTRATE ... ENDSUBSTRATE section in a    */
/* previously opened file                                      */
/***************************************************************/
LayeredSubstrate::LayeredSubstrate(FILE *f, int *pLineNum)
{
  Initialize(f,0,pLineNum);
}

/***************************************************************/
/* main body of constructor.                                   */
/* if the ErrMsg field of the class instance is nonzero on     */
/* return, something went wrong.                               */
/***************************************************************/
void LayeredSubstrate::Initialize(FILE *f, const char *FileName, int *pLineNum)
{
  /*--------------------------------------------------------------*/
  /*- initialize class fields ------------------------------------*/
  /*--------------------------------------------------------------*/
  NumInterfaces=0;
  NumLayers=1;
  MPLayer=(MatProp **)mallocEC(NumLayers*sizeof(MatProp *));
  MPLayer[0]=new MatProp("VACUUM");
  zInterface=0;
  zGP=-1.0*HUGE_VAL;
  
  // keep track of whether we're in a .substrate file or in a
  // SUBSTRATE...ENDSUBSTRATE section of a .scuffgeo file
  bool InSubstrateFile = (FileName!=0);

  /*--------------------------------------------------------------*/
  /*- read and parse lines from the file one at a time -----------*/
  /*--------------------------------------------------------------*/
  #define MAXSTR 1000
  char Line[MAXSTR];
  int LineNum = (pLineNum ? *pLineNum : 0);
  bool GotEndSubstrate=false;
  char FileLine[100]="";
  while( fgets(Line,MAXSTR,f) )
   { 
     /*--------------------------------------------------------------*/
     /*- skip blank lines and constants -----------------------------*/
     /*--------------------------------------------------------------*/
     LineNum++;
     int NumTokens;
     char *Tokens[2];
     int Length=strlen(Line);
     if (Length==0) continue;
     Line[(Length-1)]=0; // remove trailing carriage return
     NumTokens=Tokenize(Line, Tokens, 2);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

     /*--------------------------------------------------------------*/
     /*- all lines must be of one of the following forms   ----------*/
     /*-   zValue  MaterialName          ----------------------------*/
     /*-   zValue  GROUNDPLANE           ----------------------------*/
     /*-   MEDIUM  MaterialName          ----------------------------*/
     /*-   ENDSUBSTRATE                  ----------------------------*/
     /*--------------------------------------------------------------*/
     if ( !strcasecmp(Tokens[0],"ENDSUBSTRATE") )
      { GotEndSubstrate=true;
        break;
      };

     if (InSubstrateFile) snprintf(FileLine,100,"%s:%i: ",FileName,LineNum);

     if ( NumTokens!=2 )
      { ErrMsg=vstrdup("%ssyntax error",FileLine);
        return;
      };

     if ( !strcasecmp(Tokens[0],"MEDIUM") )
      { if (!InSubstrateFile)
         { ErrMsg=vstrdup("MEDIUM keyword forbidden in SUBSTRATE...ENDSUBSTRATE sections");
           return;
         };
        MPLayer[0] = new MatProp(Tokens[1]);
        if (MPLayer[0]->ErrMsg)
          { ErrMsg=vstrdup("%s%s",FileLine,MPLayer[0]->ErrMsg);
            return;
          }
        Log("Setting upper half-space medium to %s.",MPLayer[0]->Name);
        continue;
      };

     double z;
     if ( 1!=sscanf(Tokens[0],"%le",&z) )
      { ErrMsg=vstrdup("%sbad z-value %s",FileLine,Tokens[0]);
        return;
      };

     if ( !strcasecmp(Tokens[1],"GROUNDPLANE") )
      { zGP = z;
        Log(" Ground plane at z=%e.",zGP);
      }
     else
      { 
        if (NumInterfaces>0 && z>zInterface[NumInterfaces-1])
         { ErrMsg=vstrdup("%sz coordinate lies above previous layer",FileLine);
           return;
         };

        MatProp *MP = new MatProp(Tokens[1]);
        if (MP->ErrMsg)
         { ErrMsg=vstrdup("%s%s",FileLine,MP->ErrMsg);
           return;
         };
        NumInterfaces++;
        NumLayers++;
        MPLayer=(MatProp **)reallocEC(MPLayer,NumLayers*sizeof(MatProp *));
        zInterface=(double  *)reallocEC(zInterface, NumInterfaces*sizeof(double));
        MPLayer[NumLayers-1]=MP;
        zInterface[NumInterfaces-1]=z;
        Log(" Layer #%i: %s at z=%e.",NumInterfaces,MP->Name,z);
      };
   };

  /*--------------------------------------------------------------*/
  /* for the constructor entry point in which we are reading a    */
  /* SUBSTRATE...ENDSUBSTRATE section in a .scuffgeo file, check  */
  /* for the closing keyword and error out if absent              */
  /*--------------------------------------------------------------*/
  if (InSubstrateFile && GotEndSubstrate)
   Warn("%s:%i: ENDSUBSTRATE is not needed in .substrate files",FileName,LineNum);
  else if (!InSubstrateFile && !GotEndSubstrate)
   { ErrMsg=vstrdup("expected ENDSUBSTRATE before end of file");
     return;
   };
  if (pLineNum) 
   *pLineNum = LineNum;
  if (InSubstrateFile) snprintf(FileLine,100,"%s: ",FileName);

  /*--------------------------------------------------------------*/
  /*- sanity check                                               -*/
  /*--------------------------------------------------------------*/
  if (fabs(zGP)!=HUGE_VAL && zGP>zInterface[NumInterfaces-1])
   { ErrMsg=vstrdup("%sground plane must lie below all dielectric layers",FileLine);
     return;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  EpsLayer   = (cdouble *)mallocEC(NumLayers*sizeof(cdouble));
  MuLayer    = (cdouble *)mallocEC(NumLayers*sizeof(cdouble));
  OmegaCache = -1.0;

  qMaxEval  = 10000;
  qAbsTol   = 1.0e-8;
  qRelTol   = 1.0e-4;
  PPIOrder  = 9;
  PhiEOrder = 9;
  char *s;
  if ((s=getenv("SCUFF_SUBSTRATE_QMAXEVAL")))
   sscanf(s,"%i",&qMaxEval);
  if ((s=getenv("SCUFF_SUBSTRATE_QABSTOL")))
   sscanf(s,"%le",&qAbsTol);
  if ((s=getenv("SCUFF_SUBSTRATE_QRELTOL")))
   sscanf(s,"%le",&qRelTol);
  if ((s=getenv("SCUFF_SUBSTRATE_PPIORDER")))
   sscanf(s,"%i",&PPIOrder);
  if ((s=getenv("SCUFF_SUBSTRATE_PHIEORDER")))
   sscanf(s,"%i",&PhiEOrder);
  I1D=0;
  I1DRhoMin=HUGE_VAL;
  I1DRhoMax=0;

  ForceFreeSpace=false;
  StaticLimit=false;
  LayerOnly=-1;

  ErrMsg=0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
LayeredSubstrate::~LayeredSubstrate()
{
  for(int n=0; n<=NumInterfaces; n++)
   delete MPLayer[n];
  free(MPLayer);
  free(EpsLayer);
  free(MuLayer);
  free(zInterface);
  if (I1D)
   delete I1D;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::Describe(FILE *f)
{ 
  if (f==0) f=stdout;
  double zAbove=1.0/0.0, zBelow;
  int n=0;
  for(; n<NumInterfaces; n++)
   { zBelow = zInterface[n];
     printf("Region %2i (%-20s): [%10g < z < %-10g]\n",n,MPLayer[n]->Name,zAbove,zBelow);
     zAbove=zBelow;
   };
  printf("Region %2i (%-20s): [%10g < z < %-10g]\n",n,MPLayer[n]->Name,zAbove,zGP);
  if ( fabs(zGP) != HUGE_VAL) 
   printf("Ground plane at z=%e.\n\n",zGP);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::UpdateCachedEpsMu(cdouble Omega)
{
  if ( EqualFloat(OmegaCache, Omega) )
   return;
  OmegaCache=Omega;
  for(int n=0; n<=NumInterfaces; n++)
   MPLayer[n]->GetEpsMu(Omega, EpsLayer+n, MuLayer+n);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int LayeredSubstrate::GetRegionIndex(double z)
{ for(int ni=0; ni<NumInterfaces; ni++)
   if (z>zInterface[ni]) return ni;
  return NumInterfaces;
}
