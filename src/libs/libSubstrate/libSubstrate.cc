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

//using namespace scuff;

/***************************************************************/
/* if the ErrMsg field of the class instance is nonzero on     */
/* return, something went wrong                                */
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

  NumLayers=0;
  MPMedium=0;
  MPLayer=0;
  zLayer=0;
  zGP=HUGE_VAL;

#define MAXSTR 1000
  char Line[MAXSTR];
  int LineNum=0;
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
     /*- all lines must be of the form   ----------------------------*/
     /*-   zValue  MaterialName          ----------------------------*/
     /*- or                              ----------------------------*/
     /*-   MEDIUM  MaterialName          ----------------------------*/
     /*- or                              ----------------------------*/
     /*-   zValue  GROUNDPLANE           ----------------------------*/
     /*--------------------------------------------------------------*/
     if ( NumTokens!=2 )
      { ErrMsg=vstrdup("%s:%i syntax error",FileName,LineNum);
        return;
      };

     if ( !strcasecmp(Tokens[0],"MEDIUM") )
      { MPMedium = new MatProp(Tokens[1]);
        if (MPMedium->ErrMsg)
          { ErrMsg=vstrdup("%s:%i: %s",FileName,LineNum,MPMedium->ErrMsg);
            return;
          }
        Log("Setting upper half-space medium to %s.",MPMedium->Name);
        continue;
      };

     double z;
     if ( 1!=sscanf(Tokens[0],"%le",&z) )
      { ErrMsg=vstrdup("%s:%i bad z-value %s",FileName,LineNum,Tokens[0]);
        return;
      };

     if ( !strcasecmp(Tokens[1],"GROUNDPLANE") )
      { zGP = z;
        Log(" Ground plane at z=%e.",zGP);
      }
     else
      { MatProp *MP = new MatProp(Tokens[1]);
        if (MP->ErrMsg)
         { ErrMsg=vstrdup("%s:%i: %s",FileName,LineNum,MP->ErrMsg);
           return;
         };
/*
        if (NumLayers==MAXLAYER)
         { ErrMsg=vstrdup("%s:%i: too many layers",FileName,LineNum);
           return;
         };
*/
        NumLayers++;
        MPLayer=(MatProp **)reallocEC(MPLayer,NumLayers*sizeof(MatProp *));
         zLayer=(double  *)reallocEC(zLayer, NumLayers*sizeof(double));
        MPLayer[NumLayers-1]=MP;
         zLayer[NumLayers-1]=z;
        Log(" Layer #%i: %s at z=%e.",NumLayers,MP->Name,z);
      };
   };
  fclose(f);

  /*--------------------------------------------------------------*/
  /*- sanity check that ground plane lies below all layers       -*/
  /*--------------------------------------------------------------*/
  if (zGP!=HUGE_VAL)
   for(int n=0; n<NumLayers; n++)
    if ( zLayer[n] < zGP )
     { ErrMsg=vstrdup("%s: ground plane must lie below all dielectric layers",FileName);
       return;
     };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  EpsMedium  = 1.0;
  EpsLayer   = (cdouble *)mallocEC(NumLayers*sizeof(cdouble));
  MuLayer    = (cdouble *)mallocEC(NumLayers*sizeof(cdouble));
  OmegaCache = -1.0;

  qMaxEval  = 10000;
  qAbsTol   = 1.0e-12;
  qRelTol   = 1.0e-6;
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
  ErrMsg=0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
LayeredSubstrate::~LayeredSubstrate()
{
  if (MPMedium)
   delete MPMedium;
  for(int n=0; n<NumLayers; n++)
   delete MPLayer[n];
  free(MPLayer);
  free(EpsLayer);
  free(MuLayer);
  free(zLayer);
  if (I1D)
   delete I1D;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::UpdateCachedEpsMu(cdouble Omega)
{
  if ( EqualFloat(OmegaCache, Omega) )
   return;
  OmegaCache=Omega;
  if (MPMedium)
   EpsMedium = MPMedium->GetEps(Omega);
  for(int n=0; n<NumLayers; n++)
   MPLayer[n]->GetEpsMu(Omega, EpsLayer+n, MuLayer+n);
}
