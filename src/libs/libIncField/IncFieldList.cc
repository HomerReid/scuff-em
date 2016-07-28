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
* IncFieldList.cc -- new functionality for parsing a user-specified
*                 -- text file to obtain a list of IncField structures
*
* homer reid      -- 5/2016
*/

#include <string.h>
#include <stdlib.h>

#include "libIncField.h"
#include "libhrutil.h"

#define IFTYPE_PW  0
#define IFTYPE_PS  1
#define IFTYPE_MPS 2
#define IFTYPE_GB  3

#define MAXSTR 512
#define MAXTOK 12

/***************************************************************/
/***************************************************************/
/***************************************************************/
IncField *ProcessIFTokens(char **Tokens, int NumTokens, char **pErrMsg)
{
 int Type=0;
 int NumParameters=0;
 if ( !strcasecmp(Tokens[0], "PW" ) )
  { Type=IFTYPE_PW;
    NumParameters=6;
  }
 else if (    !strcasecmp(Tokens[0], "PS" )
           || !strcasecmp(Tokens[0], "EPS" )
         )
  { Type=IFTYPE_PS;
    NumParameters=6;
  }
 else if ( !strcasecmp(Tokens[0], "MPS" ) )
  { Type=IFTYPE_MPS;
    NumParameters=6;
  }
 else if ( !strcasecmp(Tokens[0], "GB" ) )
  { Type=IFTYPE_GB;
    NumParameters=10;
  }
 else
  { *pErrMsg = vstrdup("unknown field type %s",Tokens[0]);
    return 0;
  };

 if ( NumTokens != NumParameters+1 )
  { *pErrMsg = vstrdup("incorrect number of parameters "
                       "(need %i for %s)",NumTokens,NumParameters,Tokens[0]);
     return 0;
  };

 cdouble Parameters[10];
 for(int nt=1; nt<NumTokens; nt++)
  if ( S2CD(Tokens[nt],&(Parameters[nt-1])) )
   { *pErrMsg = vstrdup("invalid parameter value %s",Tokens[nt]);
     return 0;
   };

 if (Type==IFTYPE_PW)
  { double nHat[3];
    cdouble E0[3];
    nHat[0] = real(Parameters[0]);
    nHat[1] = real(Parameters[1]);
    nHat[2] = real(Parameters[2]);
    E0[0]   = Parameters[3];
    E0[1]   = Parameters[4];
    E0[2]   = Parameters[5];
    return new PlaneWave(E0, nHat);
  }
 else if (Type==IFTYPE_PS || Type==IFTYPE_MPS)
  { double X0[3];
    cdouble P[3];
    X0[0]  = real(Parameters[0]);
    X0[1]  = real(Parameters[1]);
    X0[2]  = real(Parameters[2]);
    P[0]   = Parameters[3];
    P[1]   = Parameters[4];
    P[2]   = Parameters[5];
    int SourceType = LIF_ELECTRIC_DIPOLE;
    if (Type==IFTYPE_MPS) SourceType = LIF_MAGNETIC_DIPOLE;
    return new PointSource(X0, P, SourceType);
  }
 else if (Type==IFTYPE_GB)
  { double X0[3];
    double nHat[3];
    cdouble E0[3];
    X0[0]    = real(Parameters[0]);
    X0[1]    = real(Parameters[1]);
    X0[2]    = real(Parameters[2]);
    nHat[0]  = real(Parameters[3]);
    nHat[1]  = real(Parameters[4]);
    nHat[2]  = real(Parameters[5]);
    E0[0]    = Parameters[6];
    E0[1]    = Parameters[7];
    E0[2]    = Parameters[8];
    double W = real(Parameters[9]);
    return new GaussianBeam(X0, nHat, E0, W);
  };

 return 0; // never get here
 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
IncFieldList *AddIncFieldToList(IncField *IF, char *Label, IncFieldList *IFList)
{
  if (IFList==0)
   { IFList = (IncFieldList *)mallocEC(sizeof(IncFieldList));
     IFList->IFs    = 0;
     IFList->NumIFs = 0;
     IFList->Labels = 0;
   };
 
  int N = IFList->NumIFs;
  IFList->IFs = (IncField **)realloc(IFList->IFs, (N+1)*sizeof(IncField *));
  IFList->Labels = (char **)realloc(IFList->Labels, (N+1)*sizeof(char *));
  IFList->IFs[N]    = IF;
  IFList->Labels[N] = strdup(Label ? Label : "default");
  IFList->NumIFs   += 1;

  return IFList;

}

/***************************************************************/
/* quick hack to convert e.g. cos(30) to 0.5 to make it a      */
/* little easier to write IFList files describing e.g.         */
/* angular sweeps of polarization or incident direction        */
/***************************************************************/
void ProcessSinCos(char **Tokens, int NumTokens)
{
  for(int nt=0; nt<NumTokens; nt++)
   { 
     int L = strlen(Tokens[nt]);
     if (L<5) continue;
     double ThetaDegrees;
     if ( (    !strncasecmp(Tokens[nt],"cos(",4)
            || !strncasecmp(Tokens[nt],"sin(",4)
          )
          &&  1==sscanf(Tokens[nt]+4,"%le",&ThetaDegrees)
        )
      snprintf(Tokens[nt],L,"%e", tolower(Tokens[nt][0])=='c' 
                                  ? cos(ThetaDegrees*M_PI/180.0)
                                  : sin(ThetaDegrees*M_PI/180.0) 
              );
   };
}


/***************************************************************/
/***************************************************************/
/***************************************************************/
IncFieldList *ReadIncFieldList(char *FileName)
{
  IncFieldList *IFList=0;

  FILE *f=fopen(FileName,"r");
  if (f==0)
   ErrExit("could not open file %s",FileName);

  char Line[MAXSTR];
  int NumTokens;
  char *Tokens[MAXTOK];
  char *ErrMsg=0;
  int LineNum=0;
  IncField *CompoundIF=0;
  char *CompoundIFLabel=0;
  while(fgets(Line, MAXSTR, f))
   { 
     // read line, break it up into tokens, skip blank lines and comments
     LineNum++;
     NumTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

     // convert strings like sin(45.0) to 0.7071....
     ProcessSinCos(Tokens, NumTokens);
 
     // switch off depending on whether or not we are in the middle
     // of a compound field declaration
     if (CompoundIFLabel)
      { 
        if ( !StrCaseCmp(Tokens[0], "END") )
         {
           if (NumTokens!=1)
            ErrExit("%s:%i: junk on line following keyword END",FileName,LineNum);
           if (CompoundIF==0)
            ErrExit("%s:%i: empty compound field definition",FileName,LineNum);
 
           IFList=AddIncFieldToList(CompoundIF, CompoundIFLabel, IFList);
           free(CompoundIFLabel);
           CompoundIFLabel=0;
           CompoundIF=0;
         }
        else
         { IncField *IFTerm=ProcessIFTokens(Tokens, NumTokens, &ErrMsg);
           if (IFTerm==0)
            ErrExit("%s:%i: %s",FileName,LineNum,ErrMsg ? ErrMsg : "syntax error");
           IFTerm->Next=CompoundIF;
           CompoundIF=IFTerm;
         };
      }
     else
      { 
        if (NumTokens==1)
         { if (CompoundIFLabel)
            ErrExit("%s:%i: expected END",FileName,LineNum);
           CompoundIFLabel = strdup(Tokens[0]);
         }
        else
         { IncField *IFTerm=ProcessIFTokens(Tokens+1, NumTokens-1, &ErrMsg);
           if (IFTerm==0)
            ErrExit("%s:%i: %s",FileName,LineNum,ErrMsg ? ErrMsg : "syntax error");
           IFList=AddIncFieldToList(IFTerm, Tokens[0], IFList);
         };
      };
         
   }; // while(fgets(Line, MAXSTR, f))
  fclose(f);
 
  if (CompoundIFLabel)
   ErrExit("%s:%i: expected END",FileName,LineNum);
 
  Log("Read %i incident-field descriptions from file %s.",IFList->NumIFs,FileName);
  return IFList;
 
}
