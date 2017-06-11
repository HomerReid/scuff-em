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
 * ExternalFields.cc -- some types of external electrostatic
 *                   -- field for use in defining excitations
 *                   -- (external stimulus) for SCUFF-EM
 *                   -- electrostatics calculations
 *
 * homer reid        -- 5/2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>
#include <cmatheval/cmatheval.h>
#include <libSpherical.h>
#include <libTriInt.h>

#include "libscuff.h"
#include "SSSolver.h"
#include "StaticSubstrate.h"

namespace scuff {

/***************************************************************/
/***************************************************************/
/***************************************************************/
void UserStaticField(double *x, void *UserData, double PhiE[4])
{
  UserSFData *Data = (UserSFData *)UserData;

  static const char *VariableNames[3] = { "x", "y", "z" };
  cdouble VariableValues[3];
  VariableValues[0] = cdouble(x[0], 0.0 );
  VariableValues[1] = cdouble(x[1], 0.0 );
  VariableValues[2] = cdouble(x[2], 0.0 );

  memset(PhiE, 0, 4*sizeof(double));
  if ( Data==0 || Data->PhiEvaluator==0 ) return;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  PhiE[0] = real( cevaluator_evaluate(Data->PhiEvaluator, 3, 
                                      const_cast<char **>(VariableNames), 
                                      VariableValues) );

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int Mu=0; Mu<3; Mu++)
   { if (Data->EEvaluator[Mu])
      { PhiE[1+Mu] = real( cevaluator_evaluate(Data->EEvaluator[Mu], 3, 
                                               const_cast<char **>(VariableNames), 
                                               VariableValues) 
                         );
      }
     else
      { double Delta = 1.0e-4*abs(VariableValues[Mu]);
        if (Delta==0.0) Delta=1.0e-4;
        VariableValues[Mu] += Delta;
        double PhiPlus = real( cevaluator_evaluate(Data->PhiEvaluator, 3,
                                                   const_cast<char **>(VariableNames), 
                                                   VariableValues) 
                             );
        VariableValues[Mu] -= 2.0*Delta;
        double PhiMinus= real( cevaluator_evaluate(Data->PhiEvaluator, 3,
                                                   const_cast<char **>(VariableNames), 
                                                   VariableValues) 
                             );
        VariableValues[Mu] += Delta;
        PhiE[1+Mu] = (PhiMinus - PhiPlus) / (2.0*Delta);
      };
   };
   
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ConstantStaticField(double *x, void *UserData, double PhiE[4])
{ 
  ConstantSFData *Data = (ConstantSFData *)UserData;
  double *E0 = Data->E0;

  memset(PhiE, 0, 4.0*sizeof(double));
  PhiE[0] = -(E0[0]*x[0] + E0[1]*x[1] + E0[2]*x[2]);
  PhiE[1] = E0[0];
  PhiE[2] = E0[1];
  PhiE[3] = E0[2];
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void MonopoleStaticField(double *x, void *UserData, double PhiE[4])
{ 
  MonopoleSFData *Data = (MonopoleSFData *)UserData;

  double *x0 = Data->x0;
  double Q   = Data->Q;

  double R[3];
  R[0]=x[0]-x0[0];
  R[1]=x[1]-x0[1];
  R[2]=x[2]-x0[2];
  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  if (r2==0.0)
   { memset(PhiE, 0, 4*sizeof(double));
     return;
   };
  double r=sqrt(r2);

  PhiE[0] = Q / (4.0*M_PI*r);
  PhiE[1] = PhiE[0] * R[0]/r2;
  PhiE[2] = PhiE[1] * R[1]/r2;
  PhiE[3] = PhiE[2] * R[2]/r2;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void DipoleStaticField(double *x, void *UserData, double PhiE[4])
{ 
  DipoleSFData *Data = (DipoleSFData *)UserData;

  double *x0 = Data->x0;
  double *P  = Data->P;

  double R[3];
  R[0]=x[0]-x0[0];
  R[1]=x[1]-x0[1];
  R[2]=x[2]-x0[2];
  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  if (r2==0.0)
   { memset(PhiE, 0, 4*sizeof(double));
     return;
   };
  double r=sqrt(r2);
  double r3=r*r2, r5=r3*r2;

  double PdR = P[0]*R[0] + P[1]*R[1] + P[2]*R[2];

  PhiE[0] = PdR / (4.0*M_PI*r3);
  PhiE[1] = (3.0*PdR*R[0] - P[0]*r2)/(4.0*M_PI*r5);
  PhiE[2] = (3.0*PdR*R[1] - P[1]*r2)/(4.0*M_PI*r5);
  PhiE[3] = (3.0*PdR*R[2] - P[2]*r2)/(4.0*M_PI*r5);
}

/***************************************************************/
/* (real-valued) spherical harmonic incident field *************/
/***************************************************************/
void SphericalStaticField(double *x, void *UserData, double PhiE[4])
{
  SphericalSFData *Data= (SphericalSFData *)UserData;
  int l = Data->l;
  int m = Data->m;

  double r, Theta, Phi;
  CoordinateC2S(x, &r, &Theta, &Phi);

  PhiE[0] = pow(r,l)*GetRealYlm(l,m,Theta,Phi);

  // FIXME
  PhiE[1] = PhiE[2] = PhiE[3] = 0.0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddStaticField(StaticExcitation *SE, StaticField SF, void *SFData)
{
  int NumSFs         = SE->NumSFs++;
  SE->SFs            = (StaticField *)realloc(SE->SFs, (NumSFs+1)*sizeof(StaticField));
  SE->SFData         = (void **)realloc(SE->SFData, (NumSFs+1)*sizeof(void *));
  SE->SFs[NumSFs]    = SF;
  SE->SFData[NumSFs] = SFData;
}


/***************************************************************/
/***************************************************************/
/***************************************************************/
void EvalStaticField(StaticExcitation *SE, double *x, double *PhiE)
{
  memset(PhiE, 0, 4*sizeof(double));
  for(int n=0; n<SE->NumSFs; n++)
   { 
     StaticField SF = SE->SFs[n];
     void *SFData = SE->SFData[n];
     double DeltaPhiE[4];
     SF(x, SFData, DeltaPhiE);
     VecPlusEquals(PhiE, 1.0, DeltaPhiE, 4);
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXSTR 1000
#define MAXTOK 10
StaticExcitation *ParseExcitationSection(SSSolver *SSS,
                                         FILE *f, char *Label,
                                         char *FileName, int *LineNum)
{
  RWGGeometry *G=SSS->G;
  int NS=G->NumSurfaces;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  StaticExcitation *SE = (StaticExcitation *)mallocEC(sizeof(*SE));
  SE->Label=strdup(Label);
  SE->Potentials=(double *)mallocEC(NS * sizeof(double));
  memset(SE->Potentials, 0, NS*sizeof(double));
  SE->SFs=0;
  SE->SFData=0;
  SE->NumSFs=0;

  /***************************************************************/
  /* read and parse lines from the file one by one ***************/
  /***************************************************************/
  char Line[MAXSTR];
  while( fgets(Line, MAXSTR, f) )
   { 
     (*LineNum)++;
     char *Tokens[MAXTOK];
     int NumTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
     if (!strcasecmp(Tokens[0],"ENDEXCITATION"))
      {
        return SE;
      }
    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
     else if (!strcasecmp(Tokens[0],"CONSTANT_FIELD"))
      { ConstantSFData *Data = (ConstantSFData *)mallocEC(sizeof *Data);
        if (NumTokens!=4)
         ErrExit("%s:%i: wrong number of tokens",FileName,*LineNum);
        if (    1!=sscanf(Tokens[1],"%le",Data->E0+0)
             || 1!=sscanf(Tokens[2],"%le",Data->E0+1)
             || 1!=sscanf(Tokens[3],"%le",Data->E0+2)
           )
         ErrExit("%s:%i: invalid E0 specification",FileName,*LineNum);
        AddStaticField(SE, ConstantStaticField, (void *)Data);
        Log("Excitation %s: added constant field (E0={%e,%e,%e})",Label,Data->E0[0],Data->E0[1],Data->E0[2]);
      } 
    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
     else if (!strcasecmp(Tokens[0],"MONOPOLE"))
      { MonopoleSFData *Data = (MonopoleSFData *)mallocEC(sizeof *Data);
        if (NumTokens!=5)
         ErrExit("%s:%i: wrong number of tokens",FileName,*LineNum);
        if (    1!=sscanf(Tokens[1],"%le",Data->x0+0)
             || 1!=sscanf(Tokens[2],"%le",Data->x0+1)
             || 1!=sscanf(Tokens[3],"%le",Data->x0+2)
             || 1!=sscanf(Tokens[4],"%le",&(Data->Q))
           )
         ErrExit("%s:%i: invalid monopole specification",FileName,*LineNum);
        AddStaticField(SE, MonopoleStaticField, (void *)Data);
        Log("Excitation %s: added monopole field (x0={%e,%e,%e}, Q=%e)",Label,Data->x0[0],Data->x0[1],Data->x0[2],Data->Q);
      } 
    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
     else if (!strcasecmp(Tokens[0],"DIPOLE"))
      { DipoleSFData *Data = (DipoleSFData *)mallocEC(sizeof *Data);
        if (NumTokens!=7)
         ErrExit("%s:%i: wrong number of tokens",FileName,*LineNum);
        if (    1!=sscanf(Tokens[1],"%le",Data->x0+0)
             || 1!=sscanf(Tokens[2],"%le",Data->x0+1)
             || 1!=sscanf(Tokens[3],"%le",Data->x0+2)
             || 1!=sscanf(Tokens[4],"%le",Data->P+0)
             || 1!=sscanf(Tokens[5],"%le",Data->P+1)
             || 1!=sscanf(Tokens[6],"%le",Data->P+2)
           )
         ErrExit("%s:%i: invalid dipole specification",FileName,*LineNum);
        AddStaticField(SE, DipoleStaticField, (void *)Data);
        Log("Excitation %s: added dipole field (x0={%e,%e,%e}, P={%e,%e,%e)",Label,Data->x0[0],Data->x0[1],Data->x0[2],Data->P[0],Data->P[1],Data->P[2]);
      } 
    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
     else if (!strcasecmp(Tokens[0],"PHI"))
      { 
        if (NumTokens<1)
         ErrExit("%s:%i: invalid PHI specification",FileName,*LineNum);

        UserSFData *Data = (UserSFData *)mallocEC(sizeof *Data);
        char str[MAXSTR]="";
        for(int nt=1; nt<NumTokens; nt++)
         snprintf(str,MAXSTR,"%s %s",str,Tokens[nt]);
        Data->PhiEvaluator=cevaluator_create(str);
        AddStaticField(SE, UserStaticField, (void *)Data);
        Log("Excitation %s: added user-defined field Phi(x,y,z)=%s",Label,str);
      } 
    /*--------------------------------------------------------------*/
    /*- try to interpret the line as CONDUCTOR_NAME VOLTAGE         */
    /*--------------------------------------------------------------*/
     else
      { double V;
        if ( NumTokens!=2 || 1!=sscanf(Tokens[1],"%lf",&V) )
         ErrExit("%s:%i: syntax error",FileName,*LineNum);

        int ns;
        RWGSurface *S=G->GetSurfaceByLabel(Tokens[0],&ns);
        if (!S)
         ErrExit("%s:%i: unknown conductor %s",FileName,*LineNum,Tokens[0]);
        if (!(S->IsPEC))
         ErrExit("%s:%i: attempt to assign potential to non-PEC surface %s",Tokens[0]);

        SE->Potentials[ns]=V;
        Log("Excitation %s: setting conductor %s to V=%e V.",Label,Tokens[0],V);
      };
   };

  //never get here
  ErrExit("%s:%i: internal error",__FILE__,__LINE__);
  return 0;

  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXSTR 1000
#define MAXTOK 10
StaticExcitation **ReadExcitationFile(SSSolver *SSS, char *FileName,
                                      int *pNumExcitations)
{
  int NumExcitations=0;
  StaticExcitation **SEList=0;

  /***************************************************************/
  /* read and parse lines from the file one by one ***************/
  /***************************************************************/
  FILE *f;
  if ( !FileName )
   ErrExit("no --ExcitationFile specified");
  if ( !(f=fopen(FileName,"r") ) )
   ErrExit("could not open file %s",FileName);
  char Line[MAXSTR];
  int LineNum=0;
  while( fgets(Line, MAXSTR, f) )
   { 
     LineNum++;
     char *Tokens[MAXTOK];
     int NumTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
     if (!strcasecmp(Tokens[0],"EXCITATION"))
      {
        if (NumTokens!=2) 
         ErrExit("%s:%i: syntax error",FileName,LineNum);
        SEList = (StaticExcitation **)realloc(SEList, (NumExcitations+1)*sizeof(*SEList));
        SEList[NumExcitations++] = ParseExcitationSection(SSS, f, Tokens[1], FileName, &LineNum);
      }
     else
      ErrExit("%s:%i: syntax error",FileName,LineNum);
   };
  fclose(f);

  *pNumExcitations=NumExcitations;
  return SEList;
}

/***************************************************************/
/* Parse a "potential" file to read a list of conductor        */
/* potentials.                                                 */
/* The file should look something like                         */
/*                                                             */
/*  BottomSurface 3.2                                          */
/*  TopSurface    1.4                                          */
/*                                                             */
/* where the first string is the name of a PEC object/surface  */
/* specified in the .scuffgeo file, and the second string is   */
/* the potential in volts at which that object will be held    */
/* in the scuff-static calculation.                            */
/***************************************************************/
void ParsePotentialFile(RWGGeometry *G, char *PotFile, double *Potentials)
{ 
  FILE *f=fopen(PotFile,"r");
  if (!f)
   ErrExit("could not open file %s",PotFile);

  int LineNum=0;
  char Line[100];
  while (fgets(Line,100,f))
   { 
     LineNum++;
     char *Tokens[3];
     int NumTokens=Tokenize(Line, Tokens, 3);
     if (NumTokens==0 || Tokens[0][0]=='#') 
      continue; // skip blank lines and comments

     if (NumTokens!=2) 
      ErrExit("%s:%i: syntax error",PotFile,LineNum);

     int ns;
     if (G->GetSurfaceByLabel(Tokens[0],&ns)==0)
      ErrExit("%s:%i: unknown surface",PotFile,LineNum,Tokens[0]);
     if ( ! (G->Surfaces[ns]->IsPEC) )
      ErrExit("%s:%i: attempt to assign potential to non-PEC surface %s",Tokens[0]);

     double V;
     if (1!=sscanf(Tokens[1],"%le",&V))
      ErrExit("%s:%i: invalid potential specification",PotFile,LineNum);

     Potentials[ns]=V;
     Log("Setting potential of surface %s to %e volts.\n",Tokens[0],V);
   }; 

  fclose(f);
}

/***************************************************************/
/* stub routine for creating a list of StaticExcitations       */
/* consisting of just one excitation.                          */
/* provided for backward compatibility with older options      */
/* for specifying potentials and external fields.              */
/***************************************************************/
StaticExcitation **CreateSimpleSEList(RWGGeometry *G,
                                      char *PotFile,
                                      char *ConstField,
                                      int nMonopoles, double *Monopoles,
                                      int nDipoles, double *Dipoles,
                                      char *PhiExt)
{
  int NS = G->NumSurfaces;

  StaticExcitation **SEList 
   = (StaticExcitation **)mallocEC(sizeof(**SEList));
  StaticExcitation *SE = SEList[0] = (StaticExcitation *)mallocEC(sizeof(*SE));
  SE->Label            = 0;
  SE->SFs              = 0;
  SE->SFData           = 0;
  SE->NumSFs           = 0;
  SE->Potentials       = (double *)mallocEC(NS*sizeof(double));
  memset(SE->Potentials, 0, NS*sizeof(double));

  /***************************************************************/  
  /***************************************************************/  
  /***************************************************************/  
  if (PotFile)
   ParsePotentialFile(G, PotFile, SE->Potentials);

  /***************************************************************/  
  /***************************************************************/  
  /***************************************************************/  
  if (ConstField)
   { int Direction = tolower(ConstField[0]) - 'x';
     if (Direction<0 || Direction>2) 
      ErrExit("invalid --ConstField specification");
     ConstantSFData *Data = (ConstantSFData *)mallocEC(sizeof(*Data));
     Data->E0[0]=Data->E0[1]=Data->E0[2]=0.0;
     Data->E0[Direction]=1.0;
     AddStaticField(SE, ConstantStaticField, (void *)Data);
   };

  /***************************************************************/  
  /***************************************************************/  
  /***************************************************************/  
  for(int n=0; n<nMonopoles; n++)
   { MonopoleSFData *Data = (MonopoleSFData *)mallocEC(sizeof(*Data));
     Data->x0[0] = Monopoles[4*n+0];
     Data->x0[1] = Monopoles[4*n+1];
     Data->x0[2] = Monopoles[4*n+2];
     Data->Q     = Monopoles[4*n+3];
     AddStaticField(SE, MonopoleStaticField, (void *)Data);
   };

  /***************************************************************/  
  /***************************************************************/  
  /***************************************************************/  
  for(int n=0; n<nDipoles; n++)
   { DipoleSFData *Data = (DipoleSFData *)mallocEC(sizeof(*Data));
     Data->x0[0] = Dipoles[6*n+0];
     Data->x0[1] = Dipoles[6*n+1];
     Data->x0[2] = Dipoles[6*n+2];
     Data->P[0]  = Dipoles[6*n+3];
     Data->P[1]  = Dipoles[6*n+4];
     Data->P[2]  = Dipoles[6*n+5];
     AddStaticField(SE, DipoleStaticField, (void *)Data);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (PhiExt)
   { 
     UserSFData *Data = (UserSFData *)mallocEC(sizeof(*Data));
     Data->PhiEvaluator=cevaluator_create(PhiExt);
     AddStaticField(SE, UserStaticField, (void *)Data);
   };
 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  return SEList;
  
}

} // namespace scuff
