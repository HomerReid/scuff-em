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
 * libMatProp.cc -- implementation of libMatProp
 *
 * homer reid    -- 12/2009
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libMDInterp.h>

#include "libMatProp.h"
#include "cmatheval.h"

#define MAXCONSTANTS 25
#define MAXSTR       100

/***************************************************************/
/* initialization of the static FreqUnit class variable ********/
/***************************************************************/
double MatProp::FreqUnit=2.99792458e14; // c / (1 micron)

/***************************************************************/
/* MPDP is an internally maintained list of pointers to        */
/* MatProp instances that have been created to date.           */
/*                                                             */
/* 20131220 Note that we now make a *copy* of the MatProp      */
/*          instance to which MP points, and store internally  */
/*          a pointer to this copy, rather than storing MP     */
/*          itself. Thanks to Owen Miller for suggesting this  */
/*          improvement.                                       */
/*                                                             */
/* 20140525 Now we also maintain a reference counter for each  */
/*          material in the MPDB, so that we know when it's    */
/*          safe to delete an entry from the table.            */
/***************************************************************/
static MatProp **MPDB=0;
static int MPDBLength=0;

void AddMPToMatPropDataBase(MatProp *MP)
{ 
  MPDB=(MatProp **)realloc(MPDB, (MPDBLength+1)*sizeof(MatProp *));
  MPDB[MPDBLength] = new MatProp(MP);
  MPDBLength++;
}

MatProp *FindMPInMatPropDataBase(const char *MaterialName)
{ 
  int nmp;
  for (nmp=0; nmp<MPDBLength; nmp++)
   if (!StrCaseCmp(MPDB[nmp]->Name, MaterialName) )
    return MPDB[nmp];
  return 0;
}

/***************************************************************/
/* constructor entry points ************************************/
/***************************************************************/
MatProp::MatProp()
 { Type=MP_PEC; 
   Zeroed=0;
   Name=strdupEC("PEC");
 }

MatProp::MatProp(int pType)
 { Type=pType;
   Zeroed=0;
   Name=strdupEC("VACUUM");
 }
 
MatProp::MatProp(const char *MaterialName)
 { InitMatProp(MaterialName, 0); }

MatProp::MatProp(const char *MaterialName, const char *MatPropFileName)
 { InitMatProp(MaterialName, MatPropFileName); }

// copy constructor
MatProp::MatProp(MatProp *MP)
 { 
   Type   = MP->Type;
   Zeroed = 0;
   Name   = strdup(MP->Name);
   Eps    = MP->Eps;
   Mu     = MP->Mu;

   /***************************************************************/
   /* 20140525 new handling of copy constructor for interpolated  */
   /*          materials                                          */
   /***************************************************************/
#if 0
   InterpReal = InterpImag = 0;
   if (MP->InterpReal || MP->InterpImag)
    ErrExit("%s:%i: MatProp copy constructor not implemented for interpolated materials",__FILE__,__LINE__);
#endif
   InterpReal = MP->InterpReal ? new Interp1D(MP->InterpReal) : 0;
   InterpImag = MP->InterpImag ? new Interp1D(MP->InterpImag) : 0;
   OwnsInterpolators = true;

   EpsExpression=MuExpression = 0;
   if (MP->EpsExpression)
    { 
      EpsExpression=cevaluator_create(cevaluator_get_string(MP->EpsExpression));
      cevaluator_set_var_index(EpsExpression, "w", 0); // w is thread-safe

      char **Vars;
      int NumVars;
      cevaluator_get_variables(MP->EpsExpression, &Vars, &NumVars);
      for(int nv=0; nv<NumVars; nv++)
       cevaluator_set_var(EpsExpression, Vars[nv], 
                          cevaluator_get_var(MP->EpsExpression, Vars[nv]) 
                          );
    };

   if (MP->MuExpression)
    { MuExpression=cevaluator_create(cevaluator_get_string(MP->MuExpression));
      cevaluator_set_var_index(MuExpression, "w", 0); // w is thread-safe

      char **Vars;
      int NumVars;
      cevaluator_get_variables(MP->MuExpression, &Vars, &NumVars);
      for(int nv=0; nv<NumVars; nv++)
       cevaluator_set_var(MuExpression, Vars[nv], 
                          cevaluator_get_var(MP->MuExpression, Vars[nv]) 
                         );
      
    };
   OwnsExpressions= true;

 }

/***************************************************************/
/* the actual body of the constructor **************************/
/***************************************************************/
void MatProp::InitMatProp(const char *MaterialName, const char *MatPropFileName)
{ 
  const char *p;

  /***************************************************************/
  /* set some defaults *******************************************/
  /***************************************************************/
  Eps=1.0;
  Mu=1.0;
  Zeroed=0;
  EpsExpression = MuExpression = NULL;
  InterpReal = InterpImag = NULL;
  OwnsExpressions = OwnsInterpolators = false;

  /* assume things will go OK */
  ErrMsg=0;

  /* set *p = first non-whitespace character in MaterialName */
  p=MaterialName;
  while ( p && *p && isspace(*p) )
   p++;

  /***************************************************************/
  /* now switch off based on the type of material specified.     */
  /***************************************************************/

  /*--------------------------------------------------------------*/
  /*- PEC keyword ------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if ( p==0 || *p==0 || !StrCaseCmp(p,"PEC") )
   { Type=MP_PEC;
     Name=strdupEC("PEC");
     return;
   };

  /*--------------------------------------------------------------*/
  /*- VACUUM keyword ---------------------------------------------*/
  /*--------------------------------------------------------------*/
  if ( !StrCaseCmp(p,"VACUUM") )
   { Type=MP_VACUUM;
     Name=strdupEC(p);
     return;
   };

  /*--------------------------------------------------------------*/
  /*- CONST_EPS_xxx                                              -*/
  /*- CONST_MU_xxx                                               -*/
  /*- CONST_EPS_xxx_MU_xxx                                       -*/
  /*- CONST_EPS_xxx_MU_XXX_SIGMA_xxx                             -*/
  /*-  ... etc                                                   -*/
  /*--------------------------------------------------------------*/
  if ( !strncasecmp(p,"CONST_",6) )
   { 
     char ConstString[100];
     char *Tokens[6];
     int nt, NumTokens, nConv;
     double ER, EI;
     char c1, c2;
     double Sigma=0.0;
     Name=strdupEC(p);

     strncpy(ConstString, p+6, 100);
     NumTokens=Tokenize(ConstString, Tokens, 6, "_");
     for(nt=0; nt<NumTokens; nt++)
      { 
        if ( !StrCaseCmp(Tokens[nt],"MU") )
         { if ( ++nt == NumTokens ) 
            { ErrMsg=strdupEC("no value specified for MU");
              return;
            };
           nConv=sscanf(Tokens[nt],"%le%c%le%c",&ER,&c1,&EI,&c2);
           if ( nConv == 1 )
            Mu = cdouble(ER,0.0);
           else if ( nConv==4 && c1=='+' && tolower(c2)=='i' )
            Mu = cdouble(ER,EI);
           else if ( nConv==4 && c1=='-' && tolower(c2)=='i' )
            Mu = cdouble(ER, -1.0*EI);
           else
            { ErrMsg=strdupEC("bad constant value specified for MU");
              return;
            };
         }
        else if ( !StrCaseCmp(Tokens[nt],"EPS") )
         { if ( ++nt == NumTokens )
            { ErrMsg=strdupEC("no value specified for EPS");
              return;
            };
           nConv=sscanf(Tokens[nt],"%le%c%le%c",&ER,&c1,&EI,&c2);
           if ( nConv == 1 )
            Eps = cdouble(ER,0.0);
           else if ( nConv==4 && c1=='+' && tolower(c2)=='i' )
            Eps = cdouble(ER,EI);
           else if ( nConv==4 && c1=='-' && tolower(c2)=='i' )
            Eps = cdouble(ER, -1.0*EI);
           else
            { ErrMsg=strdupEC("bad constant value specified for EPS");
              return;
            };
         }
        else if ( !StrCaseCmp(Tokens[nt],"SIGMA") )
         { if ( ++nt == NumTokens )
            { ErrMsg=strdupEC("no value specified for SIGMA");
              return;
            };
           nConv=sscanf(Tokens[nt],"%le",&Sigma);
           if ( nConv != 1 )
            { ErrMsg=strdupEC("bad constant value specified for Sigma");
              return;
            };
         }
        else
         { ErrMsg=strdupEC("invalid token");
           return;
         };
      };

     if (Sigma==0.0)
      Type=MP_CONSTANT;
     else
      { Type=MP_PARSED;
        char Expr[100];
        snprintf(Expr,100,"%e + %e*i*2.99792458e14/w",real(Eps),Sigma);
        EpsExpression = cevaluator_create(Expr);
	if (!EpsExpression)
         { ErrMsg=strdupEC("invalid SIGMA value");
	   return;
	 }
	cevaluator_set_var_index(EpsExpression, "w", 0);
        if (Mu!=1.0)
         { snprintf(Expr,100,"%e+%ei",real(Mu),imag(Mu));
           MuExpression = cevaluator_create(Expr);
	   if (!MuExpression)
            { ErrMsg=strdupEC("invalid MU value");
	      return;
	    };
	   cevaluator_set_var_index(MuExpression, "w", 0);
         };
	OwnsExpressions=true;
      };
     
     return;

   }; // if ( !strncasecmp(p,"CONST_",6) )

  /*--------------------------------------------------------------*/
  /*- FILE_xxx keyword -------------------------------------------*/
  /*--------------------------------------------------------------*/
  if ( !strncasecmp(p,"FILE_",5) )
   { 
     Type=MP_INTERP;
     Name=strdupEC(p);

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     MatProp *MP;
     if ( (MP=FindMPInMatPropDataBase(p)) )
      { InterpReal=MP->InterpReal; 
	InterpImag=MP->InterpImag; 
	OwnsInterpolators=false;
      }
     else
      { 
        ReadInterpolationTable(p+5);
	OwnsInterpolators=true;
        AddMPToMatPropDataBase(this);
      };
     return;
   };

  /*--------------------------------------------------------------*/
  /*- if none of the above, try to look up the given material     */
  /*- in a matprop.dat file. we look in several places:           */
  /*-  a) the MatPropFileName specified by the caller of this     */
  /*-     routine, if any                                         */
  /*-  b) the file specified by the SCUFF_MATPROPFILE environment */
  /*-     variable                                                */
  /*-  c) matprop.dat in the current working directory            */
  /*-  d) matprop.dat in the directory specified by the           */
  /*-     SCUFF_MATPROP_PATH environment variable                 */
  /*-  e) .matprop.dat in the user's home directory               */
  /*--------------------------------------------------------------*/
    Type=MP_PARSED;
    Name=strdupEC(p);
    MatProp *MP=FindMPInMatPropDataBase(p);
    if ( MP )
     { 
       EpsExpression = MP->EpsExpression;
       MuExpression = MP->MuExpression;
       OwnsExpressions = false;
     }
    else
     { 
       CreateUserDefinedMaterial(MatPropFileName, p);
       OwnsExpressions = true;
       AddMPToMatPropDataBase(this);
     };
}

/***************************************************************/
/* attempt to read (epsilon,mu) vs frequency data from a file  */
/*                                                             */
/* each line in the file (except blank lines and lines that    */
/* start with '#', which are skipped) must contain either two  */
/* or three complex numbers a+bi, which are interpreted as     */
/* follows:                                                    */
/*                                                             */
/* FREQ   EPS(w)   MU(w)                                       */
/*                                                             */
/* where the MU column is optional (defaults to 1 if absent).  */
/* Each frequency value must be either purely real or purely   */
/* imaginary (but the same table can contain both).            */
/***************************************************************/
void MatProp::ReadInterpolationTable(const char *FileName)
{
  /***************************************************************/
  /* attempt to read in the data file as a big matrix            */
  /***************************************************************/
  char *Dir=0;
  FILE *f=fopenPath(getenv("SCUFF_MATPROP_PATH"), FileName, "r", &Dir);
  if (!f) 
   ErrExit("could not open file %s",FileName);
  fclose(f);
  char FullFileName[MAXSTR];
  snprintf(FullFileName,MAXSTR,"%s/%s",Dir,FileName);
  Log("Found material data file %s in directory %s.",FileName,Dir);
  HMatrix *Data = new HMatrix(FullFileName, LHM_TEXT, "--strict");
  if (Data->ErrMsg)
   ErrExit(Data->ErrMsg);

  if (Data->NC != 2 && Data->NC != 3)
    ErrExit("interpolation data must have either two or three columns");

  Data->Sort(0); // Interp1D expects data to be sorted

  /***************************************************************/
  /* extract from the matrix the XPoints and YPoints arrays      */
  /* passed to the Interp1D constructor                          */
  /***************************************************************/
  int nr, NR=Data->NR;
  double *XPoints=(double *)mallocEC(NR*sizeof(double));
  double *YPoints=(double *)mallocEC(4*NR*sizeof(double));
  if (XPoints==0 || YPoints==0)
   ErrExit("insufficient memory to read interpolation data file %s",FileName);
  
  /* count number of purely real frequencies */
  int nreal = 0;
  for(nr=0; nr<NR; nr++)
    nreal += imag(Data->GetEntry(nr, 0)) == 0.0;

  /* Put the real frequencies at the beginning of the arrays
     and the imaginary frequencies at the end. */
  int ir = 0, ii = 0; // count of real and imaginary freqs so far
  for(nr=0; nr<NR; nr++) {
    cdouble w = Data->GetEntry(nr, 0);
    int r=0;
    if (imag(w) == 0.0)
      XPoints[r = ir++] = real(w);
    else if (real(w) == 0.0)
      XPoints[r = nreal + ii++] = imag(w);
    else
      ErrExit("frequencies in %s must be purely real or imaginary\n",FileName);

    YPoints[4*r + 0] = real(Data->GetEntry(nr, 1));
    YPoints[4*r + 1] = imag(Data->GetEntry(nr, 1));
    if (Data->NC == 3) {
      YPoints[4*r + 2] = real(Data->GetEntry(nr, 2));
      YPoints[4*r + 3] = imag(Data->GetEntry(nr, 2));
    }
    else { // default Mu = 1
      YPoints[4*r + 2] = 1.0; 
      YPoints[4*r + 3] = 0.0;
    }
  }

  /***************************************************************/
  /* create the interpolators ************************************/
  /***************************************************************/
  InterpReal = (nreal == 0) ? NULL : new Interp1D(XPoints, YPoints, nreal, 4);

  // (strictly speaking, we should be able to omit the imaginary parts
  //  along the imaginary-frequency axis since they should always be zero)
  InterpImag = (nreal == NR) ? NULL : new Interp1D(XPoints + nreal, YPoints + 4*nreal, NR-nreal, 4);

  /***************************************************************/
  /* deallocate temporary storage ********************************/
  /***************************************************************/
  free(XPoints);
  free(YPoints);
  delete Data;
  
} 

/***************************************************************/
/* destructor **************************************************/
/***************************************************************/
MatProp::~MatProp()
{
  free(Name);

  if (Type==MP_INTERP && OwnsInterpolators)
   { if (InterpReal) delete InterpReal;
     if (InterpImag) delete InterpImag;
   }
  else if (Type==MP_PARSED && OwnsExpressions)
   { 
      if (EpsExpression) cevaluator_destroy(EpsExpression);
      if (MuExpression) cevaluator_destroy(MuExpression);
   };

}  

/***************************************************************/
/***************************************************************/
/***************************************************************/
void MatProp::SetFreqUnit(double NewFreqUnit)
{ FreqUnit = NewFreqUnit; }

void MatProp::SetLengthUnit(double NewLengthUnit)
{ FreqUnit = 299792458.0 / NewLengthUnit; }

/***************************************************************/
/* return 1 if the material is a perfect electrical conductor  */
/***************************************************************/
int MatProp::IsPEC()
{ 
  if (Type==MP_PEC) 
   return 1;
  return 0;
}

/***************************************************************/
/* set constant eps and mu                                     */
/***************************************************************/
void MatProp::SetEpsMu(cdouble NewEps, cdouble NewMu)
{
  Type=MP_CONSTANT;
  Eps=NewEps;
  Mu=NewMu;
}

/***************************************************************/
/* get eps and mu at a given frequency *************************/
/***************************************************************/
void MatProp::GetEpsMu(cdouble Omega, cdouble *pEps, cdouble *pMu)
{ 
  cdouble EpsRV; // 'epsilon return value' 
  cdouble MuRV;   // 'mu return value' 

  if ( Zeroed )
   { EpsRV=0.0;
     MuRV=0.0;
   }
  else if ( Type==MP_VACUUM ) 
   { 
     EpsRV=1.0;
     MuRV=1.0; 
   }
  else if ( Type==MP_PEC ) // this should never happen in practice
   {
     EpsRV = -1.0e9;
     MuRV=1.0;
   }
  else if ( Type==MP_CONSTANT )
   { EpsRV=Eps;
     MuRV=Mu;
   }
  else if ( Type==MP_INTERP )
   { 
     double Data[4];
     if(imag(Omega)==0.0)
      { 
        if (InterpReal==0)
         ErrExit("no real-frequency data given in material data file %s",Name+5);
	if(InterpReal->Evaluate(real(Omega)*FreqUnit, Data)) {
	  EpsRV=cdouble(Data[0], Data[1]);
	  MuRV=cdouble(Data[2], Data[3]);
	} else {
	  EpsRV = MuRV = cdouble(NAN, NAN);
	  Warn("Setting Eps, Mu to NaN: frequency %s outside the interpolation range in data file %s.", z2s(Omega), Name+5);
	}
      }
     else if (real(Omega)==0.0)
      {
        if (InterpImag==0)
         ErrExit("no imaginary-frequency data given in material data file %s); setting Eps, Mu to NaN.",Name+5);
	if(InterpImag->Evaluate(imag(Omega)*FreqUnit, Data)) {
	  EpsRV=cdouble(Data[0], 0.0);
	  MuRV=cdouble(Data[2], 0.0);
	} else {
	  EpsRV = MuRV = cdouble(NAN, NAN);
	  Warn("Setting Eps, Mu to NaN: frequency %s outside the interpolation range in data file %s.", z2s(Omega), Name+5);
	}
      }
     else
      ErrExit("interpolated frequencies must lie on the real or imaginary axis");
   }
  else if ( Type==MP_PARSED ) 
   { 
       GetEpsMu_Parsed(Omega, &EpsRV, &MuRV);
   } // if (Type==... )
  else
    ErrExit("Invalid material property type");

  if (pEps) *pEps=EpsRV;
  if (pMu) *pMu=MuRV;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/

cdouble MatProp::GetEps(cdouble Omega) {
  cdouble EpsRV;
  GetEpsMu(Omega, &EpsRV, 0);
  return EpsRV;
}

cdouble MatProp::GetMu(cdouble Omega) {
  cdouble MuRV;
  GetEpsMu(Omega, 0, &MuRV);
  return MuRV;
}

cdouble MatProp::GetRefractiveIndex(cdouble Omega, cdouble *pZRel)
{
  cdouble MyEps, MyMu;
  GetEpsMu(Omega, &MyEps, &MyMu);
  if (pZRel)
   *pZRel = (MyEps==0.0) ? 0.0 : sqrt(MyMu/MyEps);
  return sqrt(MyEps*MyMu);
}
