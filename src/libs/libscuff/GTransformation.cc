/*
 * GTransformation.cc -- a very simple mechanism for handling geometric 
 *                    -- transformations (displacements and rotations)
 *                        
 * homer reid         -- 11/2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>

#include "GTransformation.h"

#define MAXSTR 1000
#define MAXTOK 50

namespace scuff {

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void ConstructRotationMatrix(double *ZHat, double Theta, double M[9]);

/***************************************************************/
/* create the identity transformation                          */
/***************************************************************/
GTransformation *CreateGTransformation()
 { 
   GTransformation *NGT=(GTransformation *)malloc(sizeof(*NGT));
   memset(NGT->DX, 0, 3*sizeof(double));
   memset(NGT->M, 0, 9*sizeof(double));
   NGT->M[0]=NGT->M[4]=NGT->M[8]=1.0;
   NGT->RotationIsNonTrivial=0;
   return NGT;
 }

/***************************************************************/
/* CreateOrAugmentGTransformation is a routine with multiple   */
/* entry points that takes an existing GTransformation (which  */
/* may be NULL) and augments it by tacking on the effect of a  */
/* second transformation. the return value is a pointer to     */
/* a GTransformation representing the composite transformation.*/
/***************************************************************/

/***************************************************************/
/* a transformation that displaces through vector DX           */
/***************************************************************/
GTransformation *CreateOrAugmentGTransformation(GTransformation *GT, double *DX)
{
  if (GT==0)
   GT=CreateGTransformation();

  GT->DX[0] += DX[0];
  GT->DX[1] += DX[1];
  GT->DX[2] += DX[2];

  return GT;
}

/***************************************************************/
/* a transformation that rotates through Theta degrees         */
/* (DEGREES, NOT RADIANS) about an axis that passes through    */
/* the origin and the point with coordinates ZHat[0..2]        */
/***************************************************************/
GTransformation *CreateOrAugmentGTransformation(GTransformation *GT, 
                                                double *ZHat, double Theta)
{
  if (GT==0) 
   GT=CreateGTransformation();

  double MP[9], NewDX[3], NewM[9];
  int i, j, k;

  ConstructRotationMatrix(ZHat, Theta, MP);

  for(i=0; i<3; i++)
   for(NewDX[i]=0.0, j=0; j<3; j++)
    NewDX[i] += MP[ 3*i + j ]*(GT->DX[j]);
  memcpy(GT->DX, NewDX, 3*sizeof(double));

  for(i=0; i<3; i++)
   for(j=0; j<3; j++)
    for(NewM[ 3*i + j ]=0.0, k=0; k<3; k++)
     NewM[ 3*i + j ] += MP[ 3*i + k ]*GT->M[ 3*k + j ];

  memcpy(GT->M,NewM,9*sizeof(double));

  GT->RotationIsNonTrivial=1;
  return GT;

}

/***************************************************************/
/* a transformation described by a character string.           */
/*                                                             */
/* the string should be of the form                            */
/*                                                             */
/*  DISPLACED xx yy zz                                         */
/*                                                             */
/* or                                                          */
/*                                                             */
/*  ROTATED tt ABOUT xx yy zz                                  */
/*                                                             */
/* in this case, the caller should look at the value returned  */
/* for the parameter ErrMsg; if ErrMsg==0 on return, the       */
/* string was successfully parsed and interpreted as a         */
/* GTransformation, and otherwise ErrMsg points to an error    */
/* message.                                                    */
/***************************************************************/
GTransformation *CreateOrAugmentGTransformation(GTransformation *GT, 
                                                char *TransformString,
                                                char **ErrMsg)
{
  char Line[MAXSTR];
  char *Tokens[MAXTOK];
  int NumTokens;
  
  strncpy(Line,TransformString,MAXSTR);
  NumTokens=Tokenize(Line, Tokens, MAXTOK);
  return CreateOrAugmentGTransformation(GT, Tokens, NumTokens, ErrMsg, 0);
}

GTransformation *CreateOrAugmentGTransformation(GTransformation *GT, 
                                                char **Tokens, int NumTokens, 
                                                char **ErrMsg, int *TokensConsumed)
{
 
  double Theta, ZHat[3], DX[3];

  if (NumTokens==0)
   { if (*ErrMsg) *ErrMsg = strdup("no tranformation specified");
     if (TokensConsumed) *TokensConsumed=0;
     return GT;
   };
   
  if ( !strcasecmp(Tokens[0],"DISP") || !strcasecmp(Tokens[0],"DISPLACED") )
   { 
     /*--------------------------------------------------------------*/
     /*-- DISPLACED xx yy zz ----------------------------------------*/
     /*--------------------------------------------------------------*/
     if ( NumTokens<4 )
      { if (ErrMsg) *ErrMsg=vstrdup("too few values specified for %s ",Tokens[0]);
        return 0;
      };

     if (    1!=sscanf(Tokens[1],"%le",DX+0)
          || 1!=sscanf(Tokens[2],"%le",DX+1)
          || 1!=sscanf(Tokens[3],"%le",DX+2)
        )
      { if (ErrMsg) *ErrMsg=vstrdup("bad value specified for %s",Tokens[0]);
        return 0;
      };

     if (ErrMsg) *ErrMsg=0;
     if (TokensConsumed) *TokensConsumed=4;
     return CreateOrAugmentGTransformation(GT, DX);
   }
  else if ( !strcasecmp(Tokens[0],"ROT") || !strcasecmp(Tokens[0],"ROTATED") )
   { 
     /*--------------------------------------------------------------*/
     /*-- ROTATED tt ABOUT xx yy zz ---------------------------------*/
     /*--------------------------------------------------------------*/
     if ( NumTokens<6 )
      { if (ErrMsg) *ErrMsg=vstrdup("too few values specified for %s ",Tokens[0]);
        return 0;
      };

     if ( strcasecmp(Tokens[2],"ABOUT") )
      {  if (ErrMsg) *ErrMsg=vstrdup("invalid syntax for %s statement",Tokens[0]);
         return 0;
      };

     if (    1!=sscanf(Tokens[1],"%le",&Theta)
          || 1!=sscanf(Tokens[3],"%le",ZHat+0)
          || 1!=sscanf(Tokens[4],"%le",ZHat+1)
          || 1!=sscanf(Tokens[5],"%le",ZHat+2)
        )
      { if (ErrMsg) *ErrMsg=vstrdup("bad value specified for %s",Tokens[0]);
        return 0;
      };

     if (ErrMsg) *ErrMsg=0;
     if (TokensConsumed) *TokensConsumed=6;
     return CreateOrAugmentGTransformation(GT, ZHat, Theta);
   }
  else
   { 
     if (ErrMsg) *ErrMsg=vstrdup("unknown keyword %s",Tokens[0]); 
     if (TokensConsumed) *TokensConsumed=0;
     return 0;
   };

} 

/***************************************************************/
/* GT -> DeltaGT*GT ********************************************/
/***************************************************************/
GTransformation *CreateOrAugmentGTransformation(GTransformation *GT, 
                                                GTransformation *DeltaGT)
{ 
  int i, j, k;
  double NewM[9], NewDX[3];

  if (GT==0)
   GT=CreateGTransformation();
 
  if ( DeltaGT->RotationIsNonTrivial )
   { 
     for(i=0; i<3; i++)
      for(j=0; j<3; j++)
       for(NewM[ 3*i + j ]=0.0, k=0; k<3; k++)
        NewM[ 3*i + j ] += DeltaGT->M[ 3*i + k ] * GT->M[ 3*k + j ];

     memcpy(GT->M, NewM, 9*sizeof(double));

     for(i=0; i<3; i++)
      for(NewDX[i]=0.0, j=0; j<3; j++)
       NewDX[i] += DeltaGT->M[ 3*i + j ]*GT->DX[j];

     memcpy(GT->DX, NewDX, 3*sizeof(double));

     GT->RotationIsNonTrivial=1;

   };

  GT->DX[0] += DeltaGT->DX[0];
  GT->DX[1] += DeltaGT->DX[1];
  GT->DX[2] += DeltaGT->DX[2];

  return GT;
   
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ResetGTransformation(GTransformation *GT)
{ 
  memset(GT->DX, 0, 3*sizeof(double));
  memset(GT->M,  0, 9*sizeof(double));
  GT->M[0]=GT->M[4]=GT->M[8]=1.0;
  GT->RotationIsNonTrivial = 0;

}

/***************************************************************/
/* apply the transformation (in-place) to a list of NX points  */
/* X[3*nx+0, 3*nx+1, 3*nx+2] = cartesian coords of point #nx   */
/***************************************************************/
void ApplyGTransformation(GTransformation *GT, double *X, int NX)
{ 
  if (GT==0) 
   return;

  int nx; 
  int i, j;
  double XP[3];

  if ( GT->RotationIsNonTrivial )
   { 
     // in this case we do the translation and rotation together
     for(nx=0; nx<NX; nx++)
      { memcpy(XP, GT->DX, 3*sizeof(double));
        for(i=0; i<3; i++) 
         for(j=0; j<3; j++)      
          XP[i] += GT->M[ 3*i + j ] * X[3*nx+j];
        memcpy(X+3*nx, XP, 3*sizeof(double));
      };
   }
  else
   { for(nx=0; nx<NX; nx++)
      { X[3*nx+0] += GT->DX[0];
        X[3*nx+1] += GT->DX[1];
        X[3*nx+2] += GT->DX[2];
      };
   };
}

void ApplyGTransformation(GTransformation *GT, double *X)
{ ApplyGTransformation(GT, X, 1); }

/***************************************************************/
/* like the above, but operate out-of-place.                   */
/***************************************************************/
void ApplyGTransformation(GTransformation *GT, double *X, double *XP, int NX)
{ 
  if (GT==0) 
   { memcpy(XP, X, 3*NX*sizeof(double));
     return;
   };

  int nx; 
  int i, j;
  if ( GT->RotationIsNonTrivial )
   { 
     for(nx=0; nx<NX; nx++)
      { memcpy(XP+3*nx, GT->DX, 3*sizeof(double));
        for(i=0; i<3; i++) 
         for(j=0; j<3; j++)      
          XP[3*nx+i] += GT->M[ 3*i + j ] * X[3*nx+j];
      };
   }
  else
   { for(nx=0; nx<NX; nx++)
      { XP[3*nx+0] = X[3*nx+0] + GT->DX[0];
        XP[3*nx+1] = X[3*nx+1] + GT->DX[1];
        XP[3*nx+2] = X[3*nx+2] + GT->DX[2];
      };
   };
}

void ApplyGTransformation(GTransformation *GT, double *X, double *XP)
 { ApplyGTransformation(GT, X, XP, 1); }

/***************************************************************/
/***************************************************************/
/***************************************************************/
void UnApplyGTransformation(GTransformation *GT, double *X, int NX)
{
  int i, j, nx;
  double XP[3];

  for(nx=0; nx<NX; nx++)
   { X[3*nx+0] -= GT->DX[0];
     X[3*nx+1] -= GT->DX[1];
     X[3*nx+2] -= GT->DX[2];
   };

  if ( GT->RotationIsNonTrivial )
   { for(nx=0; nx<NX; nx++)
      { 
        for(i=0; i<3; i++)
         for(XP[i]=0.0, j=0; j<3; j++)
          XP[i] += GT->M[ 3*j + i ]*X[3*nx+j];  // note matrix transpose

        memcpy(X+3*nx, XP, 3*sizeof(double));

      };
   };
}

/***************************************************************/
/* Construct the 3x3 matrix that represents a rotation of      */
/* Theta degrees (note we interpret Theta in DEGREES, NOT      */
/* RADIANS!) about the axis specified by cartesian coordinates */
/* ZHat.                                                       */
/* Algorithm:                                                  */
/*  1. Construct matrix M1 that rotates Z axis into alignment  */
/*     with ZHat.                                              */
/*  2. Construct matrix M2 that rotates through Theta about    */ 
/*     Z axis.                                                 */
/*  3. Construct matrix M=M1^{-1}*M2*M1=M1^T*M2*M1.            */  
/* Matrices are stored in row-major order.                     */  
/***************************************************************/
static void ConstructRotationMatrix(double *ZHat, double Theta, double M[9])
{ 
  int Mu, Nu, Rho;
  double ct, st, cp, sp, CT, ST;
  double M2M1[9], M2[9], M1[9];

  // first normalize ZHat
  double nZHat=sqrt(ZHat[0]*ZHat[0] + ZHat[1]*ZHat[1] + ZHat[2]*ZHat[2]);
  double NZHat[3];
  NZHat[0]=ZHat[0] / nZHat;
  NZHat[1]=ZHat[1] / nZHat;
  NZHat[2]=ZHat[2] / nZHat;
 
  /* construct M1 */
  ct=NZHat[2];
  st=sqrt(1.0-ct*ct);
  cp= ( st < 1.0e-8 ) ? 1.0 : NZHat[0] / st;
  sp= ( st < 1.0e-8 ) ? 0.0 : NZHat[1] / st;
  M1[ 3*0 + 0 ]=ct*cp;  M1[ 3*0 + 1 ]=ct*sp;   M1[ 3*0 + 2 ]=-st;
  M1[ 3*1 + 0 ]=-sp;    M1[ 3*1 + 1 ]=cp;      M1[ 3*1 + 2 ]=0.0;
  M1[ 3*2 + 0 ]=st*cp;  M1[ 3*2 + 1 ]=st*sp;   M1[ 3*2 + 2 ]=ct;

  /* construct M2 */
  CT=cos(Theta*M_PI/180.0);
  ST=sin(Theta*M_PI/180.0);
  M2[ 3*0 + 0 ]=CT;     M2[ 3*0 + 1 ]=-ST;     M2[ 3*0 + 2 ]=0.0;
  M2[ 3*1 + 0 ]=ST;     M2[ 3*1 + 1 ]=CT;      M2[ 3*1 + 2 ]=0.0;
  M2[ 3*2 + 0 ]=0.0;    M2[ 3*2 + 1 ]=0.0;     M2[ 3*2 + 2 ]=1.0;

  /* M2M1 <- M2*M1 */
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    for(M2M1[ 3*Mu + Nu ]=0.0, Rho=0; Rho<3; Rho++)
     M2M1[ 3*Mu + Nu ] += M2[ 3*Mu + Rho ] * M1[ 3*Rho + Nu ];

  /* M <- M1^T * M2M1 */
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    for(M[ 3*Mu + Nu ]=0.0, Rho=0; Rho<3; Rho++)
     M[ 3*Mu + Nu ] += M1[ 3*Rho + Mu ]*M2M1[ 3*Rho + Nu ];
}

/********************************************************************/
/* this is a helper function for the ReadTransFile() routine        */
/* below; it attempts to parse a TRANSFORMATION...ENDTRANSFORMATION */
/* section in a scuff-EM .trans file. (the file read pointer is     */
/* assumed to point to the line following TRANSFORMATION ...)       */
/* if successful, 0 is returned and *pGTC points on return to a     */
/* newly allocated GTComplex for the specified transformation.      */
/* if unsuccessful, an error message is returned.                   */
/********************************************************************/
char *ParseTRANSFORMATIONSection(char *Tag, FILE *f, int *pLineNum, GTComplex **pGTC)
{
  /***************************************************************/
  /* initialize a bare GTComplex *********************************/
  /***************************************************************/
  GTComplex *GTC=(GTComplex *)malloc(sizeof(GTComplex));
  GTC->Tag=strdup(Tag);
  GTC->NumObjectsAffected=0;
  GTC->ObjectLabel=0;
  GTC->GT=0;

  GTransformation *CurrentGT=0;

  /***************************************************************/
  /* parse each line in the TRANSFORMATION...ENDTRANSFORMATION   */
  /* section.                                                    */
  /***************************************************************/
  char Line[MAXSTR];
  int NumTokens, TokensConsumed;
  char *Tokens[MAXTOK];
  char *ErrMsg;
  while(fgets(Line, MAXSTR, f))
   { 
     // read line, break it up into tokens, skip blank lines and comments
     (*pLineNum)++;
     NumTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

     // switch off to handle the various tokens that may be encountered
     // in a TRANSFORMATION...ENDTRANSFORMATION section
     if ( !strcasecmp(Tokens[0],"OBJECT") )
      { 
        /*--------------------------------------------------------------*/
        /*-- OBJECT MyObject--------------------------------------------*/
        /*--------------------------------------------------------------*/
        if (NumTokens!=2) 
         return strdup("OBJECT keyword requires one argument");

        // increment the list of objects that will be affected by this
        // complex of transformations; save the label of the affected  
        // object and initialize the corresponding GTransformation to the
        // identity transformation
        int noa=GTC->NumObjectsAffected;
        GTC->ObjectLabel=(char **)realloc(GTC->ObjectLabel, (noa+1)*sizeof(char *));
        GTC->ObjectLabel[noa]=strdup(Tokens[1]);
        GTC->GT=(GTransformation **)realloc(GTC->GT, (noa+1)*sizeof(GTransformation *));
        GTC->GT[noa]=CreateGTransformation();
        GTC->NumObjectsAffected=noa+1;
        
        // subsequent DISPLACEMENTs / ROTATIONs will augment this GTransformation
        CurrentGT=GTC->GT[noa];
      }
     else if ( !strcasecmp(Tokens[0],"ENDTRANSFORMATION") )
      { 
        /*--------------------------------------------------------------*/
        /*-- ENDTRANSFORMATION  ----------------------------------------*/
        /*--------------------------------------------------------------*/
        *pGTC=GTC;
        return 0;
      }
     else 
      { /*--------------------------------------------------------------*/
        /*-- try to process the line as DISPLACED ... or ROTATED ... ---*/
        /*--------------------------------------------------------------*/
        CurrentGT=CreateOrAugmentGTransformation(CurrentGT,
                                                 Tokens, NumTokens,
                                                 &ErrMsg, &TokensConsumed);
        if (ErrMsg)
         return ErrMsg;
        if (TokensConsumed!=NumTokens)
         return strdup("junk at end of line");
      };

   }; // while(fgets(Line, MAXSTR, f))

  // if we made it here, the file ended before the TRANSFORMATION
  // section was properly terminated 
  return strdup("unexpected end of file");

}

/***************************************************************/
/* this routine reads a scuff-EM transformation (.trans) file  */
/* and returns an array of GTComplex structures.               */
/***************************************************************/
GTComplex **ReadTransFile(char *FileName, int *NumGTComplices)
{
  FILE *f=fopen(FileName,"r");
  if (f==0)
   ErrExit("could not open file %s",FileName);

  GTComplex *GTC, **GTCArray=0;
  int NumGTCs=0;

  char Line[MAXSTR];
  int LineNum;
  int NumTokens;
  char *Tokens[MAXTOK];
  char *ErrMsg;
  while(fgets(Line, MAXSTR, f))
   { 
     // read line, break it up into tokens, skip blank lines and comments
     LineNum++;
     NumTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

     // separately handle the two possible ways to specify complices:
     //  (a) single-line transformations (TAG ... )
     //  (b) TRANSFORMATION ... ENDTRANSFORMATION sections
     if ( !strcasecmp(Tokens[0], "TAG") )
      { 
        // ErrMsg=ProcessTAGLine(Tokens, NumTokens, &GTC);
        printf("howdatage\n");
      }
     else if ( !strcasecmp(Tokens[0], "TRANSFORMATION") )
      { if (NumTokens!=2) 
         ErrExit("%s:%i: syntax error (no name specified for transformation",FileName,LineNum);
        ErrMsg=ParseTRANSFORMATIONSection(Tokens[1], f, &LineNum, &GTC);
        if (ErrMsg) 
         ErrExit("%s:%i: %s",FileName,LineNum,ErrMsg);
      }
     else
      ErrExit("%s:%i: syntax error",FileName,LineNum,Tokens[0]);

     if (ErrMsg)
      ErrExit("%s:%i: %s",FileName,LineNum,ErrMsg);

     // if that was successful, add the new GTComplex to our 
     // array of GTComplex structures
     GTCArray=(GTComplex **)realloc(GTCArray, (NumGTCs+1)*sizeof(GTCArray[0]));
     GTCArray[NumGTCs]=GTC;
     NumGTCs++;

   }; // while(fgets(Line, MAXSTR, f))

  fclose(f);
  *NumGTComplices=NumGTCs; 
  Log("Read %i geometrical transformations from file %s.\n",NumGTCs,FileName); 
  return GTCArray;

}

} // namespace scuff
