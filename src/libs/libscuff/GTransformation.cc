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
#include "libscuff.h" // for Vec* functions

#define MAXSTR 1000
#define MAXTOK 50

namespace scuff {

/***************************************************************/
// Identity transformation:
void GTransformation::Reset()
{ 
  DX[0] = DX[1] = DX[2] = 
    M[0][1] = M[0][2] = M[1][0] = M[1][2] = M[2][0] = M[2][1] = 0.0;
  M[0][0] = M[1][1] = M[2][2] = 1.0;
}

/***************************************************************/
/***************************************************************/
// Constructors:

GTransformation::GTransformation() { Reset(); }

GTransformation::GTransformation(const double dx[3]) {
  Reset();
  Displace(dx);
}

GTransformation::GTransformation(const double ZHat[3], double ThetaDegrees) {
  Reset();
  Rotate(ZHat, ThetaDegrees);
}

GTransformation::GTransformation(const char *TransformString, char **ErrMsg) {
  Reset();
  Parse(TransformString, ErrMsg);
}

GTransformation::GTransformation(char **Tokens, int NumTokens, char **ErrMsg,
				 int *TokensConsumed) {
  Reset();
  Parse(Tokens, NumTokens, ErrMsg, TokensConsumed);
}

GTransformation::GTransformation(const GTransformation &G) {
  VecCopy(G.DX, DX);
  memcpy(&M[0][0], &G.M[0][0], sizeof(double)*9);
}

GTransformation::GTransformation(const  GTransformation *G) {
  VecCopy(G->DX, DX);
  memcpy(&M[0][0], &G->M[0][0], sizeof(double)*9);
}

void GTransformation::operator=(const GTransformation &G) {
  VecCopy(G.DX, DX);
  memcpy(&M[0][0], &G.M[0][0], sizeof(double)*9);
}

/***************************************************************/
/* a transformation that displaces through vector DX           */
/***************************************************************/
void GTransformation::Displace(const double dx[3]) {
  VecAdd(DX, dx, DX);
}

/***************************************************************/
/* a transformation that rotates through Theta degrees         */
/* (DEGREES, NOT RADIANS) about an axis that passes through    */
/* the origin and the point with coordinates ZHat[0..2]        */
/***************************************************************/

// rotate a 3-vector by Theta (given its cos & sin) around the UNIT vector U
static void RotateVec(const double U[3], double cosTheta, double sinTheta, 
		      double V[3]) {
  double Vpar[3], Vperp[3]; // components of V parallel/perpendicular to U
  double Vcross[3]; // U x V

  VecScale(VecCopy(U, Vpar), VecDot(U, V)); // compute Vpar
  VecSub(V, Vpar, Vperp); // compute Vperp
  VecCross(U, V, Vcross); // compute Vcross = U x V
  VecAdd(Vpar, VecAdd(VecScale(Vperp,cosTheta), VecScale(Vcross,sinTheta), V),
	 V); // V = Vpar + Vperp*cos + Vcross*sin
}

// same as RotateVec, but operates on V[0],V[stride],V[2*stride]
static void RotateVec(const double U[3], double cosTheta, double sinTheta,
                      double *V, int stride) {
  double Vtmp[3];
  Vtmp[0] = V[0]; Vtmp[1] = V[stride]; Vtmp[2] = V[2*stride]; 
  RotateVec(U, cosTheta, sinTheta, Vtmp);
  V[0] = Vtmp[0]; V[stride] = Vtmp[1]; V[2*stride] = Vtmp[2]; 
}

void GTransformation::Rotate(const double ZHat[3], double ThetaDegrees) {
  double U[3], cosTheta, sinTheta;
  VecNormalize(VecCopy(ZHat, U));

  ThetaDegrees *= (3.14159265358979323846 / 180.0); // convert to radians
  cosTheta = cos(ThetaDegrees); sinTheta = sin(ThetaDegrees);

  // apply rotation to DX:
  RotateVec(U, cosTheta, sinTheta, DX);

  // apply rotation to each column of M
  RotateVec(U, cosTheta, sinTheta, &M[0][0], 3);
  RotateVec(U, cosTheta, sinTheta, &M[0][1], 3);
  RotateVec(U, cosTheta, sinTheta, &M[0][2], 3);
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
/*                                                             */
/* Returns true on success, false on failure.                  */
/***************************************************************/

bool GTransformation::Parse(const char *TransformString, char **ErrMsg) {
  char Line[MAXSTR];
  char *Tokens[MAXTOK];
  int NumTokens;
  
  strncpy(Line,TransformString,MAXSTR);
  NumTokens=Tokenize(Line, Tokens, MAXTOK);
  return Parse(Tokens, NumTokens, ErrMsg, NULL);
}

bool GTransformation::Parse(char **Tokens, int NumTokens, 
			    char **ErrMsg, int *TokensConsumed)
{
 
  double Theta, ZHat[3], dx[3];

  if (NumTokens==0) {
    if (ErrMsg) *ErrMsg = strdup("no tranformation specified");
    if (TokensConsumed) *TokensConsumed=0;
    return false;
  }
  else if ( !strcasecmp(Tokens[0],"DISP") 
	    || !strcasecmp(Tokens[0],"DISPLACED") )
   { 
     /*--------------------------------------------------------------*/
     /*-- DISPLACED xx yy zz ----------------------------------------*/
     /*--------------------------------------------------------------*/
     if ( NumTokens<4 )
      { if (ErrMsg) *ErrMsg=vstrdup("too few values specified for %s ",Tokens[0]);
        return false;
      };

     if (    1!=sscanf(Tokens[1],"%le",dx+0)
          || 1!=sscanf(Tokens[2],"%le",dx+1)
          || 1!=sscanf(Tokens[3],"%le",dx+2)
        )
      { if (ErrMsg) *ErrMsg=vstrdup("bad value specified for %s",Tokens[0]);
        return false;
      };

     if (ErrMsg) *ErrMsg=0;
     if (TokensConsumed) *TokensConsumed=4;
     Displace(dx);
   }
  else if ( !strcasecmp(Tokens[0],"ROT") || !strcasecmp(Tokens[0],"ROTATED") )
   { 
     /*--------------------------------------------------------------*/
     /*-- ROTATED tt ABOUT xx yy zz ---------------------------------*/
     /*--------------------------------------------------------------*/
     if ( NumTokens<6 )
      { if (ErrMsg) *ErrMsg=vstrdup("too few values specified for %s ",Tokens[0]);
        return false;
      };

     if ( strcasecmp(Tokens[2],"ABOUT") )
      {  if (ErrMsg) *ErrMsg=vstrdup("invalid syntax for %s statement",Tokens[0]);
         return false;
      };

     if (    1!=sscanf(Tokens[1],"%le",&Theta)
          || 1!=sscanf(Tokens[3],"%le",ZHat+0)
          || 1!=sscanf(Tokens[4],"%le",ZHat+1)
          || 1!=sscanf(Tokens[5],"%le",ZHat+2)
        )
      { if (ErrMsg) *ErrMsg=vstrdup("bad value specified for %s",Tokens[0]);
        return false;
      };

     if (ErrMsg) *ErrMsg=0;
     if (TokensConsumed) *TokensConsumed=6;
     Rotate(ZHat, Theta);
   }
  else
   { 
     if (ErrMsg) *ErrMsg=vstrdup("unknown keyword %s",Tokens[0]); 
     if (TokensConsumed) *TokensConsumed=0;
     return false;
   }
  return true;
} 

/***************************************************************/
// Inversion: if Apply(X) computes Y = M*X+DX, then
//            inverse = UnApply(Y) computes X = M'*(Y-DX) = M'*Y - M'*DX

GTransformation GTransformation::Inverse() const {
  GTransformation C;

  for (int i = 0; i < 3; ++i) {
    C.DX[i] = 0.0;
    for (int j = 0; j < 3; ++j) {
      C.DX[i] -= M[j][i] * DX[j]; // -M' * DX
      C.M[i][j] = M[j][i]; // M'
    }
  }

  return C;
}

void GTransformation::Invert() {
  *this = Inverse();
}

/***************************************************************/
/* Composition: G1 + G2 or G1 - G2 form the composition of G1  */
/* with G2 or inv(G2).  The composition G1 + G2 is equivalent  */
/* to FIRST performing G2 and THEN performing G1 (= this).     */
/***************************************************************/

void GTransformation::Transform(const GTransformation *G) { // G + this
  // new Apply(X) = G.Apply(this.Apply(X)) = G.Apply(M*X + DX)
  //              = (G.M*M)*X + G.M*DX + G.DX

  G->Apply(DX); // new DX = G.M*DX + G.DX

  // new M = G.M * M
  double newM[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      newM[i][j] = 0;
      for (int k = 0; k < 3; ++k)
	newM[i][j] += G->M[i][k] * M[k][j];
    }
  memcpy(&M[0][0], &newM[0][0], sizeof(double)*9);
}

GTransformation GTransformation::operator+(const GTransformation &G2) const {
  GTransformation C(G2);
  C.Transform(this);
  return C;
}

GTransformation GTransformation::operator-(const GTransformation &G2) const {
  GTransformation C(-G2);
  C.Transform(this);
  return C;
}

/***************************************************************/
/* Apply the transformation (out-of-place) to a list of NX points  */
/* X[3*nx+0, 3*nx+1, 3*nx+2] = cartesian coords of point #nx,  */
/* storing the result in XP.                                   */
/* NOTE: we allow XP == X (i.e., can operate in-place).        */
/***************************************************************/
void GTransformation::Apply(const double *X, double *XP, int NX) const {
  for (int nx = 0; nx < NX; ++nx) {
    double x0 = X[3*nx], x1 = X[3*nx+1], x2 = X[3*nx+2];
    XP[3*nx+0] = M[0][0]*x0 + M[0][1]*x1 + M[0][2]*x2 + DX[0];
    XP[3*nx+1] = M[1][0]*x0 + M[1][1]*x1 + M[1][2]*x2 + DX[1];
    XP[3*nx+2] = M[2][0]*x0 + M[2][1]*x1 + M[2][2]*x2 + DX[2];
  }
}

// same, but apply the inverse transformation:
void GTransformation::UnApply(const double *X, double *XP, int NX) const {
  for (int nx = 0; nx < NX; ++nx) {
    double x0 = X[3*nx]-DX[0], x1 = X[3*nx+1]-DX[1], x2 = X[3*nx+2]-DX[2];
    XP[3*nx+0] = M[0][0]*x0 + M[1][0]*x1 + M[2][0]*x2;
    XP[3*nx+1] = M[0][1]*x0 + M[1][1]*x1 + M[2][1]*x2;
    XP[3*nx+2] = M[0][2]*x0 + M[1][2]*x1 + M[2][2]*x2;
  }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/

/********************************************************************/
/* this is a helper function for the ReadTransFile() routine below  */
/* that attempts to parse the shorthand syntax for specifying       */
/* transformations in a .trans file (i.e. a single line beginning   */
/* with the token TRANS).                                           */
/* on entry we have already verified that Tokens[0]==TRANS so that  */
/* does not need to be checked.                                     */
/********************************************************************/
char *ParseTRANSLine(char **Tokens, int NumTokens, GTComplex **pGTC)
{
  if ( NumTokens==1 )
   return strdup("empty TRANS specification");

  /***************************************************************/
  /* initialize a bare GTComplex *********************************/
  /***************************************************************/
  GTComplex *GTC=(GTComplex *)mallocEC(sizeof(GTComplex));
  GTC->Tag=strdup(Tokens[1]);
  GTC->NumObjectsAffected=0;
  GTC->ObjectLabel=0;
  GTC->GT=0;

  GTransformation *CurrentGT=0;

  /***************************************************************/
  /* parse sections of the token string one at at time.          */
  /* each section is treated as though it were on a separate     */
  /* line of a TRANSFORMATION...ENDTRANSFORMATION section.       */
  /***************************************************************/
  int nt=2; // start parsing the rest of the line at token #2
  int TokensConsumed;
  int noa;
  char *ErrMsg;
  while( nt<NumTokens )
   { 
     if ( !strcasecmp(Tokens[nt], "OBJECT") )
      { 
        if ( nt+1 == NumTokens )
         return strdup("syntax error");

        noa=GTC->NumObjectsAffected;
        GTC->ObjectLabel=(char **)realloc(GTC->ObjectLabel, (noa+1)*sizeof(char *));
        GTC->ObjectLabel[noa]=strdup(Tokens[nt+1]);
        GTC->GT=(GTransformation **)realloc(GTC->GT, (noa+1)*sizeof(GTransformation *));
        GTC->GT[noa]=new GTransformation();
        GTC->NumObjectsAffected=noa+1;
        CurrentGT=GTC->GT[noa];

        nt+=2;
      }
     else if (    !strcasecmp(Tokens[nt], "DISP") || !strcasecmp(Tokens[nt], "DISPLACED")
               || !strcasecmp(Tokens[nt], "ROT")  || !strcasecmp(Tokens[nt], "ROTATED")
             )
      { 
        if ( CurrentGT==0 )
         return vstrdup("no OBJECT specified for %s",Tokens[nt]);

	CurrentGT->Parse(Tokens+nt, NumTokens-nt,
			 &ErrMsg, &TokensConsumed);

        if (ErrMsg)
         return ErrMsg;

        nt += TokensConsumed;
      }
     else
      return vstrdup("unknown keyword %s %i",Tokens[nt],nt);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  *pGTC=GTC;
  return 0;

  
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
  GTComplex *GTC=(GTComplex *)mallocEC(sizeof(GTComplex));
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
        GTC->GT[noa]=new GTransformation();
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
        if (CurrentGT==0)
         return vstrdup("no OBJECT specified for %s",Tokens[0]);

	CurrentGT->Parse(Tokens, NumTokens, &ErrMsg, &TokensConsumed);
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
/***************************************************************/
/***************************************************************/
GTComplex *CreateDefaultGTComplex()
{
  GTComplex *GTC=(GTComplex *)mallocEC(sizeof(GTComplex));
  GTC->Tag=strdup("DEFAULT");
  GTC->NumObjectsAffected=0;
  GTC->ObjectLabel=0;
  GTC->GT=0;
  return GTC;
}

/***************************************************************/
/* this routine reads a scuff-EM transformation (.trans) file  */
/* and returns an array of GTComplex structures.               */
/***************************************************************/
GTComplex **ReadTransFile(char *FileName, int *NumGTComplices)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (FileName==0)
   { 
     *NumGTComplices=1;
     GTComplex **GTCList = (GTComplex **)mallocEC(sizeof(GTCList[0]));
     GTCList[0] = CreateDefaultGTComplex();
     Log("Using a single (empty) geometrical transformation with label %s.",GTCList[0]->Tag);
     return GTCList;
   };

  FILE *f=fopen(FileName,"r");
  if (f==0)
   ErrExit("could not open file %s",FileName);

  GTComplex *GTC, **GTCArray=0;
  int NumGTCs=0;

  char Line[MAXSTR];
  int LineNum;
  int NumTokens;
  char *Tokens[MAXTOK];
  char *ErrMsg=0;
  while(fgets(Line, MAXSTR, f))
   { 
     // read line, break it up into tokens, skip blank lines and comments
     LineNum++;
     NumTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

     // separately handle the two possible ways to specify complices:
     //  (a) TRANSFORMATION ... ENDTRANSFORMATION sections
     //  (b) single-line transformations (TRANS ... )
     if ( !strcasecmp(Tokens[0], "TRANSFORMATION") )
      { if (NumTokens!=2) 
         ErrExit("%s:%i: syntax error (no name specified for transformation",FileName,LineNum);
        ErrMsg=ParseTRANSFORMATIONSection(Tokens[1], f, &LineNum, &GTC);
      }
     else if ( !strcasecmp(Tokens[0], "TRANS") )
      { 
        ErrMsg=ParseTRANSLine(Tokens, NumTokens, &GTC);
      }
     else
      { 
        ErrMsg=vstrdup("unknown token %s",Tokens[0]);
      };

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
  Log("Read %i geometrical transformations from file %s.",NumGTCs,FileName); 
  return GTCArray;

}

} // namespace scuff
