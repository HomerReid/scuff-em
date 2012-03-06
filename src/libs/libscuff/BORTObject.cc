/*
 * BORTObject.cc
 *
 * homer reid    -- 11/2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "libscuffInternals.h"

#include <libTriInt.h>

namespace scuff {

#define MAXREFPTS 10

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct StraddlingEdge 
 { double QP[3], V1[3], V2[3], QM[3];
 } StraddlingEdge;

class BORTObject
 { 
   RWGObject *BaseObject;
   GTransformation *GT;
   int NumPieces;
   double *TVertices;                 // storage for transformed vertices
   StraddlingEdge *StraddlingEdges;
   int NumStraddlingEdges;

 } BORTObject;

/***************************************************************/
/* Helper function for the BORTObject constructor that parses  */
/* the relevant section of a .scuffgeo file.                   */
/***************************************************************/
void BORTObject::ParseBORTSectionInFile(FILE *f, char *FileName, int *LineNum)
{ 
  int nt, NumTokens;
  char Line[1000], *Tokens[50];
  char BaseObjectName[100], MaterialName[100], Label[100];
  
  /***************************************************************/
  /* read lines one at a time from the file and parse by keyword */
  /***************************************************************/
  GT=0;
  NumPieces=0;
  BaseObjectName[0]=MaterialName[0]=Label[0]=0;
  while ( fgets(Line,1000,f) )
   { 
     /*--------------------------------------------------------------*/
     /*- skip blank lines and comments ------------------------------*/
     /*--------------------------------------------------------------*/
     (*LineNum)++;
     NumTokens=Tokenize(Line,Tokens,50);
     if (NumTokens==0 || Tokens[0][0]=='#') 
      continue;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     if (!strcasetmp(Tokens[0],"BASE"))
      { 
        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        if (NumTokens!=2)
         ErrExit("%s:%i: syntax error",FileName,LineNum);
        strncpy(BaseObjectName, 100, Tokens[1]);
      }
     else if (!strcasetmp(Tokens[0],"MATERIAL"))
      { 
        if (NumTokens!=2)
         ErrExit("%s:%i: syntax error",FileName,LineNum);
        strncpy(MaterialName,100,Tokens[1]);
      }
     else if (!strcasetmp(Tokens[0],"LABEL"))
      { 
        if (NumTokens!=2)
         ErrExit("%s:%i: syntax error",FileName,LineNum);
        strncpy(Label,100,Tokens[1]);
      }
     else if (!strcasetmp(Tokens[0],"DISP"))
      { 
        if (GT!=0)
         ErrExit("%s:%i: multiple transformations not allowed",FileName,LineNum);
        if (NumTokens!=4)
         ErrExit("%s:%i: syntax error",FileName,LineNum);

        double DX[3];
        if (    1!=sscanf(Tokens[1],"%le",DX+0)
             || 1!=sscanf(Tokens[2],"%le",DX+1)
             || 1!=sscanf(Tokens[3],"%le",DX+2) 
           ) ErrExit("%s:%i: invalid displacement",FileName,LineNum);

        GT=CreateGTransformation(DX);
      }
     else if (!strcasecmp(Tokens[0],"ROT"))
      { 
        if (GT!=0)
         ErrExit("%s:%i: multiple transformations not allowed",FileName,LineNum);
        if (NumTokens!=5)
         ErrExit("%s:%i: syntax error",FileName,LineNum);

        double ZHat[3], Angle;
        if (    1!=sscanf(Tokens[1],"%le",ZHat+0)
             || 1!=sscanf(Tokens[2],"%le",ZHat+1)
             || 1!=sscanf(Tokens[3],"%le",ZHat+2) 
             || 1!=sscanf(Tokens[4],"%le",&Angle) 
           ) ErrExit("%s:%i: invalid rotation",FileName,LineNum);

        GT=CreateGTransformation(ZHat, Angle);
      }
     else if (!strcasecmp(Tokens[0],"PIECES"))
      { 
        if ( NumTokens!=2 || 1!=sscanf(Tokens[1],"%i",&NumPieces))
         ErrExit("%s:%i: syntax error",FileName,LineNum);
      }
     else if (!strcasetmp(Tokens[0],"ENDBORT"))
      break;
     else
      ErrExit("%s: %i: unknown keyword %s",FileName,LineNum,Tokens[0]);

   }; // while ( fgets(Line,1000,f) )

   /***************************************************************/
   /* make sure required fields were specified, and fill in       */
   /* default values for other fields                             */
   /***************************************************************/
   if (BaseObjectName[0]==0)
    ErrExit("%s: base object must be specified");
   if (NumPieces<=0)
    ErrExit("%s: invalid number of PIECES");
   if (GT==0)
    ErrExit("%s: no geometrical transformation specified");
   if (MaterialName[0]==0)
    sprintf(MaterialName,"PEC");
   static int NumBortObjects=1;
   if (Label[0]==0)
    sprintf(Label,"BORTObject_%i",NumBortObjects++);

   /***************************************************************/
   /* attempt to create the base object ***************************/
   /***************************************************************/
   BaseObject=new RWGObject(BaseObjectName, Label, MaterialName, 0, 0);

}

/***************************************************************/
/* *************************************************************/
/***************************************************************/
#define ECTHRESHOLD 1.0e-6
int EdgesCoincide(int iV1, int iV2,   double Vertices,
                  int iV1P, int iV2P, double VerticesP, double Length)
{ 
  if (    VecDistance(Vertices+3*iV1, VerticesP+3*iV1P) < ECTHRESHOLD*Length
       && VecDistance(Vertices+3*iV2, VerticesP+3*iV2P) < ECTHRESHOLD*Length 
     ) return 1;

  if (    VecDistance(Vertices+3*iV1, VerticesP+3*iV2P) < ECTHRESHOLD*Length
       && VecDistance(Vertices+3*iV2, VerticesP+3*iV1P) < ECTHRESHOLD*Length 
     ) return 1;

  return 0;
} 

/***************************************************************/
/* Helper function for the BORT constructor that identifies    */
/* RWG basis functions that straddle the boundary between the  */
/* base object and its image under the geometrical transform.  */
/*                                                             */
/* algorithm:                                                  */
/*  for each exterior edge E1 in the base object, look for a   */
/*  second exterior edge E2 such that E2' = E1, where E2' is   */
/*  the image of E2 under the transformation.                  */
/*  if we find such an pair, create a new 'straddling edge'    */
/*  structure.                                                 */
/***************************************************************/
BORTObject::IdentifyStraddlingEdges()
{
  int NumVertices=BaseObject->NumVertices;

  /*--------------------------------------------------------------*/
  /*- set OVertices = original vertices             --------------*/
  /*-     TVertices = vertices under transformation --------------*/
  /*--------------------------------------------------------------*/
  double *OVertices=BaseObject->Vertices;
  ApplyGTransformation(GT, Vertices, TVertices, NumVertices);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  StraddlingEdges=0;
  NumStraddlingEdges=0;

  RWGEdge *EE, *EEP;
  StraddlingEdge *SE;
  int nee, neep;
  int NEE=BaseObject->NumExteriorEdges;
  int AccountedFor[NEE];
  memset(AccountedFor,0,NEE*sizeof(int));
  for(nee=0; nee<NEE; nee++)
   { 
     /*--------------------------------------------------------------*/
     /*- loop over all exterior edges EE-----------------------------*/
     /*--------------------------------------------------------------*/
     if (AccountedFor[nee])
      continue;
     EE=BaseObject->ExteriorEdges[nee];

     for(neep=0; neep<NEE; neep++)
      { 
        /*--------------------------------------------------------------*/
        /* look over all exterior edges EEP, looking for EEP' = EE      */
        /*--------------------------------------------------------------*/
        if (AccountedFor[neep])
         continue;
        EEP=BaseObject->ExteriorEdges[neep];

        /*--------------------------------------------------------------*/
        /*- if we found EEP'=EE, create a new StraddlingEdge structure -*/
        /*--------------------------------------------------------------*/
        if ( EdgesCoincide( EE->iV1,  EE->iV2,  OVertices,
                            EEP->iV1, EEP->iV2, TVertices, EE->Length)
           ) { 
               SE=(StraddlingEdge *)mallocEC(sizeof(StraddlingEdge));
               memcpy(SE->QP, OVertices + 3*(EE->iQP), 3*sizeof(double));
               memcpy(SE->V1, OVertices + 3*(EE->iV1), 3*sizeof(double));
               memcpy(SE->V2, OVertices + 3*(EE->iV2), 3*sizeof(double));
               memcpy(SE->QM, TVertices + 3*(EEP->iQP), 3*sizeof(double));

               StraddlingEdges=realloc(StraddlingEdges, (NumStraddlingEdges+1)*sizeof(StraddlingEdge));
               StraddlingEdges[NumStraddlingEdges++]=SE;

               AccountedFor[nee]=AccountedFor[neep]=1;
               break;
             }; // if (EdgesCoincide ...)

      }; // for(neep=0; ...

   }; for(nee=0; ... 
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- BORTObject class constructor.                               */
/*-                                                             */
/*- On entry, f is an open file whose read pointer points to    */
/*- the line immediate following the BORT keyword in a .scuffgeo*/
/*- file. FileName is the name of this file, and *LineNum is    */
/*- the line number, which is updated on return.                */
/*-                                                             */
/*- If the function returns, then the BORTObject has been       */
/*- successfully created, and the file read pointer points to   */
/*- the line immediately after the ENDBORT keyword.             */
/*--------------------------------------------------------------*/
BORTObject::BORTObject(FILE *f, char *FileName, int *LineNum)
{
   /***************************************************************/
   /* parse the file and create the base object *******************/
   /***************************************************************/
   ParseBORTSectionInFile(f, FileName, LineNum);

   /***************************************************************/
   /* create TVertices array to store transformed vertices        */
   /***************************************************************/
   int NumVertices=BaseObject->NumVertices;
   TVertices=(double *)mallocEC(3*BaseObject->NumVertices*sizeof(double));

   /***************************************************************/
   /* identify straddling edges                                   */
   /***************************************************************/
   IdentifyStraddlingEdges();

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ABOMBData
 {
   BORTObject *BORT;
   double *TVertices;
 } ABOMBData;

void *AssembleBORTMatrixBlocksThread(void *Data)
{
  /***************************************************************/
  /* unpack fields from thread data structure ********************/
  /***************************************************************/
  ABOMBData *D=(ABOMBData *)Data;
  BORTObject *BORT=D->BORT;

  RWGObject *BO=BORT->BaseObject;
  int NumVertices

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memcpy(
  for(np=0; np<NumPieces-1; np++)
   { 
     ApplyGTransformation(BORT->GT, BORT->BaseObject->Vertices,
                          TVertices, BORT->BaseObject->
   };


}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AssembleBORTMatrixBlocks(
{
  for(np=0; np<NumPieces-1; np++)
   T
}

} // namespace scuff
