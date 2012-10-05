/*
 * TDRTGeometry .cc -- implementation of the TDRTGeometry class 
 *
 * homer reid     -- 11/2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include "libhrutil.h"
#include "libTDRT.h"

#define MAXOBJECTS 100
#define MAXTOKENS 10 
#define MAXSTR 1000

int TDRTGeometry::LogLevel=1;
Interp1D *TDRTGeometry::G1G2Interp=0;
 
/*****************************************************************/
/* TDRTGeometry class constructor: Create an instance of an      */
/* TDRTGeometry from a .tdgeo file.                              */
/*                                                               */
/* This file should have the following syntax:                   */
/*                                                               */
/* OBJECT Object1.msh [LABEL label1] [MATERIAL Material1]        */
/*    ...                                                        */
/* OBJECT ObjectN.msh [LABEL labelN] [MATERIAL MaterialN]        */
/* [MEDIUM MATERIAL ExternalMaterial]                            */
/*                                                               */
/* where                                                         */
/*                                                               */
/*  -- the LABEL keyword allows you to define an optional label  */
/*     for each object                                           */
/*                                                               */
/*  -- the MATERIAL keyword allows you to specify an optional    */
/*     libMatProp material name for each object                  */
/*                                                               */
/*  -- the MEDIUM MATERIAL line allows you to specify an         */
/*     optional material name for the external medium            */
/*****************************************************************/
TDRTGeometry::TDRTGeometry(const char *pGeoFileName)
{ 

  FILE *f;
  TDRTObject *O, *ObjectArray[MAXOBJECTS];
  char Line[MAXSTR], MeshFileName[MAXSTR], Label[MAXSTR], Material[MAXSTR];
  char *Tokens[50];
  int i, j, no, nt, nTokens, LineNum;

  /***************************************************************/
  /* initialize simple fields ************************************/
  /***************************************************************/
  NumObjects=TotalBFs=0;
  GeoFileName=strdup(pGeoFileName);
  LogLevel=1;
  MP=0;

  /***************************************************************/
  /* NOTE: i am not sure where to put this. put it here for now. */
  /***************************************************************/
  MatProp::SetLengthUnit(1.0e-6);

  /***************************************************************/
  /* try to open input file **************************************/
  /***************************************************************/
  f=fopen(GeoFileName,"r");
  if (!f)
   ErrExit("could not open %s",GeoFileName);

  /***************************************************************/
  /* read and process lines from input file one at a time        */
  /***************************************************************/
  LineNum=0; 
  while( fgets(Line,MAXSTR,f) )
   { 
     LineNum++;

     /*--------------------------------------------------------------*/
     /*- break up line into tokens; skip blank lines and comments ---*/
     /*--------------------------------------------------------------*/
     nTokens=Tokenize(Line, Tokens, 50);;
     if ( nTokens==0 || Tokens[0][0]=='#' )
      continue; 
    
     /*--------------------------------------------------------------*/
     /*- switch off based on first token ----------------------------*/
     /*--------------------------------------------------------------*/
     if ( !StrCaseCmp(Tokens[0],"MEDIUM") )
      { 
        /*--------------------------------------------------------------*/
        /*- MEDIUM keyword: set material properties of external medium -*/
        /*--------------------------------------------------------------*/
        if ( nTokens!=3 || StrCaseCmp(Tokens[1],"MATERIAL") )
         ErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
        MP=new MatProp(Tokens[2]);
        if (MP->ErrMsg)
         ErrExit("file %s:%i: error in MATERIAL value: %s",GeoFileName,LineNum,MP->ErrMsg);

      }
     else if ( !StrCaseCmp(Tokens[0],"OBJECT") )
      { 
        /*--------------------------------------------------------------*/
        /*- OBJECT keyword: add a new object to the geometry           -*/
        /*--------------------------------------------------------------*/

        if ( nTokens<2 )
         ErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
        strncpy(MeshFileName,Tokens[1], MAXSTR);

        Material[0]=0;
        sprintf(Label,"Object%i",NumObjects+1); // default label
        for(nt=2; nt<nTokens; nt++)
         { 
           if ( !StrCaseCmp(Tokens[nt],"LABEL") )
            { if ( nt+2 > nTokens )
               ErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
              strncpy(Label,Tokens[++nt],MAXSTR);
            }
           else if ( !StrCaseCmp(Tokens[nt],"MATERIAL") )
            { if ( nt+2 > nTokens )
               ErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
              strncpy(Material,Tokens[++nt],MAXSTR);
            }
           else
            ErrExit("file %s:%i: unknown keyword %s",GeoFileName,LineNum,Tokens[nt]);
         };

        if (NumObjects==MAXOBJECTS)
         ErrExit("%s:%i: too many objects \n",GeoFileName,LineNum);

        O=new TDRTObject(MeshFileName,Label,Material);
        if (O->MP->ErrMsg)
         ErrExit("file %s:%i: error in MATERIAL value: %s",GeoFileName,LineNum,O->MP->ErrMsg);

        TotalBFs+=O->NumBFs;
        ObjectArray[NumObjects++]=O;
      }
     else 
      { /*--------------------------------------------------------------*/
        /*- unknown keyword.                                            */
        /*--------------------------------------------------------------*/
        ErrExit("file %s:%i: unknown keyword %s",GeoFileName,LineNum,Tokens[0]);
      };

   }; // while( fgets(Line,MAXSTR,f) )

  /*******************************************************************/
  /* if no material was specified, set the external medium to vacuum */
  /*******************************************************************/
  if (MP==0)
   MP=new MatProp(MP_VACUUM);

  /*******************************************************************/
  /* copy array of objects */
  /*******************************************************************/
  Objects=(TDRTObject **)malloc( (NumObjects+1)*sizeof(TDRTObject *) );
  memcpy(Objects,ObjectArray,NumObjects * sizeof(TDRTObject *)); 
  Objects[NumObjects]=0;

  /*******************************************************************/
  /* set the AllPEC flag based on whether or not all material objects*/
  /* are PEC bodies                                                  */
  /*******************************************************************/
  AllPEC=1;
  for(no=0; no<NumObjects && AllPEC; no++)
   if ( !(Objects[no]->MP->IsPEC()) )
    AllPEC=0;

  /*******************************************************************/
  /* initialize BFIndexOffset array */
  /*******************************************************************/
  BFIndexOffset=(int *)malloc(NumObjects*sizeof(int));
  BFIndexOffset[0]=0;
  for(no=1; no<NumObjects; no++)
   BFIndexOffset[no]=BFIndexOffset[no-1]+Objects[no-1]->NumBFs;

  /*******************************************************************/
  /* initialize Mate[] array.                                        */
  /*                                                                 */
  /* how it works:                                                   */
  /*                                                                 */
  /* (1) two objects are considered identical if                     */
  /*     (a) they have the same mesh file, and                       */
  /*     (b) they have the same material properties (i.e. they       */
  /*         were given identical values for the MATERIAL            */
  /*         keyword in the .tdgeo file.)                            */
  /*                                                                 */
  /* (2) Mate[] array: If objects i, j, k, ... are identical and     */
  /*                   i<j<k<..., then we set                        */
  /*                   Mate[i] = -1                                  */
  /*                   Mate[j] = i                                   */
  /*                   Mate[k] = i                                   */
  /*                   etc.                                          */
  /*******************************************************************/
  Mate=(int *)malloc(NumObjects*sizeof(int));
  Mate[0]=-1;
  for(i=1; i<NumObjects; i++)
   { Mate[i]=-1;
     for(j=0; j<i && Mate[i]==-1; j++)
      if (    !strcmp(Objects[i]->MeshFileName, Objects[j]->MeshFileName)
           && !strcmp(Objects[i]->MP->Name    , Objects[j]->MP->Name)
         ) 
       Mate[i]=j;
   };

  /*******************************************************************/
  /* allocate ObjectMoved array */
  /*******************************************************************/
  ObjectMoved=(int *)malloc(NumObjects*sizeof(int));
  memset(ObjectMoved,0,NumObjects*sizeof(int));
  ObjectRotated=0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  InitHRBesselK();
  if (G1G2Interp==0)
   InitG1G2Interp();

}

/***************************************************************/
/* TDRTGeometry class destructor *******************************/
/***************************************************************/
TDRTGeometry::~TDRTGeometry()
{
  int no;

  for(no=0; no<NumObjects; no++)
   delete Objects[no];

  free(Objects);
  free(Mate);
  free(GeoFileName);
  free(BFIndexOffset);          

}

/***************************************************************/
/* portable strcasestr replacement******************************/
/***************************************************************/

static char *my_strcasestr(const char *haystack, const char *needle)
{
  while (*haystack) {
    const char *p = haystack;
    const char *n = needle;
    while (*n && *p && tolower(*p) == tolower(*n)) {
      ++n; ++p;
    }
    if (!*n) return const_cast<char*>(haystack);
    haystack++;
  }
  return NULL;
}

/***************************************************************/
/* Transform: Process a 'transformation' line that may contain */
/* separate displacements for one or more objects in the       */
/* geometry.                                                   */
/*                                                             */
/*                                                             */
/* The "transformation" line is a string containing one or     */
/* more of the following sections:                             */
/*  TAG    mystr                    (description of transform) */
/*  LABEL  obj1                     (object to transform)      */
/*  DISP   x1 y1                    (displacement)             */
/*  ROT    Theta1                   (rotation)                 */
/*  LABEL  obj2                     (object to transform)      */
/*  DISP   x2 y2                    (displacement)             */
/* etc.                                                        */
/*                                                             */
/* If any error is incurred in the processing of the           */
/* transformation line, a nonzero value is returned and an     */
/* error message is written to ErrMsg.                         */
/* Otherwise (i.e. if successful) zero is returned.            */
/***************************************************************/
int TDRTGeometry::Transform(char *TransLine, char *Tag, char *ErrMsg)
{ 
  int no, nRead, nConv;
  TDRTObject *O;
  char *p, *pp, Token[MAXSTR], ObjectLabel[MAXSTR], TagBuf[MAXSTR];

  if (ErrMsg) ErrMsg[0]=0;

  /* skip blank lines and comments */
  p=TransLine; 
  while( isspace(*p) )
   p++;
  if ( *p==0 || *p=='#' )
   return 0;

  /* assume no objects will be moved */
  memset(ObjectMoved,0,NumObjects*sizeof(int));

  /*
   * parse transformation line
   */
  sprintf(TagBuf,"notag"); 
  while( p && *p )
   {
     /*--------------------------------------------------------------*/
     /*- read next token off of line --------------------------------*/
     /*--------------------------------------------------------------*/
     nConv=sscanf(p,"%s%n",Token,&nRead);  
     p+=nRead;  
     if ( nConv<=0 || Token[0]=='\n' ) 
      break;
     
     /*--------------------------------------------------------------*/
     /* parse TAG element */
     /*--------------------------------------------------------------*/
     if ( !StrCaseCmp(Token,"TAG") )
      {
        sscanf(p,"%s%n",TagBuf,&nRead);
        p+=nRead;
        if ( nConv!=1 )
         { if (ErrMsg) sprintf(ErrMsg,"syntax error");
           return 1;
         };
      }
     /*--------------------------------------------------------------*/
     /* parse LABEL element.                                         */
     /*--------------------------------------------------------------*/
     else if ( !StrCaseCmp(Token,"LABEL") )
      {  
        /* read object label */
        nConv=sscanf(p,"%s%n",ObjectLabel,&nRead);
        p+=nRead;
        if ( nConv!=1 )
         { if (ErrMsg) sprintf(ErrMsg,"syntax error");
           return 1;
         };

        /* find matching object */
        for(no=0; no<NumObjects; no++)
         if ( !StrCaseCmp(Objects[no]->Label,ObjectLabel) )
          break;
        if(no==NumObjects) 
         { if (ErrMsg) sprintf(ErrMsg,"unknown object %s",ObjectLabel);
           return 1;
         };
        O=Objects[no];
        ObjectMoved[no]=1;

        /* now basically send everything between here and the  */
        /* next instance of the 'LABEL' keyword to object O to */
        /* process as a transformation                         */ 
        pp=my_strcasestr(p,"LABEL");
        if ( !pp )
         { if ( O->Transform(p) )
            { if (ErrMsg) sprintf(ErrMsg,"invalid transformation");
              return 1;
            };
           p=0;
         }
        else
         { *pp=0; // replace the 'L' in the LABEL keyword with end-of-string 
           if ( O->Transform(p) )
            { if (ErrMsg) sprintf(ErrMsg,"invalid transformation");
              return 1;
            };
           *pp='L'; // put back the 'L' in the LABEL keyword
           p=pp;
         };

        if (O->ObjectWasRotated)
         ObjectRotated=1;

      }   // else if ( !StrCaseCmp(Token,"LABEL") )
     else 
      { if (ErrMsg) sprintf(ErrMsg,"syntax error");
        return 1;
      };

   }; // while( p && *p ) 
  if (Tag) strcpy(Tag,TagBuf);

  return 0;

}

/***************************************************************/
/* Undo transformations. ***************************************/
/***************************************************************/
void TDRTGeometry::UnTransform()
{ int no;
  for(no=0; no<NumObjects; no++)
   Objects[no]->UnTransform();
}

/***************************************************************/
/* initialize the internal interpolation table for the G1, G2  */
/* functions                                                   */
/***************************************************************/
#define NUMFUNCTIONS 2
#define NUMPOINTS 5000
extern double *TValues;
extern double *G1G2Values;
void TDRTGeometry::InitG1G2Interp()
{
  G1G2Interp = new Interp1D(TValues, G1G2Values, NUMPOINTS, NUMFUNCTIONS);
}
