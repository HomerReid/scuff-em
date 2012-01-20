/*
 * RWGGeometry.cc -- implementation of some methods in the RWGGeometry
 *                 -- class 
 *
 * homer reid      -- 3/2007 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>

#include "libscuff.h"

#include "StaticPPI.h"

#define MAXOBJECTS 100
#define MAXSTR 1000

/***********************************************************************/
/* RWGGeometry class constructor: Create an instance of an             */
/* RWGGeometry from a .rwggeo file.                                    */
/*                                                                     */
/* This file should have the following syntax:                         */
/*                                                                     */
/*  OBJECT Object1.msh [additional options]                            */
/*   ...                                                               */
/*  OBJECT ObjectN.msh [additional options]                            */
/*   ...                                                               */
/*  MEDIUM MATERIAL ExtMaterialName                                    */
/*                                                                     */
/* where [additional options] can any assortment of the following:     */
/*                                                                     */
/*  LABEL    MyLabel                                                   */
/*  INSIDE   SomeOtherObjectLabel                                      */
/*  MATERIAL MyMaterialName                                            */
/***********************************************************************/
RWGGeometry::RWGGeometry(const char *pGeoFileName)
{ 
  FILE *f;
  RWGObject *O, *ObjectArray[MAXOBJECTS];
  char Line[MAXSTR];
  char MeshFileName[MAXSTR];
  char Token[MAXSTR];
  char Label[MAXSTR]; 
  char RotMat[MAXSTR];
  char ContainingObjectLabel[MAXSTR];
  char Material[MAXSTR];
  char *p, *Tokens[50];
  double DX[3];
  int i, j, no, np, nRead, nConv, nt, nTokens, LineNum;

  /***************************************************************/
  /* initialize simple fields ************************************/
  /***************************************************************/
  NumObjects=TotalBFs=TotalPanels=0;
  GeoFileName=strdup(pGeoFileName);
  ExteriorMP=0;

  /***************************************************************/
  /* NOTE: i am not sure where to put this. put it here for now. */
  /***************************************************************/
  MatProp::SetLengthUnit(1.0e-6);

  /***************************************************************/
  /* try to open input file **************************************/
  /***************************************************************/
  f=fopen(GeoFileName,"r");
  if (!f)
   RWGErrExit("could not open %s",GeoFileName);

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
     nTokens=Tokenize(Line, Tokens, 50);
     if ( nTokens==0 || Tokens[0][0]=='#' )
      continue; 
    
     /*--------------------------------------------------------------*/
     /*- switch off based on first token ----------------------------*/
     /*--------------------------------------------------------------*/
     if ( !strcasecmp(Tokens[0],"MEDIUM") )
      { if ( nTokens!=3 || strcasecmp(Tokens[1],"MATERIAL") )
         RWGErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
        ExteriorMP=new MatProp(Tokens[2]);
        if (ExteriorMP->ErrMsg)
         RWGErrExit("file %s:%i: error in MATERIAL value: %s",GeoFileName,LineNum,ExteriorMP->ErrMsg);
      }
     else if ( !strcasecmp(Tokens[0],"OBJECT") )
      { 
        if ( nTokens<2 )
         RWGErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
        strncpy(MeshFileName,Tokens[1], MAXSTR);

        ContainingObjectLabel[0]=Material[0]=RotMat[0]=0;
        memset(DX,0,3*sizeof(double));
        sprintf(Label,"Object%i",NumObjects+1); // default label
        for(nt=2; nt<nTokens; nt++)
         { 
           if ( !strcasecmp(Tokens[nt],"LABEL") )
            { if ( nt+2 > nTokens )
               RWGErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
              strncpy(Label,Tokens[++nt],MAXSTR);
            }
           else if ( !strcasecmp(Tokens[nt],"INSIDE") )
            { if ( nt+2 > nTokens )
               RWGErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
              strncpy(ContainingObjectLabel,Tokens[++nt],MAXSTR);
            }
           else if ( !strcasecmp(Tokens[nt],"MATERIAL") )
            { if ( nt+2 > nTokens )
               RWGErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
              strncpy(Material,Tokens[++nt],MAXSTR);
            }
           else if ( !strcasecmp(Tokens[nt],"ROTMAT") )
            { if ( nt+2 > nTokens )
               RWGErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
              strncpy(RotMat,Tokens[++nt],MAXSTR);
            }
           else if ( !strcasecmp(Tokens[nt],"DISP") )
            { if ( nt+4>nTokens )
               RWGErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
              if (    1!=sscanf(Tokens[++nt],"%le",DX)
                   || 1!=sscanf(Tokens[++nt],"%le",DX+1)
                   || 1!=sscanf(Tokens[++nt],"%le",DX+2) 
                 )
               RWGErrExit("file %s:%i: invalid DISP specification",GeoFileName,LineNum);
            }
           else
            RWGErrExit("file %s:%i: unknown keyword %s",GeoFileName,LineNum,Tokens[nt]);
         };

        if (NumObjects==MAXOBJECTS)
         RWGErrExit("%s:%i: too many objects \n",GeoFileName,LineNum);

        /*--------------------------------------------------------------*/
        /*- create the new RWG object ----------------------------------*/
        /*--------------------------------------------------------------*/
        O=new RWGObject(MeshFileName,Label,Material,RotMat,DX);
        if (O->MP->ErrMsg)
         RWGErrExit("file %s:%i: error in MATERIAL value: %s",
                     GeoFileName,LineNum,O->MP->ErrMsg);

        /*--------------------------------------------------------------*/
        /*- attempt to find the containing object if one was specified -*/
        /*--------------------------------------------------------------*/
        if ( ContainingObjectLabel[0]!=0 )
         { 
           for(no=0; no<NumObjects && O->ContainingObject==0; no++)
            if ( !strcasecmp(ObjectArray[no]->Label,ContainingObjectLabel) )
             O->ContainingObject=ObjectArray[no];  

           if ( O->ContainingObject==0 )
            RWGErrExit("file %s:%i: unknown object (%s) specified for INSIDE",
                        GeoFileName,LineNum,ContainingObjectLabel);
         };

        /*--------------------------------------------------------------*/
        /*- augment tallies to account for the new object  -------------*/
        /*--------------------------------------------------------------*/
        TotalBFs+=O->NumBFs;
        TotalPanels+=O->NumPanels;
        O->Index=NumObjects;
        ObjectArray[NumObjects++]=O;
      }
     else 
      RWGErrExit("file %s:%i: unknown keyword %s",GeoFileName,LineNum,Tokens[0]);

   }; // while( fgets(Line,MAXSTR,f) )

  /*******************************************************************/
  /* if no material was specified, set the external medium to vacuum */
  /*******************************************************************/
  if (ExteriorMP==0)
   ExteriorMP=new MatProp(MP_VACUUM);

  /***************************************************************/
  /* copy array of objects                                       */
  /***************************************************************/
  Objects=(RWGObject **)RWGMalloc( NumObjects*sizeof(RWGObject *) );
  memcpy(Objects,ObjectArray,NumObjects * sizeof(RWGObject *));

  /*******************************************************************/
  /* compute average panel area for statistical bookkeeping purposes */
  /*******************************************************************/
  AveragePanelArea=0.0; 
  for(no=0; no<NumObjects; no++)
   for(np=0; np<Objects[no]->NumPanels; np++)
    AveragePanelArea+=Objects[no]->Panels[np]->Area;
  AveragePanelArea/=((double) TotalPanels);

  /*******************************************************************/
  /* set the AllPEC flag based on whether or not all material objects*/
  /* are PEC bodies                                                  */
  /*******************************************************************/
  AllPEC=1;
  for(no=0; no<NumObjects && AllPEC; no++)
   if ( !(Objects[no]->MP->IsPEC()) )
    AllPEC=0;

  /***************************************************************/
  /* initialize arrays of basis-function and panel index offsets */
  /***************************************************************/
  BFIndexOffset=(int *)RWGMalloc(NumObjects*sizeof(int) );
  PanelIndexOffset=(int *)RWGMalloc(NumObjects*sizeof(int) );
  BFIndexOffset[0]=PanelIndexOffset[0]=0;
  for(no=1; no<NumObjects; no++)
   { BFIndexOffset[no]=BFIndexOffset[no-1] + Objects[no-1]->NumBFs;
     PanelIndexOffset[no]=PanelIndexOffset[no-1] + Objects[no-1]->NumPanels;
   };

  /***************************************************************/
  /* initialize Identical[][] and Mate[] arrays.                 */
  /*                                                             */
  /* how it works:                                               */
  /*                                                             */
  /* (1) two objects are considered identical if                 */
  /*     (a) they have the same mesh file, and                   */
  /*     (b) they have the same material properties (i.e. they   */
  /*         were given identical values for the MATERIAL        */
  /*         keyword in the .rwggeo file.)                       */
  /*                                                             */
  /* (2) Identical[][] array: We set Identical[i][j] = 1 if      */
  /*                          objects i and j are identical, =0  */
  /*                          otherwise.                         */
  /*                                                             */
  /* (3) Mate[] array: If objects i, j, k, ... are identical and */
  /*                   i<j<k<..., then we set                    */
  /*                   Mate[i] = -1                              */
  /*                   Mate[j] = i                               */
  /*                   Mate[k] = i                               */
  /***************************************************************/
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

  /***************************************************************/
  /* initialize ObjectMoved[] array                              */
  /***************************************************************/
  ObjectMoved=(int *)RWGMalloc(NumObjects*sizeof(int));

}

/***************************************************************/
/* RWGGeometry class destructor *******************************/
/***************************************************************/
RWGGeometry::~RWGGeometry()
{
  int no;

  for(no=0; no<NumObjects; no++)
   delete Objects[no];

  free(Objects);

  free(BFIndexOffset);
  free(PanelIndexOffset);
  free(Mate);
  free(ObjectMoved);
  free(GeoFileName);

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
/*  DISP   x1 y1 z1                 (displacement)             */
/*  ROT    Ax Ay Az Theta           (rotation)                 */
/*  LABEL  obj2                     (object to transform)      */
/*  DISP   x2 y2 z2                 (displacement)             */
/* etc.                                                        */
/*                                                             */
/* If any error is incurred in the processing of the           */
/* transformation line, a nonzero value is returned and an     */
/* error message is written to ErrMsg.                         */
/* Otherwise (i.e. if successful) zero is returned.            */
/***************************************************************/
int RWGGeometry::Transform(char *TransLine, char *Tag, char *ErrMsg)
{ 
  int no, nRead, nConv;
  RWGObject *O;
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
     if ( !strcasecmp(Token,"TAG") )
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
     else if ( !strcasecmp(Token,"LABEL") )
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
         if ( !strcasecmp(Objects[no]->Label,ObjectLabel) )
          break;
        if(no==NumObjects)
         { if (ErrMsg) sprintf(ErrMsg,"unknown object %s",ObjectLabel);
           return 1;
         };
        O=Objects[no];

        ObjectMoved[no]=1;

        /* now send everything between here and the  */
        /* next instance of the 'LABEL' keyword to object O to */
        /* process as a transformation                         */ 
        pp=strcasestr(p,"LABEL");
        if ( !pp )
         { if ( O->Transform(p) )
            { if (ErrMsg) sprintf(ErrMsg,"invalid transformation");
              return 1;
            };
           p=0;
         }
        else
         { *pp=0;
           if ( O->Transform(p) )
            { if (ErrMsg) sprintf(ErrMsg,"invalid transformation");
              return 1;
            };
           *pp='L';  // put back the 'L' in the LABEL keyword
           p=pp;
         };
      }   // else if ( !strcasecmp(Token,"LABEL") )
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
void RWGGeometry::UnTransform()
{ int no;
  for(no=0; no<NumObjects; no++)
   Objects[no]->UnTransform();
}

/***************************************************************/
/* Return the dimension of the linear system. ******************/
/***************************************************************/
int RWGGeometry::GetDimension()
{ return TotalBFs; }

/***************************************************************/
/* Given an index ei into the overall list of edges, figure    */
/* out which object edge #ei belongs to and get its index      */
/* within the list of edges for that object. no error checking */
/* to determine if ei is a valid edge index.                   */
/***************************************************************/
#if 0
int RWGGeometry::GetObjectAndEdgeIndex(int ei, RWGObject **pO)
{ 
  RWGObject *O;
  int no;

  for (no=0; no<(NumObjects-1); no++)
   if ( ei<EdgeIndexOffset[no+1] ) 
    break;
  if (pO) *pO=Objects[no];
  return ei-EdgeIndexOffset[no];
} 
#endif
