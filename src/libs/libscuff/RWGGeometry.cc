/*
 * RWGGeometry.cc -- implementation of some methods in the RWGGeometry class
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

#define MAXSTR 1000
#define MAXTOK 50

/***********************************************************************/
/* subroutine to parse the MEDIUM...ENDMEDIUM section in a .scuffgeo   */
/* file. (currently the only keyword supported for this section is     */
/* MATERIAL xx).                                                       */
/***********************************************************************/
void ProcessMediumSectionInFile(FILE *f, char *FileName, int *LineNum, 
                                char *ExteriorMPName)
{
  char Line[MAXSTR];
  int NumTokens;
  char *p, *Tokens[MAXTOK];
  ExteriorMPName[0]=0;
  while( fgets(Line,MAXSTR,f) )
   { 
     (*LineNum)++;
     NumTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

     if ( !strcasecmp(Tokens[0],"MATERIAL") )
      {
        if (NumTokens!=2)
         ErrExit("%s:%i: syntax error",FileName,*LineNum);
        strncpy(ExteriorMPName,Tokens[1],MAXSTR);
      }
     else if ( !strcasecmp(Tokens[0],"ENDMEDIUM") )
      { 
        return;
      }
     else
      {
        ErrExit("%s:%i: unknown keyword %s",FileName,*LineNum,Tokens[0]);
      };
     
   };

  ErrExit("%s: unexpected end of file",FileName);

}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
RWGGeometry::RWGGeometry(const char *pGeoFileName)
{ 
  /***************************************************************/
  /* NOTE: i am not sure where to put this. put it here for now. */
  /***************************************************************/
  MatProp::SetLengthUnit(1.0e-6);
   
  /***************************************************************/
  /* storage for material properties defined on-the-fly in the   */
  /* .scuffgeo file.                                             */
  /* minor garbage-collection issue: MPs allocated in this       */
  /* routine are never freed.                                    */
  /***************************************************************/
  MatProp *MP, **MPs=0;
  int NumMPs;

  /***************************************************************/
  /* initialize simple fields ************************************/
  /***************************************************************/
  NumObjects=TotalBFs=TotalPanels=0;
  GeoFileName=strdup(pGeoFileName);
  ExteriorMP=0;
  Objects=0;

  /***************************************************************/
  /* try to open input file **************************************/
  /***************************************************************/
  FILE *f=fopen(GeoFileName,"r");
  if (!f)
   ErrExit("could not open %s",GeoFileName);

  /***************************************************************/
  /* read and process lines from input file one at a time        */
  /***************************************************************/
  RWGObject *O;
  char Line[MAXSTR], Label[MAXSTR];
  int LineNum=0; 
  int nt, nTokens;
  char *p, *Tokens[MAXTOK];
  char ExteriorMPName[MAXSTR];
  while( fgets(Line,MAXSTR,f) )
   { 
     LineNum++;

     /*--------------------------------------------------------------*/
     /*- break up line into tokens; skip blank lines and comments ---*/
     /*--------------------------------------------------------------*/
     nTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( nTokens==0 || Tokens[0][0]=='#' )
      continue; 
    
     /*--------------------------------------------------------------*/
     /*- switch off based on first token on the line ----------------*/
     /*--------------------------------------------------------------*/
     if ( !strcasecmp(Tokens[0],"MEDIUM") )
      { 
        ProcessMediumSectionInFile(f,GeoFileName,&LineNum,ExteriorMPName);
      }
     else if ( !strcasecmp(Tokens[0],"MATERIAL") )
      {
        /*--------------------------------------------------------------*/
        /* hand off to MatProp class constructor to parse this section  */
        /*--------------------------------------------------------------*/
        if ( nTokens==1 )
         ErrExit("%s:%i: no name given for MATERIAL ",GeoFileName,LineNum);
        else if ( nTokens>2 )
         ErrExit("%s:%i: syntax error",GeoFileName,LineNum);
         
        //MP = new MatProp(f,Label);
        if (MP->ErrMsg)
         ErrExit("%s:%i: %s",GeoFileName,LineNum,MP->ErrMsg); 

        NumMPs++;
        MPs=(MatProp **)realloc(MPs, NumMPs*sizeof(MatProp *));
        MPs[NumMPs-1]=MP;
           
      }
     else if ( !strcasecmp(Tokens[0],"OBJECT") )
      { 
        /*--------------------------------------------------------------*/
        /* hand off to RWGObject class constructor to parse this section*/
        /*--------------------------------------------------------------*/
        if ( nTokens>2 )
         ErrExit("%s:%i: syntax error",GeoFileName,LineNum);
        else if ( nTokens==2 )
         O=new RWGObject(f,Tokens[1],&LineNum);
        else if ( nTokens==1 )
         { snprintf(Label,MAXSTR,"Object_%i",NumObjects+1);
           O=new RWGObject(f,Label,&LineNum);
         };

        if (O->ErrMsg)
         ErrExit("%s:%i: %s",GeoFileName,LineNum,O->ErrMsg); 

        NumObjects++;
        Objects=(RWGObject **)realloc(Objects, NumObjects*sizeof(RWGObject *) );
        Objects[NumObjects-1]=O;

        TotalBFs+=O->NumBFs;
        TotalPanels+=O->NumPanels;
        O->Index=NumObjects-1;
      }
     else 
      { 
        /*--------------------------------------------------------------*/
        /* unknown keyword                                              */
        /*--------------------------------------------------------------*/
        ErrExit("%s:%i: syntax error",GeoFileName,LineNum);
      };

   }; // while( fgets(Line,MAXSTR,f) )

  /***************************************************************/
  /* process material properties of exterior medium             */
  /***************************************************************/
  int nmp;
  int no, nop;
  if ( strlen(ExteriorMPName)>0 )
   {
     for(nmp=0; ExteriorMP==0 && nmp<NumMPs; nmp++)
      if ( !strcasecmp(ExteriorMPName, MPs[nmp]->Name) )
       ExteriorMP=MPs[nmp];

     if (ExteriorMP==0)
      { ExteriorMP = new MatProp(O->MPName);
        if (ExteriorMP->ErrMsg)
        ErrExit("%s: medium: error in MATERIAL value: %s",GeoFileName,O->MP->ErrMsg);
      };
   }
  else
   ExteriorMP=new MatProp(MP_VACUUM);

  /***************************************************************/
  /* process material properties of objects                      */
  /***************************************************************/
  for(no=0; no<NumObjects; no++)
   { 
     O=Objects[no];

     if ( O->MPName )
      { 
        // deallocate the default MP created by the RWGObject constructor
        delete O->MP;

        /* look first to see if the requested material was one of */
        /* the materials defined on-the-fly in the .scuffgeo file */
        for(nmp=0; O->MP==0 && nmp<NumMPs; nmp++)
         if ( !strcasecmp(O->MPName, MPs[nmp]->Name) )
          O->MP=MPs[nmp];

        /* otherwise ... */ 
        if (O->MP==0)
         { O->MP = new MatProp(O->MPName);
           if (O->MP->ErrMsg)
            ErrExit("%s: object %s (%s): error in MATERIAL value: %s",
                     GeoFileName,O->Label,O->MeshFileName,O->MP->ErrMsg); 
         };

        free(O->MPName);
      }
     else
      O->MP = new MatProp(MP_PEC);
   };

  /***************************************************************/
  /* process object nesting relationships                        */
  /***************************************************************/
  for(no=0; no<NumObjects; no++)
   { 
     O=Objects[no];

     if ( O->ContainingObjectLabel )
      { 
        /* look for an object appearing earlier in the .scuffgeo file */
        /* whose label matches the requested object label             */
        for(nop=0; O->ContainingObject==0 && nop<no; nop++)
         if ( !strcasecmp(O->ContainingObjectLabel, Objects[nop]->Label ) )
          O->ContainingObject=Objects[nop];

        /* if no matching object was found, it's an error */
        if (O->ContainingObject==0)
         { 
           for(nop=no+1; nop<NumObjects; nop++)
            if ( !strcasecmp(O->ContainingObjectLabel, Objects[nop]->Label ) )
             ErrExit("%s: object %s: containing object %s must appear earlier in file",
                     GeoFileName, O->Label, O->ContainingObjectLabel);

           ErrExit("%s: object %s: containing object %s not found",
                    GeoFileName, O->Label, O->ContainingObjectLabel);
         };

        free( O->ContainingObjectLabel );

      };
   };
 
  /*******************************************************************/
  /* compute average panel area for statistical bookkeeping purposes */
  /*******************************************************************/
  AveragePanelArea=0.0; 
  int np;
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
  for(no=1; no<NumObjects; no++)
   { Mate[no]=-1;
     for(nop=0; nop<no && Mate[no]==-1; nop++)
      if (    !strcmp(Objects[no]->MeshFileName, Objects[nop]->MeshFileName)
           && !strcmp(Objects[no]->MP->Name    , Objects[nop]->MP->Name)
         ) 
       Mate[no]=nop;
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
