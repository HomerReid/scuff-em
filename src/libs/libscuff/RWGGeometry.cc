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

namespace scuff {

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
  char *Tokens[MAXTOK];
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
RWGGeometry::RWGGeometry(const char *pGeoFileName, int pLogLevel)
{ 
  /***************************************************************/
  /* NOTE: i am not sure where to put this. put it here for now. */
  /***************************************************************/
  MatProp::SetLengthUnit(1.0e-6);
   
  /***************************************************************/
  /* initialize simple fields ************************************/
  /***************************************************************/
  LogLevel=pLogLevel;
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
  RWGObject *O=NULL;
  char Line[MAXSTR], Label[MAXSTR];
  int LineNum=0; 
  int nTokens;
  char *Tokens[MAXTOK];
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
        ExteriorMPName[0]=0;
        ProcessMediumSectionInFile(f,GeoFileName,&LineNum,ExteriorMPName);
        if ( ExteriorMPName[0]!=0 )
         { ExteriorMP=new MatProp(ExteriorMPName);
           if ( ExteriorMP->ErrMsg )
            ErrExit("%s:%i: %s",GeoFileName,LineNum,ExteriorMP->ErrMsg);
         };
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
         
        char *ErrMsg=AddMaterialToMatPropDataBase(f, GeoFileName, Tokens[1], &LineNum);
        if (ErrMsg)
         ErrExit("%s:%i: %s",GeoFileName,LineNum,ErrMsg); 

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

  if (!ExteriorMP)
   ExteriorMP=new MatProp(MP_VACUUM);

  // Autodetect nesting relationships & topologically sort
  // (so that if A contains B, then B comes after A)
  for (int no = 0; no < NumObjects; ++no)
    Objects[no]->InitkdPanels(false, LogLevel);
  for (int no = 1; no < NumObjects; ++no) {
    int noi = no;
    O = Objects[no];
    for (int nop = no-1; nop >= 0; --nop) { // innermost to outermost
      if (Objects[nop]->Contains(O)) {
	O->ContainingObject = Objects[nop];
	break;
      }
      else if (O->Contains(Objects[nop])) {
	noi = nop; // O must go at noi (or earlier) to be in topo. order
      }
    }
    // insert object O at noi:
    for (int nop = no; nop > noi; --nop) {
	Objects[nop] = Objects[nop-1];
	if (Objects[nop]->ContainingObject == O->ContainingObject
	    && O->Contains(Objects[nop]))
	  Objects[nop]->ContainingObject = O;
    }
    Objects[noi] = O;
  }
  // check user ContainingObjectLabels (if any) for consistency:
  for (int no = 0; no < NumObjects; ++no) {
    O=Objects[no];
    if (O->ContainingObjectLabel && 
	(!O->ContainingObject 
	 || strcasecmp(O->ContainingObjectLabel, O->ContainingObject->Label)))
      ErrExit("%s: object %s: containing object is %s, not %s",
	      GeoFileName, O->Label,
	      O->ContainingObject ? O->ContainingObject->Label : "(none)",
	      O->ContainingObjectLabel);
    if (LogLevel >= SCUFF_VERBOSELOGGING && O->ContainingObject)
      Log("%s: object %s contained in object %s", GeoFileName,
	  O->Label, O->ContainingObject->Label);
    free(O->ContainingObjectLabel);
    O->ContainingObjectLabel = NULL;
  }
 
  /*******************************************************************/
  /* compute average panel area for statistical bookkeeping purposes */
  /*******************************************************************/
  AveragePanelArea=0.0; 
  for(int no=0; no<NumObjects; no++)
   for(int np=0; np<Objects[no]->NumPanels; np++)
    AveragePanelArea+=Objects[no]->Panels[np]->Area;
  AveragePanelArea/=((double) TotalPanels);

  /*******************************************************************/
  /* set the AllPEC flag based on whether or not all material objects*/
  /* are PEC bodies                                                  */
  /*******************************************************************/
  AllPEC=1;
  for(int no=0; no<NumObjects && AllPEC; no++)
   if ( !(Objects[no]->MP->IsPEC()) )
    AllPEC=0;

  /***************************************************************/
  /* initialize arrays of basis-function and panel index offsets */
  /***************************************************************/
  BFIndexOffset=(int *)mallocEC(NumObjects*sizeof(int) );
  PanelIndexOffset=(int *)mallocEC(NumObjects*sizeof(int) );
  BFIndexOffset[0]=PanelIndexOffset[0]=0;
  for(int no=1; no<NumObjects; no++)
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
  Mate=(int *)mallocEC(NumObjects*sizeof(int));
  Mate[0]=-1;
  for(int no=1; no<NumObjects; no++)
   { Mate[no]=-1;
     for(int nop=0; nop<no && Mate[no]==-1; nop++)
      if (    !strcmp(Objects[no]->MeshFileName, Objects[nop]->MeshFileName)
           && !strcmp(Objects[no]->MP->Name    , Objects[nop]->MP->Name)
         ) 
       Mate[no]=nop;
   };

  /***************************************************************/
  /* initialize ObjectMoved[] array.                             */
  /* the values of this array are only defined after             */
  /* a call to the RWGGeometry::Transform() function, when we    */
  /* have ObjectMoved[i]=1 if the ith object was modified by     */
  /* the transformation.                                         */
  /***************************************************************/
  ObjectMoved=(int *)mallocEC(NumObjects*sizeof(int));

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
/* return the object whose label is Label. if pno is non-NULL  */
/* on entry, then on return it is set to the index of the      */
/* object.                                                     */
/*                                                             */
/* note: if the return value is NULL, then there are three     */
/* possibilities:                                              */
/*  (a) the Label string was NULL on entry                     */
/*  (b) the Label string was "EXTERIOR" or "MEDIUM"            */
/*  (c) the Label string was some non-empty string that did    */
/*      not match the label of any object in the geometry.     */
/*                                                             */
/* you can tell which happened by looking at pno: in cases (a) */
/* and (b) (medium was specified) then *pno==-1, whereas in    */
/* case (c) (invalid object label) *pno==-2.                   */
/***************************************************************/
RWGObject *RWGGeometry::GetObjectByLabel(const char *Label, int *pno)
{
  if (pno) *pno=-2;
 
  if (Label==0)
   { if (pno) *pno=-1;
     return NULL;
   };
  
  for(int no=0; no<NumObjects; no++)
   if ( !strcasecmp(Label, Objects[no]->Label) ) 
    { if (pno) *pno = no;
      return Objects[no];
    }

  if (pno && (!strcasecmp(Label,"EXTERIOR") || !strcasecmp(Label,"MEDIUM")))
    *pno = -1;

  return NULL;
}

/***************************************************************/
/* Apply the specified GTComplex to transform the geometry.    */
/* (Note that a 'GTComplex' is a list of GTransformations, each*/
/* of which is applied to one specific object in the geometry.)*/
/***************************************************************/
void RWGGeometry::Transform(GTComplex *GTC)
{ 
  int noa, WhichObject;
  RWGObject *O;

  // assume that no objects will be modified by this operation
  memset(ObjectMoved, 0, NumObjects*sizeof(int));

  // loop over the individual transformations in the complex
  for(noa=0; noa<GTC->NumObjectsAffected; noa++)
   { 
     // find the object corresponding to the label for this transformation
     O=GetObjectByLabel(GTC->ObjectLabel[noa], &WhichObject);

     // apply the transformation to that object
     if (O) 
      { O->Transform(GTC->GT + noa);
        ObjectMoved[WhichObject]=1;
      };
        
   };

}

/***************************************************************/
/* Undo transformations. ***************************************/
/***************************************************************/
void RWGGeometry::UnTransform()
{ 
  int no;
  for(no=0; no<NumObjects; no++)
   Objects[no]->UnTransform();
}

/***************************************************************/
/* Quick sanity check to make sure that a given list of        */
/* GTComplex structures actually makes sense for the given     */
/* geometry, which is to say that it doesn't request           */
/* transformations on any objects that don't exist in the      */
/* geometry.                                                   */
/* Returns 0 if the check passed, or an error message if not.  */
/***************************************************************/
char *RWGGeometry::CheckGTCList(GTComplex **GTCList, int NumGTCs)
{
  int ngtc, noa;
  
  for(ngtc=0; ngtc<NumGTCs; ngtc++)
   for (noa=0; noa<GTCList[ngtc]->NumObjectsAffected; noa++)
    if (!GetObjectByLabel(GTCList[ngtc]->ObjectLabel[noa]))
     return vstrdup("transformation requested for unknown object %s",
                     GTCList[ngtc]->ObjectLabel[noa]);

  return 0;
}

/***************************************************************/
/* Return the dimension of the linear system. ******************/
/***************************************************************/
int RWGGeometry::GetDimension()
{ return TotalBFs; }

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::SetLogLevel(int NewLogLevel)
 { LogLevel=NewLogLevel; }

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::SetEpsMu(const char *Label, cdouble Eps, cdouble Mu)
{ 
  if ( Label==NULL || !strcasecmp(Label,"EXTERIOR") )
   ExteriorMP->SetEpsMu(Eps, Mu); 
  else
   { RWGObject *O=GetObjectByLabel(Label);
     if (O)
      O->MP->SetEpsMu(Eps, Mu); 
     else
      Warn("unknown object %s specified in SetEpsMu() (ignoring)",Label);
   };
} 

void RWGGeometry::SetEpsMu(cdouble Eps, cdouble Mu)
{ SetEpsMu(0, Eps, Mu); }

void RWGGeometry::SetEps(const char *Label, cdouble Eps)
{ SetEpsMu(Label, Eps, 1.0); }

void RWGGeometry::SetEps(cdouble Eps)
{ SetEpsMu(0, Eps, 1.0); }


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

/***************************************************************/
/***************************************************************/
/***************************************************************/
double RWGGeometry::SWPPITol = 1.0e-5;

} // namespace scuff
