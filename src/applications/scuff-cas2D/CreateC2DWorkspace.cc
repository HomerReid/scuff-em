/*
 * CreateC2DWorkspace.cc  -- create and initialize a Casimir2D workspace
 *
 * homer reid  -- 10/2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctype.h>

#include <libhrutil.h>

#include "libTDRT.h"
#include "scuff-cas2D.h"

#define MAXTRANS 1000
#define MAXSTR   1000

/***************************************************************/
/* allocate, initialize, and return a pointer to a new         */
/* C2DWorkspace structure for subsequent passage to XQIntegrand*/
/***************************************************************/
C2DWorkspace *CreateC2DWorkspace(TDRTGeometry *G, char *TransListName,
                                 int WhichQuantities, double *Rectangle, 
                                 int NumThreads, int TETM, 
                                 int GroundPlane,
                                 int WriteHDF5, int IntCache, 
				 int VisualizeOnly)
{ 
  C2DWorkspace *W;

  char buffer[MAXSTR], TransLine[MAXSTR], Tag[MAXSTR], ErrMsg[MAXSTR];
  char *GeoFileBase;
  char PPFileName[1000];
  FILE *TransFile;
  char *p, *TransLines[MAXTRANS], *Tags[MAXTRANS];
  int nc, no, nop, NO, NBF, NBFp, N, N1, nt, nq, NTNQ;
  int LineNum;

  /***************************************************************/
  /* 1. allocate memory for the workspace structure and          */
  /* initialize the simple fields.                               */
  /***************************************************************/
  W=(C2DWorkspace *)malloc(sizeof(C2DWorkspace));
  W->G=G;
  W->TransListName=strdup(TransListName);
  N=W->N=G->TotalBFs;
  N1=G->Objects[0]->NumBFs;
  W->WhichQuantities=WhichQuantities;
  W->FixedXi=W->FixedQ=-1.0;
  W->TETM=TETM;
  W->GroundPlane=GroundPlane;
  W->WriteHDF5=WriteHDF5;
  W->IntCache=IntCache;
  W->NumContributingIVs=0;
  W->ProfileFileName=0;

  if (NumThreads==0)
   NumThreads=GetNumThreads();
  W->NumThreads=NumThreads;

  //snprintf(buffer,1000,"%s.byXQ",GetFileBase(G->GeoFileName));
  //f=CreateUniqueFile(buffer,0,buffer);
  //fclose(f);
  snprintf(buffer,1000,"%s.byXQ",GetFileBase(G->GeoFileName));
  W->ByXQFileName=strdup(buffer);
  snprintf(buffer,1000,"%s.byXi",GetFileBase(G->GeoFileName));
  W->ByXiFileName=strdup(buffer);

  /***************************************************************/
  /* 2. count the number of quantities requested. ****************/
  /***************************************************************/
  nq=0;
  if ( WhichQuantities & QUANTITY_ENERGY ) nq++;
  if ( WhichQuantities & QUANTITY_XFORCE ) nq++;
  if ( WhichQuantities & QUANTITY_YFORCE ) nq++;
  W->NumQuantities=nq;
  if (TETM) W->NumQuantities *= 2;

  /***************************************************************/
  /* 3. create image objects if a ground plane is present ********/
  /***************************************************************/
  if (W->GroundPlane)
   { W->ImageObjects=(TDRTObject **)malloc(G->NumObjects*sizeof(TDRTObject *));
     for(no=0; no<G->NumObjects; no++)
      W->ImageObjects[no]=CreateImageObject(G->Objects[no]);
   };

  /***************************************************************/
  /* 4. process transformation file to extract transformation    */
  /* lines and their associated tags                             */
  /***************************************************************/
  if ( !(TransFile=fopen(TransListName,"r")) )
   ErrExit("cannot open file %s",TransListName);

#ifdef _WIN32
  mkdir("gpmesh");
  mkdir("ppmesh");
#else
  mkdir("gpmesh",0755);
  mkdir("ppmesh",0755);
#endif

  GeoFileBase=GetFileBase(G->GeoFileName);

  sprintf(PPFileName,"ppmesh/%s.pp",GeoFileBase);
  unlink(PPFileName);
  
  nt=0;
  LineNum=0;
  while( fgets(TransLine,1000,TransFile) )
   { 
     LineNum++;

     /*--------------------------------------------------------------*/
     /* skip blank lines and comments                                */
     /*--------------------------------------------------------------*/
     p=TransLine; 
     while( isspace(*p) )
      p++;
     if ( *p==0 || *p=='#' )
      continue;

     /*--------------------------------------------------------------*/
     /* transform geometry                                           */
     /*--------------------------------------------------------------*/
     if ( G->Transform(TransLine,Tag,ErrMsg) )
      ErrExit("%s:%i: %s",TransListName,LineNum,ErrMsg);

     /*--------------------------------------------------------------*/
     /* transform images if present                                  */
     /*--------------------------------------------------------------*/
     if ( W->GroundPlane ) 
      { 
        if ( G->ObjectRotated )
         ErrExit("%s:%i: transforms may not include rotations when --GroundPlane is present");
        for(no=0; no<G->NumObjects; no++)
         { W->ImageObjects[no]->UnDisplace();
           W->ImageObjects[no]->Displace( G->Objects[no]->Displacement[0],
                                          -G->Objects[no]->Displacement[1]);
         };
      };

     /*--------------------------------------------------------------*/
     /*- write visualization files ----------------------------------*/
     /*--------------------------------------------------------------*/
 //    G->WriteGPMesh("gpmesh/%s.%s",GeoFileBase,Tag);
     G->WriteGPMeshPlus("gpmesh/%s.%s",GeoFileBase,Tag);
     G->WritePPMesh(PPFileName,Tag);

     if ( W->GroundPlane ) 
      WriteImageObjectPPMeshes(W,PPFileName,Tag);

     /*--------------------------------------------------------------*/
     /*- undo transformations ---------------------------------------*/
     /*--------------------------------------------------------------*/
     G->UnTransform();

     if ( W->GroundPlane )
      for(no=0; no<G->NumObjects; no++)
       W->ImageObjects[no]->UnDisplace();

     /*--------------------------------------------------------------*/
     /*- copy transformation tag ------------------------------------*/
     /*--------------------------------------------------------------*/
     TransLines[nt]=strdup(TransLine);
     Tags[nt]=strdup(Tag);

     /* prune any excess white space from the end of the tag */
     nc=strlen(Tags[nt]) - 1;
     while (nc>=0 && isspace(Tags[nt][nc]))
      Tags[nt][nc--]=0;
   
     nt++;
     if ( nt==MAXTRANS )
      ErrExit("too many transforms");

   };
  if ( nt==0 ) 
   ErrExit("no valid transformations found in file %s",TransListName);

  printf("Found %i valid transformations in %s.\n",nt,TransListName);

  W->NumTransforms=nt;

  W->TransLines=(char **)malloc(nt*sizeof(char *));
  memcpy(W->TransLines,TransLines,nt*sizeof(char *));

  W->Tags=(char **)malloc(nt*sizeof(char *));
  memcpy(W->Tags,Tags,nt*sizeof(char *));

  NTNQ=W->NTNQ=W->NumTransforms * W->NumQuantities;
  W->Converged=(int *)malloc(NTNQ*sizeof(int));
  memset(W->Converged,0,NTNQ*sizeof(int));

  /***************************************************************/
  /* 4. allocate matrices and vectors ****************************/ 
  /***************************************************************/
  /* T block matrices */
  W->T=(HMatrix **)malloc(G->NumObjects*sizeof(HMatrix *));
  for(no=0; no<G->NumObjects; no++)
   { if ( (nop=G->Mate[no]) > -1 )
      W->T[no]=W->T[nop];
     else
      { NBF=G->Objects[no]->NumBFs;
        W->T[no]=new HMatrix(NBF, NBF, LHM_COMPLEX, LHM_HERMITIAN);
      };
   };

  /* TI block matrices */
  if ( W->GroundPlane )
   {
     W->TI=(HMatrix **)malloc(G->NumObjects*sizeof(HMatrix *));
     for(no=0; no<G->NumObjects; no++)
      { NBF=G->Objects[no]->NumBFs;
        W->TI[no]=new HMatrix(NBF,NBF,LHM_COMPLEX,LHM_HERMITIAN);
      };
   };

  /* U block matrices. Note Uab is only allocated for b>a. */
  W->Uab=(HMatrix ***)malloc(G->NumObjects*sizeof(HMatrix **));
  for(no=0; no<G->NumObjects; no++)
   W->Uab[no]=(HMatrix **)malloc(G->NumObjects*sizeof(HMatrix *));
  for(no=0; no<G->NumObjects; no++)
   for(nop=no+1; nop<G->NumObjects; nop++)
    { NBF=G->Objects[no]->NumBFs;
      NBFp=G->Objects[nop]->NumBFs;
      W->Uab[no][nop]=new HMatrix(NBF, NBFp, LHM_COMPLEX);
    };

  // the BEM matrix is actually hermitian, but we do not take 
  // advantage of this by using lapack's compressed-storage
  // format for hermitian matrices, because the linear     
  // algebra routines run much more slowly in that case.   
  W->M=new HMatrix(N,N,LHM_COMPLEX,LHM_NORMAL);

  W->dU0bdX=(HMatrix **)malloc(G->NumObjects*sizeof(HMatrix *));
  memset(W->dU0bdX,0,G->NumObjects*sizeof(cdouble *));

  W->dU0bdY=(HMatrix **)malloc(G->NumObjects*sizeof(HMatrix *));
  memset(W->dU0bdY,0,G->NumObjects*sizeof(cdouble *));

  /*--------------------------------------------------------*/
  /*- dT0dY is the derivative of the T0 block of the matrix */
  /*- wrt y-displacements of object 0. this only exists if  */
  /*- GroundPlane==1 _and_ we are computing the y force.    */
  /*--------------------------------------------------------*/
  W->dT0dY=0;

  /*--------------------------------------------------------*/
  /* 4b. items that are only needed for energy calculation  */
  /*--------------------------------------------------------*/
  if ( WhichQuantities & QUANTITY_ENERGY )
   { 
     W->ipiv=(int *)malloc(N*sizeof(double));
     W->DRMInf=(double *)malloc(N*sizeof(double));
   };

  /*--------------------------------------------------------*/
  /* 4c. items that are only needed for force calculations  */
  /*--------------------------------------------------------*/
  if ( WhichQuantities & QUANTITY_ANYFORCE )
   { 
     NBF=G->Objects[0]->NumBFs;
 
     W->dM=new HMatrix(N,NBF,LHM_COMPLEX);

     if ( W->GroundPlane && (WhichQuantities & QUANTITY_YFORCE ) )
      W->dT0dY=new HMatrix(NBF, NBF, LHM_COMPLEX, LHM_HERMITIAN);

     for(nop=1; nop<G->NumObjects; nop++)
      { NBFp=G->Objects[nop]->NumBFs;
        if ( WhichQuantities & QUANTITY_XFORCE ) 
         W->dU0bdX[nop]=new HMatrix(NBF,NBFp,LHM_COMPLEX);
        if ( WhichQuantities & QUANTITY_YFORCE ) 
         W->dU0bdY[nop]=new HMatrix(NBF,NBFp,LHM_COMPLEX);
      };
   };

  /***************************************************************/
  /* 5a. create tables of static segment-segment integral data   */
  /*     for each object.                                        */
  /*                                                             */
  /* TSSSIDataTables[no]                                         */
  /*  = static segment-segment data table for object #no         */
  /***************************************************************/
  NO=G->NumObjects;
  if (!VisualizeOnly)
   { W->TSSSIDataTables=(StaticSSIDataTable **)malloc(NO*sizeof(StaticSSIDataTable *));
     for(no=0; no<NO; no++)
      { if (G->Mate[no]!=-1)
         W->TSSSIDataTables[no]=W->TSSSIDataTables[G->Mate[no]];
        else
         { if (TDRTGeometry::LogLevel>=2)
            Log("Creating static SSI data table for T[%i]...",no);
           W->TSSSIDataTables[no]=CreateStaticSSIDataTable(G->Objects[no], NumThreads);
         };
      };
   };

  /***************************************************************/
  /* 5b. create tables of static segment-segment integral data   */
  /*     for each pair of objects under each geometry transform. */
  /*                                                             */
  /* USSSIDataTables[nt][no + nop*NO]                            */
  /*    = static segment-segment data table for objects #no,#nop */
  /*      (no<nop) under geometry transform #nt                  */
  /***************************************************************/
  if (!VisualizeOnly)
   { W->USSSIDataTables=(StaticSSIDataTable ***)malloc(W->NumTransforms*sizeof(StaticSSIDataTable **));
     for(nt=0; nt<W->NumTransforms; nt++)
      { 
        G->Transform( W->TransLines[nt] );
   
        W->USSSIDataTables[nt]=(StaticSSIDataTable **)malloc(NO*NO*sizeof(StaticSSIDataTable **));
        for(no=0; no<NO; no++)
         for(nop=no+1; nop<NO; nop++)
          { 
            if ( nt>0 && G->ObjectMoved[no]==0 && G->ObjectMoved[nop]==0)
             continue;
   
            if (TDRTGeometry::LogLevel>=2)
             Log("Creating static SSI data table for U[%i,%i] at tag %s...",no,nop,W->Tags[nt]);
            W->USSSIDataTables[nt][no + nop*NO]
             =CreateStaticSSIDataTable(G->Objects[no],G->Objects[nop],NumThreads);
          };
   
        G->UnTransform();
      };
   };

  /***************************************************************/
  /* 6. handle user's request to retain contributions from basis */
  /*    functions associated with vertices lying inside a certain*/
  /*    rectangle.                                               */
  /***************************************************************/
  if (Rectangle)
   { 
     double XMin=Rectangle[0]; 
     double YMin=Rectangle[1]; 
     double XMax=Rectangle[2];
     double YMax=Rectangle[3];

     memcpy(W->Rectangle,Rectangle,4*sizeof(double));

     /*--------------------------------------------------------------*/    
     /*- first pass to count vertices in rectangle ------------------*/    
     /*--------------------------------------------------------------*/
     TDRTObject *O=G->Objects[0];
     double X, Y;
     int niv;
     for(niv=0; niv<O->NumIVs; niv++)
      { X=O->Vertices[ 2*O->IVs[niv]   ]; 
        Y=O->Vertices[ 2*O->IVs[niv]+1 ];
        if ( XMin<=X && X<=XMax && YMin<=Y && Y<=YMax )
         W->NumContributingIVs++;
      };  

     if (W->NumContributingIVs==0)
      ErrExit("specified rectangle contains no vertices");
     W->ContributingIVIndices=(int *)malloc(W->NumContributingIVs*sizeof(int));

     /*--------------------------------------------------------------*/
     /*- second pass to save indices of contributing IVs            -*/
     /*--------------------------------------------------------------*/
     W->NumContributingIVs=0;
     for(niv=0; niv<O->NumIVs; niv++)
      { 
        X=O->Vertices[ 2*O->IVs[niv]   ]; 
        Y=O->Vertices[ 2*O->IVs[niv]+1 ];
        if ( XMin<=X && X<=XMax && YMin<=Y && Y<=YMax )
         W->ContributingIVIndices[W->NumContributingIVs++]=niv;
      };

   };
   
  /***************************************************************/
  /* and that's it ***********************************************/
  /***************************************************************/
  return W;
    

}
