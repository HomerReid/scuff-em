/*
 * StaticPPI.cc -- libscuff routines for dealing with static 
 *              -- panel-panel integrals 
 *
 * homer reid   -- 5/2009 -- 12/2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <pthread.h>

#include <libhrutil.h>
#include <libTriInt.h>

#include "libscuff.h"

#include "StaticPPI.h"

/*****************************************************************/
/*****************************************************************/
/* Compute static panel-panel integrals using fixed-order        */
/* numerical cubature for both panels.                           */
/*****************************************************************/
/*****************************************************************/
void ComputeStaticPPIData_Fixed(double *V1, double *V2, double *V3,
                                double *V1P, double *V2P, double *V3P,
                                int Order, StaticPPIData *SPPID)
                             
{ 
  int np, nqr, npp, nqrp, NumPts;
  int i, rp, nh;
  double u, v, w, up, vp, wp;
  double A[3], B[3], AP[3], BP[3];
  double X[3], XP[3];
  double *TCR;
  double r, rPow, XmXP[3], XxXP[3];

  double XmXPoR3Int_Inner[3];
  double XxXPoR3Int_Inner[3];
  double hrnInt_Inner[RPMAX+2][NUMHS];

  /***************************************************************/
  /* preliminary setup *******************************************/
  /***************************************************************/
  VecSub(V2, V1, A);
  VecSub(V3, V1, B);

  VecSub(V2P, V1P, AP);
  VecSub(V3P, V1P, BP);

  /***************************************************************/
  /* choose order of quadrature scheme to use ********************/
  /***************************************************************/
  TCR=GetTCR(Order, &NumPts);
/*
  switch(Order)
   { case 4:  TCR=TCR4;   NumPts=6;   break;
     case 7:  TCR=TCR7;   NumPts=15;  break;
     case 9:  TCR=TCR9;   NumPts=21;  break;
     case 16: TCR=TCR16;  NumPts=55;  break;
     case 20: TCR=TCR20;  NumPts=78;  break;
     case 25: TCR=TCR25;  NumPts=120; break;
     default: fprintf(stderr,"whoops!\n\n"); exit(1); };
*/

  /***************************************************************/
  /* outer loop **************************************************/
  /***************************************************************/
  memset(SPPID->XmXPoR3Int,0,3*sizeof(double));
  memset(SPPID->XxXPoR3Int,0,3*sizeof(double));
  for(rp=-1; rp<=RPMAX; rp++)
   memset(SPPID->hrnInt[rp+1],0,NUMHS*sizeof(double));

  for(np=nqr=0; np<NumPts; np++) 
   { 
     u=TCR[nqr++]; v=TCR[nqr++]; w=TCR[nqr++];

     /***************************************************************/
     /* set X         ***********************************************/
     /***************************************************************/
     for(i=0; i<3; i++)
      X[i]=V1[i] + u*A[i] + v*B[i];

     /***************************************************************/
     /* inner loop to calculate value of inner integrand ************/
     /***************************************************************/
     memset(XmXPoR3Int_Inner,0,3*sizeof(double));
     memset(XxXPoR3Int_Inner,0,3*sizeof(double));
     for(rp=-1; rp<=RPMAX; rp++)
      memset(hrnInt_Inner[rp+1],0,NUMHS*sizeof(double));

     for(npp=nqrp=0; npp<NumPts; npp++)
      { 
        up=TCR[nqrp++]; vp=TCR[nqrp++]; wp=TCR[nqrp++];

        /***************************************************************/ 
        /* set XP           ********************************************/
        /***************************************************************/
        for(i=0; i<3; i++)
         XP[i]=V1P[i] + up*AP[i] + vp*BP[i];

        /***************************************************************/
        /* inner integrands ********************************************/
        /***************************************************************/
        VecSub(X,XP,XmXP);
        VecCross(X,XP,XxXP);

        r=VecNorm(XmXP);

        if (r<=0.0)
         rPow=0.0;
        else
         rPow=1.0/(r*r*r);
  
        XmXPoR3Int_Inner[0] += wp * rPow * XmXP[0];
        XmXPoR3Int_Inner[1] += wp * rPow * XmXP[1];
        XmXPoR3Int_Inner[2] += wp * rPow * XmXP[2];
        XxXPoR3Int_Inner[0] += wp * rPow * XxXP[0];
        XxXPoR3Int_Inner[1] += wp * rPow * XxXP[1];
        XxXPoR3Int_Inner[2] += wp * rPow * XxXP[2];

        rPow*=r*r;
        for(rp=-1; rp<=RPMAX; rp++)
         { hrnInt_Inner[rp+1][0] += wp * rPow * 1.0;
           hrnInt_Inner[rp+1][1] += wp * rPow * (up+vp);
           hrnInt_Inner[rp+1][2] += wp * rPow * vp;
           hrnInt_Inner[rp+1][3] += wp * rPow * (u+v);
           hrnInt_Inner[rp+1][4] += wp * rPow * (u+v)*(up+vp); 
           hrnInt_Inner[rp+1][5] += wp * rPow * (u+v)*vp;
           hrnInt_Inner[rp+1][6] += wp * rPow * v;
           hrnInt_Inner[rp+1][7] += wp * rPow * v*(up+vp);
           hrnInt_Inner[rp+1][8] += wp * rPow * v*vp;
           rPow*=r;
         };

      }; /* end of inner cubature loop */

     /* accumulate contributions to outer integrals */
     for(i=0; i<3; i++)
      { SPPID->XmXPoR3Int[i] += w*XmXPoR3Int_Inner[i];
        SPPID->XxXPoR3Int[i] += w*XxXPoR3Int_Inner[i];
      };
     for(rp=-1; rp<=RPMAX; rp++)
      for(nh=0; nh<NUMHS; nh++)
       SPPID->hrnInt[rp+1][nh] += w*hrnInt_Inner[rp+1][nh];

   }; /* end of outer cubature loop */

}

/*****************************************************************/
/* compute static panel-panel integrals for a single pair of     */
/* panels. this is a switchboard routine that chooses between    */
/* one of two methods for computing static panel-panel integrals.*/
/*****************************************************************/
void ComputeStaticPPIData(RWGObject *O, int np, RWGObject *OP, int npp,
                          StaticPPIData *SPPID)
{
  int nCommon, Index[3], IndexP[3];
  RWGPanel *P, *PP;
  int iV1, iV2, iV3, iV1P, iV2P, iV3P;
  double *V1, *V2, *V3, *V1P, *V2P, *V3P;

  P=O->Panels[np];
  PP=OP->Panels[npp];

  /***************************************************************/
  /* count number of vertices in common between the two panels   */
  /***************************************************************/
  if ( O==OP )
   nCommon=O->CountCommonVertices(np,npp,Index,IndexP);
  else
   nCommon=0;

  /***************************************************************/
  /* if panels have 1 or 2 common vertices, compute static PPIs  */
  /* using taylor's method.                                      */ 
  /* Otherwise compute static PPIs using 20-point numerical      */
  /* cubature on both panels.                                    */
  /* Note we never compute static PPIs for panel pairs with 3    */
  /* common vertices because for that case we use taylor's       */
  /* schemes directly to get the full panel integrals.           */
  /***************************************************************/
  if ( nCommon==1 )
   { 
     iV1=Index[0];    V1=O->Vertices + 3*P->VI[iV1];
     iV2=(iV1+1)%3;   V2=O->Vertices + 3*P->VI[iV2];
     iV3=(iV1+2)%3;   V3=O->Vertices + 3*P->VI[iV3];

     iV1P=IndexP[0];
     iV2P=(iV1P+1)%3; V2P=O->Vertices + 3*PP->VI[iV2P];
     iV3P=(iV1P+2)%3; V3P=O->Vertices + 3*PP->VI[iV3P];

     ComputeStaticPPIData_Taylor(V1, V2, V3, V2P, V3P, SPPID);
   }
  else if ( nCommon==2 )
   { 
     iV1=Index[0];     V1=O->Vertices + 3*P->VI[iV1];
     iV2=Index[1];     V2=O->Vertices + 3*P->VI[iV2];
     iV3=3-iV1-iV2;    V3=O->Vertices + 3*P->VI[iV3];

     iV1P=IndexP[0];
     iV2P=IndexP[1];
     iV3P=3-iV1P-iV2P; V3P=O->Vertices + 3*PP->VI[iV3P];

     ComputeStaticPPIData_Taylor(V1, V2, V3, V2, V3P, SPPID);
   }
  else
   { 
     V1  =  O->Vertices +  3*P->VI[0];
     V2  =  O->Vertices +  3*P->VI[1];
     V3  =  O->Vertices +  3*P->VI[2];
     V1P = OP->Vertices + 3*PP->VI[0];
     V2P = OP->Vertices + 3*PP->VI[1];
     V3P = OP->Vertices + 3*PP->VI[2];

     ComputeStaticPPIData_Fixed(V1, V2, V3, V1P, V2P, V3P, 20, SPPID);

   };

}

/*******************************************************************/
/* thread data structure and thread routine used to multithreadize */
/* the second pass loop of CreateStaticPPIData()                   */
/*******************************************************************/
typedef struct ThreadData
 { RWGObject *O1, *O2;
   StaticPPIDataTable *SPPIDT;
   int nt, nThread;
 } ThreadData;

void *CreateStaticPPIDataTable_Thread(void *pTD)
{ 
  ThreadData *TD=(ThreadData *)pTD;
  RWGObject *O1, *O2;
  RWGPanel *P, *PP;
  StaticPPIDataTable *SPPIDT;
  StaticPPIDataMap::const_iterator it;
  int np, npp, nt;
  int nPairs, TotalPairs;

  O1=TD->O1;
  O2=TD->O2;
  SPPIDT=TD->SPPIDT;

  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/
#ifdef _GNU_SOURCE
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(TD->nt,&cpuset);
  pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#endif
  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/

  /*--------------------------------------------------------------*/
  /* second pass to fill in the actual static panel-panel data.   */
  /*--------------------------------------------------------------*/
  nt=0;
  nPairs=0;
  TotalPairs = (O1==O2) ? (O1->NumPanels)*(O1->NumPanels + 1)/2 : (O1->NumPanels*O2->NumPanels);
  for(np=0, P=O1->Panels[0]; np<O1->NumPanels; P=O1->Panels[++np])
   for(npp=(O1==O2 ? np : 0), PP=O2->Panels[npp]; npp<O2->NumPanels; PP=O2->Panels[++npp])
    { 
       nPairs++;

       nt++;
       if (nt==TD->nThread) nt=0;
       if (nt!=TD->nt)
        continue;

       if( nPairs == (TotalPairs/5) )
        Log(" ... 20%% ... ");
       if( nPairs == (2*TotalPairs/5) )
        Log(" ... 40%% ... ");
       if( nPairs == (3*TotalPairs/5) )
        Log(" ... 60%% ... ");
       if( nPairs == (4*TotalPairs/5) )
        Log(" ... 80%% ... ");

       it=SPPIDT->Maps[np]->find(npp);
       if ( it==(SPPIDT->Maps[np]->end()) )
        continue;

       ComputeStaticPPIData(O1,np,O2,npp,SPPIDT->Buffer + it->second);
    };

  pthread_exit(0);
  return 0;
}

/*****************************************************************/
/* CreateStaticPPIDataTable(): create a table of static          */
/* panel-panel integral data for a given pair of objects.        */
/*                                                               */
/* static panel-panel integral data are stored for each pair of  */
/* panels that are 'nearby' each other. 'nearby' is defined as   */
/* follows: if the distance between the panel centroids is less  */
/* than DESINGULARIZATION_RADIUS times the larger of the two     */
/* panel radii, then the panels are 'nearby' each other. (the    */
/* panel radius is the largest distance from panel centroid to   */
/* any panel vertex.)                                            */
/*                                                               */
/* all static panel-panel integral data for a given pair of      */
/* panels are contained in a structure called 'StaticPPIData.'   */
/*                                                               */  
/* a 'StaticPPIDataTable' is defined for a pair of objects and   */  
/* contains two components:                                      */  
/*  1) a big buffer containing storage for all the StaticPPIData */
/*     structures we will need for all pairs of nearby panels    */
/*  2) an array of map structures, one for each panel on         */
/*     the first object (O1), that map panel indices on the      */
/*     second object (O2) into indices into the big buffer.      */
/*                                                               */
/* thus, suppose panel #7 on O1 is nearby panel #9 on object O2, */
/* but not nearby panel #251 on object O2. then, if the call     */
/*  SPPITable->SPPIDataMaps[7]->find(9)                          */
/* returns 13, then that means that the static panel-panel data  */
/* for panel pair ( {01,7}, {O2,9} ) is stored in slot 13 of the */
/* big buffer, while the call                                    */
/*  SPPITable->SPPIDataMaps[7]->find(251)                        */
/* will fail.                                                    */
/*                                                               */
/* if O1 and O2 are the SAME object, then, to avoid redundancy,  */
/* the SPPIDataMaps table for panel #np will only contain entries*/
/* for panels #npp such that npp>=np. thus, if you are looking in*/
/* a table to see if static panel-panel data exist for a pair    */
/* of panels (np1,np2), make sure np1<=np2 and swap if not.      */
/*****************************************************************/
StaticPPIDataTable *CreateStaticPPIDataTable(RWGObject *O1, RWGObject *O2, int nThread)
{ 
  StaticPPIDataTable *SPPIDT;
  RWGPanel *P, *PP;
  int np, npp;
  int NumPairs;
  double rRel, DSR=DESINGULARIZATION_RADIUS;
 
  pthread_t Threads[nThread];
  ThreadData TDS[nThread], *TD;
  int nt;

  char FileName[500];
  FILE *f;

  /*--------------------------------------------------------------*/
  /* allocate memory for the structure itself and the array of    */
  /* data maps                                                    */ 
  /*--------------------------------------------------------------*/
  SPPIDT=(StaticPPIDataTable *)RWGMalloc( sizeof(SPPIDT[0]) );
  SPPIDT->SameObject= (O1==O2) ? 1 : 0; 
  SPPIDT->Maps=(StaticPPIDataMap **)RWGMalloc(O1->NumPanels*sizeof(SPPIDT->Maps[0]));
  for(np=0; np<O1->NumPanels; np++)
   { SPPIDT->Maps[np]=new StaticPPIDataMap;
     SPPIDT->Maps[np]->set_empty_key(-1);
     SPPIDT->Maps[np]->set_deleted_key(-2);
   };
  SPPIDT->NumMaps=O1->NumPanels;

  /*--------------------------------------------------------------*/
  /* first pass to compute the total number of pairs in the table */
  /* and fill in the data maps                                    */
  /*--------------------------------------------------------------*/
  NumPairs=0;
  for(np=0, P=O1->Panels[0]; np<O1->NumPanels; P=O1->Panels[++np])
   for(npp=(O1==O2 ? np : 0), PP=O2->Panels[npp]; npp<O2->NumPanels; PP=O2->Panels[++npp])
    if (VecDistance(P->Centroid,PP->Centroid) < DSR*fmax(P->Radius,PP->Radius) )
     SPPIDT->Maps[np]->insert( std::pair<int,int>(npp,NumPairs++) );
  SPPIDT->NumPairs=NumPairs;

  /*--------------------------------------------------------------*/
  /*- now allocate space for the actual panel-panel data records  */
  /*--------------------------------------------------------------*/
  SPPIDT->Buffer=(StaticPPIData *)RWGMalloc(NumPairs*sizeof(StaticPPIData));

  /*--------------------------------------------------------------*/
  /*- first attempt to read data from file.                       */
  /*-                                                             */
  /*- static PPI data for an object named Sphere456.msh are stored*/
  /*- in a binary data file called Sphere456.SPPID.               */
  /*-                                                             */
  /*- we look in a couple of places for this file:                */
  /*-  a. current working directory                               */
  /*-  b. $(HOME)/geomsh/SPPIData                                 */
  /*-                                                             */
  /*- the binary file format is:                                  */
  /*-  NumMaps         (sizeof(int))                              */
  /*-  NumPairs        (sizeof(int))                              */
  /*-  the actual data (NumPairs*sizeof(StaticPPIData))           */
  /*--------------------------------------------------------------*/
  if (O1==O2)
   { 
     int NumMaps, NumPairs, SizeRead, KeepTrying;

     KeepTrying=1;

     snprintf(FileName,500,"%s.SPPID",GetFileBase(O1->MeshFileName));
     f=fopen(FileName,"r");
     if (!f)
      { snprintf(FileName,500,"%s/geomsh/SPPIData/%s.SPPID", getenv("HOME"),
                               GetFileBase(O1->MeshFileName));
        f=fopen(FileName,"r");
      };

     if (!f)
      { KeepTrying=0;
        f=0;
      };
     Log("%s: attempting to read static PPI data from file %s...",
          O1->MeshFileName,FileName);
   
     if ( KeepTrying )
      { if (    fread(&NumMaps, sizeof(int), 1, f) != 1
             || fread(&NumPairs, sizeof(int), 1, f) != 1 )
         { Log("...failed! invalid file format",FileName);
           KeepTrying=0;
         };
      }

     if ( KeepTrying )
      { if ( NumMaps != SPPIDT->NumMaps || NumPairs != SPPIDT->NumPairs )
         { Log("...failed! incorrect data table size",FileName);
           KeepTrying=0;
         };
      };

     if ( KeepTrying )
      { SizeRead=fread(SPPIDT->Buffer,sizeof(StaticPPIData),NumPairs,f);
        if (SizeRead!=NumPairs)
         { Log("...failed! incorrect data table size",FileName);
           KeepTrying=0;
         };
      };
  
     if (f) fclose(f);

     if (KeepTrying)
      { Log("...success! (%i pairs, %i MB of RAM used)",
             NumPairs,(NumPairs*sizeof(StaticPPIData))/1048576);
        return SPPIDT;
      };
         
   };
   
  /*--------------------------------------------------------------*/
  /* if file import failed, do a second pass to compute the       */
  /* static panel-panel data.                                     */
  /*--------------------------------------------------------------*/
  if (O1==O2)
   Log("creating static PPI data table for object %s",O1->MeshFileName);
  else
   Log("creating static PPI data table for objects %s and %s",
        O1->MeshFileName, O2->MeshFileName);
 
  for(nt=0; nt<nThread; nt++)
   {
     TD=&(TDS[nt]);
     TD->nt=nt;
     TD->nThread=nThread;
     TD->O1=O1;
     TD->O2=O2;
     TD->SPPIDT=SPPIDT;
     pthread_create( &(Threads[nt]), 0, CreateStaticPPIDataTable_Thread, (void *)TD);
   };

  for(nt=0; nt<nThread; nt++)
   pthread_join(Threads[nt],0);
  Log("...done!");

  /*--------------------------------------------------------------*/
  /*- write the table to disk (FIXME for object-object tables)    */
  /*--------------------------------------------------------------*/
  snprintf(FileName,500,"%s.SPPID",GetFileBase(O1->MeshFileName));
  WriteStaticPPIDataTable(SPPIDT, FileName);

  return SPPIDT;

}

/*****************************************************************/
/* look up and return a pointer to the StaticPPIData structure   */
/* for a given pair of panels.                                   */
/*                                                               */
/* SPPID must point to a user-allocated StaticPPIData structure. */
/* On return, if static PPI data for the given pair of panels    */
/* was found, then this structure is filled in the data and      */
/* SPPID itself is returned. If no static PPI data were found,   */
/* the return value is NULL and the structure to which SPPID     */
/* points is untouched.                                          */
/*                                                               */
/* Note that, if the StaticPPIData structure is storing data     */ 
/* for pairs of panels on the same object, then, to avoid        */
/* redundancy, we only compute and store data for panel pairs    */
/* (np1,np2) in which np1<=np2. Thus, if this routine is called  */
/* with np1<np2, we need to flip np1 and np2, and we also need   */
/* to flip primed and unprimed quantities in the resulting SPPID */
/* structure (if found).                                         */
/*****************************************************************/
#define SWAP(A,B,C) { (C)=(A); (A)=(B); (B)=(C); }
StaticPPIData *GetStaticPPIData(StaticPPIDataTable *SPPIDT, int np1, int np2,
                                StaticPPIData *SPPID)
{ 
  StaticPPIDataMap::const_iterator it;
  int Flipped, iTemp, rp;
  double dTemp;
 
  Flipped=0;
  if ( SPPIDT->SameObject && np1>np2 )
   { SWAP(np1, np2, iTemp); Flipped=1;
   };
  
  it=SPPIDT->Maps[np1]->find(np2);
   if ( it == (SPPIDT->Maps[np1]->end()) )
    return 0;

  memcpy(SPPID, SPPIDT->Buffer + it->second, sizeof(StaticPPIData) );

  if ( Flipped )
   { 
     VecScale(SPPID->XmXPoR3Int,-1.0);
     VecScale(SPPID->XxXPoR3Int,-1.0);
     for(rp=-1; rp<=RPMAX; rp++)
      { 
        SWAP(SPPID->hrnInt[rp+1][1], SPPID->hrnInt[rp+1][3], dTemp );
        SWAP(SPPID->hrnInt[rp+1][2], SPPID->hrnInt[rp+1][6], dTemp );
        SWAP(SPPID->hrnInt[rp+1][5], SPPID->hrnInt[rp+1][7], dTemp );
      };
   }; 

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//if( (np1==30 && np2==38) || (np1==38 && np2==30) || (np1==45 && np2==46) || (np1==46 && np2==45))
// printf(" PPI(%i,%i): (%e, %e, %e)\n",np1,np2,SPPID->
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  return SPPID;

}

/*****************************************************************/
/* write a StaticPPIDataTable to disk.                           */
/*  (binary) file format:                                        */
/*   a. first sizeof(int) bytes = NumMaps                        */
/*   b. next sizeof(int) bytes = NumPairs                        */
/*   c. remainder of file = contents of Buffer                   */
/*****************************************************************/
void WriteStaticPPIDataTable(StaticPPIDataTable *SPPIDT, const char *FileName)
{ 
  FILE *f;
  
  if ( !(f=fopen(FileName,"w") ) )
   { fprintf(stderr,"\n*\n* WARNING: could not write to file %s\n*\n",FileName);
     return;
   };

  fwrite(&SPPIDT->NumMaps,sizeof(int),1,f);
  fwrite(&SPPIDT->NumPairs,sizeof(int),1,f);
  fwrite(SPPIDT->Buffer,sizeof(StaticPPIData),SPPIDT->NumPairs,f);

  fclose(f);

}
