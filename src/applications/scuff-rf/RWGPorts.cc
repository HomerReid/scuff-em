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
 * RWGPorts.cc -- libRWG module for handling of 'ports' in 
 *             -- RF/microwave structures
 *                                         
 * homer reid  -- 9/2011                                        
 *                                         
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>
#include <libscuffInternals.h>
#include <libSGJC.h>

#include "RWGPorts.h"

#define ABSTOL 1.0e-8
#define RELTOL 1.0e-4
#define FREQ2OMEGA (2.0*M_PI/300.0)
#define II cdouble(0.0,1.0)

using namespace scuff;

/***************************************************************/
/* slightly tricky: the panels associated with what we call the*/
/* 'positive edge' of an RWGPort are actually the *negative*   */
/* of the two panels that define an RWG basis function         */
/* straddling the port (in the sense that on these panels the  */
/* current flows *toward* the Q vertex and is sunk there, as   */
/* opposed to being sourced from there). this confounds the    */
/* usual sense of positive and negative weights for the        */
/* contributions of the two panels that make up the RWG basis  */
/* function.                                                   */  
/*                                                             */  
/* to be precise: the panels indexed by the 'PPanelIndices'    */  
/* array in the Port structure are the panels on the positive  */  
/* side of the port, which makes them the negative of the two  */  
/* panels in the RWG basis function to which they belong.      */  
/***************************************************************/
RWGPort *CreatePort(RWGObject *PObject, int NumPEdges, int *PEdgeIndices, 
                    RWGObject *MObject, int NumMEdges, int *MEdgeIndices)
{
  RWGPort *P=(RWGPort *)malloc(sizeof *P);
  RWGEdge *E;

  P->PObject       = PObject;
  P->NumPEdges     = NumPEdges;
  P->PPanelIndices = (int *)malloc(NumPEdges*sizeof(int));
  P->PPaneliQs     = (int *)malloc(NumPEdges*sizeof(int));
  P->PLengths      = (double *)malloc(NumPEdges*sizeof(double));

  P->MObject       = MObject;
  P->NumMEdges     = NumMEdges;
  P->MPanelIndices = (int *)malloc(NumMEdges*sizeof(int));
  P->MPaneliQs     = (int *)malloc(NumMEdges*sizeof(int));
  P->MLengths      = (double *)malloc(NumMEdges*sizeof(double));
  
  int ne, EI;

  /***************************************************************/
  /* gather info on the RWG edges that define the 'positive' side*/
  /* of the port                                                 */
  /***************************************************************/
  P->PPerimeter=0.0;
  for(ne=0; ne<NumPEdges; ne++)
   { 
     EI = PEdgeIndices[ne];
     if ( EI<0 || EI>PObject->NumExteriorEdges )
      ErrExit("object %s(%s): invalid exterior edge index %i",
               PObject->Label,PObject->MeshFileName,EI);
     E=PObject->ExteriorEdges[EI];
     P->PPanelIndices[ne]=E->iPPanel;
     P->PPaneliQs[ne]=E->PIndex;
     P->PLengths[ne]=E->Length;
     P->PPerimeter+=E->Length;

     if (ne==0) 
      memcpy(P->PRefPoint, E->Centroid, 3*sizeof(double));
   };

  /***************************************************************/
  /* gather info on the RWG edges that define the 'negative' side*/
  /* of the port                                                 */
  /***************************************************************/
  P->MPerimeter=0.0;
  for(ne=0; ne<NumMEdges; ne++)
   {    
     EI=MEdgeIndices[ne]; // 'minus edge index'
     if ( EI<0 || EI>MObject->NumExteriorEdges )
      ErrExit("object %s(%s): invalid exterior edge index %i",
               MObject->Label,MObject->MeshFileName,EI);
     E=MObject->ExteriorEdges[EI];
     P->MPanelIndices[ne]=E->iPPanel;
     P->MPaneliQs[ne]=E->PIndex;
     P->MLengths[ne]=E->Length;
     P->MPerimeter+=E->Length;

     if (ne==0) 
      memcpy(P->MRefPoint, E->Centroid, 3*sizeof(double));
   };

  return P;

} 

/***************************************************************/
/* port file syntax example                                    */
/*  PORT                                                       */
/*   PEDGES pe1 pe2 ... peN                                    */
/*   MEDGES me1 me2 ... meN                                    */
/*   PREFPOINT x1 x2 x3 <optional>                             */
/*   MREFPOINT x1 x2 x3 <optional>                             */
/*   POBJECT   FirstObject                                     */
/*   MOBJECT   SecondObject                                    */
/*  ENDPORT                                                    */
/***************************************************************/
#define MAXEDGES 100 // max # of edges in (each half of) a port 
RWGPort **ParsePortFile(RWGGeometry *G, 
                        const char *PortFileName, 
                        int *pNumPorts)
{
  /***************************************************************/
  /* try to open the file ****************************************/
  /***************************************************************/
  FILE *f=fopen(PortFileName,"r");
  if (f==0)
   ErrExit("could not open file %s",PortFileName);
  
  /***************************************************************/
  /* read through the file and parse lines one-at-a-time *********/
  /***************************************************************/
  RWGPort **PortArray=0;
  int NumPorts=0;
  int NumPEdges=0, NumMEdges=0;
  int PEIndices[MAXEDGES], MEIndices[MAXEDGES];
  double PRefPoint[3], MRefPoint[3];
  int no;
  RWGObject *PObject=0, *MObject=0;

  char buffer[1000];
  char *Tokens[MAXEDGES+2];
  int nt, NumTokens;
  int LineNum=0;
  int InPortSection=0;
  int PRefPointSpecified=0;
  int MRefPointSpecified=0;
  while( fgets(buffer, 1000, f) )
   { 
     /*--------------------------------------------------------------*/
     /*- skip blank lines and tokens --------------------------------*/
     /*--------------------------------------------------------------*/
     LineNum++;
     NumTokens=Tokenize(buffer, Tokens, MAXEDGES+2); 
     if ( NumTokens==0 || Tokens[0][0]=='#')
      continue;

     /*--------------------------------------------------------------*/
     /*- parse lines separately depending on whether or not we are  -*/
     /*- in a 'PORT...ENDPORT' section                              -*/
     /*--------------------------------------------------------------*/
     if ( !InPortSection ) 
      { if ( !strcasecmp(Tokens[0],"PORT") )
         { InPortSection=1;
           NumPEdges=NumMEdges=0;
           PRefPointSpecified=MRefPointSpecified=0;
           PObject=MObject=G->Objects[0];
         }
        else
         ErrExit("%s:%i: syntax error",PortFileName,LineNum);
      }
     else // !InPortSection
      {
        if ( !strcasecmp(Tokens[0],"ENDPORT") )
         { 
           PortArray = (RWGPort **)realloc( PortArray, (NumPorts+1)*sizeof(PortArray[0]) );
           PortArray[NumPorts] = CreatePort(PObject, NumPEdges, PEIndices,
                                            MObject, NumMEdges, MEIndices);
           if (PRefPointSpecified) 
            memcpy(PortArray[NumPorts]->PRefPoint, PRefPoint, 3*sizeof(double));
           if (MRefPointSpecified) 
            memcpy(PortArray[NumPorts]->MRefPoint, MRefPoint, 3*sizeof(double));

           NumPorts++;
           InPortSection=0;
         }
        else if ( !strcasecmp(Tokens[0],"PEDGES") )
         { if (NumPEdges!=0)
            ErrExit("%s:%i: multiple PEDGES specifications",PortFileName,LineNum);
           for(nt=1; nt<NumTokens; nt++)
            if (1!=sscanf(Tokens[nt],"%i",PEIndices+(nt-1)))
             ErrExit("%s:%i: syntax error %s",PortFileName,LineNum,Tokens[nt]);
           NumPEdges=NumTokens-1;
           if (NumPEdges>=MAXEDGES)
            ErrExit("%s:%i: too many edges",PortFileName,LineNum);
         }
        else if ( !strcasecmp(Tokens[0],"MEDGES") )
         { if (NumMEdges!=0)
            ErrExit("%s:%i: multiple MEDGES specifications",PortFileName,LineNum);
           for(nt=1; nt<NumTokens; nt++)
            if (1!=sscanf(Tokens[nt],"%i",MEIndices+(nt-1)))
              ErrExit("%s:%i: syntax error",PortFileName,LineNum);
           NumMEdges=NumTokens-1;
           if (NumMEdges>=MAXEDGES)
            ErrExit("%s:%i: too many edges",PortFileName,LineNum);
         }
        else if ( !strcasecmp(Tokens[0],"PREFPOINT") )
         {
           if(    NumTokens!=4
               || 1!=sscanf(Tokens[1],"%le",PRefPoint+0)
               || 1!=sscanf(Tokens[2],"%le",PRefPoint+1)
               || 1!=sscanf(Tokens[3],"%le",PRefPoint+2)
             )
            ErrExit("%s:%i: syntax error",PortFileName,LineNum);
           PRefPointSpecified=1;
         }
        else if ( !strcasecmp(Tokens[0],"MREFPOINT") )
         {
           if(    NumTokens!=4
               || 1!=sscanf(Tokens[1],"%le",MRefPoint+0)
               || 1!=sscanf(Tokens[2],"%le",MRefPoint+1)
               || 1!=sscanf(Tokens[3],"%le",MRefPoint+2)
             ) 
            ErrExit("%s:%i: syntax error",PortFileName,LineNum);
           MRefPointSpecified=1;
         }
        else if ( !strcasecmp(Tokens[0],"POBJECT") )
         {
           if( NumTokens!=2 )
            ErrExit("%s:%i: syntax error",PortFileName,LineNum);
           for(no=0; no<G->NumObjects; no++)
            if (!strcasecmp(G->Objects[no]->Label,Tokens[1]))
             break;
           if( no==G->NumObjects )
            ErrExit("%s:%i: could not find object %s in geometry %s", PortFileName,LineNum,Tokens[1],G->GeoFileName);
           PObject=G->Objects[no]; 
         }
        else if ( !strcasecmp(Tokens[0],"MOBJECT") )
         {
           if( NumTokens!=2 )
            ErrExit("%s:%i: syntax error",PortFileName,LineNum);
           for(no=0; no<G->NumObjects; no++)
            if (!strcasecmp(G->Objects[no]->Label,Tokens[1]))
             break;
           if( no==G->NumObjects )
            ErrExit("%s:%i: could not find object %s in geometry %s", PortFileName,LineNum,Tokens[1],G->GeoFileName);
           MObject=G->Objects[no];
         }
       else
        ErrExit("%s:%i: unknown keyword %s",PortFileName,LineNum,Tokens[0]);

      }; // if (InPortSection ...

   };
  fclose(f);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  *pNumPorts=NumPorts;
  return PortArray;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddPortContributionsToRHS(RWGGeometry *G, 
                               RWGPort **Ports, int NumPorts, cdouble *PortCurrents, 
                               cdouble Omega, HVector *KN)
{
  int no, ne, npe, BFIndex, nPort;
  RWGPort *Port;
  RWGEdge *E;
  cdouble IK=II*Omega;
 
  RWGObject *SourceObject;
  int SourcePanelIndex;
  int SourcePaneliQ;
  double SourceLength;
  cdouble PortCurrent, Weight;

  RWGObject *DestObject;
  int DestPPanelIndex, DestMPanelIndex;
  int DestPPaneliQ, DestMPaneliQ;
  double DestLength;

  cdouble EFieldIntegral;

  GetPPIArgStruct GPPIArgsBuffer, *GPPIArgs=&GPPIArgsBuffer;
  InitGetPPIArgs(GPPIArgs);
  GPPIArgs->k=Omega; // this assumes the exterior medium is vacuum

  /***************************************************************/
  /* insert a check to make sure all objects are PEC             */
  /***************************************************************/

  /***************************************************************/
  /* fill in (actually augment) entries of the RHS vector one-by-*/
  /* one                                                         */
  /***************************************************************/
  for(BFIndex=0, no=0; no<G->NumObjects; no++)
   { 
     DestObject=G->Objects[no];
     for(ne=0; ne<DestObject->NumEdges; ne++, BFIndex++)
      { 
        E=DestObject->Edges[ne];

        DestPPanelIndex  = E->iPPanel;
        DestPPaneliQ     = E->PIndex;

        DestMPanelIndex  = E->iMPanel;
        DestMPaneliQ     = E->MIndex;

        DestLength       = E->Length;

        /***************************************************************/
        /* get contribution of all ports to this element of the vector.*/
        /* note: at the conclusion of this loop, we have that          */
        /* EFieldIntegral is the inner product of basis function #ne   */
        /* on object O with the total E-field produced by all driven   */
        /* ports. the corresponding contribution to the RHS vector is  */
        /* then simply -EFieldIntegral.                                */
        /***************************************************************/
        for(EFieldIntegral=0.0, nPort=0; nPort<NumPorts; nPort++)
         { 
           PortCurrent=PortCurrents[nPort];
           if (PortCurrent==0.0) continue;

           Port=Ports[nPort];

           /***************************************************************/
           /* get contributions of panels on the positive side of the port*/
           /***************************************************************/
           SourceObject=Port->PObject;
           Weight=PortCurrent/(Port->PPerimeter);
           for(npe=0; npe<Port->NumPEdges; npe++) // npe = 'num port edge'
            { 
              SourcePanelIndex  = Port->PPanelIndices[npe];
              SourcePaneliQ     = Port->PPaneliQs[npe];
              SourceLength      = Port->PLengths[npe];

              GPPIArgs->Oa  = SourceObject;
              GPPIArgs->npa = SourcePanelIndex;
              GPPIArgs->iQa = SourcePaneliQ;

              GPPIArgs->Ob  = DestObject;
              GPPIArgs->npb = DestPPanelIndex;
              GPPIArgs->iQb = DestPPaneliQ;
              GetPanelPanelInteractions(GPPIArgs);
              EFieldIntegral -= Weight*SourceLength*DestLength*IK*GPPIArgs->H[0];

              GPPIArgs->Ob  = DestObject;
              GPPIArgs->npb = DestMPanelIndex;
              GPPIArgs->iQb = DestMPaneliQ;
              GetPanelPanelInteractions(GPPIArgs);
              EFieldIntegral += Weight*SourceLength*DestLength*IK*GPPIArgs->H[0];

            }; // for ne=0; ne<Port->NumPEdges; ne++)

           /***************************************************************/
           /* get contributions of panels on the negative side of the port*/
           /***************************************************************/
           SourceObject=Port->MObject;
           Weight=PortCurrent/(Port->MPerimeter);
           for(npe=0; npe<Port->NumMEdges; npe++)
            { 
              SourcePanelIndex  = Port->MPanelIndices[npe];
              SourcePaneliQ     = Port->MPaneliQs[npe];
              SourceLength      = Port->MLengths[npe];

              GPPIArgs->Oa  = SourceObject;
              GPPIArgs->npa = SourcePanelIndex;
              GPPIArgs->iQa = SourcePaneliQ;

              GPPIArgs->Ob  = DestObject;
              GPPIArgs->npb = DestPPanelIndex;
              GPPIArgs->iQb = DestPPaneliQ;
              GetPanelPanelInteractions(GPPIArgs);
              EFieldIntegral += Weight*SourceLength*DestLength*IK*GPPIArgs->H[0];

              GPPIArgs->Ob  = DestObject;
              GPPIArgs->npb = DestMPanelIndex;
              GPPIArgs->iQb = DestMPaneliQ;
              GetPanelPanelInteractions(GPPIArgs);
              EFieldIntegral -= Weight*SourceLength*DestLength*IK*GPPIArgs->H[0];

            }; // for ne=0; ne<Port->NumMEdges; ne++)
         }; // for(nPort=0; nPort<NumPorts; nPort++)

        KN->AddEntry(BFIndex, -1.0*EFieldIntegral );

      }; // for(ne=0; ne<O->NumEdges; ne++)
   }; // for(no=0=; no<G->NumObjects; no++)

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PlotPorts(const char *GPFileName, RWGPort **Ports, int NumPorts)
{
  FILE *f=fopen(GPFileName,"w");

  RWGObject *O;
  RWGPort *Port;
  int nPort, nPanel, PanelIndex, PaneliQ;
  double *V1, *V2;
  for(nPort=0; nPort<NumPorts; nPort++)
   { 
     Port=Ports[nPort];

     fprintf(f,"%e %e %e \n", Port->PRefPoint[0], Port->PRefPoint[1], Port->PRefPoint[2]);
     fprintf(f,"\n\n");
     fprintf(f,"%e %e %e \n\n\n", Port->MRefPoint[0], Port->MRefPoint[1], Port->MRefPoint[2]);
     fprintf(f,"\n\n");

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     O=Port->PObject;
     for(nPanel=0; nPanel<Port->NumPEdges; nPanel++)
      { 
        PanelIndex = Port->PPanelIndices[nPanel]; 
        PaneliQ    = Port->PPaneliQs[nPanel]; 
        V1 = O->Vertices + 3*O->Panels[PanelIndex]->VI[(PaneliQ+1)%3];
        V2 = O->Vertices + 3*O->Panels[PanelIndex]->VI[(PaneliQ+2)%3];
        fprintf(f,"%e %e %e \n",V1[0],V1[1],V1[2]);
        fprintf(f,"%e %e %e \n",V2[0],V2[1],V2[2]);
        fprintf(f,"\n\n");
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     O=Port->MObject;
     for(nPanel=0; nPanel<Port->NumMEdges; nPanel++)
      { 
        PanelIndex = Port->MPanelIndices[nPanel]; 
        PaneliQ    = Port->MPaneliQs[nPanel]; 
        V1 = O->Vertices + 3*O->Panels[PanelIndex]->VI[(PaneliQ+1)%3];
        V2 = O->Vertices + 3*O->Panels[PanelIndex]->VI[(PaneliQ+2)%3];
        fprintf(f,"%e %e %e \n",V1[0],V1[1],V1[2]);
        fprintf(f,"%e %e %e \n",V2[0],V2[1],V2[2]);
        fprintf(f,"\n\n");
      };

   };

  fclose(f);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void DrawGMSHCircle(const char *PPFile, const char *Name, 
                    double *X0, double Theta, double Phi, double Radius)
{ 
  FILE *f=fopen(PPFile,"a");
  fprintf(f,"View \"%s\" {\n",Name);

  double CT, ST, CP, SP, M[3][3];
  CT=cos(Theta);
  ST=sin(Theta);
  CP=cos(Phi);
  SP=sin(Phi);
  M[0][0]=CT*CP;     M[0][1]=CT*SP;    M[0][2]=-ST;
  M[1][0]=-SP;       M[1][1]=CP;       M[1][2]=0.0;
  M[2][0]=ST*CP;     M[2][1]=ST*SP;    M[2][2]=CT;

  double Psi, X1P[3], X1[3], X2P[3], X2[3];
  int Mu, Nu;
  for(Psi=0; Psi<2.0*M_PI; Psi+=M_PI/50.0)
   { 
     X1P[0]=Radius*cos(Psi);
     X1P[1]=Radius*sin(Psi);
     X1P[2]=0.0;

     X2P[0]=Radius*cos(Psi+M_PI/50.0);
     X2P[1]=Radius*sin(Psi+M_PI/50.0);
     X2P[2]=0.0;

     memcpy(X1,X0,3*sizeof(double));
     memcpy(X2,X0,3*sizeof(double));
     for(Mu=0; Mu<3; Mu++)
      for(Nu=0; Nu<3; Nu++)
       { X1[Mu]+=M[Nu][Mu]*X1P[Nu];
         X2[Mu]+=M[Nu][Mu]*X2P[Nu];
       };
     fprintf(f,"SL(%e,%e,%e,%e,%e,%e) {%e,%e};\n",
                X1[0],X1[1],X1[2],X2[0],X2[1],X2[2], 0.0, 0.0);
    };

  fprintf(f,"};\n");
  fclose(f);
 
}
