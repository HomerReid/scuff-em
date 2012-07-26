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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libRWG.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
double GetMinDistanceToPanel(double *V1, double *V2, double *V3,
                             double *X, double *XMin, int *OnBoundary);

int PIOVerbose=0;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int PointInsideObject(double *X, RWGObject *O)
{ 

  if (O->MP->IsPEC()) 
   return 0;

  /*--------------------------------------------------------------*/
  /*-- find the vertex closest to X ------------------------------*/
  /*-- (MinD2V = 'minimum distance to vertex' --------------------*/
  /*--------------------------------------------------------------*/
  int np, nvi, nvMD=0;
  RWGPanel *P;
  double D2V, MinD2V = 1.0e9;
  for(np=0; np<O->NumPanels; np++)
   for(P=O->Panels[np], nvi=0; nvi<3; nvi++)
    { D2V=VecDistance(X, O->Vertices + 3*P->VI[nvi]);
      if (D2V<MinD2V)
       { MinD2V=D2V;
         nvMD=P->VI[nvi];
       };
    };
  double *VMD = O->Vertices + 3*nvMD; // 'vertex of minimum distance'

if (PIOVerbose)
 printf("nvMD=%i (%e)\n",nvMD,MinD2V);

  /*--------------------------------------------------------------*/
  /*- now find the point on the surface of the object that lies   */
  /*- closest to X, by looking at all panels that share the       */
  /*- minimum-distance vertex we just found and finding the point */
  /*- closest to X within each panel                              */
  /*--------------------------------------------------------------*/
  double D2P, MinD2P=1e9;   // 'minimum distance to panel'
  double XMin[3], SPMD[3];  // 'surface point of minimum distance'
  int npMD;
  int OnBoundary, MinIsOnBoundary=0;
  for(np=0; np<O->NumPanels; np++)
   { 
     P=O->Panels[np];
     if ( P->VI[0]!=nvMD && P->VI[1]!=nvMD && P->VI[2]!=nvMD ) 
      continue;

     D2P=GetMinDistanceToPanel( O->Vertices + 3*P->VI[0],
                                O->Vertices + 3*P->VI[1],
                                O->Vertices + 3*P->VI[2],
                                X, XMin, &OnBoundary);

     if ( D2P < MinD2P )
      { MinD2P=D2P;
        npMD=np;
        memcpy(SPMD, XMin, 3*sizeof(double) );
        MinIsOnBoundary=OnBoundary;
      };

if (PIOVerbose)
 printf(" md to panel %i: %e (%e, %e, %e) (onboundary=%i)\n",
          np,D2P, VecDistance(X,O->Vertices+3*P->VI[0]),
                  VecDistance(X,O->Vertices+3*P->VI[1]),
                  VecDistance(X,O->Vertices+3*P->VI[2]),OnBoundary);

if (PIOVerbose)
{
printf("%e %e %e \n", O->Vertices[3*P->VI[0] + 0],
                      O->Vertices[3*P->VI[0] + 1],
                      O->Vertices[3*P->VI[0] + 2]);
printf("%e %e %e \n", O->Vertices[3*P->VI[1] + 0],
                      O->Vertices[3*P->VI[1] + 1],
                      O->Vertices[3*P->VI[1] + 2]);
printf("%e %e %e \n", O->Vertices[3*P->VI[2] + 0],
                      O->Vertices[3*P->VI[2] + 1],
                      O->Vertices[3*P->VI[2] + 2]);
printf("%e %e %e \n", O->Vertices[3*P->VI[0] + 0],
                      O->Vertices[3*P->VI[0] + 1],
                      O->Vertices[3*P->VI[0] + 2]);
printf("\n\n");
printf("%e %e %e \n",X[0],X[1],X[2]);
printf("%e %e %e \n",XMin[0],XMin[1],XMin[2]);
printf("\n\n");
printf("%e %e %e \n",XMin[0],XMin[1],XMin[2]);
printf("%e %e %e \n",XMin[0] + P->Radius*P->ZHat[0], 
                     XMin[1] + P->Radius*P->ZHat[1], 
                     XMin[2] + P->Radius*P->ZHat[2]);
printf("\n\n");
};

   };

  /*--------------------------------------------------------------*/
  /*- if the closest point to X was in fact that closest vertex  -*/
  /*- we found earlier, then we test for inclusion by seeing if  -*/
  /*- the distance from X to the nearest reference point inside  -*/
  /*- the object is less than the distance from that closest     -*/
  /*- vertex to the nearest reference point                      -*/
  /*--------------------------------------------------------------*/
  if ( MinIsOnBoundary )
   { 
if (PIOVerbose)
 printf("PIOVerbose case 1...%e %e \n",VecNorm(X),VecNorm(SPMD));
      if (O->NumRefPts==0)
       return VecNorm(X) <= VecNorm(SPMD) ? 1 : 0;
      else
       { 
         int nrp, nrpMin;
         double DSRP;
         double MDSRP; // 'minimum distance from surface to reference point'
 
         MDSRP=VecDistance(SPMD, O->Vertices + 3*O->RefPtIndices[0]);
         nrpMin=0;
         for(nrp=1; nrp<O->NumRefPts; nrp++)
          { DSRP=VecDistance(SPMD, O->Vertices + 3*O->RefPtIndices[nrp]);
            if (DSRP<MDSRP)
             { MDSRP=DSRP;
               nrpMin=nrp;
             };
          };

        if (VecDistance(X, O->Vertices+3*O->RefPtIndices[nrpMin]) < MDSRP) 
         return 1;
        return 0;
       };
   }
  else
   { 
     /*--------------------------------------------------------------*/
     /*- if displacing the nearest surface point in the direction of */
     /*- outward-pointing surface normal increases the distance to X,*/
     /*- then we conclude that X lies inside the object              */
     /*--------------------------------------------------------------*/
     double XmC[3];
     P=O->Panels[npMD];
     VecSub(X, SPMD, XmC);
if (PIOVerbose)
printf("PIOVerbose case 2...%e\n",VecDot(XmC,P->ZHat));
     if ( VecDot(XmC,P->ZHat) < 0.0 )
      return 1;
     else
      return 0;
   };

  return 0;

}
