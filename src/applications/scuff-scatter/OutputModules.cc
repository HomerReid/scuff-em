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
 * OutputModules.cc -- various types of 'output modules' for EMScatter
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include <libSGJC.h>

#include "scuff-scatter.h"

#define MAXSTR   1000 

#define ABSTOL   0.0
#define RELTOL   5.0e-2
#define MAXEVALS 20000


#define II cdouble(0.0,1.0)

/***************************************************************/
/* compute scattered and total fields at a user-specified list */
/* of evaluation points                                        */
/***************************************************************/
void ProcessEPFile(SSData *SSD, char *EPFileName)
{ 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry *G  = SSD->G;
  IncField *IF    = SSD->IF;
  HVector  *KN    = SSD->KN;
  cdouble  Omega  = SSD->Omega;
  double *kBloch  = SSD->kBloch;

  /*--------------------------------------------------------------*/
  /*- try to read eval points from file --------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *XMatrix=new HMatrix(EPFileName,LHM_TEXT,"-ncol 3");
  if (XMatrix->ErrMsg)
   { fprintf(stderr,"Error processing EP file: %s\n",XMatrix->ErrMsg);
     delete XMatrix;
     return;
   };

  /*--------------------------------------------------------------*/
  /*- get components of scattered and incident fields            -*/
  /*--------------------------------------------------------------*/
  Log("Evaluating fields at points in file %s...",EPFileName);

  HMatrix *FMatrix1 = G->GetFields( 0, KN, Omega, kBloch, XMatrix); // scattered
  HMatrix *FMatrix2 = G->GetFields(IF,  0, Omega, kBloch, XMatrix); // incident

  /*--------------------------------------------------------------*/
  /*- create .scattered and .total output files and write fields -*/
  /*--------------------------------------------------------------*/
  char buffer[MAXSTR];
  snprintf(buffer,MAXSTR,"%s.scattered",GetFileBase(EPFileName));
  FILE *f1=CreateUniqueFile(buffer,1);
  snprintf(buffer,MAXSTR,"%s.total",GetFileBase(EPFileName));
  FILE *f2=CreateUniqueFile(buffer,1);

  int nr, nc; 
  SetDefaultCD2SFormat("%.8e %.8e ");
  for(nr=0; nr<FMatrix1->NR; nr++)
   { fprintf(f1,"%.8e %.8e %.8e ",XMatrix->GetEntryD(nr, 0),
                                  XMatrix->GetEntryD(nr, 1),
                                  XMatrix->GetEntryD(nr, 2));

     fprintf(f2,"%.8e %.8e %.8e ",XMatrix->GetEntryD(nr, 0),
                                  XMatrix->GetEntryD(nr, 1),
                                  XMatrix->GetEntryD(nr, 2));

     for(nc=0; nc<FMatrix1->NC; nc++)
      { 
        fprintf(f1,"%s ",CD2S(  FMatrix1->GetEntry(nr,nc)) );

        fprintf(f2,"%s ",CD2S(  FMatrix1->GetEntry(nr,nc)  
                               +FMatrix2->GetEntry(nr,nc)) );
      };

     fprintf(f1,"\n");
     fprintf(f2,"\n");

   };

  fclose(f1);
  fclose(f2);
  delete XMatrix;
  delete FMatrix1;
  delete FMatrix2;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
static char *FieldFuncs=const_cast<char *>(
 "|Ex|,|Ey|,|Ez|,"
 "sqrt(|Ex|^2+|Ey|^2+|Ez|^2),"
 "|Hx|,|Hy|,|Hz|,"
 "sqrt(|Hx|^2+|Hy|^2+|Hz|^2)");

static const char *FieldTitles[]=
 {"|Ex|", "|Ey|", "|Ez|", "|E|",
  "|Hx|", "|Hy|", "|Hz|", "|H|",
 };

#define NUMFIELDFUNCS 8

void CreateFluxPlot(SSData *SSD, char *MeshFileName)
{ 
  /*--------------------------------------------------------------*/
  /*- try to open output file ------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=vfopen("%s.pp","w",GetFileBase(MeshFileName));
  if (!f) 
   { fprintf(stderr,"warning: could not open output file %s.pp\n",GetFileBase(MeshFileName));
     return;
   };
  
  /*--------------------------------------------------------------*/
  /*- try to open user's mesh file -------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=new RWGSurface(MeshFileName);

  Log("Creating flux plot for surface %s...",MeshFileName);
  printf("Creating flux plot for surface %s...\n",MeshFileName);

  /*--------------------------------------------------------------*/
  /*- create an Nx3 HMatrix whose columns are the coordinates of  */
  /*- the flux mesh panel vertices                                */
  /*--------------------------------------------------------------*/
  HMatrix *XMatrix=new HMatrix(S->NumVertices, 3);
  int nv;
  for(nv=0; nv<S->NumVertices; nv++)
   { 
     XMatrix->SetEntry(nv, 0, S->Vertices[3*nv + 0]);
     XMatrix->SetEntry(nv, 1, S->Vertices[3*nv + 1]);
     XMatrix->SetEntry(nv, 2, S->Vertices[3*nv + 2]);
   };

  /*--------------------------------------------------------------*/
  /*- get the total fields at the panel vertices                 -*/
  /*--------------------------------------------------------------*/
  HMatrix *FMatrix=SSD->G->GetFields(SSD->IF, SSD->KN, SSD->Omega, SSD->kBloch, 
                                     XMatrix, 0, FieldFuncs);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGPanel *P;
  int nff, np, iV1, iV2, iV3;
  double *V1, *V2, *V3;
  for(nff=0; nff<NUMFIELDFUNCS; nff++)
   { 
     fprintf(f,"View \"%s\" {\n",FieldTitles[nff]);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     for(np=0; np<S->NumPanels; np++)
      {
        P=S->Panels[np];
        iV1 = P->VI[0];  V1 = S->Vertices + 3*iV1;
        iV2 = P->VI[1];  V2 = S->Vertices + 3*iV2;
        iV3 = P->VI[2];  V3 = S->Vertices + 3*iV3;

        fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                   V1[0], V1[1], V1[2], 
                   V2[0], V2[1], V2[2], 
                   V3[0], V3[1], V3[2], 
                   FMatrix->GetEntryD(iV1,nff),
                   FMatrix->GetEntryD(iV2,nff),
                   FMatrix->GetEntryD(iV3,nff));

      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     fprintf(f,"};\n\n");
   };

  fclose(f);
  delete FMatrix;
  delete XMatrix;
  delete S;

}

/***************************************************************/
/* evaluate the induced dipole moments                         */
/***************************************************************/
void GetMoments(SSData *SSD, char *MomentFile)
{
  /***************************************************************/
  /* open the file and write the frequency at the top of the line*/
  /***************************************************************/
  FILE *f=fopen(MomentFile,"a");
  if ( !f ) ErrExit("could not open file %s",MomentFile);
  
  /***************************************************************/
  /* get dipole moments ******************************************/
  /***************************************************************/
  RWGGeometry *G = SSD->G;
  HVector *KN    = SSD->KN;
  cdouble Omega  = SSD->Omega;

  HVector *PM  = G->GetDipoleMoments(Omega, KN);

  /***************************************************************/
  /* print to file ***********************************************/
  /***************************************************************/
  fprintf(f,"%s ",z2s(Omega));
  for (int ns=0; ns<G->NumSurfaces; ns++)
   { fprintf(f,"%s ",G->Surfaces[ns]->Label);
     for(int Mu=0; Mu<6; Mu++)
      fprintf(f,"%s ",CD2S(PM->GetEntry(6*ns + Mu),"%.8e %.8e "));
   };
  fprintf(f,"\n");
  fflush(f);

  delete PM;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WritePFTFile(SSData *SSD, char *PFTFile)
{
  FILE *f=fopen(PFTFile,"a");
  if (!f)
   return;

  Log("Computing power, force, torque at Omega=%s...",z2s(SSD->Omega));

  RWGGeometry *G=SSD->G;
  double PFT[8]; 

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     fprintf(f,"%s %s ",z2s(SSD->Omega),G->Surfaces[ns]->Label);

     G->GetPFT(SSD->KN, SSD->RHS, SSD->Omega, ns, PFT);

     // get scattered power as difference between total and absorbed power.
     // if the difference between the quantities is less than 10% of their 
     // magnitude, recompute using the more accurate but slower SGJ formula.
     double PScat=PFT[1] - PFT[0];
     if ( fabs(PScat) < 0.1*fabs(PFT[1]) )
      { 
        Log(" Computing PScat via SGJ method [(PAbs,PTot)=(%.3e,%.3e)]",PFT[0],PFT[1]);
        PScat=G->GetScatteredPower(SSD->KN, SSD->Omega, ns);
      };
     // put PScat in the 1 slot in the array 
     PFT[1] = PScat;

     for(int nq=0; nq<8; nq++)
      fprintf(f,"%e ",PFT[nq]);
     fprintf(f,"\n");
   };

  fclose(f);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WritePSDFile(SSData *SSD, char *PSDFile)
{
  FILE *f=fopen(PSDFile,"a");
  if (!f)
   return;

  Log("Computing panel source densities at Omega=%s...",z2s(SSD->Omega));

  static HMatrix *PSDMatrix=0;
  if (PSDMatrix==0)
   PSDMatrix=SSD->G->GetPanelSourceDensities(SSD->Omega, SSD->KN, 0); 
  else
   SSD->G->GetPanelSourceDensities(SSD->Omega, SSD->KN, PSDMatrix); 

  for(int nr=0; nr<PSDMatrix->NR; nr++)
   { fprintf(f,"%s ",z2s(SSD->Omega));
     for(int nc=0; nc<4; nc++)
      fprintf(f,"%e ",real(PSDMatrix->GetEntry(nr,nc)));
     for(int nc=4; nc<=11; nc++)
      fprintf(f,"%e %e ",real(PSDMatrix->GetEntry(nr,nc)),
                         imag(PSDMatrix->GetEntry(nr,nc)));
     for(int nc=12; nc<=12; nc++)
      fprintf(f,"%e ",real(PSDMatrix->GetEntry(nr,nc)));
     fprintf(f,"\n");
   };
  fclose(f);


  SSD->G->GetPanelSourceDensities2(SSD->Omega, SSD->KN, PSDMatrix); 

  f=vfopen("%s.Take0","w",PSDFile);
  for(int nr=0; nr<PSDMatrix->NR; nr++)
   { fprintf(f,"%s ",z2s(SSD->Omega));
     for(int nc=0; nc<4; nc++)
      fprintf(f,"%e ",real(PSDMatrix->GetEntry(nr,nc)));
     for(int nc=4; nc<=11; nc++)
      fprintf(f,"%e %e ",real(PSDMatrix->GetEntry(nr,nc)),
                         imag(PSDMatrix->GetEntry(nr,nc)));
     for(int nc=12; nc<=12; nc++)
      fprintf(f,"%e ",real(PSDMatrix->GetEntry(nr,nc)));
     fprintf(f,"\n");
   };
  fclose(f);

}
