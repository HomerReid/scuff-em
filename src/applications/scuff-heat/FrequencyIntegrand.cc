/*
 * FrequencyIntegrand.cc -- evaluate the spectral density of 
 *                       -- power transfer at a single frequency
 *
 * homer reid            -- 2/2012
 *
 */

#include "scuff-heat.h"

/***************************************************************/
/* generate plots of poynting flux and field-strength arrows   */
/* on a user-supplied surface mesh                             */
/***************************************************************/
void CreateFluxPlot(ScuffHeatData *SHD, cdouble Omega)
{ 
  FILE *f=vfopen("%s.%g.flux","w",GetFileBase(SHD->GeoFileName),real(Omega));


  RWGObject *O=new RWGObject(MeshFileName);
  RWGPanel *P;
  int np;

  Log("Creating flux plot for surface %s...",MeshFileName);
  printf("Creating flux plot for surface %s...\n",MeshFileName);

  /*--------------------------------------------------------------*/
  /*- because of the way the datasets need to be organized in    -*/
  /*- the .pp file, it is easiest to make multiple passes        -*/
  /*- through the list of panels on the flux surface: one pass   -*/
  /*- to get the E and H fields at each panel centroid, and then -*/
  /*- subsequent passes to write each of the various different   -*/
  /*- data sets to the .pp file.                                 -*/
  /*--------------------------------------------------------------*/
  cdouble *EHS=(cdouble *)malloc(6*O->NumPanels*sizeof(cdouble));
  cdouble *EHT=(cdouble *)malloc(6*O->NumPanels*sizeof(cdouble));
  if (EHS==0 || EHT==0) ErrExit("out of memory");
  for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
   GetTotalField(SSD, P->Centroid, -1, EHS + 6*np, EHT + 6*np);

  // maximum values of scattered and total fields, used below for 
  // normalization 
  double MaxESMag=0.0, MaxETMag=0.0, MaxHSMag=0.0, MaxHTMag=0.0;
  for(np=0; np<O->NumPanels; np++)
   { MaxESMag=fmax(MaxESMag,sqrt( norm(EHS[6*np+0]) + norm(EHS[6*np+1]) + norm(EHS[6*np+2]) ));
     MaxETMag=fmax(MaxETMag,sqrt( norm(EHT[6*np+0]) + norm(EHT[6*np+1]) + norm(EHT[6*np+2]) ));
     MaxHSMag=fmax(MaxHSMag,sqrt( norm(EHS[6*np+3]) + norm(EHS[6*np+4]) + norm(EHS[6*np+5]) ));
     MaxHTMag=fmax(MaxHTMag,sqrt( norm(EHT[6*np+3]) + norm(EHT[6*np+4]) + norm(EHT[6*np+5]) ));
   };

  /*--------------------------------------------------------------*/
  /*- try to open output file ------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=vfopen("%s.pp","w",GetFileBase(MeshFileName));
  if (!f) 
   { fprintf(stderr,"warning: could not open output file %s.pp\n",GetFileBase(MeshFileName));
     free(EHS); 
     free(EHT); 
     delete O;
     return;
   };
  
  /*--------------------------------------------------------------*/
  /*- poynting flux of scattered field ---------------------------*/
  /*--------------------------------------------------------------*/
  double *PV[3]; // panel vertices 
  double PF;    // poynting flux
  cdouble *E, *H;
  fprintf(f,"View \"Poynting Flux (Scattered)\" {\n");
  for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
   {
      PV[0]=O->Vertices + 3*P->VI[0];
      PV[1]=O->Vertices + 3*P->VI[1];
      PV[2]=O->Vertices + 3*P->VI[2];

      // poynting flux = (E \cross H^*) \dot (panel normal) / 2 
      E=EHS + 6*np;
      H=EHS + 6*np + 3;
      PF=0.5 * real(   (E[1]*conj(H[2]) - E[2]*conj(H[1])) * P->ZHat[0]
                      +(E[2]*conj(H[0]) - E[0]*conj(H[2])) * P->ZHat[1]
                      +(E[0]*conj(H[1]) - E[1]*conj(H[0])) * P->ZHat[2] 
                   );

      fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                 PV[0][0], PV[0][1], PV[0][2],
                 PV[1][0], PV[1][1], PV[1][2],
                 PV[2][0], PV[2][1], PV[2][2],
                 PF,PF,PF);
    };
  fprintf(f,"};\n\n");
  
  /*--------------------------------------------------------------*/
  /*- poynting flux of total field     ---------------------------*/
  /*--------------------------------------------------------------*/
  fprintf(f,"View \"Poynting Flux (Total)\" {\n");
  for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
   {
      PV[0]=O->Vertices + 3*P->VI[0];
      PV[1]=O->Vertices + 3*P->VI[1];
      PV[2]=O->Vertices + 3*P->VI[2];

      // poynting flux = (E \cross H^*) \dot (panel normal) / 2 
      E=EHT + 6*np;
      H=EHT + 6*np + 3;
      PF=0.5 * real(   (E[1]*conj(H[2])  - E[2]*conj(H[1])) * P->ZHat[0]
                      +(E[2]*conj(H[0]) - E[0]*conj(H[2])) * P->ZHat[1]
                      +(E[0]*conj(H[1]) - E[1]*conj(H[0])) * P->ZHat[2] 
                   );

      fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                 PV[0][0], PV[0][1], PV[0][2],
                 PV[1][0], PV[1][1], PV[1][2],
                 PV[2][0], PV[2][1], PV[2][2],
                 PF,PF,PF);
    };
  fprintf(f,"};\n\n");

  /*--------------------------------------------------------------*/
  /*- in order that the arrows in the vector plots of the E and H */
  /*- fields be scaled in a way that makes the plot legible       */
  /*- graphically, i normalize all field strengths to ensure that */
  /*- the length of the longest field-strength arrow is the       */
  /*- average radius of a panel in the flux mesh.                 */
  /*--------------------------------------------------------------*/
  double NormFac, AvgPanelRadius=0;

  for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
   AvgPanelRadius += P->Radius;
  AvgPanelRadius/=((double)(O->NumPanels));
   
  /*--------------------------------------------------------------*/
  /*- real part of scattered E field   ---------------------------*/
  /*--------------------------------------------------------------*/
  NormFac=MaxESMag / AvgPanelRadius;
  fprintf(f,"View \"E_Scattered (real part)\" {\n");
  for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
   {
      E=EHS + 6*np;
      fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                 P->Centroid[0], P->Centroid[1], P->Centroid[2],
                 real(E[0]/NormFac), real(E[1]/NormFac), real(E[2]/NormFac));
    };
  fprintf(f,"};\n\n");

  /*--------------------------------------------------------------*/
  /*- imag part of scattered E field   ---------------------------*/
  /*--------------------------------------------------------------*/
  fprintf(f,"View \"E_Scattered (imag part)\" {\n");
  for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
   {
      // poynting flux = (E \cross H^*) \dot (panel normal) / 2 
      E=EHS + 6*np;
      fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                 P->Centroid[0], P->Centroid[1], P->Centroid[2],
                 imag(E[0]/NormFac), imag(E[1]/NormFac), imag(E[2]/NormFac));
    };
  fprintf(f,"};\n\n");
   
  /*--------------------------------------------------------------*/
  /*- real part of total E field ---------------------------------*/
  /*--------------------------------------------------------------*/
  NormFac=MaxETMag / AvgPanelRadius;
  fprintf(f,"View \"E_Total(real part)\" {\n");
  for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
   {
      // poynting flux = (E \cross H^*) \dot (panel normal) / 2 
      E=EHT + 6*np;
      fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                 P->Centroid[0], P->Centroid[1], P->Centroid[2],
                 real(E[0]/NormFac), real(E[1]/NormFac), real(E[2]/NormFac));
    };
  fprintf(f,"};\n\n");

  /*--------------------------------------------------------------*/
  /*- imag part of total E fields --------------------------------*/
  /*--------------------------------------------------------------*/
  fprintf(f,"View \"E_Total (imag part)\" {\n");
  for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
   {
      E=EHT + 6*np;
      fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                 P->Centroid[0], P->Centroid[1], P->Centroid[2],
                 imag(E[0]/NormFac), imag(E[1]/NormFac), imag(E[2]/NormFac));
    };
  fprintf(f,"};\n\n");

  /*--------------------------------------------------------------*/
  /*- real part of scattered H field   ---------------------------*/
  /*--------------------------------------------------------------*/
  NormFac=MaxHSMag / AvgPanelRadius;
  fprintf(f,"View \"H_Scattered (real part)\" {\n");
  for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
   {
      H=EHS + 6*np + 3;
      fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                 P->Centroid[0], P->Centroid[1], P->Centroid[2],
                 real(H[0]/NormFac), real(H[1]/NormFac), real(H[2]/NormFac));
    };
  fprintf(f,"};\n\n");

  /*--------------------------------------------------------------*/
  /*- imag part of scattered H field   ---------------------------*/
  /*--------------------------------------------------------------*/
  fprintf(f,"View \"H_Scattered (imag part)\" {\n");
  for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
   {
      // poynting flux = (E \cross H^*) \dot (panel normal) / 2 
      H=EHS + 6*np + 3;
      fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                 P->Centroid[0], P->Centroid[1], P->Centroid[2],
                 imag(H[0]/NormFac), imag(H[1]/NormFac), imag(H[2]/NormFac));
    };
  fprintf(f,"};\n\n");

  /*--------------------------------------------------------------*/
  /*- real part of total H field   -------------------------------*/
  /*--------------------------------------------------------------*/
  NormFac=MaxHSMag / AvgPanelRadius;
  fprintf(f,"View \"H_Total (real part)\" {\n");
  for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
   {
      // poynting flux = (E \cross H^*) \dot (panel normal) / 2 
      H=EHT + 6*np + 3;
      fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                 P->Centroid[0], P->Centroid[1], P->Centroid[2],
                 real(H[0]/NormFac), real(H[1]/NormFac), real(H[2]/NormFac));
    };
  fprintf(f,"};\n\n");

  /*--------------------------------------------------------------*/
  /*- imag part of total H field ---------------------------------*/
  /*--------------------------------------------------------------*/
  fprintf(f,"View \"H_Total (imag part)\" {\n");
  for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
   {
      // poynting flux = (E \cross H^*) \dot (panel normal) / 2 
      H=EHT + 6*np + 3;
      fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                 P->Centroid[0], P->Centroid[1], P->Centroid[2],
                 imag(H[0]/NormFac), imag(H[1]/NormFac), imag(H[2]/NormFac));
    };
  fprintf(f,"};\n\n");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  fclose(f);

  free(EHS);
  free(EHT);
  delete O;

}

/***************************************************************/
/* integrand routine for GetPower_BF ***************************/
/***************************************************************/
typedef struct GPBFIData 
 {
   SSData *SSD;
   double R;

 } GPBFIData;

void GetPower_BF_Integrand(unsigned ndim, const double *x, void *params,
			   unsigned fdim, double *fval)
{
  /***************************************************************/
  /* extract fields from data structure **************************/
  /***************************************************************/
  GPBFIData *GPBFID = (GPBFIData *)params;
  SSData *SSD = GPBFID->SSD;
  double R    = GPBFID->R;

  /***************************************************************/
  /* get coordinates of evaluation point *************************/
  /***************************************************************/
  double CosTheta=x[0];
  double SinTheta=sqrt(1.0-CosTheta*CosTheta);
  double CosPhi=cos(x[1]);
  double SinPhi=sin(x[1]);

  double nHat[3], X[3];
  nHat[0]=SinTheta*CosPhi;
  nHat[1]=SinTheta*SinPhi;
  nHat[2]=CosTheta;
  X[0]=R*nHat[0];
  X[1]=R*nHat[1];
  X[2]=R*nHat[2];
  
  /***************************************************************/
  /* get total and scattered fields at evaluation point **********/
  /***************************************************************/
  cdouble EHS[3], EHT[3];
  GetTotalField(SSD, X, -1, EHS, EHT);

  /***************************************************************/
  /* get scattered and total poynting vectors ********************/
  /***************************************************************/
  cdouble *EScat = EHS;
  cdouble *HScat = EHS+3;
  cdouble *ETot  = EHT;
  cdouble *HTot  = EHT+3;
  double PVScat[3], PVTot[3];

  PVScat[0] = 0.5*real( EScat[1] * conj(HScat[2]) - EScat[2] * conj(HScat[1]) );
  PVScat[1] = 0.5*real( EScat[2] * conj(HScat[0]) - EScat[0] * conj(HScat[2]) );
  PVScat[2] = 0.5*real( EScat[0] * conj(HScat[1]) - EScat[1] * conj(HScat[0]) );

  PVTot[0]  = 0.5*real( ETot[1]  * conj(HTot[2])  - ETot[2]  * conj(HTot[1])  );
  PVTot[1]  = 0.5*real( ETot[2]  * conj(HTot[0])  - ETot[0]  * conj(HTot[2])  );
  PVTot[2]  = 0.5*real( ETot[0]  * conj(HTot[1])  - ETot[1]  * conj(HTot[0])  );

  /***************************************************************/
  /* components of integrand vector are radial components of     */
  /* PVScat and -PVTot                                           */
  /***************************************************************/
  fval[0] = R*R*(PVScat[0]*nHat[0] + PVScat[1]*nHat[1] + PVScat[2]*nHat[2]);
  fval[1] = -R*R*( PVTot[0]*nHat[0]  + PVTot[1]*nHat[1]  + PVTot[2]*nHat[2] );

}

/***************************************************************/
/* evaluate the scattered and absorbed powers by integrating   */
/* the poynting flux of the scattered and total fields over a  */
/* bounding sphere at radius r=R                               */
/***************************************************************/
void GetPower_BF(SSData *SSD, double R, double *PScat, double *PTot)
{ 
  double Val[2], Err[2];
  double Lower[2]={-1.0, 0.0};
  double Upper[2]={+1.0, 2.0*M_PI};

  GPBFIData MyGPBFID, *GPBFID = &MyGPBFID;
  GPBFID->SSD = SSD;
  GPBFID->R = R;

  adapt_integrate(2, GetPower_BF_Integrand, (void *)SSD, 2, 
		  Lower, Upper, 0, ABSTOL, RELTOL, Val, Err);

  *PScat=Val[0];
  *PTot=Val[1];

}

/***************************************************************/
/* evaluate the scattered and absorbed powers using the        */
/* concise BEM vector-matrix-vector product expressions        */
/***************************************************************/
void GetPower(SSData *SSD, double *PScat, double *PTot)
{ 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetFrequencyIntegrand(ScuffHeatData *SHD, cdouble Omega, double *FI)
{
  RWGGeometry *G = SHD->G;
  int nThread    = SHD->nThread;
  HMatrix *M0    = SHD->M0;
  HMatrix *M1    = SHD->M1;
  HMatrix *M2    = SHD->M2;

  /***************************************************************/
  /* assemble the three separate blocks that contribute to the   */
  /* BEM matrix:                                                 */
  /*  0: contributions of environment only                       */
  /*  1: contributions of object 1 only                          */
  /*  2: contributions of objects 2...N only                     */
  /***************************************************************/

  /*--------------------------------------------------------------*/
  /*- 0: contributions of environment only -----------------------*/
  /*--------------------------------------------------------------*/
  int no;
  for(no=0; no<G->NumObjects; no++)
   G->Objects[no]->MP->Zero();

  Log("Assembling M0 matrix...");
  G->AssembleBEMMatrix(Omega, nThread, M0);

  // multiply by the 'S' matrix 
  int nr, nc;
  for(nr=1; nr<M0->NR; nr+=2)
   for(nc=0; nc<M0->NC; nc++)
    M0->SetEntry(nr, nc, -1.0*M0->GetEntry(nr,nc));

  /*--------------------------------------------------------------*/
  /*- 1: contributions of object 1 only    -----------------------*/
  /*--------------------------------------------------------------*/
  G->ExteriorMP->Zero();
  G->Objects[0]->MP->UnZero();

  Log("Assembling M1 matrix...");
  G->AssembleBEMMatrix(Omega, nThread, M1);

  for(nr=1; nr<M1->NR; nr+=2)
   for(nc=0; nc<M1->NC; nc++)
    M1->SetEntry(nr, nc, -1.0*M1->GetEntry(nr,nc));

  /*--------------------------------------------------------------*/
  /*- 2: contributions of objects 2--N only.                      */
  /*-                                                             */
  /*-  NOTE: if we only have a single object then we set M2 = M0; */
  /*-        the calculation them amounts to computing the heat   */
  /*-        transfer the single body to the environment.         */
  /*--------------------------------------------------------------*/
  if (G->NumObjects==1)
   M2->Copy(M0);
  else // (G->NumObjects>1)
   { 
     G->Objects[0]->MP->Zero();
     for(no=1; no<G->NumObjects; no++)
      G->Objects[no]->MP->UnZero();

     Log("Assembling M2 matrix...");
     G->AssembleBEMMatrix(Omega, nThread, M2);

     for(nr=1; nr<M2->NR; nr+=2)
      for(nc=0; nc<M2->NC; nc++)
       M2->SetEntry(nr, nc, -1.0*M2->GetEntry(nr,nc));
   };
  
  // undo the zeroing out of the environment and object 1
  G->ExteriorMP->UnZero();
  G->Objects[0]->MP->UnZero();

  /***************************************************************/
  /* assemble and LU-factorize the full BEM matrix (for now we   */
  /* do this in-place using the M0 matrix)                       */
  /***************************************************************/
  if (G->NumObjects==1)
   { for(nr=0; nr<M0->NR; nr++)
      for(nc=0; nc<M0->NC; nc++)
       M0->AddEntry(nr, nc, M1->GetEntry(nr,nc) );
   }
  else
   { for(nr=0; nr<M0->NR; nr++)
      for(nc=0; nc<M0->NC; nc++)
       M0->AddEntry(nr, nc, M1->GetEntry(nr,nc) + M2->GetEntry(nr,nc) );
   };

  Log("LU-factorizing MTotal...");
  M0->LUFactorize();

  /***************************************************************/
  /* set M1 = sym(M1) and M2 = sym(M2)                           */
  /*  (note that the loop runs over only the upper triangle of   */
  /*  the matrices)                                              */
  /***************************************************************/
  cdouble Sym, SymT;
  for(nr=0; nr<M1->NR; nr++)
   for(nc=nr; nc<M1->NC; nc++)
    { 
      Sym = 0.5*(M1->GetEntry(nr, nc) + conj(M1->GetEntry(nc, nr)));
      M1->SetEntry(nr, nc, Sym );
      if(nc>nr) M1->SetEntry(nc, nr, conj(Sym) );
    }; 

  if (SHD->PlotFlux)
   {
     for(nr=0; nr<M2->NR; nr++)
      for(nc=nr; nc<M2->NC; nc++)
       { Sym  = conj(M2->GetEntry(nr, nc));
         SymT = conj(M2->GetEntry(nc, nr));
         M2->SetEntry(nc, nr, Sym );
         if (nc>nr) M2->SetEntry(nr, nc, Sym );
       };

     Log("LU-solving M1...");
     M0->LUSolve(M1,'N');
     Log("LU-solving M2...");
     M0->LUSolve(M2,'C');
     Log("Multipliying...");
     M1->Multiply(M2, M0);
     for(nr=0; nr<M1->NR; nr++)
      SHD->DV->SetEntry(nr, M1->GetEntry(nr,nr));
     Log("Plotting flux vector...");
     PlotFluxVector(SHD);

     *FI=0.0;
   }
  else
   {
     for(nr=0; nr<M2->NR; nr++)
      for(nc=nr; nc<M2->NC; nc++)
       { Sym = 0.5*(M2->GetEntry(nr, nc) + conj(M2->GetEntry(nc, nr)));
         M2->SetEntry(nr, nc, Sym );
         if(nc>nr) M2->SetEntry(nc, nr, conj(Sym) );
       };

     /***************************************************************/
     /* set M1 <= M^{-1'} * M1 **************************************/
     /* set M2 <= M^{-1}  * M2 **************************************/
     /***************************************************************/
     Log("LU-solving M1...");
     M0->LUSolve(M1,'C');
     Log("LU-solving M2...");
     M0->LUSolve(M2,'N');
   
     /***************************************************************/
     /* set M0 = M1*M2                                              */
     /***************************************************************/
     Log("Multiplying M1*M2...");
     M1->Multiply(M2, M0);
   
     /***************************************************************/
     /* the value of the frequency integrand is now the trace of M0 */
     /***************************************************************/
     *FI = real( M0->GetTrace() ) / 8.0 ;
     Log("...done!");

     /***************************************************************/
     /* write the result to the frequency-resolved output file ******/
     /***************************************************************/
     FILE *f=fopen(SHD->ByOmegaFile, "a");
     fprintf(f,"%s %e\n",z2s(Omega),*FI);
     fclose(f);
   };
   
}
