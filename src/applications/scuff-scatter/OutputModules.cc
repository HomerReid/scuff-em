/*
 * OutputModules.cc -- various types of 'output modules' for EMScatter
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include <libSGJC.h>

#include "scuff-scatter.h"

#define ABSTOL   1.0e-4
#define RELTOL   1.0e-2
#define MAXEVALS 10000

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetTotalField(SSData *SSD, double *X, int WhichObject, 
                   cdouble *EHS, cdouble *EHT)
{ 
  /*--------------------------------------------------------------*/
  /*- get scattered field ----------------------------------------*/
  /*--------------------------------------------------------------*/
  SSD->G->GetFields(X,WhichObject,SSD->Omega,SSD->KN,SSD->nThread,EHS);
  memcpy(EHT,EHS,6*sizeof(cdouble));

  /*--------------------------------------------------------------*/
  /*- add incident field only if we are in the external region    */
  /*--------------------------------------------------------------*/
  if (WhichObject==-1)
   { int Mu;
     cdouble EH2[6];
     EHIncField(X, SSD->opIFD, EH2);
     for(Mu=0; Mu<6; Mu++) 
      EHT[Mu]+=EH2[Mu];
   };

}

/***************************************************************/
/* compute scattered and total fields at a user-specified list */
/* of evaluation points                                        */
/***************************************************************/
void ProcessEPFile(SSData *SSD, char *EPFileName)
{ 

  HMatrix *EPMatrix=new HMatrix(EPFileName,LHM_TEXT,"-ncol 3");
  if (EPMatrix->ErrMsg)
   { fprintf(stderr,"Error processing EP file: %s\n",EPMatrix->ErrMsg);
     delete EPMatrix;
     return;
   };
 
  /***************************************************************/
  /* FIXME *******************************************************/
  /***************************************************************/
  int WhichObject=-1;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int nep;
  double X[3]; 
  cdouble EHS[6], EHT[6];
  FILE *f1, *f2;
  char buffer[200];
 
  Log("Evaluating fields at points in file %s...",EPFileName);
  printf("Evaluating fields at points in file %s...\n",EPFileName);

  sprintf(buffer,"%s.scattered",GetFileBase(EPFileName));
  f1=CreateUniqueFile(buffer,1);
  sprintf(buffer,"%s.total",GetFileBase(EPFileName));
  f2=CreateUniqueFile(buffer,1);
  SetDefaultCD2SFormat("%18.12e %18.12e");
  for(nep=0; nep<EPMatrix->NR; nep++)
   { 
     X[0]=EPMatrix->GetEntryD(nep, 0);
     X[1]=EPMatrix->GetEntryD(nep, 1);
     X[2]=EPMatrix->GetEntryD(nep, 2);

     GetTotalField(SSD, X, WhichObject, EHS, EHT); 

     fprintf(f1,"%e %e %e ",X[0],X[1],X[2]);
     fprintf(f1,"%s %s %s ",CD2S(EHS[0]),CD2S(EHS[1]),CD2S(EHS[2]));
     fprintf(f1,"%s %s %s ",CD2S(EHS[3]),CD2S(EHS[4]),CD2S(EHS[5]));
     fprintf(f1,"\n");

     fprintf(f2,"%e %e %e ",X[0],X[1],X[2]);
     fprintf(f2,"%s %s %s ",CD2S(EHT[0]),CD2S(EHT[1]),CD2S(EHT[2]));
     fprintf(f2,"%s %s %s ",CD2S(EHT[3]),CD2S(EHT[4]),CD2S(EHT[5]));
     fprintf(f2,"\n");
 
   };
  fprintf(f1,"\n\n");
  fprintf(f2,"\n\n");
}

/***************************************************************/
/* generate plots of poynting flux and field-strength arrows   */
/* on a user-supplied surface mesh                             */
/***************************************************************/
void CreateFluxPlot(SSData *SSD, char *MeshFileName)
{ 
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
  cdouble EHS[6], EHT[6];
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
  /* -PVTot and PVScat                                           */
  /***************************************************************/
  fval[0] = -R*R*( PVTot[0]*nHat[0] + PVTot[1]*nHat[1]  + PVTot[2]*nHat[2] );
  fval[1] =  R*R*(PVScat[0]*nHat[0] + PVScat[1]*nHat[1] + PVScat[2]*nHat[2]);

}

/***************************************************************/
/* evaluate the scattered and absorbed powers by integrating   */
/* the poynting flux of the scattered and total fields over a  */
/* bounding sphere at radius r=R                               */
/***************************************************************/
void GetPower_BF(SSData *SSD, double R, double *PBF, double *EBF)
{ 
  double Lower[2]={-1.0, 0.0};
  double Upper[2]={+1.0, 2.0*M_PI};

  GPBFIData MyGPBFID, *GPBFID = &MyGPBFID;
  GPBFID->SSD = SSD;
  GPBFID->R = R;

  adapt_integrate_log(2, GetPower_BF_Integrand, (void *)GPBFID, 2, 
	     	      Lower, Upper, MAXEVALS, ABSTOL, RELTOL, 
                      PBF, EBF, "SGJC.log",15);

}

/***************************************************************/
/* get the scattered and absorbed power using steven's formulas*/
/* involving vector-matrix-vector products                     */
/***************************************************************/
void GetPower_SGJ(SSData *SSD, double *PSGJ)
{
  RWGGeometry *G = SSD->G;
  HMatrix *M     = SSD->M;
  HVector *RHS   = SSD->RHS;
  HVector *KN    = SSD->KN;

  /***************************************************************/
  /* a quick sanity check ****************************************/
  /***************************************************************/
  cdouble *ZM=M->ZM, *ZRHS=RHS->ZV, *ZKN=KN->ZV;
  if (ZM==0 || ZRHS==0 || ZKN==0)
   { PSGJ[0]=PSGJ[1]=0.0;
     return;
   };

  /***************************************************************/
  /* get the M0 matrix. we overwrite the M matrix for this       */
  /* purpose, since we won't need that until the next frequency, */
  /* at which point we will have to re-assemble it anyway.       */
  /***************************************************************/
  int no;
  for(no=0; no<G->NumObjects; no++)
   G->Objects[no]->MP->Zero();

  G->AssembleBEMMatrix(SSD->Omega, SSD->nThread, M);

  for(no=0; no<G->NumObjects; no++)
   G->Objects[no]->MP->UnZero();

  /***************************************************************/
  /* compute the vector-matrix-vector and vector-vector products */
  /* that enter into the formulas for the scattered and absorbed */
  /* power                                                       */
  /***************************************************************/
  // nr runs over rows, nc over columns, ne over matrix entries
  int nr, nc, ne, N=G->TotalBFs;
  double VMV, VV;
  double Sign;
  for(VMV=VV=0.0, Sign=1.0, ne=nc=0; nc<N; nc++)
   for(nr=0; nr<N; nr++, ne++, Sign*=-1.0)
    { 
      if (nr==nc) 
       VV += Sign*real( conj(ZKN[nr]) * ZRHS[nr] );

      VMV += -Sign*real( conj(ZKN[nr]) * ZM[ne] * ZKN[nc] );
    };
    
  double PScat, PAbs;
  PAbs  =  0.25 * ZVAC * (VV+VMV);
  PScat = -0.25 * ZVAC * (VV-VMV);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PSGJ[0]=PAbs;
  PSGJ[1]=PScat;

}

/***************************************************************/
/* evaluate the scattered and absorbed powers using the        */
/* concise BEM vector-matrix-vector product expressions        */
/***************************************************************/
void GetPower(SSData *SSD, char *PowerFile)
{
  RWGGeometry *G = SSD->G;
  HVector *RHS   = SSD->RHS;
  HVector *KN    = SSD->KN;

  Log("  Computing scattered and absorbed power...");

  /***************************************************************/
  /* open the file and write the frequency at the top of the line*/
  /***************************************************************/
  FILE *f=fopen(PowerFile,"a");
  if ( !f ) ErrExit("could not open file %s",PowerFile);

  if ( imag(SSD->Omega)==0.0 )
   fprintf(f,"%15.8e ",real(SSD->Omega));
  else if ( real(SSD->Omega)==0.0 )
   fprintf(f,"%15.8ei ",imag(SSD->Omega));
  else 
   fprintf(f,"%15.8e+%15.8ei ",real(SSD->Omega),imag(SSD->Omega));

  /***************************************************************/
  /* get absorbed and scattered powers                           */
  /***************************************************************/
  double PAbs=0.0, PScat=0.0;
  double OTimes;
  int no, nea, neb, Offset;
  RWGObject *O;
  cdouble ka, na, nb, vE, vH;
  for(no=0; no<G->NumObjects; no++)
   { 
     O=G->Objects[no];
     if ( O->ContainingObject != 0 ) // skip interior surfaces of nested objects
      continue; 

     Offset=G->BFIndexOffset[no];
     if ( O->MP->IsPEC() )
      { 
        for(nea=0; nea<O->NumEdges; nea++)
         { ka = KN->GetEntry(Offset  + nea );
           vE = RHS->GetEntry(Offset + nea );
           PScat += real( conj(ka*vE) );
         }; // for(nea==...
      }
     else
      {
        for(nea=0; nea<O->NumEdges; nea++)
         { 
           ka = KN->GetEntry(Offset + 2*nea + 0 );
           na = KN->GetEntry(Offset + 2*nea + 1 );

           vE = RHS->GetEntry(Offset + 2*nea + 0 );
           vH = RHS->GetEntry(Offset + 2*nea + 1 );

           PScat+= real( conj(ka)*vE - conj(na)*vH );

           for(neb=0; neb<O->NumEdges; neb++)
            { 
             O->GetOverlap(nea, neb, &OTimes);
             if (OTimes==0.0) 
              continue;

             nb = KN->GetEntry(Offset + 2*neb + 1 );
             PAbs += real( conj(ka) * OTimes * nb );
            }; // for (neb= ... 

         }; // for(nea==...
      }; // if ( O->MP->IsPEC() )
   }; // for(no=...)
  
  PAbs *= -0.5*ZVAC;
  PScat = -PAbs + 0.5*ZVAC*PScat;
  fprintf(f,"%e %e  ",PAbs, PScat);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double PSGJ[2];
  GetPower_SGJ(SSD, PSGJ);
  fprintf(f,"%e %e  ",PSGJ[0], PSGJ[1]);

  /***************************************************************/
  /* if the user specified a nonzero PowerRadius, repeat the     */
  /* power calculation by brute-force integration of the total   */
  /* and scattered Poynting vectors over a sphere of that radius */
  /***************************************************************/
  if (SSD->PowerRadius > 0.0 )
   { double PBF[2], EBF[2]; 
     GetPower_BF(SSD, SSD->PowerRadius, PBF, EBF);
     fprintf(f,"%e %e %e %e ",PBF[0],EBF[0],PBF[1],EBF[1]);
   };

  fprintf(f,"\n");
  fclose(f);
   
}
