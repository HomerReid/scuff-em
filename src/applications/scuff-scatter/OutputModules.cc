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
/***************************************************************/
/***************************************************************/
void GetTotalField(SSData *SSD, double *X, cdouble *EHS, cdouble *EHT)
{ 
  SSD->G->GetFields(      0, SSD->KN, SSD->Omega, X, EHS, 0); // scattered
  SSD->G->GetFields(SSD->IF,       0, SSD->Omega, X, EHT, 0); // incident
 
  for (int Mu=0; Mu<6; Mu++)
   EHT[Mu]+=EHS[Mu];
}

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

  HMatrix *FMatrix1 = G->GetFields( 0, KN, Omega, XMatrix); // scattered
  HMatrix *FMatrix2 = G->GetFields(IF,  0, Omega, XMatrix); // incident

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
static char *FieldFuncs=
 "|Ex|,|Ey|,|Ez|,"
 "sqrt(|Ex|^2+|Ey|^2+|Ez|^2),"
 "|Hx|,|Hy|,|Hz|,"
 "sqrt(|Hx|^2+|Hy|^2+|Hz|^2)";

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
  RWGObject *O=new RWGObject(MeshFileName);

  Log("Creating flux plot for surface %s...",MeshFileName);
  printf("Creating flux plot for surface %s...\n",MeshFileName);

  /*--------------------------------------------------------------*/
  /*- create an Nx3 HMatrix whose columns are the coordinates of  */
  /*- the flux mesh panel vertices                                */
  /*--------------------------------------------------------------*/
  HMatrix *XMatrix=new HMatrix(O->NumVertices, 3);
  int nv;
  for(nv=0; nv<O->NumVertices; nv++)
   { 
     XMatrix->SetEntry(nv, 0, O->Vertices[3*nv + 0]);
     XMatrix->SetEntry(nv, 1, O->Vertices[3*nv + 1]);
     XMatrix->SetEntry(nv, 2, O->Vertices[3*nv + 2]);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *FMatrix=SSD->G->GetFields(SSD->IF, SSD->KN, SSD->Omega, XMatrix, 0, FieldFuncs);

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
     for(np=0; np<O->NumPanels; np++)
      {
        P=O->Panels[np];
        iV1 = P->VI[0];  V1 = O->Vertices + 3*iV1;
        iV2 = P->VI[1];  V2 = O->Vertices + 3*iV2;
        iV3 = P->VI[2];  V3 = O->Vertices + 3*iV3;

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
  delete O;

}

/***************************************************************/
/* generate plots of poynting flux and field-strength arrows   */
/* on a user-supplied surface mesh                             */
/***************************************************************/
#if 0
void CreateFluxPlot(SSData *SSD, char *MeshFileName)
{ 
  RWGObject *O=new RWGObject(MeshFileName);
  RWGPanel *P;
  int np;

  Log("Creating flux plot for surface %s...",MeshFileName);
  printf("Creating flux plot for surface %s...\n",MeshFileName);

  /*--------------------------------------------------------------*/
  /*- create an Nx3 HMatrix whose columns are the coordinates of  */
  /*- the centroids of the panels on the flux mesh                */
  /*--------------------------------------------------------------*/
  int NP=O->NumPanels;
  HMatrix *XMatrix=new HMatrix(NP, 3);
  for(np=0; np<O->NumPanels; np++)
   { 
     P=O->Panels[np]; 
     XMatrix->SetEntry(np, 0, P->Centroid[0]);
     XMatrix->SetEntry(np, 1, P->Centroid[1]);
     XMatrix->SetEntry(np, 2, P->Centroid[2]);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *FSMatrix=SSD->G->GetFields( 0, SSD->KN, SSD->Omega, XMatrix); // scattered
  HMatrix *FTMatrix=SSD->G->GetFields(SSD->IF,  0, SSD->Omega, XMatrix); // incident

  // set total = incident + scattered
  FTMatrix->AddBlock(FSMatrix, 0, 0);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble *EHS=(cdouble *)malloc(6*O->NumPanels*sizeof(cdouble));
  cdouble *EHT=(cdouble *)malloc(6*O->NumPanels*sizeof(cdouble));
  if (EHS==0 || EHT==0) 
   ErrExit("out of memory");
  for(np=0; np<NP; np++)
   { 
      EHS[6*np + 0] = FSMatrix->GetEntry(np, 0);
      EHS[6*np + 1] = FSMatrix->GetEntry(np, 1);
      EHS[6*np + 2] = FSMatrix->GetEntry(np, 2);
      EHS[6*np + 3] = FSMatrix->GetEntry(np, 3);
      EHS[6*np + 4] = FSMatrix->GetEntry(np, 4);
      EHS[6*np + 5] = FSMatrix->GetEntry(np, 5);

      EHT[6*np + 0] = FTMatrix->GetEntry(np, 0);
      EHT[6*np + 1] = FTMatrix->GetEntry(np, 1);
      EHT[6*np + 2] = FTMatrix->GetEntry(np, 2);
      EHT[6*np + 3] = FTMatrix->GetEntry(np, 3);
      EHT[6*np + 4] = FTMatrix->GetEntry(np, 4);
      EHT[6*np + 5] = FTMatrix->GetEntry(np, 5);
   };
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
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
  delete FSMatrix;
  delete FTMatrix;
  delete XMatrix;
  delete O;

}
#endif 

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
  (void) ndim;
  (void) fdim;

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
  GetTotalField(SSD, X, EHS, EHT);

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

  double PVSI[3]; // 'poynting vector, scattered--incident'
  cdouble EHI[6];
  cdouble *EInc=EHI;
  cdouble *HInc=EHI+3;
  for(int n=0; n<6; n++) 
   EHI[n] = EHT[n] - EHS[n];
  
  PVSI[0]  = 0.5*real( EScat[1] * conj(HInc[2])  - EScat[2] * conj(HInc[1]) );
  PVSI[1]  = 0.5*real( EScat[2] * conj(HInc[0])  - EScat[0] * conj(HInc[2]) );
  PVSI[2]  = 0.5*real( EScat[0] * conj(HInc[1])  - EScat[1] * conj(HInc[0]) );

  PVSI[0] += 0.5*real( EInc[1] * conj(HScat[2])  - EInc[2] * conj(HScat[1]) );
  PVSI[1] += 0.5*real( EInc[2] * conj(HScat[0])  - EInc[0] * conj(HScat[2]) );
  PVSI[2] += 0.5*real( EInc[0] * conj(HScat[1])  - EInc[1] * conj(HScat[0]) );

  /***************************************************************/
  /* components of integrand vector are radial components of     */
  /* PVSI and PVScat                                             */
  /***************************************************************/
  //fval[0] = -R*R*( PVTot[0]*nHat[0] +  PVTot[1]*nHat[1] +  PVTot[2]*nHat[2]);
  fval[0] =  R*R*(  PVSI[0]*nHat[0] +   PVSI[1]*nHat[1] +   PVSI[2]*nHat[2]);
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
  PBF[0] = -(PBF[0] + PBF[1]);
  EBF[0] =  (EBF[0] + EBF[1]);

}

/***************************************************************/
/* get the scattered and absorbed power using the SGJ formulas */
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

  G->AssembleBEMMatrix(SSD->Omega, M, SSD->nThread);

  for(no=0; no<G->NumObjects; no++)
   G->Objects[no]->MP->UnZero();

  /***************************************************************/
  /* compute the vector-matrix-vector and vector-vector products */
  /* that enter into the formulas for the scattered and absorbed */
  /* power                                                       */
  /***************************************************************/
  // nr runs over rows, nc over columns, ne over matrix entries
  int nr, nc, ne, N=G->TotalBFs;
  double PTot, PAbs, PScat, Sign;
  for(PTot=PScat=0.0, Sign=1.0, ne=nc=0; nc<N; nc++)
   for(nr=0; nr<N; nr++, ne++, Sign*=-1.0)
    { 
      if (nr==nc) 
       PTot += Sign*real( conj(ZKN[nr]) * (-1.0*ZRHS[nr]) );

      PScat -= Sign*real( conj(ZKN[nr]) * ZM[ne] * ZKN[nc] );
    };
    
  PTot  *= 0.5 * ZVAC;
  PScat *= 0.5 * ZVAC;

  PAbs  = PTot - PScat;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PSGJ[0]=PAbs;
  PSGJ[1]=PScat;

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#if 0
double PTotBH[2], PScatBQ[4];
double *WhichPTot, *WhichPScat;
memset(PTotBH, 0, 2*sizeof(double));
memset(PScatBQ, 0, 4*sizeof(double));
for(Sign=1.0, ne=nc=0; nc<N; nc++)
 for(nr=0; nr<N; nr++, ne++, Sign*=-1.0)
  { 
    if ( (nr%2)==0 && (nc%2)==0 )
     { WhichPTot  = PTotBH + 0;
       WhichPScat = PScatBQ + 0;
     }
    else if ( (nr%2)==0 && (nc%2)==1 )
     { WhichPTot  = PTotBH + 0;
       WhichPScat = PScatBQ + 1;
     }
    else if ( (nr%2)==1 && (nc%2)==0 )
     { WhichPTot  = PTotBH + 1;
       WhichPScat = PScatBQ + 2;
     }
    else if ( (nr%2)==1 && (nc%2)==1 )
     { WhichPTot  = PTotBH + 1;
       WhichPScat = PScatBQ + 3;
     };

    if ( nr==nc ) 
     *WhichPTot += Sign*real( conj(ZKN[nr]) * (-1.0*ZRHS[nr]) );

    *WhichPScat -= Sign*real( conj(ZKN[nr]) * ZM[ne] * ZKN[nc] );
  };
PTotBH[0] *= 0.5*ZVAC;
PTotBH[1] *= 0.5*ZVAC;
PScatBQ[0] *= 0.5*ZVAC;
PScatBQ[1] *= 0.5*ZVAC;
PScatBQ[2] *= 0.5*ZVAC;
PScatBQ[3] *= 0.5*ZVAC;
FILE *ff=fopen("byQuadrant.out","a");
fprintf(ff,"%e %.12e %.12e %.12e %.12e %.12e %.12e \n",
            real(SSD->Omega), PTotBH[0], PTotBH[1], 
            PScatBQ[0], PScatBQ[1], PScatBQ[2], PScatBQ[3]);
fclose(ff);
#endif
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

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

  Log("Computing scattered and absorbed power...");

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
  double PTot=0.0, PAbs=0.0, PScat;
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
           PTot += real( conj(ka*vE) );
         }; // for(nea==...
      }
     else
      {
        for(nea=0; nea<O->NumEdges; nea++)
         { 
           ka = KN->GetEntry(Offset + 2*nea + 0 );
           na = KN->GetEntry(Offset + 2*nea + 1 );

           vE = -1.0*RHS->GetEntry(Offset + 2*nea + 0 );
           vH = -1.0*RHS->GetEntry(Offset + 2*nea + 1 );

           PTot += real( conj(ka)*vE - conj(na)*vH );

           for(neb=0; neb<O->NumEdges; neb++)
            { 
              O->GetOverlap(nea, neb, &OTimes);
              if (OTimes==0.0) 
               continue;
            }; // for (neb= ... 

         }; // for(nea==...
      }; // if ( O->MP->IsPEC() )
   }; // for(no=...)
  
  PTot *= 0.5*ZVAC;
  PAbs *= 0.5*ZVAC;
  PScat = PTot - PAbs;

  /***************************************************************/
  /* the HR formula computes the scattered power as the          */
  /* difference between the total and absorbed power; if the     */
  /* total and absorbed powers agree to three decimal places,    */
  /* then the scattered power computed this way will probably be */
  /* inaccurate, so in this case we use the SGJ scattered-power  */
  /* fomulas instead                                             */
  /***************************************************************/
  if ( fabs(PScat) < 1.0e-3*fabs(PTot) ) 
   { 
     double PSGJ[2];
     GetPower_SGJ(SSD, PSGJ);
     Log("Using SGJ formulas (PTot,PAbs,PScat) = (%.5e,%.5e,%.5e) (HR) (%.5e,%.5e,%.5e) (SGJ)",
          PTot,PAbs,PScat,PSGJ[0]+PSGJ[1],PSGJ[0],PSGJ[1]);
     PScat=PSGJ[1];
   };
  fprintf(f,"%.12e %.12e  ",PAbs,PScat);

  /***************************************************************/
  /* if the user specified a nonzero PowerRadius, repeat the     */
  /* power calculation by brute-force integration of the total   */
  /* and scattered Poynting vectors over a sphere of that radius */
  /***************************************************************/
  if (SSD->PowerRadius > 0.0 )
   { double PBF[2], EBF[2]; 
     GetPower_BF(SSD, SSD->PowerRadius, PBF, EBF);

     fprintf(f,"%.12e %.12e %.12e %.12e ",PBF[0],EBF[0],PBF[1],EBF[1]);
   };

  fprintf(f,"\n");
  fclose(f);
   
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
  setlinebuf(f);
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
  int no, Mu;
  fprintf(f,"%s ",z2s(Omega));
  for (no=0; no<G->NumObjects; no++)
   { fprintf(f,"%s ",G->Objects[no]->Label);
     for(Mu=0; Mu<6; Mu++)
      fprintf(f,"%s ",CD2S(PM->GetEntry(6*no + Mu),"%.8e %.8e "));
   };
  fprintf(f,"\n");

  delete PM;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define OVERLAP_OVERLAP     0
#define OVERLAP_CROSS       1
#define OVERLAP_XBULLET     2
#define OVERLAP_XNABLANABLA 3
#define OVERLAP_XTIMESNABLA 4
#define OVERLAP_YBULLET     5
#define OVERLAP_YNABLANABLA 6
#define OVERLAP_YTIMESNABLA 7
#define OVERLAP_ZBULLET     8
#define OVERLAP_ZNABLANABLA 9
#define OVERLAP_ZTIMESNABLA 10
void GetForce(SSData *SSD, char *ForceFile)
{
  RWGGeometry *G = SSD->G;
  HVector *KN    = SSD->KN;
  cdouble Omega  = SSD->Omega;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=fopen(ForceFile,"a");
  if (!f)
   { Warn("could not open file %s for append",ForceFile);
     return;
   };
  fprintf(f,"%s ",z2s(Omega));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble Eps, Mu;
  G->ExteriorMP->GetEpsMu(Omega, &Eps, &Mu);
  cdouble Z2 = ZVAC*ZVAC*Mu/Eps;
  cdouble OOZ2 = 1.0/Z2;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Force[3];
  double Overlaps[11];
  double OiBullet, OiNablaNabla, OiTimesNabla;
  int no, nfc, neAlpha, neBeta, Offset, IsPEC;
  cdouble KAlpha, NAlpha=0.0, KBeta, NBeta=0.0;
  cdouble FK2 = 4.0*Omega*Omega; 
  double PreFac=-0.5;
  cdouble M11, M12, M21, M22;
  RWGObject *O;
  for(no=0; no<G->NumObjects; no++)
   { 
     O=G->Objects[no];
     IsPEC = O->MP->IsPEC() ? 1 : 0;
     Offset=G->BFIndexOffset[no];
     memset(Force,0,3*sizeof(double));
     for(neAlpha=0; neAlpha<O->NumEdges; neAlpha++)
      for(neBeta=0; neBeta<O->NumEdges; neBeta++)
       { 
         if (IsPEC) 
          { KAlpha = KN->GetEntry( Offset + neAlpha );
            KBeta  = KN->GetEntry( Offset + neBeta  );
          }
         else
          { KAlpha =       KN->GetEntry( Offset + 2*neAlpha + 0 );
            NAlpha = -ZVAC*KN->GetEntry( Offset + 2*neAlpha + 1 );
            KBeta  =       KN->GetEntry( Offset + 2*neBeta  + 0 );
            NBeta  = -ZVAC*KN->GetEntry( Offset + 2*neBeta  + 1 );
          };

         O->GetOverlaps(neAlpha, neBeta, Overlaps);
         if (Overlaps[0]==0.0)
          continue; 
 
         for(nfc=0; nfc<3; nfc++)
          { 
            OiBullet     = Overlaps[ 2 + (nfc*3) + 0 ];
            OiNablaNabla = Overlaps[ 2 + (nfc*3) + 1 ];
            OiTimesNabla = Overlaps[ 2 + (nfc*3) + 2 ];

            M11 = Z2*(-OiBullet + OiNablaNabla/FK2); 
            M12 = OiTimesNabla / (II*Omega);
            M21 = OiTimesNabla / (II*Omega);
            M22 = OOZ2*(-OiBullet + OiNablaNabla/FK2); 

            Force[nfc] += PreFac*real(   conj(KAlpha)*M11*KBeta 
                                       + conj(KAlpha)*M12*NBeta
                                       + conj(NAlpha)*M21*KBeta 
                                       + conj(NAlpha)*M22*NBeta );
          };

       };

     fprintf(f,"%s %e %e %e ",O->Label,Force[0],Force[1],Force[2]);

   };// for(no=...)

  fprintf(f,"\n");
  fclose(f);

}
