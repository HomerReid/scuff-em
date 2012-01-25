/*
 * OutputModules.cc -- various types of 'output modules' for EMScatter
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include <readline/readline.h>
#include <readline/history.h>

#include <libhrutil.h>
#include <libIncField.h>
#include <libhmat.h>

#include "libRWG.h"
#include "EMScatter.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetIncidentField(EMSData *EMSD, double *X, cdouble EH[6])
{ 
  if (EMSD->opPWD)
   EHPlaneWave(X,EMSD->opPWD,EH);
  else if (EMSD->opGBD)
   EHGaussianBeam(X,EMSD->opGBD,EH);
  else if (EMSD->opPSD)
   EHPlaneWave(X,EMSD->opPSD,EH);
  else if (EMSD->opMFD)
   EHMagneticFrill(X,EMSD->opMFD,EH);
  else
   ErrExit("%s:%i: internal error",__FILE__,__LINE__);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetTotalField(EMSData *EMSD, double *X, int WhichObject, cdouble *EHS, cdouble *EHT)
{ 
  /*--------------------------------------------------------------*/
  /*- get scattered field ----------------------------------------*/
  /*--------------------------------------------------------------*/
  EMSD->G->GetFields(X,WhichObject,EMSD->Frequency,EMSD->RealFreq,EMSD->KN,EMSD->nThread,EHS);
  memcpy(EHT,EHS,6*sizeof(cdouble));

  /*--------------------------------------------------------------*/
  /*- add incident field only if we are in the external region    */
  /*--------------------------------------------------------------*/
  if (WhichObject==-1)
   { int Mu;
     cdouble EH2[6];
     GetIncidentField(EMSD,X,EH2);
     for(Mu=0; Mu<6; Mu++) 
      EHT[Mu]+=EH2[Mu];
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Console(EMSData *EMSD)
{ 
  char *p;
  double X[3]; 
  cdouble EHS[6], EHT[6];
  int i;

  using_history();
  read_history(0);
  SetDefaultCD2SFormat("(%+12.5e,%+12.5e)");
  for(;;)
   { 
     p=readline("Enter x y z: "); 
     add_history(p);
     write_history(0);
     sscanf(p,"%le %le %le ",X+0,X+1,X+2);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     GetTotalField(EMSD, X, -1, EHS,EHT); 

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     printf("%25s %25s %25s\n",
            "          ESCAT          ",
            "          EINC           ",
            "          ETOT           ");
     for(i=0; i<3; i++)
      printf("%s %s %s\n", CD2S(EHS[i]), CD2S(EHT[i]-EHS[i]), CD2S(EHT[i]));

     printf("\n");
     printf("%25s %25s %25s\n",
            "          HSCAT          ",
            "          HINC           ",
            "          HTOT           ");
     for(i=3; i<6; i++)
      printf("%s %s %s\n", CD2S(EHS[i]), CD2S(EHT[i]-EHS[i]), CD2S(EHT[i]));

   };

}


/***************************************************************/
/* compute scattered and total fields at a user-specified list */
/* of evaluation points                                        */
/***************************************************************/
void ProcessEPFile(EMSData *EMSD, char *EPFileName, int WhichObject)
{ 
  int nep;
  double X[3]; 
  cdouble EHS[6], EHT[6];
  FILE *f1, *f2;
  char buffer[200];

  HMatrix *EPMatrix=new HMatrix(EPFileName,LHM_TEXT,"-ncol 3");
  if (EPMatrix->ErrMsg)
   { fprintf(stderr,"Error processing EP file: %s\n",EPMatrix->ErrMsg);
     delete EPMatrix;
     return;
   };

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

     GetTotalField(EMSD, X, WhichObject, EHS, EHT); 

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
void CreateFluxPlot(EMSData *EMSD, char *MeshFileName)
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
   GetTotalField(EMSD, P->Centroid, -1, EHS + 6*np, EHT + 6*np);

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
      PF=0.5 * real(  (E[1]*conj(H[2]) - E[2]*conj(H[1])) * P->ZHat[0]
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
      PF=0.5 * real(  (E[1]*conj(H[2]) - E[2]*conj(H[1])) * P->ZHat[0]
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
