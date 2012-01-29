/*
 * scuff-caspol.cc   -- a standalone code within the scuff-EM suite for 
 *                   -- computing casimir-polder and van-der-Waals 
 *                   -- potentials for polarizable molecules in the 
 *                   -- vicinity of material bodies
 *
 * documentation at: http://homerreid.com/scuff-em/scuff-caspol
 *
 * homer reid        -- 10/2006 -- 2/2012
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include <libhrutil.h>
#include <libIncField.h>
#include <libhmat.h>

#include "libscuff.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXSTR  1000
#define MAXXI   10

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define FILETYPE_BYXI 0
#define FILETYPE_OUT  1
void WriteFilePreamble(FILE *f, int argc, char *argv[], int FileType)
{ 

  time_t MyTime=time(0);
  struct tm *MyTm=localtime(&MyTime);
  char TimeString[200];
  strftime(TimeString,30,"%D::%T",MyTm);

  fprintf(f,"# scuff-caspol running on %s (%s)\n",getenv("HOST"),TimeString);
  fprintf(f,"#");
  fprintf(f,"# command line:\n");
  fprintf(f,"# scuff-caspol ");
  for(narg=0; narg<(argc-1); narg++)
   fprintf(f,"# %s%s",argv[narg+1],(narg%4)==3 ? "\n# " : " ");
  fprintf(f,"#\n");
  fprintf(f,"#\n");
  fprintf(f,"# data columns:\n");
  fprintf(f,"1: X \n");
  fprintf(f,"2: Y \n");
  fprintf(f,"3: Z \n");
  if (FileType==FILETYPE_BYXI)
   { fprintf(f,"4: Xi \n");
     fprintf(f,"5: casimir-polder potential per unit frequency at (X, Y, Z, Xi)\n");
   }
  else
   fprintf(f,"4: casimir-polder potential at (X, Y, Z)\n");
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  char *EPFile=0;
  char *PolFile=0;
  double XiValues;	int nXiValues;
  char *XiFile=0;
  int nThread=0;
  double Temperature=0.0;
  double AbsTol=ABSTOL, RelTol=RELTOL, XiMin=XIMIN;
  OptStruct OSArray[]=
   { {"geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,  0,          ".scuffgeo file"},
     {"PolFile",        PA_STRING,  1, 1,       (void *)&PolFile,  0,          "molecule polarizability file"},
     {"EPFile",         PA_STRING,  1, 1,       (void *)&EPFile,   0,          "list of evaluation points"},
     {"Xi",             PA_DOUBLE,  1, MAXXI,   (void *)XiValues,  &nXiValues, "imaginary frequency"},
     {"XiList",         PA_DOUBLE,  1, 1,       (void *)&XiFile,   0,          "list of imaginary frequencies"},
     {"Temperature",    PA_DOUBLE,  1, 1,       (void *)&Temperature, 0,       "temperature in kelvin"},
     {"AbsTol",         PA_DOUBLE,  0, 1,       (void *)&AbsTol,   0,          "absolute tolerance for summation/integration"},
     {"RelTol",         PA_DOUBLE,  0, 1,       (void *)&RelTol,   0,          "relative tolerance for summation/integration"},
     {"XiMin",          PA_DOUBLE,  0, 1,       (void *)&XiMin,    0,          "relative tolerance for summation/integration"},
     {"nThread",        PA_INT,     0, 1,       (void *)&nThread,  0,          "number of CPU threads to use"},
     {0,0,0,0,0,0,0}
   };

  if (GeoFile==0)
   OSUsage(argv[0], OSArray,"--geometry option is mandatory");
  if (EPFile==0)
   OSUsage(argv[0], OSArray,"--EPfile option is mandatory");
  if (nThread==0)
   nThread=GetNumProcs();

  /*******************************************************************/
  /* create the RWGGeometry, allocate BEM matrix and RHS vector, and */
  /* set up the data structure passed to the computational routines  */
  /*******************************************************************/
  SCPData MySCPData, *SCPD=&MySCPData;

  SCPD->G  = new RWGGeometry(GeoFile); 
  SCPD->M  = G->AllocateBEMMatrix(IMAG_FREQ);
  SCPD->KN = G->AllocateRHSVector(IMAG_FREQ);

  SCPD->nThread=nThread;
  SCPD->AbsTol=AbsTol;
  SCPD->RelTol=RelTol;
  SCPD->XiMin=XiMin;

  SetLogFileName("scuff-caspol.log");

  /*******************************************************************/
  /* process polarization file                                       */
  /*******************************************************************/
  ProcessPolarizationFile(PolFile, SCPD);

  /*******************************************************************/
  /* process list of evaluation points *******************************/
  /*******************************************************************/
  HMatrix *EPList=SCPD->EPList=new HMatrix(EPFile,LHM_TEXT,"--ncol 3");
  if (EPList->ErrMsg)
   ErrExit(EPList->ErrMsg);

  int nEvalPoints=EPList->NR;
  // storage for values of CP potential at eval points
  double U=(double *)malloc(nEvalPoints * sizeof(double)); 

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run simulations                         */
  /*******************************************************************/
  HVector *XiList=0;
  int NumFreqs=0;
  // first process --XiFile option if present
  if (XiFile) 
   { XiList = new HVector(XiFile,LHM_TXT);
     if (XiList->ErrMsg)
      ErrExit(ErrMsg);
     NumFreqs=XiList->N;
   };

  // now add any individually specified --Xi options
  if (nXiValues>0)
   { 
     NumFreqs += nXiValues;

     HVector *XiList0=XiList;
     XiList=new HVector(NumFreqs);

     int nXi;
     if (XiList0)
      { for(nXi=0; nXi<XiList0->N; nXi++)
         XiList->SetEntry(nXi, XiList0->GetEntryD(nFreq));
        delete XiList0;
      };

     int nxv;
     for(nxv=0; nxv<nXiValues; nxv++)
      XiList->SetEntry(nXi+nxv, XiValues[nxv]);
   };

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  char GeoFileBase[MAXSTR], ByXiFileName[MAXSTR]; 
  snprintf(GeoFileBase, MAXSTR, GetFileBase(GeoFile));
  snprintf(ByXiFileName, MAXSTR, "%s.byXi",GeoFileBase);
  SCPD->ByXiFile=CreateUniqueFile(ByXiFileName,1,ByXiFileName);
  WriteFilePreamble(SCPD->ByXiFile, argc, argv, FILETYPE_BYXI);

  /*******************************************************************/
  /* main program body:                                              */
  /*                                                                 */
  /*  (a) if the user specified a set of frequencies, then we just   */
  /*      calculate the contributions to the casimir-polder potential*/
  /*      from those frequencies.                                    */
  /*                                                                 */
  /*  (b) if the user specified a temperature, then we evaluate the  */
  /*      matsubara sum for the full casimir-polder potential at     */
  /*      that temperature.                                          */
  /*                                                                 */
  /*  (c) if the user specified neither of the above, then we        */
  /*      evaluate the integral over the imaginary-frequency axis    */
  /*      to get the full zero-temperature casimir-polder potential. */
  /*                                                                 */
  /*******************************************************************/
  if (XiList)
   {
     for(nXi=0; nXi<XiList->N; nXi++)
      { Xi=XiList->GetEntry(N);
        GetCPIntegrand(SCPD, Xi, U);
      };
   }
  else if ( Temperature != 0.0 )
   { 
     EvaluateMatsubaraSum(SCPD, Temperature, U);
   }
  else
   { 
     EvaluateFrequencyIntegral(SCPD, U);
   };
  fclose(SCPData->ByXiFile);

  /***************************************************************/
  /* write frequency-integrated or matsubara-summed results to   */
  /* output file if that was requested                           */
  /***************************************************************/
  if (!XiList)
   {
     char OutFileName[MAXSTR];
     snprintf(OutFileName, MAXSTR, "%s.out",GeoFileBase);
     FILE *OutFile=CreateUniqueFile(OutFileName,1,OutFileName);
     WriteFilePreamble(OutFileName, argc, argv, FILETYPE_OUT);
     int nep;
     for(nep=0; nep<EPList->NR; nep++)
      fprintf(OutFile,"%e %e %e %e \n", EPList->GetEntryD(nep,0),
                                        EPList->GetEntryD(nep,1),
                                        EPList->GetEntryD(nep,2),
                                        U[nep]);
     fclose(f);
   };

  /***************************************************************/
  /***************************************************************/  
  /***************************************************************/
  fprintf("Frequency-resolved data written to file %s.\n",ByXiFileName);
  if (Temperature>0.0)
   fprintf("Matsubara-summed data written to file %s.\n",OutFileName);
  else if (XiList!=0)
   fprintf("Frequency-integrated data written to file %s.\n",OutFileName);
  printf("Thank you for your support.\n");

}
