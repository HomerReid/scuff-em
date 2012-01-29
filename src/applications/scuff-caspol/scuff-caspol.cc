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
  OptStruct OSArray[]=
   { {"geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,  0,          ".scuffgeo file"},
     {"PolFile",        PA_STRING,  1, 1,       (void *)&PolFile,  0,          "molecule polarizability file"},
     {"EPFile",         PA_STRING,  1, 1,       (void *)&EPFile,   0,          "list of evaluation points"},
     {"Xi",             PA_DOUBLE,  1, MAXXI,   (void *)XiValues,  &nXiValues, "imaginary frequency"},
     {"XiList",         PA_DOUBLE,  1, 1,       (void *)&XiFile,   0,          "list of imaginary frequencies"},
     {"Temperature",    PA_DOUBLE,  1, 1,       (void *)&Temperature, 0,       "temperature in kelvin"},
     {"nThread",        PA_BOOL,    0, 1,       (void *)&nThread,  0,          "number of CPU threads to use"},
     {0,0,0,0,0,0,0}
   };

  if (GeoFile==0)
   OSUsage(argv[0], OSArray,"--geometry option is mandatory");
  if (EPFile==0)
   OSUsage(argv[0], OSArray,"--EPfile option is mandatory");
 // if (PolFile==0) 
 //  OSUsage(argv[0], OSArray,"--PolFile option is mandatory");
  if (nThread==0)
   nThread=GetNumProcs();

  /*******************************************************************/
  /* create the RWGGeometry, allocate BEM matrix and RHS vector      */
  /*******************************************************************/
  SCPData MySCPData, *SCPD=&MySCPData;

  SCPD->G  = new RWGGeometry(GeoFile); 
  SCPD->M  = G->AllocateBEMMatrix(IMAG_FREQ);
  SCPF->KN = G->AllocateRHSVector(IMAG_FREQ);

  SCPD->nThread=nThread;

  char GeoFileBase[MAXSTR];
  snprintf(GeoFile, MAXSTR, GetFileBase(GeoFile));
  SetLogFileName("scuff-caspol.log");

  /*******************************************************************/
  /* process polarization file                                       */
  /*******************************************************************/
  memset(SCPD->AlphaMP, 0, 9*sizeof(SCPD->AlphaMP[0]) );
  ProcessPolarizationFile(

  /*******************************************************************/
  /* process list of evaluation points *******************************/
  /*******************************************************************/
  HMatrix *EPList=new HMatrix(EPFile,LHM_TEXT,"--ncol 3");
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
        GetCPIntegrand(SCPData, Xi, EPList, U);
      };
   }
  else if ( Temperature != 0.0 )
   { 
     EvaluateMatsubaraSum(SCPData, EPList, U);
   }
  else
   { 
     EvaluateFrequencyIntegral(SCPData, EPList, U);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}
