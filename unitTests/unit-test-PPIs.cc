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
 * unit-test-PPIs.cc -- SCUFF-EM unit tests for panel-panel integrals 
 * 
 * homer reid        -- 11/2005 -- 10/2011
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include "libscuff.h"
#include "libscuffInternals.h"

using namespace scuff;

#define II cdouble (0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define NCLW 0 
#define NCMW 1
#define NCSW 2
#define CVLW 3 
#define CVMW 4 
#define CVSW 5 
#define CELW 6 
#define CEMW 7 
#define CESW 8 
#define CTLW 9
#define CTMW 10 
#define CTSW 11

const char *CaseStrings[]=
 { "no common vertices, long wavelength",
   "no common vertices, moderate wavelength",
   "no common vertices, short imaginary wavelength",
   "common-vertex, long wavelength",
   "common-vertex, moderate wavelength",
   "common-vertex, short imaginary wavelength",
   "common-edge, long wavelength",
   "common-edge, moderate wavelength",
   "common-edge, short imaginary wavelength",
   "common-triangle, long wavelength",
   "common-triangle, moderate wavelength",
   "common-triangle, short imaginary wavelength",
 };

cdouble ExactHValues[]=
 { //NCLW
   cdouble(-5.380869460e+00,-7.929329541e-01) ,
   cdouble(-1.203784538e-04,-1.244908468e-07) ,
   //NCMW
   cdouble(-1.097563488e-03,-1.833260738e-03) ,
   cdouble(-8.115226086e-04,+3.380858352e-04) ,
   //NCSW
   cdouble(+4.102927698e-18,+3.052119470e-17) ,
   cdouble(-5.053807314e-16,-5.115425623e-16) ,
   //CVLW
   cdouble(-3.524853903e+01,-7.957111381e-01) ,
   cdouble(+6.121515147e-04,+2.221998597e-09) ,
   //CVMW
   cdouble(-6.533461672e-03,-1.297335912e-02) ,
   cdouble(+8.426339263e-04,+2.352340213e-04) ,
   //CVSW
   cdouble(-6.821962338e-05,-2.337628513e-05) ,
   cdouble(+7.297969569e-05,+3.992918286e-05) ,
   //CELW
   cdouble(-6.239923968e+01,-7.957716374e-01) ,
   cdouble(+2.402244187e-04,+1.434069003e-15) ,
   //CEMW
   cdouble(-2.178203417e-02,-1.571388971e-02) ,
   cdouble(+2.489349158e-04,+4.229357165e-07) ,
   //CESW
   cdouble(-4.206537219e-04,-2.440499471e-05) ,
   cdouble(+1.756575807e-04,+1.876616003e-05) ,
   //CTLW
   cdouble(-1.313680741e+02,-7.957096522e-01) ,
   cdouble(+0.000000000e+00,+0.000000000e+00) ,
   //CTMW
   cdouble(-4.041769199e-02,-1.278652037e-02) ,
   cdouble(+0.000000000e+00,+0.000000000e+00) ,
   //CTSW
   cdouble(+4.906526745e-03,+1.348488281e-03) ,
   cdouble(+0.000000000e+00,+0.000000000e+00)
};

cdouble ExactGradHValues[]=
 { 
   //NCLW
   cdouble(+1.502636050e+00,+1.549165675e-03) ,
   cdouble(+6.043309045e-05,-3.999852750e-08) ,
   cdouble(+3.433439134e+00,+3.538868074e-03) ,
   cdouble(+1.497424520e-04,-7.547129141e-08) ,
   cdouble(-1.070314832120e-01,-1.109692895744e-04) ,
   cdouble(-1.857527798e-05,-8.480699031e-09) ,

   //NCMW
   cdouble(+3.974917953e-03,-1.695100547e-03) ,
   cdouble(-4.784204396e-04,-1.695686924e-03) ,
   cdouble(+9.057717451e-03,-3.869356261e-03) ,
   cdouble(-9.884719990e-04,-3.846876777e-03) ,
   cdouble(-2.975204991e-04,+1.220203417e-04) ,
   cdouble(-3.158616218e-05,+2.115827537e-04) ,

   //NCSW
   cdouble(-8.992250606e-17,-2.390013435e-16) ,
   cdouble(+5.221056603e-15,+3.235122753e-15) ,
   cdouble(-2.196127334e-16,-5.608924910e-16) ,
   cdouble(+1.188236819e-14,+7.360878996e-15) ,
   cdouble(-1.241700346e-18,+8.826165882e-18) ,
   cdouble(-6.075993693e-16,-3.945105243e-16) ,
 };

/***************************************************************/
/* returns 0 on success, 1 on failure.                         */
/***************************************************************/
int ComparePPIs(int WhichCase, GetPPIArgStruct *Args)
{
  
  int ErrorsDetected = 0;

  cdouble *ExactH = ExactHValues + 2*WhichCase;

  for(int n=0; n<2; n++)
   if ( !EqualFloat(Args->H[n], ExactH[n]) )
    { Log("%s: H[%i]: (%s,%s)\n",CaseStrings[WhichCase],n,CD2S(Args->H[n]),CD2S(ExactH[n]));
      Warn("%s: H[%i]: (%s,%s)\n",CaseStrings[WhichCase],n,CD2S(Args->H[n]),CD2S(ExactH[n]));
      ErrorsDetected++;
    };

  if ( Args->NumGradientComponents > 0 )
   { cdouble *ExactGradH = ExactGradHValues + 6*WhichCase;
     for(int n=0; n<6; n++)
      if ( !EqualFloat(Args->GradH[n], ExactGradH[n]) )
       { Log("%s: GradH[%i]: (%s,%s)\n",CaseStrings[WhichCase],n,CD2S(Args->GradH[n]),CD2S(ExactGradH[n]));
         Warn("%s: GradH[%i]: (%s,%s)\n",CaseStrings[WhichCase],n,CD2S(Args->GradH[n]),CD2S(ExactGradH[n]));
         ErrorsDetected++;
       };
   };

  return (ErrorsDetected==0) ? 0 : 1;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  SetLogFileName("scuff-unit-tests.log");
  Log("SCUFF-EM PPI unit tests running on %s",GetHostName());

  int FailedCases=0;
  RWGGeometry *G = new RWGGeometry("PECSphere_R0P75_414.scuffgeo");
  G->SetLogLevel(SCUFF_VERBOSELOGGING);

  GetPPIArgStruct MyArgs, *Args=&MyArgs;
  InitGetPPIArgs(Args);

  Args->Sa = Args->Sb = G->Surfaces[0];

  SetDefaultCD2SFormat("(%+.8e,%+.8e)");

  /*--------------------------------------------------------------*/
  /*- no common vertices, long wavelength   ----------------------*/
  /*--------------------------------------------------------------*/
  Args->npa=0;
  Args->iQa=0;
  Args->npb=223;
  Args->iQb=0;
  Args->k=0.1;
  Args->NumGradientComponents=3;
  GetPanelPanelInteractions(Args);
  FailedCases += ComparePPIs(NCLW, Args);

  /*--------------------------------------------------------------*/
  /*- no common vertices, moderate wavelength   ------------------*/
  /*--------------------------------------------------------------*/
  Args->npa=0;
  Args->iQa=0;
  Args->npb=223;
  Args->iQb=0;
  Args->k=5.0;
  Args->NumGradientComponents=3;
  GetPanelPanelInteractions(Args);
  FailedCases+=ComparePPIs(NCMW, Args);

  /*--------------------------------------------------------------*/
  /*- no common vertices, short imaginary wavelength -------------*/
  /*--------------------------------------------------------------*/
  Args->npa=0;
  Args->iQa=0;
  Args->npb=223;
  Args->iQb=0;
  Args->k = 5.0 + 20.0*II;
  Args->NumGradientComponents=3;
  GetPanelPanelInteractions(Args);
  FailedCases+=ComparePPIs(NCSW, Args);

  /*--------------------------------------------------------------*/
  /*- common-vertex, long wavelength   ---------------------------*/
  /*--------------------------------------------------------------*/
  Args->npa=0;
  Args->iQa=0;
  Args->npb=65;
  Args->iQb=0;
  Args->k=0.1;
  Args->NumGradientComponents=0;
  GetPanelPanelInteractions(Args);
  FailedCases += ComparePPIs(CVLW, Args);

  /*--------------------------------------------------------------*/
  /*- common-vertex, moderate wavelength    ----------------------*/
  /*--------------------------------------------------------------*/
  Args->npa=0;
  Args->iQa=0;
  Args->npb=65;
  Args->iQb=0;
  Args->k=5.0;
  Args->NumGradientComponents=0;
  GetPanelPanelInteractions(Args);
  FailedCases+=ComparePPIs(CVMW, Args);

  /*--------------------------------------------------------------*/
  /*- common-vertex, short imaginary wavelength ------------------*/
  /*--------------------------------------------------------------*/
  Args->npa=0;
  Args->iQa=0;
  Args->npb=65;
  Args->iQb=0;
  Args->k = 5.0 + 20.0*II;
  Args->NumGradientComponents=0;
  GetPanelPanelInteractions(Args);
  FailedCases+=ComparePPIs(CVSW, Args);

  /*--------------------------------------------------------------*/
  /*- common-edge, long wavelength   -----------------------------*/
  /*--------------------------------------------------------------*/
  Args->npa=0;
  Args->iQa=0;
  Args->npb=66;
  Args->iQb=0;
  Args->k=0.1;
  Args->NumGradientComponents=0;
  GetPanelPanelInteractions(Args);
  FailedCases += ComparePPIs(CELW, Args);

  /*--------------------------------------------------------------*/
  /*- common-edge, moderate wavelength    ------------------------*/
  /*--------------------------------------------------------------*/
  Args->npa=0;
  Args->iQa=0;
  Args->npb=66;
  Args->iQb=0;
  Args->k=5.0;
  Args->NumGradientComponents=0;
  GetPanelPanelInteractions(Args);
  FailedCases+=ComparePPIs(CEMW, Args);

  /*--------------------------------------------------------------*/
  /*- common-edge, short imaginary wavelength --------------------*/
  /*--------------------------------------------------------------*/
  Args->npa=0;
  Args->iQa=0;
  Args->npb=66;
  Args->iQb=0;
  Args->k = 5.0 + 20.0*II;
  Args->NumGradientComponents=0;
  GetPanelPanelInteractions(Args);
  FailedCases+=ComparePPIs(CESW, Args);

  /*--------------------------------------------------------------*/
  /*- common-triangle, long wavelength   -------------------------*/
  /*--------------------------------------------------------------*/
  Args->npa=0;
  Args->iQa=0;
  Args->npb=0;
  Args->iQb=0;
  Args->k=0.1;
  Args->NumGradientComponents=0;
  GetPanelPanelInteractions(Args);
  FailedCases += ComparePPIs(CTLW, Args);

  /*--------------------------------------------------------------*/
  /*- common-triangle, moderate wavelength    --------------------*/
  /*--------------------------------------------------------------*/
  Args->npa=0;
  Args->iQa=0;
  Args->npb=0;
  Args->iQb=0;
  Args->k=5.0;
  Args->NumGradientComponents=0;
  GetPanelPanelInteractions(Args);
  FailedCases+=ComparePPIs(CTMW, Args);

  /*--------------------------------------------------------------*/
  /*- common-triangle, short imaginary wavelength ----------------*/
  /*--------------------------------------------------------------*/
  Args->npa=0;
  Args->iQa=0;
  Args->npb=0;
  Args->iQb=0;
  Args->k = 5.0 + 20.0*II;
  Args->NumGradientComponents=0;
  GetPanelPanelInteractions(Args);
  FailedCases+=ComparePPIs(CTSW, Args);

  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 
  if (FailedCases>0)
   abort();

  printf("All tests successfully passed.\n");
  return 0;

}
