/*
 * IncFieldData.cc -- implementation of the IncFieldData base class,  
 *                    together with the EHIncField() non-class method
 *
 * homer reid      -- 1/2012
 */

#include <string.h>

#include "libIncField.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
void IncFieldData::SetFrequency(cdouble Omega)
{
  IncFieldData *IFD;
  for(IFD=this; IFD; IFD=IFD->Next)
   IFD->Omega=Omega;
}

void IncFieldData::SetFrequencyAndEpsMu(cdouble Omega, cdouble Eps, double Mu)
{
  IncFieldData *IFD;
  for(IFD=this; IFD; IFD=IFD->Next)
   { IFD->Omega=Omega;
     IFD->Eps=Eps;
     IFD->Mu=Mu;
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void EHIncField(double *X, void *UserData, cdouble EH[6])
{
  memset(EH, 0, 6*sizeof(cdouble));
 
  IncFieldData *IFD, *IFDList=(IncFieldData *)UserData;
  cdouble PEH[6]; 
  int nc;
  for(IFD=IFDList; IFD; IFD=IFD->Next)
   { 
     IFD->GetFields(X, PEH);
     for(nc=0; nc<6; nc++)
      EH[nc] += PEH[6];
   };
}
