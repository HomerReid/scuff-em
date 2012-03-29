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
void IncFieldData::SetFrequency(cdouble pOmega)
{
  IncFieldData *IFD;
  for(IFD=this; IFD; IFD=IFD->Next)
   IFD->Omega=pOmega;
}

void IncFieldData::SetFrequencyAndEpsMu(cdouble pOmega, cdouble pEps, cdouble pMu)
{
  IncFieldData *IFD;
  for(IFD=this; IFD; IFD=IFD->Next)
   { IFD->Omega=pOmega;
     IFD->Eps=pEps;
     IFD->Mu=pMu;
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void EHIncField(const double X[3], void *UserData, cdouble EH[6])
{
  memset(EH, 0, 6*sizeof(cdouble));
 
  IncFieldData *IFD, *IFDList=(IncFieldData *)UserData;
  cdouble PEH[6]; 
  int nc;
  for(IFD=IFDList; IFD; IFD=IFD->Next)
   { 
     IFD->GetFields(X, PEH);
     for(nc=0; nc<6; nc++)
      EH[nc] += PEH[nc];
   };
}
