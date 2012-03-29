/*
 * IncField.cc -- implementation of the IncField base class,  
 *                    together with the EHIncField() non-class method
 *
 * homer reid      -- 1/2012
 */

#include <string.h>

#include "libIncField.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
void IncField::SetFrequency(cdouble pOmega)
{
  IncField *IFD;
  for(IFD=this; IFD; IFD=IFD->Next)
   IFD->Omega=pOmega;
}

void IncField::SetFrequencyAndEpsMu(cdouble pOmega, cdouble pEps, cdouble pMu)
{
  IncField *IFD;
  for(IFD=this; IFD; IFD=IFD->Next)
   { IFD->Omega=pOmega;
     IFD->Eps=pEps;
     IFD->Mu=pMu;
   };
}

/***************************************************************/

void IncField::GetTotalFields(const double X[3], cdouble EH[6])
{
  memset(EH, 0, 6*sizeof(cdouble));
  for(IncField *IFD=this; IFD; IFD=IFD->Next)
   { 
     cdouble PEH[6]; 
     IFD->GetFields(X, PEH);
     for(int nc=0; nc<6; nc++)
      EH[nc] += PEH[nc];
   };

}

/***************************************************************/
/***************************************************************/
void EHIncField(const double X[3], void *UserData, cdouble EH[6])
{
 
  IncField *IFD=(IncField *)UserData;
  IFD->GetTotalFields(X, EH);
}
