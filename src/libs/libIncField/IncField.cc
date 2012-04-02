/*
 * IncField.cc -- implementation of the IncField base class,  
 *                    together with the EHIncField() non-class method
 *
 * homer reid      -- 1/2012
 */

#include <string.h>
#include <stdlib.h>

#include "libIncField.h"

/***************************************************************/
/***************************************************************/

IncField::~IncField() {
  free(Object);
}

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

void IncField::SetObject(const char *Label) {
  if (Object) free(Object);
  Object = Label ? strdup(Label) : 0;
  if (!Label) ObjectIndex = -1; // exterior medium is always index == -1
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

/* Function for passing to AssembleRHSVector.  The ObjectIndex
   for each IncField should be set to the medium that it is in;
   if this == exterior_index, the field is added, if it ==
   interior_index, the field is subtracted, and otherwise it is 0. */

void EHIncField2(const double X[3], void *UserData, cdouble EH[6],
		int exterior_index, int interior_index)
{
  memset(EH, 0, 6*sizeof(cdouble));
  for(IncField *IFD=(IncField *)UserData; IFD; IFD=IFD->Next) {
    cdouble PEH[6]; 
    if (IFD->ObjectIndex == exterior_index) {
      IFD->GetFields(X, PEH);
      for(int nc=0; nc<6; nc++)
	EH[nc] += PEH[nc];
    }
    else if (IFD->ObjectIndex == interior_index) {
      IFD->GetFields(X, PEH);
      for(int nc=0; nc<6; nc++)
	EH[nc] -= PEH[nc];
    }
   }
}

// unlike EHIncField, return the total field always
void EHIncField(const double X[3], void *UserData, cdouble EH[6])
{
  IncField *IFD=(IncField *)UserData;
  IFD->GetTotalFields(X, EH);
}
