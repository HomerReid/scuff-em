/*
 * IncField.cc -- implementation of the IncField base class
 *
 * homer reid      -- 1/2012
 */

#include <string.h>
#include <stdlib.h>

#include "libIncField.h"

/***************************************************************/
/* base class constructor **************************************/
/***************************************************************/
IncField::~IncField()
{ 
  Eps=cdouble(1.0,0.0);
  Mu=cdouble(1.0,0.0);

  // Omega is initialized to an absurd value to help catch cases
  // in which GetFields() is called without a prior call to 
  // SetFrequency()
  Omega=-1.0;

  Next = NULL;

  // field sources lie in the exterior medium by default
  ObjectLabel = NULL;
  ObjectIndex = -1;
}

/***************************************************************/
/* base class destructor ***************************************/
/***************************************************************/
IncField::~IncField()
{
  free(ObjectLabel);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void IncField::SetFrequency(cdouble pOmega, bool Traverse)
{
  Omega=pOmega;

  if (Traverse)
   for(IncField *IFD=this->Next; IFD; IFD=IFD->Next)
    IFD->Omega=pOmega;
}

void IncField::SetFrequencyAndEpsMu(cdouble pOmega, 
                                    cdouble pEps, cdouble pMu,
                                    bool Traverse)
{
  Omega=pOmega;
  Eps=pEps;
  Mu=pMu;
 
  if (Traverse)
   for(IncField *IFD=this->Next; IFD; IFD=IFD->Next)
    { 
      IFD->Omega=pOmega;
      IFD->Eps=pEps;
      IFD->Mu=pMu;
    };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void IncField::SetObjectLabel(const char *Label) 
{
  if (ObjectLabel) free(ObjectLabel);
  ObjectLabel = Label ? strdup(Label) : 0;
  if (!Label) ObjectIndex = -1; // exterior medium is always index == -1
}

/***************************************************************/
/***************************************************************/
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
