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
IncField::IncField()
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
