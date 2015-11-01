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
#include "libhrutil.h"

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

  // incident fields are non-periodic by default
  LDim=0;
  LBV[0]=LBV[1]=0;
  kBloch[0]=kBloch[1]=0.0;

  // field sources lie in the exterior region by default
  RegionLabel = NULL;
  RegionIndex = 0;

  Next = NULL;
}

/***************************************************************/
/* base class destructor ***************************************/
/***************************************************************/
IncField::~IncField()
{
  if (RegionLabel)
   free(RegionLabel);
}

/***************************************************************/
/* non-class-method function to delete an entire linked chain  */
/* of IncFields                                                */
/***************************************************************/
void DeleteIncFieldChain(IncField *IF)
{ 
  IncField *ThisIF=IF;
  while (ThisIF)
   { IncField *NextIF=ThisIF->Next;
     delete ThisIF;
     ThisIF=NextIF;
   };
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
void IncField::SetLattice(int NewLDim, double NewLBV[2][2], bool Traverse)
{
  if (LDim<0 || LDim>2)
   ErrExit("%s: %i: periodic lattices must have dimension {0,1,2}",
            __FILE__,__LINE__);

  LDim = NewLDim;
  if (LDim>=1)
   { LBV1[0] = NewLBV[0][0];
     LBV1[1] = NewLBV[0][1];
     LBV[0] = LBV1;
   };
  if (LDim==2)
   { LBV2[0] = NewLBV[1][0];
     LBV2[1] = NewLBV[1][1];
     LBV[1] = LBV2;
   };

  if (Traverse)
   for(IncField *IFD=this->Next; IFD; IFD=IFD->Next)
    IFD->SetLattice(NewLDim,NewLBV,false);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void IncField::SetkBloch(double *NewkBloch, bool Traverse)
{ 
  if ( LDim==0 )
   ErrExit("%s:%i: attempt to set kBloch in a non-periodic IncField");

  if (NewkBloch)
   memcpy(kBloch, NewkBloch, LDim*sizeof(double));

  if (Traverse)
   for(IncField *IFD=this->Next; IFD; IFD=IFD->Next)
    memcpy(IFD->kBloch, NewkBloch, LDim*sizeof(double));
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void IncField::SetRegionLabel(const char *Label) 
{
  if (RegionLabel) free(RegionLabel);
  RegionLabel = Label ? strdupEC(Label) : 0;
  if (!Label) RegionIndex = 0; // exterior medium is always index == 0
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

/***************************************************************/
/* get field gradients by finite-differencing; this method may */
/* be overridden by subclasses who know how to compute their   */
/* field gradients analytically                                */
/* dEH[i][0+j] = \partial_i E_j                                */
/* dEH[i][3+j] = \partial_i H_j                                */
/***************************************************************/
void IncField::GetFieldGradients(const double X[3], cdouble dEH[3][6])
{
  cdouble EH[6], EHP[6], EHM[6];
  GetFields(X, EH);
  for(int i=0; i<3; i++)
   { 
     double xTweaked[3];
     xTweaked[0]=X[0];
     xTweaked[1]=X[1];
     xTweaked[2]=X[2];

     double Delta = (X[i]==0.0 ? 1.0e-4 : 1.0e-4*fabs(X[i]));

     xTweaked[i] += Delta;
     GetFields(xTweaked, EHP);
     xTweaked[i] -= 2.0*Delta;
     GetFields(xTweaked, EHM);
     for(int j=0; j<6; j++)
      dEH[i][j] = (EHP[j]-EHM[j])/(2.0*Delta);
   };
}
