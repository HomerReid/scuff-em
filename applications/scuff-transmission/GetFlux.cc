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
 * GetFlux.cc -- compute upward-traveling and downward-traveling power flux
 *               by integrating the Poynting vector over the unit cell
 *               area at points above and below the structure
 *
 * homer reid -- 9/2012
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libIncField.h>
#include <libscuff.h>

using namespace scuff;

#define II cdouble (0.0, 1.0)
#define PWFLUX (1.0/(2.0*ZVAC)) // power flux of a plane wave in vacuum

/*******************************************************************/
/* this is a 9th-order, 17-point cubature rule for the unit square */
/* with corners {(0,0) (1,0) (1,1) (1,0)}.                         */
/* array entries:                                                  */
/*  x_0, y_0, w_0,                                                 */
/*  x_1, y_1, w_1,                                                 */
/*  ...                                                            */
/*  x_16, y_16, w_16                                               */
/* where (x_n, y_n) and w_n are the nth cubature point and weight. */
/*******************************************************************/
double SCR9[]={
  +5.0000000000000000e-01, +5.0000000000000000e-01, +1.3168724279835392e-01,
  +9.8442498318098881e-01, +8.1534005986583447e-01, +2.2219844542549678e-02,
  +9.8442498318098881e-01, +1.8465994013416559e-01, +2.2219844542549678e-02,
  +1.5575016819011134e-02, +8.1534005986583447e-01, +2.2219844542549678e-02,
  +1.5575016819011134e-02, +1.8465994013416559e-01, +2.2219844542549678e-02,
  +8.7513854998945029e-01, +9.6398082297978482e-01, +2.8024900532399120e-02,
  +8.7513854998945029e-01, +3.6019177020215176e-02, +2.8024900532399120e-02,
  +1.2486145001054971e-01, +9.6398082297978482e-01, +2.8024900532399120e-02,
  +1.2486145001054971e-01, +3.6019177020215176e-02, +2.8024900532399120e-02,
  +7.6186791010721466e-01, +7.2666991056782360e-01, +9.9570609815517519e-02,
  +7.6186791010721466e-01, +2.7333008943217640e-01, +9.9570609815517519e-02,
  +2.3813208989278534e-01, +7.2666991056782360e-01, +9.9570609815517519e-02,
  +2.3813208989278534e-01, +2.7333008943217640e-01, +9.9570609815517519e-02,
  +5.3810416409630857e-01, +9.2630786466683113e-01, +6.7262834409945196e-02,
  +5.3810416409630857e-01, +7.3692135333168873e-02, +6.7262834409945196e-02,
  +4.6189583590369143e-01, +9.2630786466683113e-01, +6.7262834409945196e-02,
  +4.6189583590369143e-01, +7.3692135333168873e-02, +6.7262834409945196e-02 
};

/*******************************************************************/
/* get the upward-traveling and downward-traveling power fluxes    */
/* through the planes at zAbove and zBelow by integrating the      */
/* poynting vector over the area of the unit cell.                 */
/*                                                                 */
/* Note: Flux values are normalized by the flux of an unit-strength*/
/* plane wave in vacuum.                                           */
/*                                                                 */
/* Return values: Flux[0] = upward-traveling flux at ZAbove        */
/*                Flux[1] = downward-traveling flux at ZBelow      */
/*                                                                 */
/* If nHatAbove/Below are non-null, the flux in those directions   */
/* is computed, i.e. nHatFlux \cdot E^* \times H                   */
/*******************************************************************/
#ifndef NUMPOLS
#define NUMPOLS 2
#endif
void GetFlux(RWGGeometry *G, IncField *IF, HVector *KN,
             cdouble Omega, double *kBloch, int NQPoints,
             double ZAbove, double ZBelow, bool FromAbove,
             PlaneWave *ReflectedPW[NUMPOLS],
             PlaneWave *TransmittedPW[2],
             double *Flux, cdouble tIntegral[NUMPOLS], cdouble rIntegral[NUMPOLS],
             double *nHatAbove, double *nHatBelow)
{
  double *SCR=SCR9;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NCP; // number of cubature points
  if (NQPoints==0)
   NCP=17; 
  else
   NCP=NQPoints*NQPoints;

  /***************************************************************/
  /* on the first invocation we allocate space for the matrices  */
  /* of evaluation points and fields.                            */
  /***************************************************************/
  static HMatrix *XMatrixAbove = 0, *XMatrixBelow = 0;
  static HMatrix *FMatrixAbove = 0, *FMatrixBelow = 0;
  if (XMatrixAbove==0)
   { XMatrixAbove = new HMatrix(NCP, 3 );
     XMatrixBelow = new HMatrix(NCP, 3 );
     FMatrixAbove = new HMatrix(NCP, 6, LHM_COMPLEX);
     FMatrixBelow = new HMatrix(NCP, 6, LHM_COMPLEX);
   }

  /***************************************************************/
  /* fill in coordinates of evaluation points.                   */
  /* the first NCP points are for the upper surface; the next    */
  /* NCP points are for the lower surface.                       */
  /***************************************************************/
  double x, y, LBV[2][3];
  for(int nd=0; nd<2; nd++)
   for(int nc=0; nc<3; nc++)
    LBV[nd][nc] = G->LBasis->GetEntryD(nc,nd);
  if (NQPoints==0)
   for(int ncp=0; ncp<NCP; ncp++)
    { 
      x=SCR[3*ncp+0];
      y=SCR[3*ncp+1];

      XMatrixAbove->SetEntry(ncp, 0, x*LBV[0][0] + y*LBV[1][0]);
      XMatrixAbove->SetEntry(ncp, 1, x*LBV[0][1] + y*LBV[1][1]);
      XMatrixAbove->SetEntry(ncp, 2, ZAbove);

      XMatrixBelow->SetEntry(ncp, 0, x*LBV[0][0] + y*LBV[1][0]);
      XMatrixBelow->SetEntry(ncp, 1, x*LBV[0][1] + y*LBV[1][1]);
      XMatrixBelow->SetEntry(ncp, 2, ZBelow);

    }
  else
   for(int nqpx=0, ncp=0; nqpx<NQPoints; nqpx++)
    for(int nqpy=0; nqpy<NQPoints; nqpy++, ncp++)
     { 
       x = ((double)nqpx + 0.5) / ((double)NQPoints);
       y = ((double)nqpy + 0.5) / ((double)NQPoints);

       XMatrixAbove->SetEntry(ncp, 0, x*LBV[0][0] + y*LBV[1][0]);
       XMatrixAbove->SetEntry(ncp, 1, x*LBV[0][1] + y*LBV[1][1]);
       XMatrixAbove->SetEntry(ncp, 2, ZAbove);

       XMatrixBelow->SetEntry(ncp, 0, x*LBV[0][0] + y*LBV[1][0]);
       XMatrixBelow->SetEntry(ncp, 1, x*LBV[0][1] + y*LBV[1][1]);
       XMatrixBelow->SetEntry(ncp, 2, ZBelow);
     }

  /***************************************************************/
  /* get scattered fields at all cubature points                 */
  /***************************************************************/
  G->GetFields(IF, KN, Omega, kBloch, XMatrixAbove, FMatrixAbove);
  G->GetFields(IF, KN, Omega, kBloch, XMatrixBelow, FMatrixBelow);

  /***************************************************************/
  /* integrate poynting vector over upper and lower surfaces.    */
  /* Note: The jacobian in this cubature is the area of the unit */
  /*       cell, so omitting that factor is equivalent to        */
  /*       dividing the integrated power by the unit-cell area,  */
  /*       which is what we want to do anyway.                   */
  /***************************************************************/
  double PAbove=0.0, PBelow=0.0;
  cdouble rDenominator[2]={0.0, 0.0};
  cdouble tDenominator[2]={0.0, 0.0};
  tIntegral[0]=tIntegral[1]=rIntegral[0]=rIntegral[1]=0.0;
  double pzHat[3]={0,0,1.0}, mzHat[3]={0,0,-1.0};
  if (nHatAbove==0) nHatAbove=pzHat;
  if (nHatBelow==0) nHatBelow=mzHat;
  for(int ncp=0; ncp<NCP; ncp++)
   {
     double w;
     if (NQPoints==0) 
      w=SCR[3*ncp+2];  // cubature weight
     else
      w=1.0/((double)(NCP));

     cdouble *E, *H;
     FMatrixAbove->GetEntries(ncp, "0:2", E);
     FMatrixAbove->GetEntries(ncp, "3:5", H);
     PAbove += 0.5*w*real(   nHatAbove[0]*(E[1]*conj(H[2]) - E[2]*conj(H[1]) )
                           + nHatAbove[1]*(E[2]*conj(H[0]) - E[0]*conj(H[2]) )
                           + nHatAbove[2]*(E[0]*conj(H[1]) - E[1]*conj(H[0]) )
                         );

     FMatrixBelow->GetEntries(ncp, "0:2", E);
     FMatrixBelow->GetEntries(ncp, "3:5", H);
     PBelow += 0.5*w*real(   nHatBelow[0]*(E[1]*conj(H[2]) - E[2]*conj(H[1]) )
                           + nHatBelow[1]*(E[2]*conj(H[0]) - E[0]*conj(H[2]) )
                           + nHatBelow[2]*(E[0]*conj(H[1]) - E[1]*conj(H[0]) )
                         );

     double XSource[3],  XDest[3];
     cdouble EHSource[6], EHDest[6];
     cdouble EHTE[6], EHTM[6];

     if (FromAbove)
      { XMatrixAbove->GetEntriesD(ncp,"0:2",XSource);
        FMatrixAbove->GetEntries(ncp,"0:5",EHSource);
        XMatrixBelow->GetEntriesD(ncp,"0:2",XDest);
        FMatrixBelow->GetEntries(ncp,"0:5",EHDest);
      }
     else
      { XMatrixBelow->GetEntriesD(ncp,"0:2",XSource);
        FMatrixBelow->GetEntries(ncp,"0:5",EHSource);
        XMatrixAbove->GetEntriesD(ncp,"0:2",XDest);
        FMatrixAbove->GetEntries(ncp,"0:5",EHDest);
      }
      
     ReflectedPW[POL_TE]->GetFields(XSource, EHTE);
     ReflectedPW[POL_TM]->GetFields(XSource, EHTM);
     rIntegral[POL_TE]    += w*VecHDot(EHSource,  EHTE, 3);
     rIntegral[POL_TM]    += w*VecHDot(EHSource,  EHTM, 3);
     rDenominator[POL_TE] += w*VecHDot(EHTE, EHTE, 3);
     rDenominator[POL_TM] += w*VecHDot(EHTM, EHTM, 3);
      
     TransmittedPW[POL_TE]->GetFields(XDest, EHTE);
     TransmittedPW[POL_TM]->GetFields(XDest, EHTM);
     tIntegral[POL_TE]    += w*VecHDot(EHDest,  EHTE, 3);
     tIntegral[POL_TM]    += w*VecHDot(EHDest,  EHTM, 3);
     tDenominator[POL_TE] += w*VecHDot(EHTE, EHTE, 3);
     tDenominator[POL_TM] += w*VecHDot(EHTM, EHTM, 3);
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Flux[0]=PAbove / PWFLUX;
  Flux[1]=PBelow / PWFLUX;
  rIntegral[POL_TE] /= rDenominator[POL_TE];
  rIntegral[POL_TM] /= rDenominator[POL_TM];
  tIntegral[POL_TE] /= tDenominator[POL_TE];
  tIntegral[POL_TM] /= tDenominator[POL_TM];

}
