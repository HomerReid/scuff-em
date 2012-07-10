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
 * PBCGeometry.h  -- header file for libscuff implementation of 
 *                -- periodic boundary conditions
 *
 * homer reid  -- 4/2011 -- 7/2012
 */

#ifndef PBCGEOMETRY_H
#define PBCGEOMETRY_H

/***************************************************************/
/* PBCGeometry is a class built atop RWGGeometry for handling  */
/* periodic boundary conditions.                               */
/***************************************************************/
class PBCGeometry
 {
   /*--------------------------------------------------------------*/
   /* public class methods                                        -*/
   /*--------------------------------------------------------------*/
public:

   // constructor / destructor 
   PBCGeometry::PBCGeometry(RWGGeometry *G, double LBV[2][2]);
   ~PBCGeometry::PBCGeometry(RWGGeometry *G, double LBV[2][2]);

   // assemble BEM matrix
  HMatrix *AssembleBEMMatrix(cdouble Omega, double P[2], HMatrix *M=0);

  // assemble RHS vector; this just devolves to a call to the usual
  // RWGGeometry routine for assembling the BEM matrix, but i will 
  // make it a PBCGeometry routine for completeness
  HVector *AssembleRHSVector(cdouble Omega, IncField *IF, HVector *RHS)
   { return G->AssembleRHSVector(Omega, IF, RHS); }

   // get scattered fields
  void GetFields(double *X, IncField *IF, cdouble *EH);

   /*--------------------------------------------------------------*/
   /*- class data (would be private if we were fastidious about    */
   /*-             such things)                                    */
   /*--------------------------------------------------------------*/
   RWGGeometry *G;       // unit cell geometry

   double LBV[2][2];     // lattice basis vectors 

   // NumStraddlers[ 2*no+i ] = the number of basis functions on object #no
   //                           that straddle unit-cell boundary #i
   int *NumStraddlers;  

   // contributions to the BEM matrix from lattice sites 
   // (n1,n2) = { (1,1), (1,-1), (1,0), (0,1), (0,0) }
   HMatrix *MPP, *MPM, *MPZ, *MZP, *MZZ;

   // this field stores the value of Omega that was passed to 
   // the most recent invocation of AssembleBEMMatrix(). we store it 
   // because, if the user makes a second call to AssembleBEMMatrix() 
   // with the same value of Omega (but presumably a different value
   // of the bloch vector P), then we can reuse the MPP...MZZ matrix 
   // blocks. 
   cdouble CurrentOmega;

   /* interpolation tables to accelerate the calculation of */
   /* the periodic green's function for the exterior medium */
   /* and for the medium interior to each object            */
   Interp3D *GBarAB9_Exterior;
   Interp3D **GBarAB9_Interior;

   /*--------------------------------------------------------------*/
   /*- class methods which would be private if we were fastidious  */
   /*- about such things; these are helper routines for the public */
   /*- class methods                                               */
   /*--------------------------------------------------------------*/
   void GetInnerContributions();
   void AddOuterContributions();

 };

/***************************************************************/
/* this is really a helper routine for the PBCGeometry class   */
/* constructor, but i am going to make it a standalone routine */
/* for now                                                     */
/***************************************************************/
void AddStraddlers(RWGObject *O, double **LBV, int NumStraddlers[2]);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GBarVDEwald(double R, cdouble k, double *P, double *LBV[2],
                 double E, int ExcludeFirst9, cdouble *GBarVD);

typedef struct GBarData 
 { 
   cdouble k;           // wavenumber 
   double P[2];         // bloch vector 
   double **LBV;        // lattice basis vectors 
   double E;            // separation parameter
   bool ExcludeInner9;  
 
 } GBarData;

void GBarVDPhi3D(double X1, double X2, double X3, 
                 void *UserData, double *PhiVD);

#endif
