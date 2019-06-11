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
 * libIncField.h -- routines for evaluating the electric and magnetic
 *               -- fields of various known types of field configurations
 *
 * homer reid    -- 11/2009
 */
#ifndef LIBINCFIELD_H
#define LIBINCFIELD_H

#include <complex>
#include <cmath>
#include "libhmat.h"

#ifndef cdouble
  typedef std::complex<double> cdouble;
#endif

#ifndef ZVAC
#define ZVAC 376.73031346177   // impedance of free space
#endif
  // Omega is initialized to an absurd value to help catch cases
  // in which GetFields() is called without a prior call to 
  // SetFrequency()
/**********************************************************************/
/* IncField is a general base class from which specific classes       */
/* for various types of field are derived.                            */
/**********************************************************************/
class IncField
 { 
 public:

   cdouble Omega;
   cdouble Eps;
   cdouble Mu;

   // lattice basis and Bloch vector for Bloch-periodic case
   // (we also store a basis for the reciprocal lattice)
   HMatrix *LBasis, *RLBasis;
   double LVolume, RLVolume;
   double kBloch[3];

   IncField *Next;

   // label and/or index for the region within which the field sources 
   // lie (NULL == "EXTERIOR")
   char *RegionLabel; 
   int RegionIndex;

   // constructor just initializes the simple fields
   IncField(); 
   IncField(const IncField &IF);
   virtual ~IncField();

   void SetFrequency(cdouble Omega, bool Traverse=true);
   void SetFrequencyAndEpsMu(cdouble Omega, cdouble Eps, cdouble Mu, bool Traverse=true);
   void SetRegionLabel(const char *Label = 0);

   void SetLattice(HMatrix *LBasis, bool Traverse=true);
   void SetkBloch(const double *NewkBloch, bool Traverse=true);

   // obsolete calling convention retained for backward compatibility
   void SetObjectLabel(const char *Label) { SetRegionLabel(Label); }

   // if true, returns the location of the source, used to
   // set the RegionIndex if !RegionLabel.
   virtual bool GetSourcePoint(double X[3]) const { (void) X; return false; }
   
   virtual void GetFields(const double X[3], cdouble EH[6]) = 0 ;
   void GetTotalFields(const double X[3], cdouble EH[6]);

   // the default implementation of this routine uses finite-differencing;
   // subclasses may override it in cases where they know how to compute
   // field gradients directly
   virtual void GetFieldGradients(const double X[3], cdouble dEH[3][6]);
 };

/**********************************************************************/
/* non-class-method function to delete an entire linked chain of IncFields */
/**********************************************************************/
void DeleteIncFieldChain(IncField *IF);

/**********************************************************************/
/* Next come the various possible types of incident field, implemented*/
/* as structs derived from IncField.                                  */
/**********************************************************************/

/**********************************************************************/
/* plane wave *********************************************************/
/**********************************************************************/
class PlaneWave : public IncField
 { 
 public:
   cdouble E0[3];         /* E-field polarization vector */
   double nHat[3];        /* unit vector in direction of propagation */

   PlaneWave(const cdouble E0[3], const double nHat[3], const char *Label = 0);
   PlaneWave(const PlaneWave &PW);
   PlaneWave();
   ~PlaneWave();

   void SetE0(cdouble pE0[3]);
   void SetnHat(double nHat[3]);

   void GetFields(const double X[3], cdouble EH[6]);
   void GetFieldGradients(const double X[3], cdouble dEH[3][6]);

 };

/**********************************************************************/
/* point dipole source ************************************************/
/**********************************************************************/  
#define LIF_ELECTRIC_DIPOLE 0
#define LIF_MAGNETIC_DIPOLE 1
class PointSource: public IncField
 { 
 public:
   double X0[3];         /* location */
   cdouble P[3];         /* strength */
   int Type;             /* LIF_ELECTRIC_DIPOLE or LIF_MAGNETIC_DIPOLE */

   PointSource(const double X0[3], const cdouble P[3], 
               int Type = LIF_ELECTRIC_DIPOLE, const char *Label = 0);
   PointSource(const PointSource &PS);
   PointSource();
   ~PointSource();

   void InitPointSource(const double X0[3], const cdouble P[3],
                        int Type, const char *Label);
   void SetX0(double X0[3]);
   void SetP(cdouble P[3]);
   void SetType(int pType);

   void GetFields(const double X[3], cdouble EH[6]);
   void GetFields_Periodic(const double X[3], cdouble EH[6]);
   void Get2DPeriodicFields_Fourier(const double X[3], cdouble EH[6]);

   bool GetSourcePoint(double X[3]) const;

   bool UseEwaldFields;
 };

/**********************************************************************/
/* focused gaussian beam  **********8**********************************/
/**********************************************************************/
class GaussianBeam: public IncField
 { 
 public:
   double X0[3];            /* beam center point */
   double KProp[3];         /* beam propagation vector */
   cdouble E0[3];           /* complex field-strength vector */
   double W0;               /* beam waist */

   // constructor 
   GaussianBeam(const double X0[3], const double KProp[3], 
                const cdouble E0[3], double W0, const char *Label = 0);
   ~GaussianBeam();

   void SetX0(double pX0[3]);
   void SetKProp(double pKProp[3]);
   void SetE0(cdouble pE0[3]);
   void SetW0(double pW0);

   void GetFields(const double X[3], cdouble EH[6]);

   double TotalBeamFlux();

 };

/**********************************************************************/
/* vector spherical wave **********************************************/
/* note: SW_MAGNETIC and SW_ELECTRIC labelt the functions commonly    */
/* \mathbf{M}_{\ell m} and \mathbf{N}_{\ell m} respectively; in       */
/* particular, \mathbf{M} has no radial component.                    */
/**********************************************************************/
#define SW_MAGNETIC 0
#define SW_ELECTRIC 1
class SphericalWave : public IncField
 { 
 public:
   int L, M, P;     // spherical wave indices; P=0,1 for M,N-type wave

   SphericalWave(int L=1, int M=0, int P=0);

   void SetL(int NewL);
   void SetM(int NewM);
   void SetP(int NewP);
   void SetType(int NewP); // DELETEME legacy routine for backward compat

   void GetFields(const double X[3], cdouble EH[6]);

 };


/**********************************************************************/
/* magnetic 'frill' (annulus)    **************************************/
/**********************************************************************/
#if 0
class MagneticFrill: public IncField
 {
 public:
   double X0[3];          /* center of frill */
   double Theta, Phi;     /* angles of frill axis */
   double RIn, ROut;      /* inner and outer radii */
 } MagneticFrill;

/**********************************************************************/
/* magnetic solenoid                            ***********************/
/**********************************************************************/
class MagneticSolenoid: public IncField
 { 
 public:
   double X0[3];          /* center */
   double Theta, Phi;     /* angles of solenoid axis */
   double Radius, L;      /* radius and length  */
 } MagneticSolenoid;
#endif

/**********************************************************************/
/* 20160522 new functionality for parsing an input text file to yield */
/* a list of incident fields; used to implement "advanced mode" in    */
/* scuff-scatter and buff-scatter                                     */
/**********************************************************************/
typedef struct IncFieldList
 { IncField **IFs;
   char **Labels;
   int NumIFs;
 } IncFieldList;

IncFieldList *ReadIncFieldList(char *FileName);
IncFieldList *AddIncFieldToList(IncField *IF, char *Label=0, IncFieldList *IFList=0); 

#endif
