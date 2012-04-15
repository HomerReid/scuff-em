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

#ifndef cdouble
  typedef std::complex<double> cdouble;
#endif

#ifndef ZVAC
#define ZVAC 376.73031346177   // impedance of free space
#endif

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
   IncField *Next;

   // label and/or index for the object within which the field sources 
   // lie (NULL == "EXTERIOR")
   char *ObjectLabel; 
   int ObjectIndex;

   // constructor just initializes the simple fields
   IncField(); 
   ~IncField();

   void SetFrequency(cdouble Omega, bool Traverse=true);
   void SetFrequencyAndEpsMu(cdouble Omega, cdouble Eps, cdouble Mu, bool Traverse=true);
   void SetObjectLabel(const char *Label = 0);

   // if true, returns the location of the source, used to
   // set the ObjectIndex if !ObjectLabel.
   virtual bool GetSourcePoint(double X[3]) const { (void) X; return false; }
   
   virtual void GetFields(const double X[3], cdouble EH[6]) = 0 ;
   void GetTotalFields(const double X[3], cdouble EH[6]);
 };

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

   void GetFields(const double X[3], cdouble EH[6]);

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

   void GetFields(const double X[3], cdouble EH[6]);
   bool GetSourcePoint(double X[3]);
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

   void GetFields(const double X[3], cdouble EH[6]);

   double TotalBeamFlux();

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

#endif
