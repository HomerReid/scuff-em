/*
 * libIncField.h -- routines for evaluating the electric and magnetic
 *               -- fields of various known types of field configurations
 *
 *               -- documentation at
 *               --  http://homerreid.com/scuff-EM/libIncField
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
/* IncFieldData is a general base class from which specific classes   */
/* for various types of field are derived.                            */
/**********************************************************************/
typedef class IncFieldData
 { 
 public:
   cdouble Omega;
   cdouble Eps;
   double Mu;
   IncFieldData *Next;

   // constructor just initializes the simple fields 
   IncFieldData() : Eps(1.0,0.0), Mu(1.0), Next(0) {}

   void SetFrequency(cdouble Omega);
   void SetFrequencyAndEpsMu(cdouble Omega, cdouble Eps, double Mu);
   
   virtual void GetFields(double *X, cdouble *EH) = 0 ;

 } IncFieldData;


/**********************************************************************/
/* this is the function with the proper prototype for passage to the  */
/* AssembleRHSVector() routine in libscuff. its implementation within */
/* libIncField expects the UserData field to point to the head of a   */
/* linked list of structures derived from IncFieldData. the field     */
/* returned is the sum of the fields associated with each structure in*/
/* the list.                                                          */
/**********************************************************************/
void EHIncField(double *X, void *UserData, cdouble EH[6]);

/**********************************************************************/
/* Next come the various possible types of incident field, implemented*/
/* as structs derived from IncFieldData.                              */
/**********************************************************************/

/**********************************************************************/
/* plane wave *********************************************************/
/**********************************************************************/
struct PlaneWaveData : public IncFieldData
 { 
   cdouble E0[3];         /* E-field polarization vector */
   double nHat[3];        /* unit vector in direction of propagation */

   // constructor 
   PlaneWaveData(cdouble E0[3], double nHat[3]);
   void GetFields(double *X, cdouble *EH);

 };

/**********************************************************************/
/* point dipole source ************************************************/
/**********************************************************************/  
#if 0
#define LIF_TYPE_PSEC 0   // 'point source of electric current'
#define LIF_TYPE_PSMC 1   // 'point source of magnetic current'
struct PointSourceData: public IncFieldData
 { 
   double X0[3];         /* location */
   cdouble S[3];         /* strength */
   int Type;             /* LIF_TYPE_PSEC / PSMC for electric/magnetic */ 

   // constructor 
   PointSourceData(double X0[3], cdouble S[3], int Type);
   PointSourceData(double X0[3], cdouble S[3]); // defaults to Type=LIF_TYPE_PSEC

 };

/**********************************************************************/
/* focused gaussian beam  **********8**********************************/
/**********************************************************************/
struct GaussianBeamData: public IncFieldData
 { 
   double X0[3];            /* beam center point */
   double KProp[3];         /* beam propagation vector */
   cdouble E0[3];           /* complex field-strength vector */
   double W0;               /* beam waist */

   // constructor 
   GaussianBeamData(double X0[3], double KProp[3], cdouble E0[3], double W0);

 };
#endif

/**********************************************************************/
/* magnetic 'frill' (annulus)    **************************************/
/**********************************************************************/
#if 0
struct MagneticFrillData: public IncFieldData
 {
   double X0[3];          /* center of frill */
   double Theta, Phi;     /* angles of frill axis */
   double RIn, ROut;      /* inner and outer radii */
 } MagneticFrillData;

/**********************************************************************/
/* magnetic solenoid                            ***********************/
/**********************************************************************/
struct MagneticSolenoidData: public IncFieldData
 { 
   double X0[3];          /* center */
   double Theta, Phi;     /* angles of solenoid axis */
   double Radius, L;      /* radius and length  */
 } MagneticSolenoidData;
#endif

#endif
