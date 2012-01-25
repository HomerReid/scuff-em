/*
 * libIncField.h -- routines for evaluating the electric and magnetic
 *               -- fields of plane waves and of point sources
 *
 *               -- originally designed for use with libRWG scattering 
 *               -- codes but now split out into a separate standalone 
 *               -- library
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
/* IncFieldData is a general base class from which various specific   */
/* classes are derived.                                               */
/**********************************************************************/
typedef struct IncFieldData;
 { 
   cdouble Omega;
   cdouble Eps;
   double Mu;
   struct IncFieldData *Next;

 } IncFieldData;

void 

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/


/**********************************************************************/
/* plane wave *********************************************************/
/**********************************************************************/
typedef struct PlaneWaveData
 { double Frequency;      /* real or imaginary frequency */
   int RealFreq;          /* =1 for real frequency, 0 for imaginary frequency */
   double Eps, Mu;        /* RELATIVE constants of medium in which wave travels */
   cdouble E0[3];         /* E-field polarization vector */
   double nHat[3];        /* unit vector in direction of propagation */
 } PlaneWaveData;

void EHPlaneWave(double X[3], void *opPWD, cdouble *EH);

/**********************************************************************/
/* point source *******************************************************/
/**********************************************************************/
#define LIF_TYPE_PSEC 0
#define LIF_TYPE_PSMC 1

typedef struct PointSourceData
 {
   double Frequency;      /* real or imaginary frequency */
   int RealFreq;          /* =1 for real frequency, 0 for imaginary frequency */
   double Eps, Mu;        /* properties of medium through which wave travels */
   double X0[3];          /* location of point source */
   double nHat[3];        /* direction of point source */
   int SourceType;        /* either LIF_TYPE_PSEC or LIF_TYPE_PSMC */
 } PointSourceData;

void EHPointSource(double X[3], void *opPSD, cdouble *EH);

/**********************************************************************/
/* complex point source  **********************************************/
/**********************************************************************/
typedef struct ComplexPointSourceData
 { double Frequency;      /* real or imaginary frequency */
   int RealFreq;          /* =1 for real frequency, 0 for imaginary frequency */
   double Eps, Mu;        /* properties of medium through which wave travels */
   double X0[3];          /* location of point source */
   cdouble P[3];          /* components of source */
   int SourceType;        /* either LIF_TYPE_PSEC or LIF_TYPE_PSMC */
 } ComplexPointSourceData;

void EHComplexPointSource(double X[3], void *opCPSD, cdouble *EH);

/**********************************************************************/
/* magnetic 'frill' (annulus)    **************************************/
/**********************************************************************/
typedef struct MagneticFrillData
 { double Frequency;      /* real or imaginary frequency */
   int RealFreq;          /* =1 for real frequency, 0 for imaginary frequency */
   double Eps, Mu;        /* properties of medium through which wave travels */
   double X0[3];          /* center of frill */
   double Theta, Phi;     /* angles of frill axis */
   double RIn, ROut;      /* inner and outer radii */
 } MagneticFrillData;

void EHMagneticFrill(double X[3], void *opMFD, cdouble *EH);
void DrawMagneticFrill(void *opMFD, const char *PPFile);

/**********************************************************************/
/* magnetic solenoid                            ***********************/
/**********************************************************************/
typedef struct MagneticSolenoidData
 { double Frequency;      /* real or imaginary frequency */
   int RealFreq;          /* =1 for real frequency, 0 for imaginary frequency */
   double Eps, Mu;        /* properties of medium through which wave travels */
   double X0[3];          /* center */
   double Theta, Phi;     /* angles of solenoid axis */
   double Radius, L;      /* radius and length  */
 } MagneticSolenoidData;

void EHMagneticSolenoid(double X[3], void *opMDD, cdouble *EH);
void DrawMagneticSolenoid(void *opMDD, const char *PPFile);

/**********************************************************************/
/* focused gaussian beam  **********8**********************************/
/**********************************************************************/
typedef struct GaussianBeamData
 { double Frequency;        /* real or imaginary frequency */
   int RealFreq;            /* =1 for real frequency, 0 for imaginary frequency */
   double Eps, Mu;          /* properties of medium through which wave travels */
   double W0;               /* beam waist */
   double X0[3];            /* beam center point */
   double KProp[3];         /* beam propagation vector */
   cdouble E0[3];           /* complex field-strength vector */
 } GaussianBeamData;

void EHGaussianBeam(double X[3], void *opGBD, cdouble *EH);

#endif
