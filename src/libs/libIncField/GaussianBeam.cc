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
 * GaussianBeam.cc -- routine for computing the electric and magnetic
 *                    fields of a focused gaussian beam
 *
 * homer reid        -- 4/2011
 *
 * new approach:
 * johannes feist    -- 11/2011
 *
 * note: this calculation follows CJR Sheppard & S Saghafi, J Opt Soc Am A 16, 1381,
 * approximating a Gaussian beam as the field of the sum of an 
 * x-polarized electric and y-polarized magnetic dipole, but located at
 * the _complex_ source point X = (0,0,i z0).
 * This produces a field that looks like a Gaussian beam in paraxial approximation
 * and exactly solves the field-free Maxwell equations everywhere in (real) space.
 * 
 * note that the core of the calculation is carried out in a coordinate
 * system in which the beam propagates in the +z direction and is
 * polarized such that the E-field is (mostly) in the +x direction. after
 * carrying out the computation described in the memo to obtain the field
 * components in this coordinate system, we then rotate the field components
 * as appropriate for the user's specified propagation and polarization
 * vectors.
 */

#include <stdio.h>
#include <complex>
#include <libhrutil.h>

#include "libIncField.h"
#include "libVec.h"

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
GaussianBeam::GaussianBeam(const double pX0[3], const double pKProp[3], 
			   const cdouble pE0[3], double pW0, const char *Label)
{
  memcpy(X0, pX0, 3*sizeof(double));
  memcpy(KProp, pKProp, 3*sizeof(double));
  memcpy(E0, pE0, 3*sizeof(cdouble));
  W0=pW0;
  SetRegionLabel(Label);
}

GaussianBeam::~GaussianBeam()
{ 
  // no malloc'ed data to free
}

void GaussianBeam::SetX0(double pX0[3])       { memcpy(X0,    pX0,    3*sizeof(double)); }
void GaussianBeam::SetKProp(double pKProp[3]) { memcpy(KProp, pKProp, 3*sizeof(double)); }
void GaussianBeam::SetE0(cdouble pE0[3])      { memcpy(E0,    pE0,    3*sizeof(cdouble)); }
void GaussianBeam::SetW0(double pW0)          { W0 = pW0; }

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GaussianBeam::GetFields(const double X[3], cdouble EH[6])
{
  if ( imag(Eps) !=0.0 || imag(Mu) != 0.0 )
   ErrExit("%s:%i: gaussian beams not implemented for dispersive media");
  if ( imag(Omega) !=0.0 )
   ErrExit("%s:%i: gaussian beams not implemented for imaginary frequencies");
  if ( LBasis )
   ErrExit("%s:%i: gaussian beams not implemented for bloch-periodic geometries");

  const cdouble IU(0,1);
  double EpsR = real(Eps);
  double MuR  = real(Mu);
  double k    = sqrt(EpsR*MuR)*real(Omega); // wavenumber of medium 
  double ZR   = sqrt(MuR/EpsR);             // relative wave impedance of medium

  zVec zvE0     = E0;     // complex field-strength vector
                          // containing information on the 
                          // field strength and polarization 

  // we follow CJR Sheppard & S Saghafi, J Opt Soc Am A 16, 1381 and
  // approximate the Gaussian beam as the field of the sum of an   
  // x-polarized electric and y-polarized magnetic dipole, located at
  // the complex point X = (0,0,i z0)
  // this has the advantage that the field can be analytically calculated
  // and exactly solves Maxwell's equations everywhere in space
  double z0 = k*W0*W0/2;
  double kz0 = k*z0;
  dVec Xrel = dVec(X) - dVec(X0);

  // the field has NO cylindrical symmetry! this means that for
  // complex polarization vectors, we have to do a separate calculation
  // for the real and for the complex part

  // first, we do everything that is not direction dependent, i.e.
  // where we only need the z-coordinate and the radial distance rho
  dVec zHat = KProp; zHat.normalize();
  double z, rho;

  // this is from libVec.h
  GetLocalCylinderCoordinates(Xrel, zHat, rho, z);

  // HR 20130915 the cos, sin below can overflow if kR has large 
  // imaginary part, so in that case we use the 'rescaled' versions 
  // of f and g, defined as f,g divided by exp(kz0). x
  bool UseRescaledFG = false;

  cdouble zc = z - IU*z0;
  cdouble Rsq  = rho*rho + zc*zc, R = sqrt(Rsq), kR = k*R, kRsq = kR*kR, kR3 = kRsq*kR;
  cdouble f,g,fmgbRsq;
  // we have to be careful: R can go to zero, leading to numerical problems
  if (std::abs(kR)>1.e-4) 
   {
    cdouble coskR, sinkR;
    if ( fabs(imag(kR))>30.0 )
     { UseRescaledFG = true;
       cdouble ExpI     = exp( IU*real(kR) );
       cdouble ExpPlus  = exp( imag(kR) - kz0 );
       cdouble ExpMinus = exp( -(imag(kR) + kz0) );
       coskR = 0.5*( ExpI*ExpMinus + conj(ExpI)*ExpPlus);
       sinkR = -0.5*IU*( ExpI*ExpMinus - conj(ExpI)*ExpPlus);
     }
    else
     { coskR = cos(kR); 
       sinkR = sin(kR);
     };
    f   = -3.  *            (coskR/kRsq - sinkR/kR3);
    //g =  1.5 * (sinkR/kR + coskR/kRsq - sinkR/kR3)
    g   =  1.5 *  sinkR/kR - 0.5 * f;
    fmgbRsq = (f-g)/Rsq;
  } else {
    cdouble kR4 = kRsq*kRsq;
    // use a series expansion for small R
    // fourth order term is already at most 1e-16*3/280!
    f = kR4   /280. - kRsq/10. + 1.;
    g = kR4*3./280. - kRsq/5.  + 1.;
    // note: this is (f(kR)-g(kR))/R^2, not /kR^2 - so we get an additional k^2 term
    fmgbRsq = (kR4/5040. - kRsq/140. + 0.1) * (k*k);
  }
  cdouble i2fk = 0.5*IU*f*k;

  // now calculate the actual coordinates for having either zvE0.real() or zvE0.imag() as the x axis
  // local coordinate system given by ^z = kProp, ^x ~ Re(E0) or Im(E0), ^y ~ ^z x ^x  
  zVec E, H;

  double rnorm = norm(zvE0.real());
  if (rnorm>1e-13) {
    // calculate fields in local coordinate system
    dVec xHat = zvE0.real() / rnorm;
    dVec yHat = cross(zHat,xHat);
    double  x  = dot(xHat,Xrel);
    double  y  = dot(yHat,Xrel);
    
    cdouble Ex = g + fmgbRsq * x * x  + i2fk * zc;
    cdouble Ey =     fmgbRsq * x * y;
    cdouble Ez =     fmgbRsq * x * zc - i2fk * x;
    cdouble Hx = Ey;
    cdouble Hy = g + fmgbRsq * y * y  + i2fk * zc;
    cdouble Hz =     fmgbRsq * y * zc - i2fk * y;

    // go back to the laboratory frame
    E += cdouble(rnorm) * (Ex * zVec(xHat) + Ey * zVec(yHat) + Ez * zVec(zHat));
    H += cdouble(rnorm) * (Hx * zVec(xHat) + Hy * zVec(yHat) + Hz * zVec(zHat));
  } 
  
  double inorm = norm(zvE0.imag());
  if (inorm>1e-13) {
    // calculate fields in local coordinate system
    dVec xHat = zvE0.imag() / inorm;
    dVec yHat = cross(zHat,xHat);
    double  x  = dot(xHat,Xrel);
    double  y  = dot(yHat,Xrel);
    
    cdouble Ex = g + fmgbRsq * x * x  + i2fk * zc;
    cdouble Ey =     fmgbRsq * x * y;
    cdouble Ez =     fmgbRsq * x * zc - i2fk * x;
    cdouble Hx = Ey;
    cdouble Hy = g + fmgbRsq * y * y  + i2fk * zc;
    cdouble Hz =     fmgbRsq * y * zc - i2fk * y;
    
    // go back to the laboratory frame
    E += IU * inorm * (Ex * zVec(xHat) + Ey * zVec(yHat) + Ez * zVec(zHat));
    H += IU * inorm * (Hx * zVec(xHat) + Hy * zVec(yHat) + Hz * zVec(zHat));
  }

  // the field as calculated above is not normalized, so we get the field strength at the origin
  // (for E0 == 1)
  // this can be simplified very much by using that for x=y=z=0, R = sqrt((-i z0)**2) = i z0

  // 20130915 HR see comments above; note sinh(kz0)/exp(kz0) = 0.5(1-exp(-2*kz0)) 
  double Eorig; 
  if (UseRescaledFG)
   Eorig = 3./(2*kz0*kz0*kz0) * (kz0*(kz0-1) + 0.5*(1.0-exp(-2.0*kz0)) );
  else
   Eorig = 3./(2*kz0*kz0*kz0) * (exp(kz0)*kz0*(kz0-1) + sinh(kz0));
  
  // now scale the fields to have E(0,0,0) = E0
  E /= Eorig;
  EH[0] = E[0]; EH[1] = E[1]; EH[2] = E[2];
  H /= (Eorig*ZVAC*ZR);
  EH[3] = H[0]; EH[4] = H[1]; EH[5] = H[2];
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
double GaussianBeam::TotalBeamFlux() {
  // for analytical calculation see gaussian_beam_complexpointsource.nb
  double k  = real(sqrt(Eps*Mu)*Omega);   // wavenumber of medium
  double ZR = real(sqrt(Mu/Eps));         // relative wave impedance of medium
  double z0 = k*W0*W0/2, kz0 = k*z0;
  double Eorig = 3./(2*kz0*kz0*kz0) * (exp(kz0)*kz0*(kz0-1) + sinh(kz0));
  double flux  = 9*M_PI / (32*kz0*kz0*kz0*k*k);
  flux *= (-1 - 2*kz0 * (kz0-1) + (1+2*kz0*(2*kz0-1))*cosh(2*kz0) + 2*kz0*(2*kz0-1)*sinh(2*kz0));
  zVec zvE0     = E0;
  return flux * norm2(zvE0) / (Eorig*Eorig*ZVAC*ZR);
}
