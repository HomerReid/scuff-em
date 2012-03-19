/*
 * EHGaussianBeam.cc -- routine for computing the electric and magnetic
 *                      fields of a focused gaussian beam
 *
 * homer reid        -- 4/2011
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
#include <libVec.h>

#include "libIncField.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
void EHGaussianBeam(double X[3], void *parms, cdouble EH[6]) {
  const cdouble IU(0,1);

  /***************************************************************/
  /* unpack parameters from structure ****************************/
  /***************************************************************/
  GaussianBeamData *GBD=(GaussianBeamData *)parms;

  double Omega  = GBD->Frequency;

  zVec E0 = GBD->E0;     // complex field-strength vector
                         // containing information on the 
                         // field strength and polarization 
  dVec KProp = GBD->KProp;  // propagation vector
  double W0     = GBD->W0;     // beam waist parameter
  double EpsR   = GBD->Eps;    // permittivity of medium is assumed real
  double MuR    = GBD->Mu;
  double k =sqrt(EpsR*MuR)*Omega; // wavenumber of medium 
  double ZR=sqrt(MuR/EpsR);       // relative wave impedance of medium

  // we follow CJR Sheppard & S Saghafi, J Opt Soc Am A 16, 1381 and
  // approximate the Gaussian beam as the field of the sum of an 
  // x-polarized electric and y-polarized magnetic dipole, located at
  // the complex point X = (0,0,i z0)
  // this has the advantage that the field can be analytically calculated
  // and exactly solves Maxwell's equations everywhere in space
  double z0 = k*W0*W0/2;
  dVec Xrel = dVec(X) - GBD->X0;

  // the field has NO cylindrical symmetry! this means that for
  // complex polarization vectors, we have to do a separate calculation
  // for the real and for the complex part

  // first, we do everything that is not direction dependent, i.e.
  // where we only need the z-coordinate and the radial distance rho
  dVec zHat = KProp; zHat.normalize();
  double z, rho;
  // this is from libVec.h
  GetLocalCylinderCoordinates(Xrel, zHat, rho, z);

  cdouble zc = z - IU*z0;
  cdouble Rsq  = rho*rho + zc*zc, R = sqrt(Rsq), kR = k*R, kRsq = kR*kR, kR3 = kRsq*kR;
  cdouble f,g,fmgbRsq;
  // we have to be careful: R can go to zero, leading to numerical problems
  if (std::abs(kR)>1.e-4) {
    cdouble coskR = cos(kR), sinkR = sin(kR);
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

  // now calculate the actual coordinates for having either E0.real() or E0.imag() as the x axis
  // local coordinate system given by ^z = kProp, ^x ~ Re(E0) or Im(E0), ^y ~ ^z x ^x  
  zVec E, H;

  double rnorm = norm(E0.real());
  if (rnorm>1e-13) {
    // calculate fields in local coordinate system
    dVec xHat = E0.real() / rnorm;
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
    E += rnorm * (Ex * zVec(xHat) + Ey * zVec(yHat) + Ez * zVec(zHat));
    H += rnorm * (Hx * zVec(xHat) + Hy * zVec(yHat) + Hz * zVec(zHat));
  } 
  
  double inorm = norm(E0.imag());
  if (inorm>1e-13) {
    // calculate fields in local coordinate system
    dVec xHat = E0.imag() / inorm;
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
  double kz0 = k*z0;
  double Eorig = 3./(2*kz0*kz0*kz0) * (exp(kz0)*kz0*(kz0-1) + sinh(kz0));
  
  // now scale the fields to have E(0,0,0) = E0
  E /= Eorig;
  EH[0] = E[0]; EH[1] = E[1]; EH[2] = E[2];
  H /= (Eorig*ZVAC*ZR);
  EH[3] = H[0]; EH[4] = H[1]; EH[5] = H[2];
}


void EHGaussianBeam_old(double R[3], void *parms, cdouble EH[6])
{
  /***************************************************************/
  /* unpack parameters from structure ****************************/
  /***************************************************************/
  GaussianBeamData *GBD=(GaussianBeamData *)parms;

  double Omega  = GBD->Frequency;

  cdouble *E0   = GBD->E0;     // complex field-strength vector
                               // containing information on the 
                               // field strength and polarization 

  double *KProp = GBD->KProp;  // propagation vector

  double W0     = GBD->W0;     // beam waist parameter

  double EpsR   = GBD->Eps;    // permittivity of medium is assumed real
  double MuR    = GBD->Mu;

  double k=sqrt(EpsR*MuR)*Omega; // wavenumber of medium 
  double ZR=sqrt(MuR/EpsR);      // relative wave impedance of medium

  /***************************************************************/
  /* extract (Theta,Phi) angles from user's specified propagation*/
  /* vector.                                                     */
  /* note: (Theta, Phi) are the angular coordinates of the point */
  /* on a unit sphere, centered at the beam center point, through*/
  /* which the beam propagation axis passes.                     */
  /***************************************************************/
  double Theta, Phi;
  Theta=atan2( sqrt(KProp[0]*KProp[0] + KProp[1]*KProp[1]), KProp[2] );
  Phi=atan2( KProp[1], KProp[0] );

  /***************************************************************/
  /* convert coordinates of observation point into the beam      */
  /* coordinate system (defined as a coordinate system with      */
  /* origin at the beam center and in which the beam travels in  */
  /* the positive z direction)                                   */
  /***************************************************************/
  double RPrime[3], RTemp[3];
  double M[3][3];
  double CT, ST, CP, SP;
  int Mu, Nu;
  RPrime[0]=R[0] - GBD->X0[0];
  RPrime[1]=R[1] - GBD->X0[1];
  RPrime[2]=R[2] - GBD->X0[2];
  CT=cos(Theta);
  ST=sin(Theta);
  CP=cos(Phi);
  SP=sin(Phi);
  M[0][0]=CT*CP;     M[0][1]=CT*SP;    M[0][2]=-ST;
  M[1][0]=-SP;       M[1][1]=CP;       M[1][2]=0.0;
  M[2][0]=ST*CP;     M[2][1]=ST*SP;    M[2][2]=CT;
  memset(RTemp,0,3*sizeof(double));
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    RTemp[Mu] += M[Mu][Nu] * RPrime[Nu];
  memcpy(RPrime,RTemp,3*sizeof(double));
  
  /***************************************************************/
  /* evaluate the field of the beam in the beam coordinate system*/
  /* assuming it is linearly polarized with E-field in the       */
  /* positive X direction and has unit field-strength prefactor  */
  /***************************************************************/
  double x=RPrime[0], y=RPrime[1], z=RPrime[2];
  double s, s2, s3, Xi, Xi2, Eta, Zeta, Rho2, Rho4;  
  cdouble Q, Q2, Q3, Q4, ExpFac, PreFac;
  cdouble Psi0, Psi2;
  cdouble dPsi0dXi, dPsi0dEta, dPsi0dZeta;
  cdouble d2Psi0dXi2, d2Psi0dXidEta, d2Psi0dXidZeta;
  cdouble dPsi2dXi, dPsi2dEta, dPsi2dZeta; 

  const cdouble IU(0,1);

  s=1/(k*W0);
  s2=s*s;
  s3=s2*s;
  Xi=x/W0;
  Eta=y/W0;
  Zeta=s*z/W0;
  Rho2=Xi*Xi + Eta*Eta;
  Rho4=Rho2*Rho2;
  Q=1.0/(2.0*Zeta - IU);
  Q2=Q*Q;
  Q3=Q2*Q;
  Q4=Q2*Q2;

  /*--------------------------------------------------------------*/
  /*- compute Psi0, Psi2, and all needed derivatives -------------*/
  /*--------------------------------------------------------------*/
  ExpFac=exp(IU*Q*Rho2);
  Psi0=IU*Q*ExpFac;
  Psi2=-1.0*(2.0*IU*Q + IU*Q3*Rho4)*Psi0;

  dPsi0dXi=2.0*IU*Q*Xi*Psi0;
  dPsi0dEta=2.0*IU*Q*Eta*Psi0;
  dPsi0dZeta= (-2.0*Q - 2.0*IU*Q2*Rho2)*Psi0;
  d2Psi0dXi2=(2.0*IU*Q - 4.0*Q2*Xi*Xi)*Psi0;
  d2Psi0dXidEta=-4.0*Q2*Xi*Eta*Psi0;
  d2Psi0dXidZeta=-1.0*4.0*IU*Q*Q*Xi*Psi0 + 2.0*IU*Q*Xi*dPsi0dZeta;

  dPsi2dXi=-1.0*4.0*IU*Q3*Rho2*Xi*Psi0 + 2.0*IU*Q*Xi*Psi2;
  dPsi2dEta=-1.0*4.0*IU*Q3*Rho2*Eta*Psi0 + 2.0*IU*Q*Eta*Psi2;
  dPsi2dZeta =  (4.0*IU*Q2 + 6.0*IU*Q4*Rho4)*Psi0
               -(2.0*IU*Q + IU*Q3*Rho4)*dPsi0dZeta;

  /*--------------------------------------------------------------*/
  /*- assemble contributions to field components                 -*/
  /*--------------------------------------------------------------*/
  PreFac=exp(IU*k*z);
  EH[0]=PreFac*(Psi0 + s2*(Psi2 + d2Psi0dXi2));
  EH[1]=1.0*PreFac*s2*d2Psi0dXidEta;
  EH[2]=1.0*PreFac*( IU*s*dPsi0dXi +
                      s3*(IU*dPsi2dXi + d2Psi0dXidZeta));
   
  PreFac/=(ZVAC*ZR);

  EH[3]=0.0;
  EH[4]=PreFac*(Psi0 + s2*(Psi2 - IU*dPsi0dZeta) ); 
  EH[5]=1.0*PreFac*(IU*s*dPsi0dEta + IU*s3*dPsi2dEta);

  /***************************************************************/
  /* the above computation assumed that the field was linearly   */
  /* polarized with E-field pointing in the X direction; we now  */
  /* form linear combinations of the field components as         */
  /* appropriate given the complex polarization vector specified */
  /* by the user. this is equations 4a, 4b of the memo.          */
  /* note: E0BS[0], E0BS[1] are what we called E_{0x}, E_{0y} in */
  /* the memo; these are the x and y components of the field-    */
  /* strength vector E0 in the beam coordinate system, which we  */
  /* compute by rotating the user's vector E0 into the beam      */
  /* system using the M matrix computed above.                   */
  /* if there is a non-negligible Z component of E0BS, then      */
  /* the user screwed up by specifying nonorthogonal             */
  /* polarization and propagation vectors.                       */
  /***************************************************************/
  cdouble E0BS[3]; // "E0 in the beam system"
  memset(E0BS,0,3*sizeof(cdouble));
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    E0BS[Mu] += M[Mu][Nu] * E0[Nu];
  if ( abs(E0BS[2]) > 1.0e-8 ) 
   fprintf(stderr,"** warning: nonorthogonal propagation and polarization vectors in EHGaussianBeam\n");

  cdouble EHPrime[6];
  EHPrime[0] = E0BS[0]*EH[0] - E0BS[1]*EH[1]; 
  EHPrime[1] = E0BS[0]*EH[1] + E0BS[1]*EH[0];
  EHPrime[2] = E0BS[0]*EH[2] + E0BS[1]*EH[2]; 

  EHPrime[3] =           0.0 - E0BS[1]*EH[4];
  EHPrime[4] = E0BS[0]*EH[4] + 0.0;
  EHPrime[5] = E0BS[0]*EH[5] + E0BS[1]*EH[5]; 

  /***************************************************************/
  /* finally, rotate the field components in the beam coordinate */
  /* system into the user's coordinate system. (note that this   */
  /* involves the transpose of M).                               */
  /***************************************************************/
  memset(EH,0,6*sizeof(cdouble));
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    { EH[Mu]   += M[Nu][Mu] * EHPrime[Nu];
      EH[Mu+3] += M[Nu][Mu] * EHPrime[Nu+3];
    };

} 
