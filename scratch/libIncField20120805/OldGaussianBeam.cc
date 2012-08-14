// this is the old gaussian beam routine as originally implemented by 
// homer; the calculation follows L. W. Davis, PRA _19_ 1177 (1979), 
// with the exception that our time-dependence convention has the opposite  
// sign, i.e. we assume all quantities vary in time like e^{-i\omega t}
// instead of e^{+i\omega t}.
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
