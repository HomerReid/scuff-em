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
 * AnalyticalDGFs.cc -- calculation of dyadic Green's functions 
 *                   -- in some analytically tractable situations
 *
 * Homer Reid 7/2015
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libIncField.h>
#include <libhmat.h>
#include <libTriInt.h>
#include <libSpherical.h>
#include <libscuff.h>

#define II cdouble (0.0,1.0)

using namespace scuff;

/***************************************************************/
/* data structure containing everything needed to evaluate the */
/* integrand of the half-space DGF integral                    */
/***************************************************************/
typedef struct HalfSpaceData
 {
   HMatrix *XMatrix;
   cdouble Omega;     // angular frequency
   cdouble Epsilon;   // half-space permittivity
   cdouble Mu;        // half-space permeability

   double *kBloch;    // point in Brillouin zone
   bool Accumulate; 
   
   double Workspace[12];
   double qrOffset;
   bool uqTransform;
   bool Polar;

   int nCalls;

 } HalfSpaceData;

/***************************************************************/
/* Integrand[ 18*nx + 0*9 + 3*Mu + Nu ] = G^{E}_{Mu,Nu}        */
/* Integrand[ 18*nx + 1*9 + 3*Mu + Nu ] = G^{M}_{Mu,Nu}        */
/***************************************************************/
void HalfSpaceDGFIntegrand(const double *q, HalfSpaceData *Data,
                           cdouble *Integrand)
{
  cdouble EpsRel    = Data->Epsilon;
  cdouble MuRel     = Data->Mu;
  cdouble k0        = Data->Omega;
  HMatrix *XMatrix  = Data->XMatrix;
  bool Polar        = Data->Polar;

  Data->nCalls++;

  int IDim = 18*XMatrix->NR;
  if (Data->Accumulate == false)
   memset(Integrand, 0, IDim*sizeof(cdouble));

  // Polar = true --> we have already integrated out 
  //                  q_Theta to yield Bessel functions,
  //                  and what we are evaluating here 
  //                  is just the integrand of the 1-dimensional
  //                  q_r integral 
  //
  // Polar = false--> we are evaluating the 2-dimensional
  //                  (qx,qy) integral
  //
  double q2, qMag;
  cdouble One, Cos, Sin, Cos2, Sin2, CosSin;
  if (Polar)
   { q2        = q[0]*q[0];
     qMag      = q[0];
   }
  else
   { q2       = q[0]*q[0] + q[1]*q[1];
     qMag     = sqrt(q2);
     double qxHat = (qMag==0.0) ? 1.0 : q[0] / qMag;
     double qyHat = (qMag==0.0) ? 0.0 : q[1] / qMag;
     One      = 1.0;
     Cos      = qxHat;
     Sin      = qyHat;
     Cos2     = qxHat*qxHat;
     CosSin   = qxHat*qyHat;
     Sin2     = qyHat*qyHat;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble k02 = k0*k0;
  cdouble qz2 = q2 - k0*k0;
  if (qz2==0.0)
   return;
  cdouble qz = sqrt(k0*k0 - q2);
  cdouble qzPrime = sqrt(EpsRel*MuRel*k0*k0 - q2);
  if ( imag(qz)<0.0 )
   qz*=-1.0;
  if ( imag(qzPrime)<0.0 )
   qzPrime*=-1.0;

  cdouble rTE, rTM;
  if (EpsRel==0.0 && MuRel==0.0) // PEC case
   { rTE = -1.0;
     rTM = +1.0;
   }
  else
   { rTE = (MuRel*qz - qzPrime) / (MuRel*qz + qzPrime);
     rTM = (EpsRel*qz - qzPrime) / (EpsRel*qz + qzPrime);
   };

  cdouble MTE[3][3], MTM[3][3];
  MTE[0][2] = MTE[1][2] = MTE[2][0] = MTE[2][1] = MTE[2][2] = 0.0;

  bool TwoPointDGF = (XMatrix->NC>=6);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nx=0; nx<XMatrix->NR; nx++)
   { 
     double XDest[3], XSourceBuffer[3];
     double *XSource = (TwoPointDGF) ? XSourceBuffer : XDest;
     XMatrix->GetEntriesD(nx,"0:2",XDest);
     if (TwoPointDGF)
      XMatrix->GetEntriesD(nx,"3:5",XSource);

     if ( abs(imag(qz*(XSource[2] + XDest[2])) > 40.0 ) )
      continue;

     double R[3];
     VecSub(XDest, XSource, R);
     double Rho=sqrt( R[0]*R[0] + R[1]*R[1] );
     double xHat = (Rho==0.0) ? 1.0 : R[0] / Rho;
     double yHat = (Rho==0.0) ? 0.0 : R[1] / Rho;

     double qDotRho=0.0;
     if (Polar)
      { cdouble J[3];
        double qRho = qMag*Rho;
        double TPQ = 2.0*M_PI*qMag;
        AmosBessel('J', qRho, 0.0, 3, false, J, Data->Workspace);
        cdouble J1oqRho = (qRho==0.0 ? 0.0 : J[1]/qRho);
        cdouble Bracket = (J[0] - 2.0*J1oqRho - J[2]);
        One     = TPQ*J[0];
        Cos     = II*TPQ*J[1]*xHat;
        Sin     = II*TPQ*J[1]*yHat;
        Cos2    = TPQ*(0.5*Bracket*xHat*xHat + J1oqRho);
        CosSin  = TPQ*0.5*Bracket*xHat*yHat;
        Sin2    = TPQ*(0.5*Bracket*yHat*yHat + J1oqRho);
      }
     else
      qDotRho = q[0]*R[0] + q[1]*R[1];
     
     MTE[0][0] = Sin2;
     MTE[1][1] = Cos2;
     MTE[0][1] = MTE[1][0] = -1.0*CosSin;

     MTM[0][0] = qz2*Cos2 / k02;
     MTM[1][1] = qz2*Sin2 / k02;
     MTM[2][2] = q2*One  / k02;
     MTM[0][1] = MTM[1][0] = qz2*CosSin / k02;
     MTM[2][0] =      qMag*qz*Cos / k02;
     MTM[0][2] = -1.0*MTM[2][0];
     MTM[2][1] =      qMag*qz*Sin / k02;
     MTM[1][2] = -1.0*MTM[2][1];

     cdouble ExpArg = II*( qDotRho + qz*(XSource[2]+XDest[2]) );
     cdouble Factor = II*exp(ExpArg) / (8.0*M_PI*M_PI*qz);
  
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       { Integrand[18*nx + 0*9 + 3*Mu + Nu] 
          += Factor * (rTE*MTE[Mu][Nu] + rTM*MTM[Mu][Nu]);
         Integrand[18*nx + 1*9 + 3*Mu + Nu] 
          += Factor * (rTM*MTE[Mu][Nu] + rTE*MTM[Mu][Nu]);
       };
#if 0
     if (LogFile)
     { fprintf(LogFile,"%e %s %s ",q,CD2S(qz),CD2S(qzPrime));
       fprintf(LogFile,"%e %e %e ",R[0],R[1],XA[2]+XB[2]);
       fprintf(LogFile,"%s ",CD2S(ExpFac));
       for(int n=0; n<9; n++)
        fprintf(LogFile,"%s ",CD2S(zf[18*nx + n]));
       fprintf(LogFile,"\n");
     };
#endif

   };
}

/***************************************************************/
/* summand function passed to GetLatticeSum() to evaluate the  */
/* reciprocal-lattice sum for the BZ integrand at a single     */
/* kBloch point                                                */
/***************************************************************/
void HalfSpaceDGFSummand(double *Gamma, void *UserData, double *Sum)
{
  HalfSpaceData *Data = (HalfSpaceData *)UserData;
  double *kBloch    = Data->kBloch;
  
  double q[2];
  q[0] = kBloch[0] + Gamma[0];
  q[1] = kBloch[1] + Gamma[1];
  
  HalfSpaceDGFIntegrand(q, (HalfSpaceData *)UserData, (cdouble *)Sum);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetHalfSpaceDGFs_BZ(HMatrix *XMatrix,
                         cdouble Omega, double kBloch[2],
                         HMatrix *RLBasis, double BZVolume, 
                         MatProp *MP,
                         double RelTolSum, double AbsTolSum, 
                         int MaxCells,
                         HMatrix *GMatrix)
{ 
  int IDim = 18*XMatrix->NR;
  static int IDimSave=0;
  static cdouble *Sum=0;
  if (IDimSave<IDim)
   { IDimSave=IDim;
     Sum = (cdouble *)reallocEC(Sum, IDim*sizeof(cdouble));
   };

  cdouble Epsilon=0.0, Mu=0.0;
  if (MP && MP->IsPEC() == false )
   MP->GetEpsMu(Omega, &Epsilon, &Mu);

  HalfSpaceData MyData, *Data=&MyData;
  Data->XMatrix     = XMatrix;
  Data->Omega       = Omega;
  Data->Epsilon     = Epsilon;
  Data->Mu          = Mu;
  Data->kBloch      = kBloch;
  Data->Polar       = false;
  Data->uqTransform = false;
  Data->nCalls      = 0;
  Data->Accumulate  = true;
 
  Log("Evaluating BZ sum for DGF integrand at kBloch=(%e,%e)...", kBloch[0], kBloch[1]);
  GetLatticeSum(HalfSpaceDGFSummand, (void *)Data, 2*IDim, RLBasis,
                (double *)Sum, AbsTolSum, RelTolSum, MaxCells);
  Log("...%i lattice cells summed",Data->nCalls);

  int NX = XMatrix->NR;
  for(int nx=0; nx<NX; nx++)
   for(int ng=0; ng<18; ng++)
    GMatrix->SetEntry(nx, ng, BZVolume*Sum[18*nx + ng]);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int HalfSpaceDGFIntegrand_Polar(unsigned ndim, const double *u, void *UserData, unsigned fdim, double *fval)
{ 
  (void) ndim;
  (void) fdim;

  HalfSpaceData *Data = (HalfSpaceData *)UserData;
  bool uqTransform = Data->uqTransform;
  int qrOffset     = Data->qrOffset;
  int IDim         = 18*(Data->XMatrix->NR);

  cdouble *Integrand = (cdouble *)fval;
  memset(Integrand, 0, IDim*sizeof(cdouble));

  double qr, Jacobian;
  if (uqTransform)
   { double Denom = 1.0 - fabs(u[0]);
     if (Denom==0.0)
      return 0;
     qr       = qrOffset + u[0] / Denom;
     Jacobian = 1.0/(Denom*Denom);
   }
  else
   { qr=u[0];
     Jacobian=1.0;
   };

  HalfSpaceDGFIntegrand(&qr, Data, Integrand);
  VecScale(Integrand, Jacobian, IDim);

  return 0;

}

void GetHalfSpaceDGFs_Polar(HMatrix *XMatrix, cdouble Omega,
                            MatProp *MP,
                            double RelTol, double AbsTol,
                            int MaxEvals, HMatrix *GMatrix)
{ 
 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NX = XMatrix->NR;
  int IDim = 18*NX;
  static int IDimSave=0;
  static cdouble *Integral1=0, *Integral2=0, *Error=0;
  if (IDimSave!=IDim)
   { IDimSave = IDim;
     Integral1 = (cdouble *)reallocEC(Integral1, IDim * sizeof(cdouble));
     Integral2 = (cdouble *)reallocEC(Integral2, IDim * sizeof(cdouble));
     Error     = (cdouble *)reallocEC(Error, IDim * sizeof(cdouble));
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble Epsilon=0.0, Mu=0.0;
  if ( MP && MP->IsPEC() == false )
   MP->GetEpsMu(Omega, &Epsilon, &Mu);

  Log("Evaluating qr integral for DGFs at %i points...",NX);

  double Lower, Upper;
  HalfSpaceData MyData, *Data=&MyData;
  Data->XMatrix    = XMatrix;
  Data->Omega      = Omega; 
  Data->Epsilon    = Epsilon;
  Data->Mu         = Mu;     
  Data->Polar      = true;
  Data->Accumulate = false;

  Lower=0.0;
  Upper=abs(Omega);
  Data->qrOffset    = 0.0;
  Data->uqTransform = false;
  Data->nCalls      = 0;
  hcubature(2*IDim, HalfSpaceDGFIntegrand_Polar, (void *)Data, 1,
            &Lower, &Upper, MaxEvals, AbsTol, RelTol, 
            ERROR_INDIVIDUAL, (double *)Integral1, (double *)Error);
  Log(" small-q integral: %i calls",Data->nCalls);

  Lower=0.0;
  Upper=1.0;
  Data->qrOffset    = abs(Omega);
  Data->uqTransform = true;
  Data->nCalls      = 0;
  hcubature(2*IDim, HalfSpaceDGFIntegrand_Polar, (void *)Data, 1, 
            &Lower, &Upper, MaxEvals, AbsTol, RelTol, 
            ERROR_INDIVIDUAL, (double *)Integral2, (double *)Error);
  Log(" large-q integral: %i calls",Data->nCalls);

 
  for(int nx=0; nx<NX; nx++)
   for(int ng=0; ng<18; ng++)
    GMatrix->SetEntry(nx, ng, Integral1[18*nx + ng] + Integral2[18*nx + ng]);
  
}

/***************************************************************/
/* Ground-plane DGFs computed by the image-source method       */
/***************************************************************/
void GetGroundPlaneDGFs(HMatrix *XMatrix,
                        cdouble Omega, double *kBloch, HMatrix *LBasis, 
                        HMatrix *GMatrix)
{
  /***************************************************************/
  /* create an IncidentField structure describing the field of a */
  /* point source or a periodic array of point sources at the    */
  /* location of the image of XSource                            */
  /***************************************************************/
  PointSource PS;
  PS.SetFrequency(Omega);
  if (LBasis)
   { PS.SetLattice(LBasis);
     PS.SetkBloch(kBloch);
   };

  bool TwoPointDGF = (XMatrix->NC >= 6);

  int NX = XMatrix->NR;
  for(int nx=0; nx<NX; nx++)
   { 
     double XDest[3], XSource[3];
     XMatrix->GetEntriesD(nx, "0:2", XDest);
     if (TwoPointDGF)
      XMatrix->GetEntriesD(nx, "3:5", XSource);
     else
      XMatrix->GetEntriesD(nx, "0:2", XSource);
     XSource[2]*=-1.0;
     PS.SetX0(XSource);

     /***************************************************************/
     /* columns of DGF are fields of image source at eval point    **/
     /***************************************************************/
     for(int Nu=0; Nu<3; Nu++)
      { 
        cdouble P[3], EH[6];

        memset(P, 0, 3*sizeof(cdouble));
        P[Nu] = (Nu==2) ? 1.0 : -1.0;
        PS.SetP(P);
        PS.SetType(LIF_ELECTRIC_DIPOLE);
        PS.GetFields(XDest, EH);
        for(int Mu=0; Mu<3; Mu++)
         GMatrix->SetEntry(nx, 0*9 + 3*Mu + Nu, EH[Mu] / (Omega*Omega));

        P[Nu] *= -1.0;
        PS.SetP(P);
        PS.SetType(LIF_MAGNETIC_DIPOLE);
        PS.GetFields(XDest, EH);
        for(int Mu=0; Mu<3; Mu++)
         GMatrix->SetEntry(nx, 1*9 + 3*Mu + Nu, EH[3+Mu] / (Omega*Omega));
      };
   };

}

/***************************************************************/
/* Vacuum DGFs computed by the plane-wave decomposition        */
/***************************************************************/
#if 0
int GetVacuumDGFs(double *XSource, double *XDest,
                  cdouble Omega, double kBloch[2],
                  HMatrix *RLBasis, double BZVolume,
                  double RelTol, double AbsTol, int MaxCells,
                  cdouble GE[3][3], cdouble GH[3][3])
{ 
 
  HalfSpaceData MyHalfSpaceData, *Data=&MyHalfSpaceData;
  Data->XSource = XSource;
  Data->XDest   = XDest;
  Data->Omega   = Omega;
  Data->kBloch  = kBloch;
  Data->Epsilon = 0.0;

  cdouble Sum[18];
  int NumCells=GetLatticeSum(DGFSummand, (void *)Data, 36, RBasis,
                            (double *)Sum, AbsTol, RelTol, MaxCells);

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { GE[Mu][Nu] = BZVolume*Sum[0 + 3*Mu + Nu]/(4.0*M_PI*M_PI);
      GH[Mu][Nu] = BZVolume*Sum[9 + 3*Mu + Nu]/(4.0*M_PI*M_PI);
    };

  return NumCells;
}
#endif

/***************************************************************/
/* Compute half-space DGFs by direct methods, i.e. bypassing   */
/* the Brillouin-zone sum.                                     */
/***************************************************************/
void ProcessHalfSpaceDGFs(HVector *OmegaPoints,
                          char **EPFiles, int nEPFiles,
                          char *HalfSpace,
                          double RelTol, double AbsTol, int MaxEvals)
{
  MatProp::SetLengthUnit(1.0e-6);
  MatProp *MP = HalfSpace ? new MatProp(HalfSpace) : 0;
  char *Label = HalfSpace ? HalfSpace : const_cast<char *>("GroundPlane");

  for(int nep=0; nep<nEPFiles; nep++)
   { 
     HMatrix *XMatrix = new HMatrix(EPFiles[nep]);
     int NX = XMatrix->NR;
     bool TwoPointDGF = (XMatrix->NC==6);
     HMatrix *GMatrix = new HMatrix(NX, 18, LHM_COMPLEX);
     
     FILE *f=vfopen("%s.%s.DGFs","w",GetFileBase(EPFiles[nep]),Label);
     fprintf(f,"# scuff-ldos ran on %s (%s)\n",GetHostName(),GetTimeString());
     fprintf(f,"# data file columns: \n");
     fprintf(f,"#  1, 2, 3 {x,y,z} (source point)\n");
     int nc=4;
     if (TwoPointDGF)
      { fprintf(f,"#  %i, %i, %i {x,y,z} (evaluation point)\n", nc, nc+1, nc+2);
        nc+=3;
      };
     fprintf(f,"# %i, %i    re,im Omega\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    electric, magnetic LDOS\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GExx\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GExy\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GExz\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GEyx\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GEyy\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GEyz\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GEzx\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GEzy\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GEzz\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GMxx\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GMxy\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GMxz\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GMyx\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GMyy\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GMyz\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GMzx\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GMzy\n",nc,nc+1); nc+=2;
     fprintf(f,"# %i, %i    re,im GMzz\n",nc,nc+1); nc+=2;

     for(int nOmega=0; nOmega<OmegaPoints->N; nOmega++)
      { cdouble Omega = OmegaPoints->GetEntry(nOmega);
        if (HalfSpace)
         GetHalfSpaceDGFs_Polar(XMatrix, Omega, MP,
                                RelTol, AbsTol, MaxEvals, GMatrix);
        else
         GetGroundPlaneDGFs(XMatrix, Omega, 0, 0, GMatrix);

        double PreFac = abs(Omega)/M_PI;
        for(int nx=0; nx<NX; nx++)
         { for(int Mu=0; Mu<(TwoPointDGF ? 6 : 3); Mu++)
            fprintf(f,"%e ",XMatrix->GetEntryD(nx,Mu));
           fprintf(f,"%e %e ",real(Omega),imag(Omega));
           
           fprintf(f,"%e  ",PreFac*imag( GMatrix->GetEntry(nx, 0 + 0)
                                        +GMatrix->GetEntry(nx, 0 + 4)
                                        +GMatrix->GetEntry(nx, 0 + 8)
                                       ));
           fprintf(f,"%e  ",PreFac*imag( GMatrix->GetEntry(nx, 9 + 0)
                                        +GMatrix->GetEntry(nx, 9 + 4)
                                        +GMatrix->GetEntry(nx, 9 + 8)
                                       ));
           for(int m=0; m<18; m++)
            fprintf(f,"%e %e ",real(GMatrix->GetEntry(nx, m)),
                               imag(GMatrix->GetEntry(nx, m)));
           fprintf(f,"\n");
         };
        fprintf(f,"\n\n");
        fflush(f);
      };

     delete XMatrix; 
     delete GMatrix; 
     fclose(f);

   };

}

/***************************************************************/
/* data structure containing everything needed to evaluate the */
/* integrand of the cylinder DGF integral                      */
/***************************************************************/
typedef struct CylinderData
 {
   HMatrix *XMatrix;
   cdouble Omega;          // angular frequency
   cdouble EpsOut, MuOut;  // outer material properties
   cdouble EpsIn,  MuIn;   // outer material properties

   double kx;              // (one-dimensional) bloch vector
   bool Accumulate; 
   
   bool hTransform;
   double hOffset;

   int nCalls;

 } CylinderData;

/***************************************************************/
/* Integrand[ 2*nx + 0 ] = G^E_{transverse}                    */
/* Integrand[ 2*nx + 1 ] = G^E_{axial direction}               */
/***************************************************************/
void CylinderDGFIntegrand(const double h, CylinderData *Data,
                          cdouble *Integrand)
{
  HMatrix *XMatrix = Data->XMatrix;
  cdouble Omega    = Data->Omega;
  cdouble EpsOut   = Data->EpsOut;
  cdouble EpsIn    = Data->EpsIn;
  cdouble MuOut    = Data->MuOut; 
  cdouble MuIn     = Data->MuIn;   
  double kx        = Data->kx;
  bool Accumulate  = Data->Accumulate;

  Data->nCalls++;

  int IDim = 2*XMatrix->NR;
  if (Data->Accumulate == false)
   memset(Integrand, 0, IDim*sizeof(cdouble));



}
