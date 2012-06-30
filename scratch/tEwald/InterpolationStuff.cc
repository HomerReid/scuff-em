
/***************************************************************/
/* an interface to GBarVDEwald that has the right prototype to */
/* be passed to my Interp3D class constructor                  */
/***************************************************************/
typedef struct GBarVDEwaldData
 {
   cdouble Beta;
   double *K;
   double *L1;
   double *L2;
   double D;

 } GBarVDEwaldData;

void GBarVDEwaldWrapper3D(double X1, double X2, double X3, 
                          void *UserData, double *PhiVD)
{  
  GBarVDEwaldData *GD=(GBarVDEwaldData *)UserData;

  double R[3]={X1, X2, X3};
  cdouble GBarVD[8];
  GBarVDEwald(GD->Beta, GD->K, GD->L1, GD->L2, R, GBarVD);

  int nd;
  for(nd=0; nd<8; nd++)
   { PhiVD[0 + nd] = real(GBarVD[nd]);
     PhiVD[8 + nd] = imag(GBarVD[nd]);
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ReInitializeGBarInterp3D(cdouble Beta, double *K, double *L1, double *L2, 
                              int nThread, Interp3D *GBarInterp3D)
{ 
  GBarVDEwaldData MyData, *GD=&MyData;
  
  GD->Beta=Beta;
  GD->K=K;
  GD->L1=L1;
  GD->L2=L2;

  Log(" Initializing 3D interpolation table (%i grid points)",
        GBarInterp3D->N1 * GBarInterp3D->N2 * GBarInterp3D->N3);
  GBarInterp3D->ReInitialize(nThread, GBarVDEwaldWrapper3D, (void *)GD);

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GBarVDEwaldWrapper2D(double X1, double X2, void *UserData, double *PhiVD)
{  
  GBarVDEwaldData *GD=(GBarVDEwaldData *)UserData;

  cdouble GBarVD[8];
  double R[3]={X1, X2, GD->D};
  GBarVDEwald(GD->Beta, GD->K, GD->L1, GD->L2, R, GBarVD);

// f[0] = real(Phi)
// f[1] = imag(Phi)
// f[2] = real(dPhidZ)
// f[3] = imag(dPhidZ)

  PhiVD[0 ] = real(GBarVD[0]);  // Phi
  PhiVD[1 ] = real(GBarVD[1]);  // dPhidX
  PhiVD[2 ] = real(GBarVD[2]);  // dPhidY
  PhiVD[3 ] = real(GBarVD[4]);  // dPhidXdY

  PhiVD[4 ] = imag(GBarVD[0]);
  PhiVD[5 ] = imag(GBarVD[1]); 
  PhiVD[6 ] = imag(GBarVD[2]);
  PhiVD[7 ] = imag(GBarVD[4]);

  PhiVD[8 ] = real(GBarVD[3]);  // dPhidZ
  PhiVD[9 ] = real(GBarVD[5]);  // d2PhidXdZ
  PhiVD[10] = real(GBarVD[6]);  // d2PhidYdZ
  PhiVD[11] = real(GBarVD[7]);  // d3PhidXdYdZ

  PhiVD[12] = imag(GBarVD[3]);
  PhiVD[13] = imag(GBarVD[5]);
  PhiVD[14] = imag(GBarVD[6]);
  PhiVD[15] = imag(GBarVD[7]);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ReInitializeGBarInterp2D(cdouble Beta, double *K, double *L1, double *L2, double D,
                              int nThread, Interp2D *GBarInterp2D)
{ 
  GBarVDEwaldData MyData, *GD=&MyData;
  
  GD->Beta=Beta;
  GD->K=K;
  GD->L1=L1;
  GD->L2=L2;
  GD->D=D;

  Log(" Initializing 2D interpolation table (%i grid points) at D=%g",
        GBarInterp2D->N1 * GBarInterp2D->N2, D);
  GBarInterp2D->ReInitialize(nThread, GBarVDEwaldWrapper2D, (void *)GD);

} 
