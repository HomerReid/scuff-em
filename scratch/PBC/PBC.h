#ifndef PBC_H
#define PBC_H

typedef struct PBCAccelerator
 {
   HMatrix *ImageBlocks;

   Interp3D *GBarAB9Interp;

 };

void AddStraddlers(RWGObject *O, double **LBV);

void AssembleBEMMatrix_PBC(RWGGeometry *G, double **LBV, double *P, 
                           PBCAccelerator *PBCA, HMatrix *M);

PBCAccelerator *CreatePBCAccelerator(RWGGeometry *G);
void SetPBC(PBCAccelerator PBCA, cdouble Omega, double *P);


#endif
