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


#endif
