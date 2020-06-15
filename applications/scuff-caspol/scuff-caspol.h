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
 * scuff-caspol.h    -- header file for scuff-caspol
 *
 * homer reid        -- 1/2012
 *
 */
#ifndef SCUFF_CASPOL_H
#define SCUFF_CASPOL_H

#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libMDInterp.h>
#include <libTriInt.h>
#include <BZIntegration.h>
#include <libscuff.h>

using namespace scuff;

// values for the Which field of the PolModel constructor
#define PM_BUILTIN 0 
#define PM_FILE    1 

// values for the Class field of the PolModel class
#define PM_SCALAR   0
#define PM_DIAGONAL 1
#define PM_GENERAL  2

#define FILETYPE_OUT   0
#define FILETYPE_BYXI  1
#define FILETYPE_BYXIK 2

/***************************************************************/
/* PolModel is a simple class used to describe the frequency-  */
/* dependent polarizability of atoms and molecules. it exports */
/* just a single class method, GetPolarizability(), which takes*/
/* an imaginary frequency as an input and computes values of   */
/* the polarizability tensor at that frequency.                */
/*                                                             */
/* for the time being, the only available polarizabilities     */
/* are built-in models for a few atomic species. this could    */
/* eventually be generalized to arbitrary user-specified       */
/* polarizability data.                                        */
/***************************************************************/
class PolModel
 {
public:

   // constructor entry point
   PolModel(char *Name, int Which=PM_FILE);

   // constructors for built-in and user-defined models
   void InitPolModel_BI(char *Atom);
   void InitPolModel_UD(char *FileName);

   // routine to compute the polarizability tensor
   void GetPolarizability(double Xi, HMatrix *Alpha);

   // fields used for all classes and implementations
   char *Name;
   char *ErrMsg;

   // data that are present only for MatProp models
   MatProp *MP;

   // data that are present only for interpolated models
   int NumPoints;
   double *XiPoints, *PolPoints;
   Interp1D *PolInterp;
   double LargeXiCoefficient;

 };

/***************************************************************/
/* SCPData ('scuff-caspol-data') is the basic structure passed */
/* around among the various routines; it contains everything   */
/* needed to compute the casimir-polder potential              */
/***************************************************************/
typedef struct SCPData
 {
   RWGGeometry *G;
   HMatrix *M;
   HVector *KN;

   int NumAtoms;
   PolModel **PolModels;
   HMatrix **Alphas;

   HMatrix *EPMatrix;

   double RelTol;
   double AbsTol;

   char *ByXiFileName; 
   char *ErrMsg;

   double Xi;

   // these items only used for periodic geometries
   void **ABMBCache;
   char *ByXikFileName;
   GetBZIArgStruct *BZIArgs;

 } SCPData; 

/***************************************************************/
/***************************************************************/
/***************************************************************/
SCPData *CreateSCPData(char *GeoFile,
                       char **Atoms, int NumBIAtoms,
                       char **Particles, int NumParticles,
                       char *EPFile, char *FileBase);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetCPIntegrand(void *pSCPD, cdouble Omega,
                    double *kBloch, double *U);
void GetXiIntegrand(SCPData *SCP, double Xi, double *U);
void EvaluateFrequencyIntegral(SCPData *SCP, double *U);
void EvaluateMatsubaraSum(SCPData *SCPD, double Temperature, double *U);
/***************************************************************/
void LoadPolarizability(SCPData *SCPD, double Xi);

#endif
