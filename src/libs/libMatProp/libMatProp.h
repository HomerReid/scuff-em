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
 * libMatProp.h  -- a simple C++ class for managing material
 *               -- properties (epsilon and mu as a function
 *               -- of frequency) for scuff-EM
 *
 * homer reid    -- 12/2009
 */

#ifndef LIBMATPROP_H
#define LIBMATPROP_H

#include <libMDInterp.h>

/***************************************************************/
/* constants      **********************************************/
/***************************************************************/

// values for the Type field of the MatProp class

#define MP_VACUUM   0  // free space 
#define MP_PEC      1  // perfect electrical conductor 
#define MP_CONSTANT 2  // constant frequency-independent properties
#define MP_INTERP   3  // interpolate from user-supplied data tables
#define MP_PARSED   4  // user-supplied parsed expressions for eps,mu
#define MP_OTHER    5  // everything else 

/***************************************************************/
/* MatProp class definition ************************************/
/***************************************************************/
class MatProp
 { 
  public:  

   /***************************************************************/
   /* public class methods ****************************************/
   /***************************************************************/

   /* constructor entry points */
   MatProp();
   MatProp(int Type);
   MatProp(const char *MaterialName);
   MatProp(const char *MaterialName, const char *MatPropFileName);
   MatProp(MatProp *MP);

  /* destructor */
   ~MatProp();

   /* set the angular frequency unit; the Omega parameter     */
   /* to each of the following routines will be multiplied by */
   /* parameter.                                              */
   /* SetLengthUnit(L) is equivalent to SetFreqUnit(L/c).     */
   static void SetFreqUnit(double NewFreqUnit);
   static void SetLengthUnit(double NewLengthUnit);

   /* get epsilon and mu at a given frequency */
   void GetEpsMu(cdouble Omega, cdouble *Eps, cdouble *Mu);
   cdouble GetEps(cdouble Omega);
   cdouble GetMu(cdouble Omega);

   /* get index of refraction and relative wave impedance at given freq*/
   cdouble GetRefractiveIndex(cdouble Omega, cdouble *ZRel=0);

   /* set constant eps/mu */
   void SetEpsMu(cdouble pEps, cdouble pMu);
   void SetEps(cdouble pEps) { SetEpsMu(pEps, 1.0); }

   /* return 1 if the material is a perfect electrical conductor, 0 othrws */
   int IsPEC();

   /* the following two functions implement a simple mechanism for */
   /* temporarily zeroing out a MatProp so that calls to GetEpsMu()*/
   /* return Eps=Mu=0.0 at all frequencies                         */
   void Zero() { Zeroed=1; }
   void UnZero() { Zeroed=0; }

   /* if ErrMsg is not NULL after the class constructor is invoked, there  */
   /* was an error.                                                        */
   char *ErrMsg;

 // private:

   /* the actual body of the class constructor */
   void InitMatProp(const char *MaterialName, const char *MatPropFileName);

   /* constructor helper function for tabulated-data materials */
   void ReadInterpolationTable(const char *FileName);

   /* constructor helper functions for user-defined materials */
   void CreateUserDefinedMaterial(const char *MatPropFileName, const char *MaterialName);
   int ReadMaterialFromFile(const char *FileName, const char *MaterialName);
   int ParseMaterialSectionInFile(FILE *f, const char *FileName, int *LineNum);
   void GetEpsMu_Parsed(cdouble Omega, cdouble *pEps, cdouble *pMu);

   /***************************************************************/
   /* class data **************************************************/
   /***************************************************************/
   int Type;
   int Zeroed;
   char *Name;

   // constant values for the case Type=MP_CONSTANT 
   cdouble Eps;
   cdouble  Mu;

   // interpolator for real and imaginary frequency axes, for Type=MP_INTERP
   Interp1D *InterpReal, *InterpImag;
   bool OwnsInterpolators;

   // opaque pointer to cmatheval parsed expressions
   void *EpsExpression, *MuExpression;
   bool OwnsExpressions;

   // angular frequency unit (common to all instances of MatProp)
   static double FreqUnit;

 };

/***************************************************************/
/* this routine allows user-defined materials to be defined    */
/* on-the-fly without needed to be found in a MatProp database */
/* file..                                                      */
/***************************************************************/
char *AddMaterialToMatPropDataBase(FILE *f, char *FileName,
                                   char *MaterialName, int *LineNum);

#endif
