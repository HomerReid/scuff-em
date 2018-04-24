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
 * GTransformation.h  -- a very simple data structure for 
 *                    -- handling geometric transformations
 *                        
 * homer reid         -- 11/2011
 */

#ifndef GTRANSFORMATION_H
#define GTRANSFORMATION_H

namespace scuff {

/*****************************************************************/
/* a GTransformation consists of a 3x3 matrix M and a 3-vector   */
/* DX; it maps a point with cartesian coordinates X[0..2] into   */
/* a new point with coordinates XP[0..2] such that XP = M*X + DX */
/*****************************************************************/
class GTransformation {
public:
     GTransformation();
     GTransformation(const double dx[3]); // displacement
     GTransformation(const double ZHat[3], double ThetaDegrees); // rotation
     GTransformation(const char *TransformString, char **ErrMsg = 0); // parse
     GTransformation(char **Tokens, int NumTokens, char **ErrMsg = 0,
		     int *TokensConsumed = 0); // parse

     // copy constructors
     GTransformation(const GTransformation &G);
     GTransformation(const GTransformation *G);
     void operator=(const GTransformation &G);

     ~GTransformation() {}
     // no need for destructor or assignment op since no private allocation

     void Reset(); // reset to identity transformation

     void Displace(const double dx[3], const double Scale=1.0);
     void Rotate(const double ZHat[3], double ThetaDegrees);
     //void Mirror(int Axis);

     void Transform(const GTransformation *G); // compose G * this
     void Transform(const GTransformation &G) { Transform(&G); }

     void Invert(); // invert
     GTransformation Inverse() const; // return inverse

     bool Parse(const char *TransformString, char **ErrMsg = 0);
     bool Parse(char **Tokens, int NumTokens, char **ErrMsg = 0, 
		int *TokensConsumed = 0);

     // compose this with G or inverse(G)
     GTransformation operator+(const GTransformation &G) const;
     GTransformation operator-(const GTransformation &G) const;
     GTransformation operator-() const { return Inverse(); };

     // apply out-of-place to one or NX 3-vectors,
     // where i-th point is (X[3*i+0],X[3*i+1],X[3*i+2]),
     // defined as: Apply(X) = M*X + DX.
     void Apply(const double *X, double *XP, int NX) const;
     void Apply(const double X[3], double XP[3]) const { Apply(X, XP, 1); }

     // apply in-place to one or NX 3-vectors
     void Apply(double *X, int NX) const { Apply(X, X, NX); }
     void Apply(double X[3]) const { Apply(X, X, 1); }

     // apply inverse transformation
     void UnApply(const double *X, double *XP, int NX) const;
     void UnApply(const double X[3], double XP[3]) const { UnApply(X, XP, 1); }
     void UnApply(double *X, int NX) const { UnApply(X, X, NX); }
     void UnApply(double X[3]) const { UnApply(X, X, 1); }

     // apply rotations only
     void ApplyRotation(const double X[3], double XP[3]) const;
     void ApplyRotation(double X[3]) const { ApplyRotation(X, X); }
     void UnApplyRotation(const double X[3], double XP[3]) const;
     void UnApplyRotation(double X[3]) const { UnApplyRotation(X, X); }

     bool IsIdentity(double LengthScale=1.0, double Tol=1.0e-7);
     bool IsIdentical(GTransformation GT, double LengthScale=1.0, double Tol=1.0e-7);
     
//private:
     double DX[3]; // translation
     double M[3][3]; // rotation matrix
};

typedef std::vector<GTransformation *> GTList;

/***************************************************************/
/* a GTComplex is collection of GTransformations, each of      */
/* which is carried out on a specific surface in a geometry.   */
/***************************************************************/
typedef struct GTComplex
 {
   char *Tag;          // a label for this entire complex
   GTList GTs;         // GTs[i] is the transformation applied to the ith surface
   sVec SurfaceLabels; // SurfaceLabel[i] is the label of the ith transformed surface
 } GTComplex;

// create a GTComplex that doesn't do anything
GTComplex *CreateGTComplex(const char *Tag=0);

typedef std::vector<GTComplex*> GTCList;

// this routine reads a SCUFF-EM transformation (.trans) file
// and returns an array of GTComplex structures.
GTCList ReadTransFile(const char *FileName);

void DestroyGTCList(GTCList GTCs);

} // namespace scuff

#endif // #ifndef GTRANSFORMATION_H
