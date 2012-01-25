/*
 * GTransformation.h  -- a very simple data structure for 
 *                    -- handling geometric transformations
 *                        
 * homer reid         -- 11/2011
 */

#ifndef GTRANSFORMATION_H
#define GTRANSFORMATION_H

#define GTRANSFORMATION_DISPLACEMENT 1
#define GTRANSFORMATION_ROTATION     2

/***************************************************************/
/* a GTransformation maps a point with cartesian coordinates   */
/* X[0..2] into a new point with coordinates XP[0..2] such that*/
/*  XP = M*X + DX                                              */
/***************************************************************/
typedef struct GTransformation
 { int Type;
   double DX[3];
   double M[3][3];
 } GTransformation;

/***************************************************************/
/* routines to create a new GTransformation, or to augment an  */
/* existing one                                                */
/***************************************************************/
GTransformation *CreateOrAugmentGTransformation(GTransformation *GT, double *DX);
GTransformation *CreateGTransformation(double *DX);

GTransformation *CreateOrAugmentGTransformation(GTransformation *GT,
                                                double *ZHat, double Theta);
GTransformation *CreateGTransformation(double *ZHat, double Theta);

void ResetGTransformation(GTransformation *GT);

/* apply in-place */
void ApplyGTransformation(GTransformation *GT, double *X, int NX);
void ApplyGTransformation(GTransformation *GT, double *X);

/* apply out-of-place */
void ApplyGTransformation(GTransformation *GT, double *X, double *XP, int NX);
void ApplyGTransformation(GTransformation *GT, double *X, double *XP);

void UnApplyGTransformation(GTransformation *GT, double *X, int NX);

#endif // #ifndef GTRANSFORMATION_H
