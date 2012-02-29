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

namespace scuff{

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

// create or augment using a known displacement vector
GTransformation *CreateOrAugmentGTransformation(GTransformation *GT, double *DX);
GTransformation *CreateGTransformation(double *DX);
 
// create or augment using a known rotation angle and axis 
GTransformation *CreateOrAugmentGTransformation(GTransformation *GT,
                                                double *ZHat, double Theta);
GTransformation *CreateGTransformation(double *ZHat, double Theta);

// create or augment by parsing a character string.   
GTransformation *CreateOrAugmentGTransformation(GTransformation *GT, 
                                                char *TransformString,
                                                char **ErrMsg);
GTransformation *CreateOrAugmentGTransformation(GTransformation *GT, 
                                                char **Tokens, int NumTokens,
                                                char **ErrMsg);
// identity transformation 
GTransformation *CreateGTransformation(); // identity transformation

void AugmentGTransformation(GTransformation *DeltaGT, GTransformation *GT);

void ResetGTransformation(GTransformation *GT);

/***************************************************************/
/* routines to apply a GTransformation to a single point or to */
/* a list of points.                                           */
/* coordinates are assumed to be ordered as follows:           */
/*  X[0, 1, 2] == x,y,z coords of first point                  */
/*  X[3, 4, 5] == x,y,z coords of second point                 */
/* etc.                                                        */
/***************************************************************/
/* in-place */
void ApplyGTransformation(GTransformation *GT, double *X, int NX);
void ApplyGTransformation(GTransformation *GT, double *X);

/* out-of-place */
void ApplyGTransformation(GTransformation *GT, double *X, double *XP, int NX);
void ApplyGTransformation(GTransformation *GT, double *X, double *XP);

void UnApplyGTransformation(GTransformation *GT, double *X, int NX);

/***************************************************************/
/* a GTComplex is collection of GTransformations, each of      */
/* which is carried out on a different object in a geometry.   */
/***************************************************************/
typedef struct GTComplex
 {
   char *Tag;                  // a label for this entire complex 
   int NumObjectsAffected;     // number of objects transformed by this complex
   char **ObjectLabel;         // ObjectLabel[i] is the label of the ith transformed object
   GTransformation **GT;       // GT[i] is the transformation applied to the ith object

 } GTComplex;

// this routine reads a scuff-EM transformation (.trans) file
// and returns an array of GTComplex structures.
GTComplex **ReadTransFile(char *FileName, int *NumGTComplices);

} // namespace scuff

#endif // #ifndef GTRANSFORMATION_H
