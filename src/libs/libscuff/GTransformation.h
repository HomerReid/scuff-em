/*
 * GTransformation.h  -- a very simple data structure for 
 *                    -- handling geometric transformations
 *                        
 * homer reid         -- 11/2011
 */

#ifndef GTRANSFORMATION_H
#define GTRANSFORMATION_H

namespace scuff {

/***************************************************************/
/* a GTransformation maps a point with cartesian coordinates   */
/* X[0..2] into a new point with coordinates XP[0..2] such that*/
/*  XP = M*X + DX                                              */
/***************************************************************/
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

     void Displace(const double dx[3]);
     void Rotate(const double ZHat[3], double ThetaDegrees);

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
     
private:
     double DX[3]; // translation
     double M[3][3]; // rotation matrix
};

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

// create a GTComplex that doesn't do anything
GTComplex *CreateDefaultGTComplex();

// this routine reads a scuff-EM transformation (.trans) file
// and returns an array of GTComplex structures.
GTComplex **ReadTransFile(char *FileName, int *NumGTComplices);

} // namespace scuff

#endif // #ifndef GTRANSFORMATION_H
