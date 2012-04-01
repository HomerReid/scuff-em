/*
   FieldGrid.h -- classes and functions so that we can compute
                  arbitrary functions of the fields on an arbitrary 2d grid 
		  for convenient visualization and other analyses

   SGJ, 3/2012
*/

#ifndef SCUFF_FIELDGRID_H
#define SCUFF_FIELDGRID_H 1

#include <libhmat.h>
#include "GTransformation.h"

namespace scuff {

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
// Functions of the fields (and coords, etc) used for visualization etc.

// subclass to define arbitrary functions of the fields that we will
// evaluate on 2d grids to obtain matrices for visualization
class FieldFunc {
public:
     // whether the function is purely real-valued at all points
     virtual bool IsReal() const { return false; }

     // the value of the function for each point and fields
     virtual cdouble Eval(const double X[3], // (x,y,z)
			  const double dA[3], // unit-normal area element
			  const cdouble EH[3], // E and H fields
			  cdouble Eps, cdouble Mu) const = 0; // materials at X

     virtual ~FieldFunc() {} // need virtual destructor for subclasses to work
};

// Parse arbitrary user-defined cmatheval string in terms of the variables:
//             Ex,Ey,Ez/Hx,Hy,Hz -- fields in Cartesian coords, SI units
//             x,y,z -- Cartesian coordinates
//             dAx,dAy,dAz -- surface-normal area of grid cell
//             dA -- |dA|
//             nx,ny,nz -- unit normal to grid surface = dA / |dA|
//             Eps,Mu -- relative permittivity/permeability
//             eps,mu -- absolute permittivity/permeability
//             eps0,m0,c,Z0 -- vacuum permittivity/permeability/speed/impedance
// (This should be plenty fast in comparison to the cost of GetField
//  and should be convenient enough for most common field computations
//  that most users will not need or want any other subclasses.)
class ParsedFieldFunc : public FieldFunc {
private:
     void *Expr; // parsed expression from cmatheval
     bool isReal; // whether expression is real-valued
public:
     ParsedFieldFunc(const char *ExprString);
     ~ParsedFieldFunc();

     cdouble Eval(const double X[3], // (x,y,z)                         
		  const double dA[3], // unit-normal area element      
		  const cdouble EH[3], // E and H fields                
		  cdouble Eps, cdouble Mu) const;

     bool IsReal() const { return isReal; }

     // return (simplified) expression string ... caller must free
     char *String() const; 
};

// A few useful non-trivial expressions for the parser:
#define POYNTINGx_EXPR "0.5*real(conj(Ey)*Hz-conj(Ez)*Hy)" // x Poynting vector
#define POYNTINGy_EXPR "0.5*real(conj(Ez)*Hx-conj(Ex)*Hz)" // y Poynting vector
#define POYNTINGz_EXPR "0.5*real(conj(Ex)*Hy-conj(Ey)*Hx)" // z Poynting vector
#define POYNTING_EXPR "0.5*real(nx*(conj(Ey)*Hz-conj(Ez)*Hy) + ny*(conj(Ez)*Hx-conj(Ex)*Hz) + nz*(conj(Ex)*Hy-conj(Ey)*Hx))" // Poynting vector in grid-normal direction
#define ENERGYDENSITY_EXPR "0.25*real(eps*(|Ex|^2+|Ey|^2+|Ez|^2) + mu*(|Hx|^2+|Hy|^2+|Hz|^2))" // electromagnetic energy density per unit volue (neglecting dispersion)
#define E_ENERGYDENSITY_EXPR "0.25*real(eps*(|Ex|^2+|Ey|^2+|Ez|^2))" // electric energy density per unit volume (neglecting dispersion)
#define H_ENERGYDENSITY_EXPR "0.25*real(mu*(|Hx|^2+|Hy|^2+|Hz|^2))" // magnetic energy density per unit volume (neglecting dispersion)
#define ABSORPTION_EXPR "0.5*imag(eps*(|Ex|^2+|Ey|^2+|Ez|^2) + mu*(|Hx|^2+|Hy|^2+|Hz|^2))" // absorption loss per unit time & volume, divided by omega

// Maxwell stress tensor (problematic in non-vacuum materials)
#define FORCEx_DENSITY_EXPR "0.5*real(eps*conj(Ex)*(Ex*nx+Ey*ny+Ez*nz) + mu*conj(Hx)*(Hx*nx+Hy*ny+Hz*nz) - 0.5*nx*(eps*(|Ex|^2+|Ey|^2+|Ez|^2) + mu*(|Hx|^2+|Hy|^2+|Hz|^2)))"
#define FORCEy_DENSITY_EXPR "0.5*real(eps*conj(Ey)*(Ex*nx+Ey*ny+Ez*nz) + mu*conj(Hy)*(Hx*nx+Hy*ny+Hz*nz) - 0.5*ny*(eps*(|Ex|^2+|Ey|^2+|Ez|^2) + mu*(|Hx|^2+|Hy|^2+|Hz|^2)))"
#define FORCEz_DENSITY_EXPR "0.5*real(eps*conj(Ez)*(Ex*nx+Ey*ny+Ez*nz) + mu*conj(Hz)*(Hx*nx+Hy*ny+Hz*nz) - 0.5*nz*(eps*(|Ex|^2+|Ey|^2+|Ez|^2) + mu*(|Hx|^2+|Hy|^2+|Hz|^2)))"
#define FORCE_DENSITY_EXPR "nx * " FORCEx_DENSITY_EXPR " + ny * " FORCEy_DENSITY_EXPR " + nz * " FORCEz_DENSITY_EXPR // force normal to grid surface

// (It is convenient to #define these rather than using const char *,
//  since this way you can combine them just by putting adjacent strings.
//  e.g. the area-weighted Poynting flux is POYNTING_EXPR "* dA".)

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/* The following classes are used to describe 2d surfaces that are mapped
   onto 2d rectangular grids (matrices) */

// subclass to define arbitrary 2d grids in 3-space
class SurfaceGrid {
public:
     int N1, N2; // grid dimensions
     
     SurfaceGrid() : N1(0), N2(0) {}
     SurfaceGrid(int n1, int n2) : N1(n1), N2(n2) {}
     
     HMatrix *AllocateGrid(bool IsReal) const {
	  return new HMatrix(N1,N2, IsReal ? LHM_REAL : LHM_COMPLEX);
     }
     
     // return coordinates X for a point (n1,n2) in [0,N1-1] x [0,N2-1]
     // along with the unit normal area element dA
     virtual void GetPoint(int n1, int n2, double X[3], double dA[3]) const=0;

     virtual void Transform(const GTransformation &G) = 0;
};

/***************************************************************************/
// Useful predefined SurfaceGrid subclasses.

enum CartesianDirection { XDIR = 0, YDIR, ZDIR };

// planar surface
class PlaneGrid : public SurfaceGrid {
public:
     double X0[3]; // (0,0) corner
     double S1[3]; // (N1,0) corner is X0 + S1*N1
     double S2[3]; // (0,N2) corner is X0 + S2*N2
     double dA0[3]; // area element, in normal direction

     // create a planar grid with center c0 and corners
     // (for infinite resolution) at c0+/-(s1/2)+/-(s2/2)
     PlaneGrid(int n1, int n2, const double c0[3],
	       const double s1[3], const double s2[3]);

     // constructor for common case of rectangular s1xs2 grid normal to x/y/z
     PlaneGrid(int n1, int n2, const double c0[3],
	       double s1, double s2, CartesianDirection normal);
     
     void GetPoint(int n1, int n2, double X[3], double dA[3]) const;

     void Transform(const GTransformation &G) {
	  G.Apply(X0); G.ApplyRotation(S1); G.ApplyRotation(S2);
	  G.ApplyRotation(dA0);
     }
};

#if 0 // TODO

// spherical surface (latitude/longitude grid)
class SphereGrid : public SurfaceGrid {
public:
     double X0[3]; // origin
     double R; // radius

     SphereGrid(int n1, int n2, const double x0[3], double r);
     
     void GetPoint(int n1, int n2, double X[3], double dA[3]) const;

     void Transform(const GTransformation &G) { G.Apply(X0); }
};

// cylindrical surface (angle/z grid)
class CylinderGrid : public SurfaceGrid {
public:
     double X0[3]; // center of cylinder bottom
     double S2[3]; // X0+S2 is center of cylinder top
     double R; // radius

     CylinderGrid(int n1, int n2, const double x0[3], 
		  const double s2[3], double r);

     // common case of a cylinder along the x/y/z axes
     CylinderGrid(int n1, int n2, const double x0[3], 
		  double s2, CartesianDirection axis);
     
     void GetPoint(int n1, int n2, double X[3], double dA[3]) const;

     void Transform(const GTransformation &G) {
	  G.Apply(X0); G.ApplyRotation(S2);
     }
};

#endif

/***************************************************************************/

} // namespace scuff

#endif // SCUFF_FIELDGRID_H
