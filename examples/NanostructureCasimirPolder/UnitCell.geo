//////////////////////////////////////////////////
// GMSH geometry description for the unit cell of
// the structure shown in Figure 1(a) of the 
// paper ''Trapped atoms in one-dimensional
// photonic crystals'' by Hung et al, New Journal
// of Physics *15* (2013) 083036
//
// This geometry may also be used to make meshes for
// a 2D-lattice version of the same geometry, or 
// for a non-periodic geometry consisting of just 
// one isolated cell, by specifying the "LDim"
// ("lattice dimension") parameter on the GMSH 
// command line:
// 
//  gmsh -setnumber LDim 0
//  gmsh -setnumber LDim 1 (default)
//  gmsh -setnumber LDim 2
//
// Homer Reid 11/2015
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// user-specified parameters
//////////////////////////////////////////////////
DefineConstant[ a  = 0.367 ];
DefineConstant[ W  = 0.845 ];
DefineConstant[ T  = 0.825 ];
DefineConstant[ Hx = 0.246 ];
DefineConstant[ Hy = 0.745 ];

DefineConstant[ LDim = 1 ];

//////////////////////////////////////////////////
// this sets the meshing resolution (higher is finer).
// N is roughly the number of triangle edges per lattice
// constant
//////////////////////////////////////////////////
DefineConstant[ N = 4 ];

//////////////////////////////////////////////////
// end of user-specified parameters
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// derived parameters 
//////////////////////////////////////////////////
L=a/N;
NX=N;
NY=Ceil[W/L];
NZ=Ceil[T/L];

DeltaX = 0.5*(a-Hx);
DeltaY = 0.5*(W-Hy);

X1 = DeltaX;
X2 = 0.5*a;
X3 = a-DeltaX;

YMin = 0; // -W/2
YMax = W; // +W/2

Y1 = YMin + DeltaY;
Y2 = YMin + DeltaY + 0.5*Hx;
Y3 = YMax - DeltaY - 0.5*Hx;
Y4 = YMax - DeltaY;

Z1 = -T/2;

//////////////////////////////////////////////////
// points and lines defining the cross-sectional 
// profile in the xy plane
//////////////////////////////////////////////////
Point(100)  = {X2, Y2, Z1, L};
Point(101)  = {X2, Y1, Z1, L};
Point(102)  = {X1, Y2, Z1, L};
Point(103)  = {X1, Y3, Z1, L};
Point(104)  = {X2, Y4, Z1, L};
Point(105)  = {X3, Y3, Z1, L};
Point(106)  = {X3, Y2, Z1, L};
Point(110)  = {X2, Y3, Z1, L};

Circle(101) = {101,100,102};
Line(102)   = {102,103};
Circle(103) = {103,110,104};
Circle(104) = {104,110,105};
Line(105)   = {105,106};
Circle(106) = {106,100,101};
Line Loop(100) = {101, 102, 103, 104, 105, 106};

Point(201)  = {0,  0.0, Z1, L};
Point(202)  = {0,    W, Z1, L};
Point(203)  = {a,    W, Z1, L};
Point(204)  = {a,  0.0, Z1, L};
Line(201)   = {201, 202};
Line(202)   = {202, 203};
Line(203)   = {203, 204};
Line(204)   = {204, 201};
Line Loop(200) = {201, 202, 203, 204};

Ruled Surface(300)={200,100};

//////////////////////////////////////////////////
// extrude cross section in the T direction to define
// surfaces
//////////////////////////////////////////////////
Extrude{ 0, 0, T } { Surface{300}; Layers{NZ}; }

//////////////////////////////////////////////////
// if we are making a mesh for a compact (non-periodic)
// geometry, we need to include the side walls in the 
// X and Y directions.
// 
// if we are making a mesh for the unit cell of a 
// 1D periodic geometry, then we don't include the 
// X sidewalls, but we do include the Y sidewalls.
// 
// if we are making a mesh for the unit cell of a 
// 2D periodic geometry, then we exclude both
// the X and Y sidewalls.
//////////////////////////////////////////////////
//////////////////////////////////////////////////
If (LDim==0)
  Physical Surface(1) = { 300, 315, 319, 323, 327, 331, 335, 339, 343, 347, 351, 352 };
EndIf

If (LDim==1)
  Physical Surface(1) = { 300, 319, 327, 331, 335, 339, 343, 347, 351, 352 };
EndIf

If (LDim==2)
  Physical Surface(1) = { 300, 331, 335, 339, 343, 347, 351, 352 };
EndIf
