//#################################################
// GMSH geometry file for a coplanar waveguide
// 
// Usage: gmsh -2 CPW.geo -o OutFile.msh [ options ]
// 
// Options:
//  --setnumber L  10.00       # length of waveguide in Y direction
//  --setnumber XA  0.25       # X coordinate of right edge of center conductor
//  --setnumber XB  0.50       # X coordinate of left edge of right conductor
//  --setnumber XC  1.50       # X coordinate of right edge of right conductor
//  --setnumber N   4          # triangles per unit length (mesh fineness) 
// 
// Reference: Figure 2 (a) in the following paper:
//  -- Chen and Chou, "Characteristics of CTLs on Multilayer Substrates,"
//     *IEEE Transactions on Microwave Theory and Techniques* **45** 6 (1997)
//#################################################

DefineConstant [ L  = 10.0 ];
DefineConstant [ XA = 0.25 ];
DefineConstant [ XB = 0.50 ];
DefineConstant [ XC = 1.50 ];
DefineConstant [ N  = 4    ];

Point(1) = { -XA, -0.5*L, 0.0 };
Point(2) = {  XA, -0.5*L, 0.0 };
Point(3) = {  XB, -0.5*L, 0.0 };
Point(4) = {  XC, -0.5*L, 0.0 };

Line(1) = {1,2};
Transfinite Line{1} = Ceil[2*XA*N] + 1;

Line(2) = {3,4};
Transfinite Line{2} = Ceil[(XC-XB)*N] + 1;

Extrude { 0, L, 0  } { Line{1,2}; Layers{N*L}; }

Physical Surface(1) = { 6  };
Physical Surface(2) = { 10 };
