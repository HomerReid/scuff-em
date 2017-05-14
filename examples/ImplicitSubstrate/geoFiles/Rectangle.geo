//#################################################
// GMSH geometry file for a rectangle
// 
// Usage: gmsh -2 Rectangle.geo -o OutFile.msh [ options ]
// 
// Options:
//  --setnumber LX  0.50       # width
//  --setnumber LY 10.00       # length
//  --setnumber N   4          # triangles per unit length (mesh fineness) 
//#################################################
DefineConstant [ LX = 0.50 ];
DefineConstant [ LY = 10.0 ];
DefineConstant [ N  = 4    ];

NX = Ceil[LX*N];
NY = Ceil[LY*N];

Point(1) = { -0.5*LX, -0.5*LY, 0.0 };
Point(2) = { +0.5*LX, -0.5*LY, 0.0 };
Line(1) = {1,2};
Transfinite Line{1} = Ceil[NX] + 1;

Extrude { 0, LY, 0  } { Line{1}; Layers{NY}; }
