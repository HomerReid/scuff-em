//////////////////////////////////////////////////
// GMSH geometry file for a structured rectangle
//
// Usage:
//  gmsh [options] -2 Rectangle.geo -o OutputMesh.msh
//
// Options:
//  -setnumber W  1   // width
//  -setnumber L  10  // length
//  -setnumber N  2   // edges per unit area
//  -setnumber FR 3   // meshing fineness ratio, left/right edges
//////////////////////////////////////////////////

// user-adjustable constants
DefineConstant[  W =  1.0 ];
DefineConstant[  L = 10.0 ];
DefineConstant[  N =    2 ];
DefineConstant[ FR =    1 ];

// derived parameters
NW = Ceil[N*W];
NL = Ceil[N*L];

//////////////////////////////////////////////////
// geometry description
//////////////////////////////////////////////////
If (FR==1)
  Point(1)  = {-0.5*W, 0.0, 0.0};
  Point(2)  = {+0.5*W, 0.0, 0.0};
  Line(1)   = {1,2}; 
  Transfinite Line{1} = NW+1;
EndIf

If (FR!=1)
  lF = 1.0/N;
  lC = FR*lF;
  Point(1)  = {-0.5*W, 0.0, 0.0, lC };
  Point(2)  = {+0.5*W, 0.0, 0.0, lF };
  Line(1)   = {1,2}; 
EndIf

Extrude{ 0, L, 0} { Line{1}; Layers{NL}; }
