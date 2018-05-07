//////////////////////////////////////////////////// 
// GMSH .geo file for the patch antenna analyzed by Wu et al,
// "Feeding Structure Contribution to Radiation by Patch Antennas 
// with Rectangular Boundaries", *IEEE Trans. Ant. Prop.* **40** 1245–1249 (Oct. 1992.)
//
// Usage: gmsh -2 EdgeFedAntenna.geo -o OutputMesh.msh [options]
//
// Available options:
//
//  -setnumber W     12.448  // width of patch in X direction
//  -setnumber L     12.448  // width of patch in Y direction
//
//  -setnumber WFeed  2.334  // width of feed line
//  -setnumber LFeed  8.0    // length of feed line
//  -setnumber OFeed  8.169  // offset of feed from right edge of patch
//
//  -setnumber N      2      // meshing fineness (#triangle edges per unit length)
//////////////////////////////////////////////////// 

// user-tweakable constants

DefineConstant [ W     = 12.448 ]; // width of patch (X direction)
DefineConstant [ L     = 16.0   ]; // length of patch (Y direction)

DefineConstant [ WFeed = 2.334  ]; // width of feed
DefineConstant [ LFeed = 8.0    ]; // length of feed (may be zero)
DefineConstant [ OFeed = 8.169  ]; // offset of feed from right side of patch

DefineConstant [ N     = 2      ]; // meshing fineness (triangle edges per unit length)

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
X1 = -0.5*W;
X2 = +0.5*W - OFeed - WFeed;
X3 = +0.5*W - OFeed;
X4 = +0.5*W;

// place origin at center of feed line
X0 = 0.5*(X2+X3);

//////////////////////////////////////////////////
// patch 
//////////////////////////////////////////////////
Point(101) = {X1-X0, LFeed, 0.0};
Point(102) = {X2-X0, LFeed, 0.0};
Point(103) = {X3-X0, LFeed, 0.0};
Point(104) = {X4-X0, LFeed, 0.0};
Line(101)  = {101,102};
Line(102)  = {102,103};
Line(103)  = {103,104};
Transfinite Line{101} = 1 + Ceil[N*(X2-X1)];
Transfinite Line{102} = 1 + Floor[N*(X3-X2)];
Transfinite Line{103} = 1 + Floor[N*(X4-X3)];
Extrude { 0, L, 0  } { Line{101,102,103}; Layers{Floor[N*L]}; }

//////////////////////////////////////////////////
// feed line
//////////////////////////////////////////////////
If (LFeed>0.0)
  Extrude { 0, -LFeed, 0 } { Line{102}; Layers{Ceil[N*LFeed]}; }
EndIf
