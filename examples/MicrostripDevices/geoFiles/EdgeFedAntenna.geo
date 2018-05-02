//////////////////////////////////////////////////// 
// GMSH .geo file for the patch antenna analyzed by Wu et al,
// "Feeding Structure Contribution to Radiation by Patch Antennas 
// with Rectangular Boundaries", *IEEE Trans. Ant. Prop.* **40** 1245–1249 (Oct. 1992.)
//////////////////////////////////////////////////// 

DefineConstant [ W     = 12.448 ]; // width of patch (X direction)
DefineConstant [ L     = 16.0   ]; // length of patch (Y direction)

DefineConstant [ WFeed = 2.334  ]; // width of feed
DefineConstant [ LFeed = 8.0    ]; // length of feed (may be zero)
DefineConstant [ OFeed = 8.169  ]; // offset of feed from right side of patch

DefineConstant [ NFeed = WFeed  ];

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
X1 = -0.5*W;
X2 = +0.5*W - OFeed - WFeed;
X3 = +0.5*W - OFeed;
X4 = +0.5*W;

N = Ceil[NFeed/ WFeed];

//////////////////////////////////////////////////
// patch 
//////////////////////////////////////////////////
Point(101) = {X1, LFeed, 0.0};
Point(102) = {X2, LFeed, 0.0};
Point(103) = {X3, LFeed, 0.0};
Point(104) = {X4, LFeed, 0.0};
Line(101)  = {101,102};
Line(102)  = {102,103};
Line(103)  = {103,104};
Transfinite Line{101} = 1 + Ceil[0.5*N*(X2-X1)];
Transfinite Line{102} = 1 + NFeed;
Transfinite Line{103} = 1 + Ceil[0.5*N*(X4-X3)];
Extrude { 0, L, 0  } { Line{101,102,103}; Layers{Ceil[0.5*N*L]}; }

//////////////////////////////////////////////////
// feed line
//////////////////////////////////////////////////
If (LFeed>0.0)
  Extrude { 0, -LFeed, 0 } { Line{102}; Layers{Ceil[N*LFeed]}; }
EndIf
