//////////////////////////////////////////////////// 
// GMSH .geo file for the point-fed patch antenna analyzed in Figure 2
// of Jedlicka et al, "Measured Mutual Coupling Between Microstrip Antennas,"
// *IEEE Trans. Ant. Prop.* **29** 147 (January 1981).
// Note that the feed point (which lies in the interior of the
// patch, displaced from the center) is the origin of coordinates.
////////////////////////////////////////////////////

DefineConstant [ L      = 105.7  ]; // width of patch (X direction)
DefineConstant [ W      = 65.5   ]; // length of patch (Y direction)
DefineConstant [ Offset = 34.8   ]; // offset of feed point from edge
DefineConstant [ N      = 0.5    ]; // meshing fineness (panel edges per unit area)( (panel edges per unit area)( (panel edges per unit area)( (panel edges per unit area)

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
X1 = -0.5*L;
X2 = +0.5*L-Offset;
X3 = +0.5*L;

Point(1) = { X1-X2, 0.0, 0.0};
Point(2) = { 0.0,   0.0, 0.0};
Point(3) = { X3-X2, 0.0, 0.0};
Line(1)  = {1,2};
Line(2)  = {3,2};
 
// use structured meshes if N==0
If (N>0.0)
  Transfinite Line{1} = 1 + Ceil[N*(X2-X1)];
  Transfinite Line{2} = 1 + Ceil[N*(X3-X2)];
  Extrude{ 0,  0.5*W, 0 } { Line{1,2}; Layers{ Ceil[0.5*N*W] }; }
  Extrude{ 0, -0.5*W, 0 } { Line{1,2}; Layers{ Ceil[0.5*N*W] }; }
EndIf

// otherwise use unstructured meshes
If (N==0.0)
  Extrude{ 0,  0.5*W, 0 } { Line{1}; }
  Extrude{ 0, -0.5*W, 0 } { Line{1}; }
EndIf
