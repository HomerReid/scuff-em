//////////////////////////////////////////////////
// GMSH surface-mesh geometry for a square to be 
// discretized into a structured mesh 
//
// homer reid 4/2013
// 
// options:
//  -setnumber N   10         (#triangle edges per unit length)
//  -setnumber L   1.0        (side length)
//  -setnumber Z0  1.0        (height above xy plane)
//////////////////////////////////////////////////

DefineConstant[ N  = 10   ];
DefineConstant[ L  = 1.0  ];
DefineConstant[ Z0 = 1.0  ];

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
// center 
X0 = 0;
Y0 = 0;

NTotal=N*L;

////////////////////////////////////////////////////
// end of user-specifiable parameters 
////////////////////////////////////////////////////
Point(1) = { 0, 0, Z0 };
Point(2) = { L, 0, Z0 };

Line(1)  = { 1, 2 };

Transfinite Line{1} = NTotal+1;
Extrude { 0, L, 0  } { Line{1}; Layers{NTotal}; }
