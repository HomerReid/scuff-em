//////////////////////////////////////////////////
// GMSH surface-mesh geometry for a square to be 
// discretized into a structured mesh 
//
// homer reid 4/2013
// 
// options:
//  -setnumber N   10         (#triangle edges per unit length)
//  -setnumber L   1.0        (side length)
//////////////////////////////////////////////////

DefineConstant[ N             = 10   ];
DefineConstant[ L             = 1.0  ];

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
Point(1) = { 0, 0, 0, lc};
Point(2) = { L, 0, 0, lc};
Point(3) = { L, L, 0, lc};
Point(4) = { 0, L, 0, lc};

Line(1)  = { 1, 2 };
//Line(2)  = { 2, 3 };
//Line(3)  = { 3, 4 };
//Line(4)  = { 4, 1 };

Transfinite Line{1} = NTotal+1;
Extrude { 0, L, 0  } { Line{1}; Layers{NTotal}; }
