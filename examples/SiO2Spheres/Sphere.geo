//
// gmsh geometry specification for a sphere of radius R=1
// 

//************************************************************
//* input parameters      
//************************************************************
DefineConstant[ R          = 1.0  ];
DefineConstant[ lNorthPole = 0.35 ];
DefineConstant[ lEquator   = 0.35 ];
DefineConstant[ lSouthPole = 0.35 ];
DefineConstant[ Hemispheres = 0   ];

Point(0)  = {  0,    0,   0,  lEquator*R};
Point(1)  = {  0,    0,   R,  lNorthPole*R};
Point(2)  = {  R,    0,   0,  lEquator*R};
Point(3)  = {  0,    0,  -R,  lSouthPole*R};

If(Hemispheres==0)
  Circle(1) = {1,0,2};
  Circle(2) = {2,0,3};
  result = Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Line{1,2}; };
  result = Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Line{result[0],result[3]}; };
  result = Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Line{result[0],result[3]}; };
EndIf

If(Hemispheres==1)
  Circle(1) = {1,0,2};
  Circle(2) = {2,0,3};
  Line(3)   = {0, 2};
  result = Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Line{1,2,3}; };
  result = Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Line{result[0],result[3],result[6]}; };
  result = Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Line{result[0],result[3],result[6]}; };
  Physical Surface(1) = { 6, 15, 24}; // northern hemisphere
  Physical Surface(2) = { 9, 18, 27}; // southern hemisphere
  Physical Surface(3) = {12, 21, 30}; // equatorial plane   
EndIf
