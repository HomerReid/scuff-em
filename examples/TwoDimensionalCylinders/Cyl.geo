//
// gmsh geometry specification for unit circle with extra 
// refinement near the point (1,0)
//
// homer reid 
//

//Mesh.CharacteristicLengthFactor=0.08;

//************************************************************
//* input parameters      
//************************************************************
R = 1.0;    // radius

//************************************************************
//* meshing finenesses ***************************************
//************************************************************
lcoarse=0.3;
lmiddle=0.2;
lfine=0.1;

//************************************************************
//************************************************************
Point(1) = {  0 ,    0,   0,  lmiddle};
Point(2) = {  R ,    0,   0,  lfine};
Point(3) = {  0 ,    R,   0,  lmiddle};
Point(4) = { -R ,    0,   0,  lcoarse};
Point(5) = {  0 ,   -R,   0,  lmiddle};
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Line Loop(5)={1,2,3,4};
Physical Line(1)={5};
