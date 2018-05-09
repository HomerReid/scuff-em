//
// GMSH geometry file for a thin round plate with a circular hole
// 
// homer reid 06/2012
//

//
// user-specifiable parameters 
//

ROuter = 10;      // outer plate radius
RInner =  1;      // inner plate radius
T      =  1;      // thickness

lOuter = ROuter/4; // meshing fineness near outer circumference
lInner = RInner/4; // meshing fineness near hole


//
// geometry description
//

Point(100) = {       0,         0,  T/2,  lInner};
Point(101) = {  RInner,         0,  T/2,  lInner};
Point(102) = {       0,    RInner,  T/2,  lInner};
Point(103) = { -RInner,         0,  T/2,  lInner};
Point(104) = {       0,   -RInner,  T/2,  lInner};

Circle (100) = { 101, 100, 102 };
Circle (101) = { 102, 100, 103 };
Circle (102) = { 103, 100, 104 };
Circle (103) = { 104, 100, 101 };
Line Loop(100) = { 100, 101, 102, 103};

Point(201) = {  ROuter,         0,  T/2,  lOuter};
Point(202) = {       0,    ROuter,  T/2,  lOuter};
Point(203) = { -ROuter,         0,  T/2,  lOuter};
Point(204) = {       0,   -ROuter,  T/2,  lOuter};

Circle (200) = { 201, 100, 202 };
Circle (201) = { 202, 100, 203 };
Circle (202) = { 203, 100, 204 };
Circle (203) = { 204, 100, 201 };

Line Loop(200) = { 200, 201, 202, 203};

Plane Surface(1) = { 200, 100 }; 

Extrude { 0, 0, -T } { Surface{1}; }
