R=5.0;   // disc radius
T=5.0;   // disc thickness
C=1.0;   // sidewall curvature minor axis (=0 for flat side walls)

lc = R/2;

Point(100) = {       0,   0,    0, lc };
Point(101) = {       R,   0,    0, lc };
Point(102) = {       0,   R,    0, lc };
Point(103) = {      -R,   0,    0, lc };
Point(104) = {       0,  -R,    0, lc };

Circle(101)    = { 101, 100, 102 };
Circle(102)    = { 102, 100, 103 };
Circle(103)    = { 103, 100, 104 };
Circle(104)    = { 104, 100, 101 };

Extrude{0,0,-T} { Line{101, 102, 103, 104}; }

Line Loop(1) = { 101, 102, 103, 104 };
Line Loop(2) = { 105, 109, 113, 117 };

Ruled Surface(1) = { 1 };
Ruled Surface(2) = { 2 };

Field[1] = Attractor;
Field[1].NodesList = {100};

// We then define a Threshold field, which uses the return value of
// the Attractor Field[1] in order to define a simple change in
// element size around the attractors (i.e., around point 5 and line
// 1)
//
// LcMax -                         /------------------
//                               /
//                             /
//                           /
// LcMin -o----------------/
//        |                |       |
//     Attractor       DistMin   DistMax
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc / 10;
Field[2].LcMax = lc;
Field[2].DistMin = 0.5;
Field[2].DistMax = 5;

Background Field=2;
