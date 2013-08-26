// GMSH geometry file for a finite-thickness disc with 
// radial sinusoidal corrugations on the upper surface
// homer reid 8/2013

R=10.0;  // disc radius
T=3.0;   // disc thickness
C=T;     // sidewall curvature radius

//H=0.75;  // corrugation amplitude
H=0.00;  // corrugation amplitude
SF=3;    // "spatial frequency" 
         //  (SF = number of ripples between 0 and R)

lf = 10/5;  // fine   meshing fineness
lc = 10/2;  // coarse meshing fineness
lm = (lf+lc)/2.0;

//////////////////////////////////////////////////
// corrugated upper surface
//////////////////////////////////////////////////
Delta=1/10;
Point(100) = {         0,   0,   0,                             0.5*lf };
Point(101) = { 1*Delta*R,   0,   H*(Cos(2*Pi*SF*1*Delta) - 1),  lf };
Point(102) = { 2*Delta*R,   0,   H*(Cos(2*Pi*SF*2*Delta) - 1),  lf };
Point(103) = { 3*Delta*R,   0,   H*(Cos(2*Pi*SF*3*Delta) - 1),  lf };
Point(104) = { 4*Delta*R,   0,   H*(Cos(2*Pi*SF*4*Delta) - 1),  lf };
Point(105) = { 5*Delta*R,   0,   H*(Cos(2*Pi*SF*5*Delta) - 1),  lf };
Point(106) = { 6*Delta*R,   0,   H*(Cos(2*Pi*SF*6*Delta) - 1),  lf };
Point(107) = { 7*Delta*R,   0,   H*(Cos(2*Pi*SF*7*Delta) - 1),  lf };
Point(108) = { 8*Delta*R,   0,   H*(Cos(2*Pi*SF*8*Delta) - 1),  lf };
Point(109) = { 9*Delta*R,   0,   H*(Cos(2*Pi*SF*9*Delta) - 1),  lf };
Point(110) = { 10*Delta*R,  0,   H*(Cos(2*Pi*SF*10*Delta) - 1),  lf };

Spline(100) = { 100, 101, 102, 103, 104, 105, 
                106, 107, 108, 109, 110 };

Extrude{ {0, 0, 1}, {0, 0, 0}, 2*Pi/3 } { Line{100}; }
Extrude{ {0, 0, 1}, {0, 0, 0}, 2*Pi/3 } { Line{101}; }
Extrude{ {0, 0, 1}, {0, 0, 0}, 2*Pi/3 } { Line{104}; }

//////////////////////////////////////////////////
// curved sidewall and lower surface 
/////////////////////////////////////////////////
Point(200) = {       0,   0,   -T, lc };
Point(201) = {       R,   0,   -T, lc };
Point(202) = {       R,   0, -T/2, lm };
Point(203) = {   R+C/2,   0, -T/2, lm };
Point(204) = {       R,   0,    0, lf };
Point(205) = {       0,   0,    0, lf };

Line(200)    = { 200, 201 };
Ellipse(201) = { 201, 202, 203, 203 };
Ellipse(202) = { 203, 202, 203, 204 };
Line(203)    = { 204, 205 };

Extrude{ {0, 0, 1}, {0, 0, 0}, 2*Pi/3 } { Line{200,201,202,203}; }
Extrude{ {0, 0, 1}, {0, 0, 0}, 2*Pi/3 } { Line{204,207,211,215}; }
Extrude{ {0, 0, 1}, {0, 0, 0}, 2*Pi/3 } { Line{218,221,225,229}; }

//////////////////////////////////////////////////
// list of all surfaces that we want to be meshed
//////////////////////////////////////////////////
Physical Surface(1) = { 103, 106, 109, 
                        206, 210, 214,
                        220, 224, 228,
                        234, 238, 242 };
