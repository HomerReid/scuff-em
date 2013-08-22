R=5.0;  // disc radius
T=1.0;  // disc thickness
C=1.0;  // sidewall curvature radius

lf = R/4;  // fine   meshing fineness
lc = R/2;  // coarse meshing fineness
lm = (lf+lc)/2.0;

//////////////////////////////////////////////////
// curved sidewall and lower surface 
/////////////////////////////////////////////////
Point(200) = {       0,   0,   0, lf };
Point(201) = {       R,   0,   0, lf };
Point(202) = {       R,   0, T/2, lm };
Point(203) = {   R+C/2,   0, T/2, lm };
Point(204) = {       R,   0,   T, lc };
Point(205) = {       0,   0,   T, lc };

Line(200)    = { 200, 201 };
Ellipse(201) = { 201, 202, 203, 203 };
Ellipse(202) = { 203, 202, 203, 204 };
Line(203)    = { 204, 205 };

Extrude{ {0, 0, 1}, {0, 0, 0}, 2*Pi/3 } { Line{200,201,202,203}; }
Extrude{ {0, 0, 1}, {0, 0, 0}, 2*Pi/3 } { Line{204,207,211,215}; }
Extrude{ {0, 0, 1}, {0, 0, 0}, 2*Pi/3 } { Line{218,221,225,229}; }

Physical Surface(1) = { 206, 210, 214, 217, 
                        220, 224, 228, 231,  
                        234, 238, 242  234};
