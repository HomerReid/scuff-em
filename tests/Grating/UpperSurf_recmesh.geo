/////////////////////////////////////////////////
//
// GMSH geometry description for the unit cell of
// a rectangular periodic grating,
// whose period is in x direction
// and translational invariant in y direction.
//
// Length Unit: um
//
// Jun Wang @LevLab 2020
/////////////////////////////////////////////////

// In order to be compatible with SCUFF-EM

Mesh.MshFileVersion = 2.2; 

Geometry.CopyMeshingMethod = 1;

//////////////////////////////////////////////////
// user-specified parameters
////////////////////////////////////////////////

// Depth of the trench
DefineConstant[ a  = 1.000 ];
// Period of the grating, along x
DefineConstant[ Px = 2.000 ];
// Width of the plateau
DefineConstant[ w  = 1.000 ];

// Mesh resolution
DefineConstant[ l  = 0.15 ];
// Number of period in x
DefineConstant[ NpX = 1 ];
// Dimension in y
DefineConstant[ Ly = 10.00 ];


//////////////////////////////////////////////////
// end of user-specified parameters
////////////////////////////////////////////////

////////////////////////////////////////////////
// Derived params
////////////////////////////////////////////////

X2 = 0.5*w;
X1 = 0.5*Px;

Y1 = -0.5*Ly;

Ny = Ceil(Ly/l*0.5);
Nz = Ceil(a/l);
Nx_p = Ceil(w/l)+1;
Nx_g = Ceil(0.5*(Px-w)/l)+1;

////////////////////////////////////////////////
// End of Derived params
////////////////////////////////////////////////

Point(1) = {-X1, Y1, -a};
Point(2) = {-X2, Y1, -a};
Point(3) = {-X2, Y1, 0};

Point(4) = {X2, Y1, 0};
Point(5) = {X2, Y1, -a};
Point(6) = {X1, Y1, -a};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};

Transfinite Line{1} = Nx_g;
Transfinite Line{2} = Nz;
Transfinite Line{3} = Nx_p;
Transfinite Line{4} = Nz;
Transfinite Line{5} = Nx_g;

For j In {1:Floor(NpX/2)}
	Translate {Px*j, 0, 0} {
	  Duplicata { Line{1}; Line{2}; Line{3}; Line{4}; Line{5}; }
	}
	Translate {-Px*j, 0, 0} {
	  Duplicata { Line{1}; Line{2}; Line{3}; Line{4}; Line{5}; }
	}
EndFor

Extrude {0, Ly, 0} {
  Line{:}; Layers{Ny};
}
