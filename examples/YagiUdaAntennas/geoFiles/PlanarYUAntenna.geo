////////////////////////////////////////////////////
// GMSH geometry file for a planar Yagi-Uda antenna
//
// usage:
//  gmsh -2 PlanarYUAntenna.geo [options]
//
// Available options:
//
//  -setnumber FCenter    2.45 (center frequency in GHz)
//
//  -setnumber LDrive     0.49 (length of driven element / wavelength)
//  -setnumber LDirector  0.4  (length of director / wavelength)
//  -setnumber LReflector 0.6  (length of reflector / wavelength)
//
//  -setnumber WDrive     0.1  (width of driven element / wavelength)
//  -setnumber WDirector  0.1  (width of director / wavelength)
//  -setnumber WReflector 0.1  (width of reflector / wavelength)
//
//  -setnumber LFeed      0.05 (length of feed gap / wavelength)
//
//  -setnumber N          10   (meshing fineness: triangle edges per wavelength)
//
//  Homer Reid 9/2017
//
////////////////////////////////////////////////////

// user-defined parameters
DefineConstant[FCenter    = 2.45 ];

DefineConstant[LDrive     = 0.49 ];
DefineConstant[LDirector  = 0.40 ];
DefineConstant[LReflector = 0.60 ];

DefineConstant[WDrive     = 0.05 ];
DefineConstant[WDirector  = 0.05 ];
DefineConstant[WReflector = 0.05 ];

DefineConstant[LFeed      = 0.01 ];

DefineConstant[N          = 30   ]/;

//////////////////////////////////////////////////
// derived parameters
//////////////////////////////////////////////////
Lambda = 3.0e2/FCenter;         // wavelength in mm 

NLDrive     = Ceil[N*LDrive/2];
NLDirector  = Ceil[N*LDirector];
NLReflector = Ceil[N*LReflector];

NWDrive     = Ceil[N*WDrive];
NWDirector  = Ceil[N*WDirector];
NWReflector = Ceil[N*WReflector];

LDrive     = LDrive     * Lambda;
LDirector  = LDirector  * Lambda;
LReflector = LReflector * Lambda;
LFeed      = LFeed      * Lambda;
   
WDrive     = WDrive     * Lambda;
WDirector  = WDirector  * Lambda;
WReflector = WReflector * Lambda;

//////////////////////////////////////////////////
// driven element
//////////////////////////////////////////////////
Point(101) = {0.0, -0.5*WDrive, LFeed};
Point(102) = {0.0,  0.5*WDrive, LFeed};
Line(101)  = {101, 102};
Transfinite Line{101} = NWDrive + 1;

Point(111) = {0.0, -0.5*WDrive, -LFeed};
Point(112) = {0.0,  0.5*WDrive, -LFeed};
Line(111)  = {111, 112};
Transfinite Line{111} = NWDrive + 1;

Extrude { 0.0, 0.0,  0.5*LDrive} { Line{101}; Layers{NLDrive}; }
Extrude { 0.0, 0.0, -0.5*LDrive} { Line{111}; Layers{NLDrive}; }

Physical Surface(1) = {115, 119};

//////////////////////////////////////////////////
// reflector
//////////////////////////////////////////////////
Point(201) = {0.0, -0.5*WReflector, -0.5*LReflector};
Point(202) = {0.0, +0.5*WReflector, -0.5*LReflector};
Line(201)  = {201, 202};
Transfinite Line{201} = NWReflector + 1;

Extrude { 0.0, 0.0, LReflector} { Line{201}; Layers{NLReflector}; }

Physical Surface(2) = {205};

//////////////////////////////////////////////////
// director
//////////////////////////////////////////////////
Point(301) = {0.0, -0.5*WDirector, -0.5*LDirector};
Point(302) = {0.0, +0.5*WDirector, -0.5*LDirector};
Line(301)  = {301, 302};
Transfinite Line{301} = NWDirector + 1;

Extrude { 0.0, 0.0, LDirector} { Line{301}; Layers{NLDirector}; }

Physical Surface(3) = {305};
