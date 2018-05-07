###############################################################
# python code for analyzing a coplanar waveguide section in scuff-em
# Homer Reid 20180507
###############################################################
scuffExample="MicrostripDevices"
execfile('setupSCUFF.py')

###################################################
# create the solver, add metal traces, define ports,
# initialize substrate
###################################################
Solver=scuff.RFSolver()

Res = "Coarse" # or "Medium" or "Fine"
Solver.AddMetalTraceMesh("CPWCenter" + Res + ".msh");
Solver.AddMetalTraceMesh("CPWFlank"  + Res + ".msh", "DISPLACED -1.75 0 0");
Solver.AddMetalTraceMesh("CPWFlank"  + Res + ".msh", "DISPLACED +1.75 0 0");

XP  = 0.6
XM1 = 0.61
XM2 = 3.0
Y1  = 0.0;
Y2  = 10.0;
ZZ  = 0.0
Port1P  = [-XP,  Y1, ZZ,   +XP,  Y1, ZZ];
Port1MA = [+XM1, Y1, ZZ,   +XM2, Y1, ZZ];
Port1MB = [-XM2, Y1, ZZ,   -XM1, Y1, ZZ];
Solver.AddPortTerminal('+',Port1P);
Solver.AddPortTerminal('-',Port1MA);
Solver.AddPortTerminal('-',Port1MB);
Solver.AddPortTerminal('+',Port2P);
Solver.AddPortTerminal('-',Port2MA);
Solver.AddPortTerminal('-',Port2MB);

Solver.SetSubstratePermittivity(4.4)
Solver.SetSubstrateThickness(1.500)

Solver.PlotGeometry("CPW.scuffpy.pp");

###################################################
# open output files, write file preamble           
###################################################
ZParmFile = open(Geometry + ".Zparms", 'w')
SParmFile = open(Geometry + ".Sparms", 'w')

###################################################
# compute Z- and S-parameters over a range of frequencies
###################################################
Freqs = np.linspace(0.5,20,40) # frequencies at which to calculate (GHz)
for Freq in Freqs:

        print("Working at F={:.2g} GHz: ".format(Freq))
        tic = time.time();
	Solver.AssembleSystemMatrix(Freq)
        toc = time.time() - tic;
        print(" Matrix assembly: {.1f}s ".format(toc))
        tic = time.time();
        ZMatrix=Solver.GetZMatrix()
        print(" ZMatrix: {.1f}s ".format(toc))

        Z11 = ZMatrix.GetEntry(0,0);
        ZParmFile.write("%e %e %e \n" % (Freq,np.real(Z11),np.imag(Z11)))
        ZParmFile.flush()

        SMatrix=Solver.Z2S(ZMatrix);
        S11 = SMatrix.GetEntry(0,0);
        SParmFile.write("%e %e %e \n" % (Freq,np.real(S11),np.imag(S11)))
        SParmFile.flush()

ZParmFile.close()
SParmFile.close()
