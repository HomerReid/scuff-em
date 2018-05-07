###############################################################
# python code for analyzing a microstrip patch antenna in scuff-em
# Homer Reid 20180422
###############################################################
from __future__ import print_function;

scuffExample="MicrostripDevices"
execfile('setupSCUFF.py')

# set meshing resolution ("Coarse" or "Fine")
Res = "Coarse"

###################################################
# create the solver, add metal traces, define ports,
# initialize substrate
###################################################
Solver=scuff.RFSolver()

Solver.AddMetalTraceMesh("EFAntenna_L8_" + Res + ".msh");

Solver.AddPort([-5, 0, 0, 5, 0, 0])

Solver.SetSubstratePermittivity(2.2)
Solver.SetSubstrateThickness(0.794)

Solver.PlotGeometry("EFAntenna_L8" + Res + ".pp");

###################################################
# open output files, write file preamble
###################################################
ZParmFile = open("EFAntenna_L8_" + Res + ".Zparms", 'w')
ZParmFile.write("# Data columns:\n");
ZParmFile.write("# 1   frequency (GHz)\n");
ZParmFile.write("# 2,3 re,im Z11 (Ohms)\n");

SParmFile = open("EFAntenna_L8_" + Res + ".Sparms", 'w')
SParmFile.write("# Data columns:\n");
SParmFile.write("# 1   frequency (GHz)\n");
SParmFile.write("# 2,3 re,im S11       \n");

###################################################
# compute Z- and S-parameters over a range of frequencies
###################################################
Freqs = np.linspace(2.5,20,36) # frequencies at which to calculate (GHz)
for Freq in Freqs:

        print("F={:.1f} GHz: ".format(Freq),end='')

        t0  = time.time();
	Solver.AssembleSystemMatrix(Freq)        # assemble SIE system
        t1  = time.time()

        ZMatrix = Solver.GetZMatrix();           # compute Z matrix
        SMatrix = Solver.Z2S(ZMatrix);           # convert to S matrix
        Z11 = ZMatrix.GetEntry(0,0);
        S11 = SMatrix.GetEntry(0,0);
        t2  = time.time()

        print(" S11={:+.2f},{:+.2f}   System: {:.1f}s  ZMat: {:.1f}s".format(S11.real,S11.imag,t1-t0,t2-t1));

        ZParmFile.write("%e %e %e \n" % (Freq,Z11.real,Z11.imag))
        ZParmFile.flush()

        SParmFile.write("%e %e %e \n" % (Freq,S11.real,S11.imag))
        SParmFile.flush()

ZParmFile.close()
SParmFile.close()
