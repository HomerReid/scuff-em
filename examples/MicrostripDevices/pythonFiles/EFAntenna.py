###############################################################
# python code for analyzing a microstrip patch antenna in scuff-em
# Homer Reid 20180422
###############################################################
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
SParmFile = open("EFAntenna_L8_" + Res + ".Sparms", 'w')

###################################################
# compute Z- and S-parameters over a range of frequencies
###################################################
Freqs = np.linspace(2.5,20,36) # frequencies at which to calculate (GHz)
for Freq in Freqs:

        print("Working at F={:.2g} GHz: ".format(Freq))
        tic = time.time();
	Solver.AssembleSystemMatrix(Freq)
        toc = time.time() - tic;
        print(" Matrix assembly: {:1f}s ".format(toc))
        tic = time.time();
        ZMatrix=Solver.GetZMatrix()
        print(" ZMatrix: {:1f}s ".format(toc))

        Z11 = ZMatrix.GetEntry(0,0);
        ZParmFile.write("%e %e %e \n" % (Freq,np.real(Z11),np.imag(Z11)))
        ZParmFile.flush()

        SMatrix = Solver.Z2S(ZMatrix);
        S11 = SMatrix.GetEntry(0,0);
        SParmFile.write("%e %e %e \n" % (Freq,np.real(S11),np.imag(S11)))
        SParmFile.flush()

ZParmFile.close()
SParmFile.close()
