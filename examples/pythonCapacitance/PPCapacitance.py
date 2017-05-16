##################################################
# python code for studying finite-size capacitors in scuff-em
# Homer Reid 20170515
##################################################
import os;
import numpy;
from scipy.optimize import curve_fit;
import subprocess;
import scuff;

###################################################
# set some environment variables so that SCUFF-EM knows where
# to look for input files
# (this assumes we are running from the share/examples/pythonCapacitance
#  directory of the scuff-em installation, so that e.g.
#  geoFiles and mshFiles are subdirectories of the current
#  working directory; this assumption is also made by the
#  gmsh shell commands below)
###################################################
os.environ["SCUFF_MESH_PATH"]="mshFiles"
os.environ["SCUFF_GEO_PATH"]="scuffgeoFiles"
os.environ["SCUFF_SUBSTRATE_PATH"]="substrateFiles"

Fineness = 3      # meshing fineness (triangle edges per unit length)

###################################################
# loop over square side lengths L
###################################################
LMin     = 10
LMax     = 20
LPoints  = 11
LVector=[]
CPUAVector=[]   # 'capacitance per unit area'
DataFile = open('PPCapacitor.CvsL','w')
DataFile.truncate();
for L in numpy.linspace(LMin, LMax, LPoints).tolist():
#
    #--------------------------------------------------
    # run gmsh to generate mesh file for square of side L
    #--------------------------------------------------
    subprocess.call(['gmsh', '-2', 'geoFiles/Rectangle.geo',
                     '-setnumber', 'LX', str(L),
                     '-setnumber', 'LY', str(L),
                     '-setnumber', 'N',  str(Fineness),
                     '-o', 'mshFiles/PPCapacitor.msh'])
#
    #--------------------------------------------------
    # use scuff-em to compute capacitance
    #--------------------------------------------------
    print "Computing capacitance at L=", format(L)
    Solver=scuff.SSSolver("PPCapacitor.scuffgeo", "SiliconGP.substrate");
    CMatrix=Solver.GetCapacitanceMatrix()
    CPUA=CMatrix[0,0] / (L*L)
    LVector.append(L)
    CPUAVector.append(CPUA)
    DataFile.write('%-10s %-10s\n' % (format(L), format(CPUA)));
    DataFile.flush()

DataFile.close()

###################################################
# fit CPUA versus L data to the form
#  C(L) = CInfinity + Beta/L^2
###################################################
def FunctionalForm(x, CInf, Beta):
    return CInf + Beta/(x*x)

CInfBeta = curve_fit(FunctionalForm, LVector, CPUAVector)[0]
CInf=CInfBeta[0]

print "\n*\n*\n"
print "Capacitance per unit area, extracted to L=infinity limit = ", format(CInf)
