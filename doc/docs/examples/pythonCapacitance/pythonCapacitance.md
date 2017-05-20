# Studying finite-size effects in metal-on-substrate capacitors using data fitting in python

In this example, we use the [[python]] interface to
<span class=SC>scuff-em</span>---specifically, to the
[<span class=SC>scuff-em</span> electrostatics module][scuffStatic]---to
study finite-size effects in capacitors formed by metal traces on 
(infinite-area) dielectric substrates with and without ground planes.

Our calculation will exploit [[scuff-em]]'s capability for
[implicit treatment of layered dielectric substrates][ImplicitSubstrate],
eliminating the need to mesh substrate layers or ground planes
and greatly reducing computational cost. 

Our goal will be to estimate the capacitance per unit length or area
for infinite-size capacitors of various shapes.
We will do this by using the 
<span class=SC>python</span> interface to [[scuff-em]]
to fit numerical data for several finite values of
the length $L$ or area $A$ to a functional form in $L$ or $A$, then 
taking the limits $L,A\to \infty$.

### Layout of input files

The input files for the calculations discussed here may be
found in the `share/scuff-em/examples/pythonCapacitance` subdirectory
of your <span class=SC>scuff-em</span> installation.
The various files here are organized into the following subdirectories:

+ `geoFiles`: [<span class=SC>gmsh</span>][gmsh] geometry files
+ `mshFiles`: <span class=SC>gmsh</span> mesh files
+ `scuffgeoFiles`: [[scuff-em]] geometry files
+ `substrateFiles`: [<span class=SC>scuff-em</span> substrate-definition files][SubstrateFile]

### Customizable <span class=SC>gmsh</span> geometry for rectangular plates

The `geoFiles` subdirectory contains a [[gmsh]] geometry file named
[`Rectangle.geo`](Rectangle.geo) which we will use to produce
surface meshes describing infinitesimally thin metallic traces of
various thicknesses. This file contains adjustable parameters
to allow the dimensions and meshing fineness to be specified on the
[[gmsh]] command line. Here are two examples:

```bash
# square of side length 10 with 4 triangles per unit area:

% gmsh -2 Rectangle.geo -o Square.msh \
       -setnumber Lx 10 -setnumber Ly 10 -setnumber N 4

# center conductor for coplanar transmission line: width 1, length 10

% gmsh -2 Rectangle.geo -o Center.msh \
       -setnumber Lx 1 -setnumber Ly 10 -setnumber N 4
```

### Implicit substrate definition files

The geometries we consider will consist of one or more metal
plates lying on the upper surface of an infinite-area dielectric
substrate, possibly with a perfectly-conducting ground
plane underneath.
The substrate, which will be handled implicitly by [[scuff-em]]
using the method described [in this memo][SubstrateMemo]
and demonstrated [in this example][ImplicitSubstrate],
is described by a [simple text file][SubstrateFile].
We will consider two substrate files: one describing a
freestanding finite-thickness layer of silicon 
(relative permittivity $\epsilon=11.7$):
```bash
0.0   CONST_EPS_11.7
-1.0  VACUUM
```
and one describing the same layer but now with a ground plane 
underneath:
```bash
0.0   CONST_EPS_11.7
-1.0  GROUNDPLANE
```
These two files are named `Silicon.substrate` and `SiliconGP.substrate`
and live in the `substrateFiles` subdirectory of the `pythonCapacitance` example folder.

## Parallel-plate capacitor

Our first calculation will be for a parallel-plate capacitor consisting
of a metal square of side length $L$ and area $A=L^2$ on a dielectric
substrate of height $h$ and relative permittivity $\epsilon$
lying atop a perfectly-conducting ground plane.
Neglecting finite-size effects, the capacitance
per unit area (in units of $1/\epsilon_0$, the permittivity of free space)
should be
$$ \lim_{A\to \infty} \frac{C(L)}{\epsilon_0 L^2} = \frac{\epsilon_r}{h}.$$

For the substrate described by the file `SiliconGP.substrate` 
($h=1, \epsilon_r=11.7$) we expect to find the numerical value
\[\lim_{L\to \infty} \frac{C(L)}{\epsilon_0 L^2} = 11.7 \qquad\textbf{(1)}\]

However, the finite size of the upper plate will cause calculated
values to deviate from this prediction by an amount which (for fixed $h$)
[we expect to scale asymptotically like $\frac{1}{L}$:][ParallelPlateCapacitorPaper1937]
\[ \text{for finite } L: \,\,
   \frac{C}{\epsilon_0 L^2} = \frac{\epsilon_r}{h} + \frac{\beta}{L}
   +\text{higher-order terms}
   \qquad \textbf{(2)}
\]
where $\beta$ is a constant.
Thus we will use [[scuff-em]] to compute the capacitance per
unit area for various finite values of $L$, fit these data
to the functional form $C_0+\beta/L$, and identify the
constant $C_0$ in the fit as the $L\to\infty$ value of 
the capacitance per unit area.

Here's the [[python]] code named [`PPCapacitance.py`](PPCapacitance.py)
that does this. For a set of
$L$ values, this script **(1)** invokes [[gmsh]] to produce
a surface mesh for an $L\times L$ square, **(2)** calls routines
in the [[scuff-em]] [[python]] module to get the capacitance 
matrix for the resulting geometry. (Since the geometry only 
has one conductor, the capacitance matrix is a $1\times 1$ matrix.) 
Then we use the `curve_fit`
routine provided by [[scipy]] to fit the data to the functional
form **(2)** and extract our estimate of the capacitance per
unit area in the $L\to\infty$ limit.

```python
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
#  C(L) = CInfinity + Beta/L
###################################################
def FunctionalForm(L, CInf, Beta):
    return CInf + Beta/L

CInfBeta = curve_fit(FunctionalForm, LVector, CPUAVector)[0]
CInf=CInfBeta[0]

print "\n*\n*\n"
print "Capacitance per unit area, extracted to L=infinity limit = ", format(CInf)

```

## Results

Running the python script from the command line produces, eventually,
the following output:

```bash
% python PPCapacitor.py
...
...
...
Capacitance per unit area, extracted to L=infinity limit =  11.764354485
```

Comparing against equation **(1)**, we see that we recover the
correct theoretical value to 3 decimal places.

Also produced is a file named `PPCapacitor.CvsL`, which
tabulates the finite-$L$ values of the capacitance. Here's
a plot of these data, together with the fitting function
and the $L\to\infty$ extrapolation:

![Plot of capacitance data and functional fit](PPCapacitorData.png)

(Here's the [<span class=SC>gnuplot</span> script](Plotter.PPCapacitance)
that produces this plot).

Take-home messages:

+ The $L\to\infty$ extrapolation recovers the correct infinite-area 
result to an accuracy of better than 1%.

+ *Without* this extrapolation, we would incur much more severe
errors. For example, if we were to approximate the $L\to\infty$
capacitance by the capacitance computed for the largest
value of $L$ considered here ($L=20$), we would incur an error
of 9.8%, some 20&times; greater than the error in the infinite-$L$
extrapolation.


--------------------------------------------------

[GMSH]:                 http://www.geuz.org/gmsh
[scuffStatic]:          ../../applications/scuff-static/scuff-static.md
[ImplicitSubstrate]:    ../ImplicitSubstrate/ImplicitSubstrate.md
[scuffInstallation]:    ../../reference/Installing
[scuffMaterials]:       ../../reference/Materials
[SubstrateMemo]:	../../tex/StaticDielectricSubstrate.pdf
[SubstrateFile]:	../ImplicitSubstrate/ImplicitSubstrate.md#SubstrateFile
[ParallelPlateCapacitorPaper1937]:	http://ieeexplore.ieee.org/document/6540485/
