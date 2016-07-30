# Double-dipole detection of unexploded ordnance in weakly conductive soil

This example demonstrates the use of <span class=SC>scuff-em</span>
to model the real-world problem of detecting unexploded ordnance (UXO),
buried below ground in conducting soil layers, by electromagnetic induction
with above-ground excitation and probe antennas.

## Toy problem: Dipole source above permeable sphere in vacuum

Although we will ultimately be interested in objects of
complex shapes buried in conducting soil layers, let's begin
by warming up with a toy problem: using
<span class=SC>scuff-em</span> to reproduce
known analytical results for the response of
of a conducting, permeable spherical particle in vacuum
to a dipole excitation. 
More specifically, we'll compute
[scattering dyadic Green's functions][LDOSMemo] (DGFs)
$\boldsymbol\mathcal{G}( \mathbf{x}_{\text{\tiny{D}}}
                        ,\mathbf{x}_{\text{\tiny{S}}})$
of the sphere---that is, the scattered fields at
destination point $\mathbf{x}_{\text{\tiny{D}}}$
due to excitation by a point dipole source at 
source point $\mathbf{x}_{\text{\tiny{S}}}$---for
various pairs of source and destination points in the
vicinity of the sphere.

### Geometry description

#### Length and frequency units

As discussed [here][UnitsInSCUFFEM], the default units for
length and angular frequency in <span class=SC>scuff-em</span>
are typically $L_0=1\,\mu$m and
$\omega_0=\frac{c}{L_0}=3\cdot 10^{14}$ rad/sec.
However, in view of the the typical lengthscales and excitation
frequencies that appear in UXO problems, it's convenient 
here to use $L_0=1$ m and $\omega_0=\frac{c}{L_0}=3\cdot 10^8$ rad/sec.
As discussed in the FAQ item referenced above, the only place we will
ever need explicitly to recognize that we have made this choice is 
in defining frequency-dependent material-property functions (see below).

#### Surface meshes

We will consider spheres of radius $R=100$ cm.
Here's a [<span class=SC>gmsh</span>][GMSH] geometry file
describing a sphere with an adjustable parameter `R` that
may be set on the <span class=SC>gmsh</span> command line 
to define the radius: [`Sphere_R.scuffgeo`](Sphere_R.scuffgeo).
To define surface meshes with radius 0.1 with two different
meshing finenesses, go like this:

````bash
 % gmsh -2 -clscale 1.0 -setnumber R 0.1 Sphere_R.geo -o Sphere_R0P1.msh
 % RenameMesh Sphere_R0P1.msh
 % gmsh -2 -clscale 0.5 -setnumber R 0.1 Sphere_R.geo -o Sphere_R0P1.msh
 % RenameMesh Sphere_R0P1.msh
````

(Here [`RenameMesh`][RenameMesh] is a little script that uses 
[<span class=SC>scuff-analyze</span>][scuff-analyze] to count
the number of interior edges in a mesh and rename the file
accordingly.) This produces mesh files
[`Sphere_R0P1_504.msh`](Sphere_R0P1_504.msh)
and
[`Sphere_R0P1_1485.msh`](Sphere_R0P1_1485.msh).

#### Material properties 

We will consider three constitutent materials for the
sphere:

+ PEC (perfect electric conductor)

+ copper (DC conductivity $\sigma = 6\cdot 10^7 \mho$/m)

+ iron (DC conductivity $\sigma = 1\cdot 10^7 \mho$/m, relative
        DC permeability $\mu$=5000)

Since we will be interested only in low-frequency behavior, we
characterize the electric response of metals entirely in terms
of the DC conductivity $\sigma$, corresponding to a frequency-dependent
relative dielectric function
%====================================================================%
$$ \epsilon^{\text{\scriptsize{rel}}}(\omega)
   = 1 + i\{sigma}{\epsilon_0 \omega}
$$
%====================================================================%
Recalling that the permittivity of vacuum is
$\epsilon_0=\frac{1}{Z_0 c}$ where $c$ is the vacuum speed of light
and $Z_0\approx 377 \Omega$ is the impedance of vacuum,
the dimensionless quantity in the second term may be written
%====================================================================%
$$  \frac{\sigma}{\epsilon_0 \omega}
   =\frac{\sigma Z_0 c }{\omega}
    _{ \frac{377 \cdot \left(\sigma \text{ in }\mho \text{ m}^{-1}\right)}
            {\left(\omega \text{ in }\mho 3\cdot 10^{8} rad/sec\right)}
     }
   =\frac{\sigma Z_0 c }{\omega}
    _{ \frac{3.77 \cdot 10^{-4} \cdot \left(\sigma \text{ in }\mho \text{ m}^{-1}\right)}
            {\left(\omega \text{ in }\mho 3\cdot 10^{14} rad/sec\right)}
     }
$$
%====================================================================%

#### [[scuff-em]] geometry files

Here are 
[<span class=SC>scuff-em</span> geometry files][Geometries] for our
PEC, copper, and iron spheres with the coarser meshing. Note
that I displace the sphere downward by one radius, so that its
north pole is at $z=0;$ this means that the $z$ coordinates
of evaluation points above the sphere
agree with the point-surface distance $h$.

+ [`PECSphere_R0P1_504.scuffgeo`](PECSphere_R0P1_504.scuffgeo)

````bash
OBJECT Sphere
	MESHFILE Sphere_R0P1_504.scuffgeo
	DISPLACED 0 0 -0.1
ENDOBJECT	

+ [`CopperSphere_R0P1_504.scuffgeo`](CopperSphere_R0P1_504.scuffgeo])

````bash
MATERIAL Copper
 # 3.77e-4 * 6e7 = 2.26e4
 SIGMA  = 2.26e4 
 Eps(w) = 1.0 + i*SIGMA/w
ENDMATERIAL

OBJECT Sphere
	MESHFILE Sphere_R0P1_504.scuffgeo
	MATERIAL Copper
	DISPLACED 0 0 -0.1
ENDOBJECT	
````bash

+ [`IronSphere_R0P1_504.scuffgeo`](IronSphere_R0P1_504.scuffgeo])

````bash
MATERIAL Iron
 # 3.77e-4 * 1e7 = 3.77e3
 SIGMA  = 3.77e3 
 Eps(w) = 1.0 + i*SIGMA/w
 Mu(w)  = 5000
ENDMATERIAL

OBJECT Sphere
	MESHFILE Sphere_R0P1_504.scuffgeo
	MATERIAL Iron
	DISPLACED 0 0 -0.1
ENDOBJECT	
````bash

### Two ways to proceed

As it happens, this problem can actually be attacked in two
distinct but equivalent ways using two distinct application
codes in the <span class=SC>scuff-em</span> suite:
[<span class=SC>scuff-scatter</span>][scuff-scatter]
or
[<span class=SC>scuff-ldos</span>][scuff-ldos].
(Of course a third way would be to write a C++ or Python
program calling the 
[<span class=SC>scuff-em</span> API][libscuff], but we
won't consider that here.) The basic differences 

+ <span class=SC>scuff-scatter</span> is a general-purpose
full-featured scattering code that allows the scattering problem
to be customized in all sorts of ways---arbitrary
[incident fields][IncidentFields] (not just fields of dipole sources),
[geometric transformation][Transformations] among the 
scattering bodies, etc.---and also offers more types of outputs
(induced moments, absorbed/scattered power and force/torque,
visualization diagrams, etc.) beyond just the scattered field.

+ <span class=SC>scuff-ldos</span> is more narrowly
targeted at the specific problem of computing dyadic Green's
functions (scattered fields due to pointlike dipole excitation),
and it uses a specialized computational approach that is 
particularly efficient for this particular problem, but 
which does not extend to tweaked versions of the problem
(for example, it wouldn't work if we wanted to excite the
sphere with plane waves instead of point sources.)
    One advantage of <span class=SC>scuff-ldos</span> over 
<span class=SC>scuff-scatter</span> that will be relevant for
our purposes here is that <span class=SC>scuff-ldos</span>
comes with built-in support for 
[Brillouin-zone integration][BZI], which we will eventually
need to compute the response of an infinitely periodic 
system to point-source excitation. 

For purposes of illustration we will show how our toy
problem may be set up and solved in both codes.

### Solution using <span class=SC>scuff-scatter</span>

In <span class=SC>scuff-scatter</span>, each calculation
of the DGFs for a given source point $\mathbf{x}_{\text{\tiny{S}}}$
corresponds to a different incident field for a scattering
problem, so in this case our first task is to write
a [list of incident fields][IFList]. For each of 21
logarithmically-spaced values of $h$ in the range
$0.1R \le h \le 10R$, I will define four incident fields,
corresponding to $x$- and $z$-directed electric and
magnetic dipoles located at $(0,0,h)$. This gives 
a total of 84 incident fields, which are specified
one per line in a file named `IFFile`:

````bash

````



### Solution using <span class=SC>scuff-ldos</span>

In <span class=SC>scuff-ldos</span>, it is understood 
that we will be computing the scattered fields 
due to point sources, so the only things we need to 
specify are the locations of the point sources
and (if they are different) the locations of the
evaluation points. Thus the separate `IFFile` and `EPFile`
inputs that we prepared for [[scuff-scatter]] are 
replaced by a single `EPFile` input---but which may 
now have 
