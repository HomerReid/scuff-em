<h1>Photonic LDOS above a dielectric half-space</h1>

In this example, we exploit [[scuff-em]]'s
[support for 2D periodic geometries][ExtendedGeometries]
by using [<span class="SC">scuff-ldos</span>][scuff-ldos]
to compute the electromagnetic local density
of states (LDOS) at evaluation points lying above an
infinite planar dielectric interface.

Because this geometry happens to be amenable to an analytical 
treatment, we will also *check* the numerical results of
[[scuff-em]] against the results of analytical calculations of
the LDOS.

The files for this example may be found in the
`share/scuff-em/examples/HalfSpaceLDOS` subdirectory
of your [[scuff-em]] installation.

Also, the computational procedure implemented by [scuff-ldos]
is described in this memo:
[Computation of Green's Functions and LDOS in <span class="SC">scuff-em</span>][LDOSMemo].

### Slight complication: The need for Brillouin-zone integration in periodic geometries

In general, the quantity in which we will be interested
is the LDOS at a single point $\bf x$ in space---basically,
the scattered field at $\bf x$ due to a point dipole
source at $\bf x$. For a compact geometry such as a 
dielectric nanoparticle, this requires just a single 
scattering calculation (technically 6 calculations
with different incident fields, but calculations with
additional incident fields come basically for free
once we've assembled and factorized the BEM matrix).

However, the situation is slightly more complicated for
periodic geometries. To compute the LDOS at a point
in the vicinity of (say) an infinitely extended
(periodically replicated) surface, we need to solve a
scattering problem in which the incident field is not 
Bloch-periodic, but is instead the field of a localized 
point source, which certainly does not obey
Bloch-periodic boundary conditions.
However, when working with periodic geometries in [[scuff-em]],
each individual scattering problem has the property
that *everything* in the problem---including the incident
field---must be Bloch-periodic with a specific Bloch vector
$\bf k_{\hbox{B}}$. Thus, when we do a single scattering
problem with a point-source incident field in a 
periodic [[scuff-em]] geometry at Bloch vector
$\bf k_{\hbox{B}}$, the incident field is really
the field of an infinite phased *array* of point sources, 
with the source in lattice cell $\bf L$ phased
by the Bloch factor $e^{i\bf k_B\cdot \bf L}.$
This is not the incident field we want to consider 
when computing LDOS.

Of course, it is still possible to get what we want---the
response of our geometry to a non-periodic point source---but
it requires summing the results of many scattering problems,
with each problem phased in such a way as to cancel out
the effects of all source in all lattice cells except the 
particular cell we want. In other words, to get the 
LDOS at a single point in a periodic geometry
requires a *Brillouin-zone integration*, in which
we perform and sum the results of many scattering calculations
at various Bloch vectors $\bf k_B$ lying in the 
Brillouin zone of the reciprocal lattice.

There are two ways to do calculations like this in [[scuff-em]]:

+ You can design your own cubature scheme for integration over
the Brillouin zone, and simply ask [[scuff-ldos]] to give you
values of the integrand at specific Bloch vectors $\bf k_B$.
In this case you will pass the ``--OmegakBlochFile`` command-line
argument to [[scuff-ldos]]

+ Alternatively, you can use the built-in Brillouin-zone integrator
in [[scuff-em]] to perform the BZ integral automatically,
resulting in reports

The general topic of Brillouin-zone integration in [[scuff-em]] codes is
discussed in more detail 
on the page 
[Brillouin-zone integration in <span class="SC">scuff-em</span>](../../reference/BrillouinZoneIntegration.md.)



--------------------------------------------------

# [[gmsh]] geometry file for unit-cell geometry 

The [[gmsh]] geometry file [`Square_N.geo`](Square_N.geo)
describes the portion of the half-space interface
that lies within the *unit cell,*
i.e. the cell that is infinitely periodically
replicated to yield the full geometry.

I call this file `Square_N.geo` to remind myself that 
it contains a parameter `N` that describes the meshing 
fineness; more specifically, `N` defines the number of 
segments per unit length. (The file also contains 
a parameter named `L` that defines the unit-cell
length in microns; here we will keep this parameter
fixed at $L$=1 micron.)

To produce a discretized surface-mesh
representation of this geometry, we run it through 
[[gmsh]], using the [[gmsh]] command-line argument
`-setnumber` to fix a value for the `N` parameter.

````bash
% gmsh -2 -setnumber N 4 Square_N.geo -o Square.msh
% RenameMesh Square.msh

% gmsh -2 -setnumber N 8 Square_N.geo -o Square.msh
% RenameMesh Square.msh
````

(Here [`RenameMesh`][RenameMesh] is a simple `bash` script
that uses [[scuff-analyze]] to count the number of interior
edges in a surface mesh and rename the mesh file accordingly.)
These commands produce the files `Square_40.msh`
and `Square_176.msh`.
These meshes may be visualized in [[gmsh]]:

````bash
% gmsh Square_40.msh
% gmsh Square_176.msh
````

Note the following:

 * For 2D periodic geometries in [[scuff-em]], the 
   lattice vectors must lie in the $xy$ plane.

 * For surfaces that straddle the unit-cell boundaries
   (as is the case here), each triangle edge that lies
   on any edge of the unit cell must have an identical
   image edge on the opposite side of the unit cell.
   An easy way to achieve this is to use *extrusions*
   in [[gmsh]], as in the `.geo` file above.

 * In this case the unit cell dimensions are 
   $L_x\times L_y$ where $L_x=L_y=1\, \mu\text{m}$.
   (More generally, $L_x$ and $L_y$ may be any arbitrary
   nonzero values, and they need not equal each other.)

--------------------------------------------------
# [[scuff-em]] geometry files

The 
[<span class="SC">scuff-em</span> geometry files][Geometries]
describing an infinite-area PEC ground plane at $z$=0
are 
[`PECPlate_40.scuffgeo`](PECPlate_40.scuffgeo)
and 
[`PECPlate_176.scuffgeo`](PECPlate_176.scuffgeo).

The 
[<span class="SC">scuff-em</span> geometry files][Geometries]
describing an infinite aluminum half-space occupying
the region $z<0$ are 
[`AlHalfSpace_40.scuffgeo`](AlHalfSpace_40.scuffgeo)
and 
[`AlHalfSpace_176.scuffgeo`](AlHalfSpace_176.scuffgeo).

Here's the content of the file `AlHalfSpace_176.scuffgeo`:

````bash
# this comes from Phys Rev B **68** 245405
MATERIAL ALUMINUM
    wp = 1.747e16; 
    gamma = 7.596e13;
    Eps(w) = 1 - wp^2 / (w * (w + i*gamma));
ENDMATERIAL

LATTICE
	VECTOR 1.0  0.0
	VECTOR 0.0  1.0 
ENDLATTICE

REGION Exterior       MATERIAL Aluminum
REGION UpperHalfSpace MATERIAL Vacuum

SURFACE Plate
	MESHFILE Square_176.msh
	REGIONS Exterior Aluminum
ENDSURFACE
````

--------------------------------------------------
# List of evaluation points

We'll compute the LDOS at two points, located 
a distance of 0.1 $\mu$m and 1 $\mu$m above the 
origin on the $z$ axis. The `EPFile` looks like
this:

````
0.0 0.0 0.1
0.0 0.0 1.0
````
--------------------------------------------------
# Launching a Bloch-vector-resolved run

Before running a full Brillouin-zone-integrated
calculation to get the full LDOS at a given
frequency $\omega$, we'll first run calculations
at a set of pre-specified individual Bloch vectors
$\mathbf k_{\small{\text{B}}}$ in the Brillouin zone.

We'll take $\omega=3\times 10^{14}$ rad/sec
(or $\omega=1$ in [[scuff]] units) and 
will consider Bloch
vectors of the form $\mathbf k_{\small{\text{B}}}=(0,k_y)$ for 
values of $k_y$ running from $0$ to $\pi/L$
(where $L$=1 $\mu$m is the lattice constant
in this case).
Thus we create a text file called [`OKBFile`](OKBFile) 
that looks like this:

````
1.0 0.0 0.00
1.0 0.0 0.10
...
1.0 0.0 3.14
````

As a sanity check, we will do two [[scuff-ldos]] runs,
one in which the LDOS is computed using a semi-analytical
approach (plane-wave decomposition) and another in which
the LDOS is computed using the [[scuff-em]] core
library. The semi-analytical approach is implemented
by [scuff-ldos] and is requested by giving all the same
command-line arguments you would give to do an ordinary
[scuff-ldos] calculation, but with the additional argument
`--HalfSpace Aluminum.` (Note that you must specify a `.scuffgeo`
file for the `--HalfSpace` calculation, even though the 
meshed geometry is not used in this calculation; the
`.scuffgeo` file is used to specify the lattice.)

#### Command-line arguments for quasi-analytical calculation using plane-wave decomposition:

````bash
  #!/bin/bash
  ARGS=""
  ARGS="${ARGS} --geometry AlHalfSpace_40.scuffgeo" 
  ARGS="${ARGS} --EPFile  EPFile"
  ARGS="${ARGS} --OmegakBlochFile OKBFile"
  ARGS="${ARGS} --HalfSpace Aluminum"
  scuff-ldos ${ARGS}
````

#### Command-line arguments for fully numerical calculation using discretized surface meshes:

````bash
  #!/bin/bash
  ARGS=""
  ARGS="${ARGS} --geometry AlHalfSpace_40.scuffgeo" 
  ARGS="${ARGS} --EPFile  EPFile"
  ARGS="${ARGS} --OmegakBlochFile OKBFile"
  scuff-ldos ${ARGS}
````

Both of these runs produce an output file 
named `AlHalfSpace_40.byOmegakBloch.` (If you
don't rename the file after the first run, 
the data from the second run will simply be 
appended to the file below the first set of data).

Here's a plot of $\mathcal{G}^E_{zz}$ (the $zz$
component of the electric dyadic Green's function)
versus the $y$-component of the Bloch vector
as computed using the analytical and [[scuff-em]]
approaches.

![aluminum LDOS data](AluminumLDOS.png)

--------------------------------------------------
# Launching a Brillouin-zone integrated run

Finally, we'll ask [[scuff-ldos]] to perform
the Brillouin-zone integrations at each frequency
to compute the full LDOS at that frequency.
To do this we simply use `--OmegaFile`` instead
of `--OmegakBlochFile`` to specify a list of 
frequencies instead of a list of 
$(\omega,\mathbf{k}_{\small{\hbox{B}}})$
points; by not specifying individual Bloch vectors
we are implicitly asking [[scuff-ldos]] to perform
the Brillouin-zone integral at each frequency.

Thus we go like this:

````bash
  #!/bin/bash
  ARGS=""
  ARGS="${ARGS} --geometry AlHalfSpace_40.scuffgeo" 
  ARGS="${ARGS} --EPFile  EPFile"
  ARGS="${ARGS} --OmegaFile OmegaFile"
  scuff-ldos ${ARGS}
````

where [`OmegaFile`](OmegaFile) is a simple text
file specifying a list of angular frequencies in
units of $3\cdot 10^{14}$ rad/sec.

[Geometries]:          ../../reference/Geometries.md
[ExtendedGeometries]:  ../../reference/Geometries.md#Extended
[RenameMesh]:          ../../examples/SiO2Spheres/RenameMesh
[scuff-ldos]:          ../../applications/scuff-ldos/scuff-ldos.md
[LDOSMemo]:            ../../tex/scuff-ldos.pdf
