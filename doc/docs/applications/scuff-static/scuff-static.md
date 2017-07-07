<h1> Solving electrostatics problems with
     <span class="SC">scuff-static</span>
</h1>

[[scuff-static]] is a tool within the
[[scuff-em]] code suite for solving
a broad class of electrostatics problems.

The various calculations that [[scuff-static]] may
be divided into two categories: 
**(1)** calculations that yield
the response of the geometry to one or more 
specific excitation stimuli ("type-1" calculations),
and 
**(2)** calculations
that yield intrinsic properties of the material geometry,
independent of any external excitation ("type-2"
calculations).

**Type-1 calculations**

For type-1 calculations, you pass command-line options
to <span class=SC>scuff-static</span> that specify a specific
electrostatic "excitation." A single excitation consists of
one or both of the following:

 + a set of electrostatic potential values at which to maintain
   the conductors in the geometry (option `--potfile`);
   conductors for which no value is specified are maintained 
   at 0V by default.
   
 + a description of an external electrostatic field,
   which may be a superposition of
      + the field of one or point charges (`--monopole`),
      + the field of one or point dipoles (`--dipole`),
      + a spatially-constant electrostatic field (`--ConstField`)
      + an spatially-varying electrostatic potential defined by a
        user-specified functional form (`--PhiExt`).

You may also use the ``--ExcitationFile`` option to
specify a *list* of excitations, each consisting of one or both
of the above ingredients. This is equivalent to running
<span class=SC>scuff-static</span> multiple times
with options specifying a single excitation 
(`--potfile`, `--monopole`, etc.) each time---but is much *faster*,
because of the structure of BEM solvers: once we have
done the set-up work needed to solve for one
excitation (namely, forming and factorizing the 
[BEM matrix][scuffStaticMemo]), solving for any number
of additional excitations is essentially "free" computationally.

Having specified one or more excitations for a type-1
calculation, you may request multiple types of output,
which will be computed separately for each excitation:

+ Values of the electrostatic potential $\phi(\mathbf x)$ and
  electric field $\mathbf E(\mathbf x)$ components
  at arbitrary user-specified **E**valuation **P**oints $\mathbf x$ 
  inside or outside material bodies, output in the form of
  numerical data in a text file (option `--EPFile`)

+ Field visualization files plotting the potential and
  fields on one or more arbitrary user-specified **F**ield-**V**isualization
  mesh(es) (option `--FVMesh`), with the field-visualization
  mesh(es)
  optionally subject to one or more 
  [geometric transformations][GeometricTransformations]
  (option `--FVMeshTransFile`)

+ Surface-charge visualization files showing
  the distribution of electrostatic surface charges
  induced by the excitation (option `--PlotFile`)

**Type-2 calculations**

In contrast to type-1 calculations, type-2 calculations
do not involve any user-specified excitation. The available
type-2 calculations are:

+ Computation of the capacitance matrix for the conductors 
  in the geometry (option `--capfile`)

+ Computation of the DC polarizabilities of the bodies 
  (both conductors and dielectrics) in the geometry (option `--polfile`)

+ Computation of the *C-matrix*, a sort of electrostatic
  version of the ["T-matrix"](../scuff-tmatrix/scuff-tmatrix.md)
  used to characterize the scattering properties
  of bodies at nonzero frequencies (options `--CMatrixFile`
  and/or `--CMatrixHDF5File`). The C-matrix
  was shown in [this paper](http://dx.doi.org/10.1103/PhysRevLett.114.151602)
  to be related to quantum-mechanical entanglement entropy.

**Geometrical transformations**

For both type-1 and type-2 calculations, you may use
the optional `--transfile` option to specify a list
of [geometric transformations][GeometricTransformations]
to be applied to the scattering geometry; you will then
get back results for each specified output quantity
under each transformation.

**Implicit dielectric substrates**

You may use the optional ``--substrate`` option to specify
that your geometry exists in the presence of
a substrate consisting
of one or more (infinite-area) stacked dielectric layers with an
optional ground plane. This substrate will be treated *implicitly*
by <span class=SC>scuff-static</span>; your `.scuffgeo` file
and the `.msh` files it specifies need only define the
(finite-area) conducting and/or dielectric objects
that lie above (or within the substrate.

**Under the hood: Implementation of <span class=SC>scuff-static</span>**

As a technical detail, we note that the implementation of
[[scuff-static]] actually differs in some significant ways
from the other codes in the [[scuff-em]] suite; in particular,
as compared to the [[scuff-em]] core library,
[[scuff-static]] uses different basis functions and a
different formulation of the boundary-element method, as
appropriate for zero-frequency problems. (More specifically,
[[scuff-static]] expands
surface electric charge densities on PEC and dielectric
surfaces using ``pulse'' basis functions, which are
constant on individual triangles and vanishing everywhere
else.) However, from the implementation standpoint, it
turns out that the calculations needed to implement the
electrostatics calculations in [[scuff-static]] are a
proper *subset* of the calculations already implemented
in [[scuff-em]]. Moreover, from the user's standpoint,
the work needed to set up a [[scuff-static]] problem
(create surface meshes, write geometry files, etc.)
is similar to the setup needed for the nonzero-frequency
codes in the [[scuff-em]] suite.
This is why it makes sense to package these codes together.

Here is a brief [technical memo][scuffStaticMemo]
discussing the implementation of [[scuff-static]],
including both the underlying BEM electrostatics formulation
and the execution of the various types of calculation
(capacitance, polarizability, etc.) that the code can do.

[//]: ###################################################
[//]: ###################################################
[//]: ###################################################

[TOC]

--------------------------------------------------

<a name="CommandLineOptions"></a>
# 1. <span class="SC">scuff-static</span> command-line options

### Options defining the geometry 

````bash
--geometry    MyGeometry.scuffgeo
--substrate   MySubstrate.substrate
--TransFile   Transformations.trans
````

The mandatory `--geometry` option specifies a
 [<span class=SC>scuff-em</span> geometry file][Geometries]
defining one or more surface meshes that comprise
your geometry.

The optional `--substrate` option specifies an 
file describing a substrate consisting of zero or
more dielectric layers stacked atop an optional
ground plane; see [here][ImplicitSubstrateExample]
for an example.

The optional `--TransFile` option specifies a list
of [geometric transformations][GeometricTransformations]
to be applied to your geometry; each output calculation
you request will be repeated once for each
geometric transformation.

### Options defining individual excitations for type-1 calculations
<a name="Type1Excitations"></a>

To run a "type-1" calculation with just a single
excitation, you may specify any combination
of the following options.

````bash
--PotFile    MyPotFile
--ConstField Ex Ey Ez
--Monopole   x y z Q
--Dipole     x y z Px Py Pz
--PhiExt     PhiExpression(x,y,z)
````

Here `MyPotFile` should be a simple text file
consisting of a list of (surface label, potential value) pairs,
where the surface label is the label assigned to a surface
in the `.scuffgeo` file. For example, if your geometry
contains conductors labeled `UpperSurface` and `LowerSurface`,
which you wish to maintain at 1.2 volts and -3.4 volts
respectively, then `MyPotFile` would look like this:

````bash
 UpperSurface  1.2
 LowerSurface -3.4
````

Any conductors that are not specified in the `--PotFile`
are maintained at 0 volts by default.

`--ConstField Ex Ey Ez` specifies a spatially-constant
    external electrostatic field with components
    **E**=(`Ex`,`Ey`,`Ez`).

`--Monopole x y z Q` specifies a point charge of strength `Q`
    at cartesian coordinates (`x`,`y`,`z`).

`--Dipole x y z px py pz` specifies a point dipole
    at cartesian coordinates (`x`,`y`,`z`) 
    with dipole moment **p**=(`px`,`py`,`pz`).

`--PhiExt MyPhiExpression(x,y,z)` specifies a spatially-varying
    electrostatic potential defined by the given function of cartesian
    coordinates `x`, `y`, `z`. Example: `--PhiExt cos(3*x)*cosh(3*y)`

### Options defining multiple excitations for type-1 calculations
<a href="ExcitationFile"></a>

````bash
 --ExcitationFile MyExcitationFile
````

Specifies that `MyExcitationFile` contains a list of
excitations. The file should consist of one or more 
clauses of the form `EXCITATION Label ... ENDEXCITATION`
where `Label` is an arbitrary label you assign to the
excitation (which will be used to tag the corresponding
output data). Each clause should contain one or more
lines, each defining either **(a)** a {ConductorLabel, PhiValue}
pair, or **(b)** an external-field specification.
(Note that all external-field specifications defined within
a single `EXCITATION` clause are present *simultaneously*
when that excitation is active.)

Here's an example of an excitation file for a geometry
containing conductors labeled `UpperSurface` and `LowerSurface`:

````bash
 EXCITATION UpperOnly
   UpperSurface   1
 ENDEXCITATION

 EXCITATION LowerOnly
   LowerSurface  -1
 ENDEXCITATION

 EXCITATION Both
   UpperSurface   1
   LowerSurface  -1
 ENDEXCITATION

 EXCITATION BothWithMonopoles
   UpperSurface   1
   LowerSurface  -1
   Monopole       0 0 +1 1
   Monopole       0 0 -1 1
 ENDEXCITATION

 EXCITATION Ex
   CONSTFIELD     1.0 0.0 0.0
 ENDEXCITATION

 EXCITATION EyAndDipole
   CONSTFIELD     0.0  1.0 0.0
   DIPOLE         -0.1 0.2 0.3 0.4 -0.5 0.2
 ENDEXCITATION
````

### Options defining outputs for type-1 calculation

````bash
--EPFile     MyEPFile
````

Requests that the electrostatic potential and field 
components be computed at each evaluation point
listed in the file `MyEPFile`, which is a simple text
file containing three numbers per line (the coordinates 
of a single evaluation point.) The resulting output file
will contain a header explaining how to interpret its contents.

````bash
--FVMesh          MyFVMesh.msh
--FVMeshTransFile MyFVMeshTransFile
````

Requests that electrostatic fields be computed
on a **F**ield-**V**isualization mesh described by the 
[<span class=SC>gmsh</span>](http://geuz.org/gmsh) mesh
file `MyFVMesh.msh.` This will produce
in an output file with extension `.pp` that can be opened directly
in [[gmsh]] to visualize the fields in any region of space.

The optional `FVMeshTransFile` defines a list of
 [geometric transformations][GeometricTransformations] to be applied
to `MyFVMesh.msh.` This allows you to obtain visualization data on a
large region foliated by multiple translated or rotated copies of a
single mesh screen.


````bash
--PlotFile   MyPlotFile.pp
````

Requests creation of a [[gmsh]] visualization file named `MyPlotFile.pp`
plotting the induced charge density on all conducting and
dielectric surfaces in the geometry.

### Options requesting outputs for type-2 calculations

````bash
--CapFile    MyCapacitanceMatrix.dat
````

If you specify a file name using `--CapFile`, [[scuff-static]]
will compute the full capacitance matrix for your geometry
and write the data to the specified file. (The file will be
overwritten if it already exists.)


````bash
--PolFile    MyPolFile.dat
````

If you specify a file name using `--PolFile`, [[scuff-static]]
will compute the DC polarizability of each object in your
geometry and write the data to the specified file. (The file
will be overwritten if it already exists.)

--------------------------------------------------
<a name="Examples"></a>
# 2. Examples of calculations using <span class="SC">scuff-static</span>

+ [Polarizability of platonic solids][PlatonicSolids]

+ [Self- and mutual-capacitance of pairs of conductors][TwoBodyCapacitors]

+ [Electrostatic fields in the vicinity of a complicated gate array][PlatonicSolids]

[PlatonicSolids]:           ../../examples/PlatonicSolids/PlatonicSolids.md
[TwoBodyCapacitors]:        ../../examples/TwoBodyCapacitors/TwoBodyCapacitors.md
[CommonOptions]:            ../GeneralReference.md#CommonOptions
[GeometricTransformations]: ../../reference/Transformations.md
[Geometries]:               ../../reference/Geometries.md
[scuffStaticMemo]:          scuff-static.pdf
[ImplicitSubstrateExample]: ../../examples/ImplicitSubstrate/ImplicitSubstrate.md
