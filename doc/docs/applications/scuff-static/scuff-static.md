<h1> Solving electrostatics problems with
     <span class="SC">scuff-static</span>
</h1>

[[scuff-static]] is a tool within the
[[scuff-em]] code suite for solving
a broad class of electrostatics problems.

The calculations that [[scuff-static]] can
perform include the following:

+ Compute the capacitance matrix (i.e. the self- and mutual-
capacitances) of a collection of conductors.

+ Compute the DC polarizability of a conducting or
dielectric body.

+ Compute the electrostatic potential and field
at arbitrary user-specified points in the vicinity
of conducting or dielectric bodies, with the
conductors maintained at arbitrary user-specified
potentials and (optionally) an arbitrary user-specified
external forcing field.

+ Compute the *C-matrix*, a sort of electrostatic
version of the
["T-matrix"](../scuff-tmatrix/scuff-tmatrix.md)
used to characterize the scattering properties
of bodies at nonzero frequencies. The C-matrix
was shown in
[this paper](http://dx.doi.org/10.1103/PhysRevLett.114.151602)
to be related to quantum-mechanical entanglement
entropy.

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

Here is a brief [technical memo](scuff-static.pdf)
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

### Common options

[[scuff-static]] recognizes the following subset of the
[list of commonly accepted options to
 <span class="SC">scuff-em</span> command-line codes][CommonOptions].

  ````
--geometry
--TransFile
--EPFile
  ````
{.toc}

### Options requesting capacitance-matrix output

  ````
--CapFile    MyCapacitanceMatrix.dat
--TransFile  MyTransFile.trans
  ````
{.toc}

If you specify a file name using `--CapFile`, [[scuff-static]]
will compute the full capacitance matrix for your geometry
and write the data to the specified file. (The file will be
overwritten if it already exists.)

The optional `--TransFile` option may be used to specify a
list of
[geometrical transformations](../../reference/Transformations.md)
to be applied to your geometry. In this case, the full capacitance
matrix will be computed for each specified transformation.

### Options requesting polarizability output

  ````
--PolFile    MyPolFile.dat
  ````
{.toc}

If you specify a file name using `--PolFile`, [[scuff-static]]
will compute the DC polarizability of each object in your
geometry and write the data to the specified file. (The file
will be overwritten if it already exists.)

### Options requesting computation of electrostatic fields

  ````
--EPFile     MyEPFile
--PlotFile   MyPlotFile.pp
--PotFile    MyPotFile
--ConstField [X|Y|Z]
--PhiExt     PhiFile
  ````
{.toc}

If you specify a list of field evaluation points using
`--EPFile MyEPFile`, then [[scuff-static]] will compute
the electrostatic potential and field at each evaluation
point and write the results to a file named `MyEPFile.out`.
(The file will be overwritten if it already exists.)
The file will contain a header explaining how to interpret
its contents.

The `--PlotFile` option may be used to request
creation of a [[gmsh]] visualization file plotting
the induced charge density on all conducting and
dielectric surfaces in the geometry.
If you say `--PlotFile MyPlotFile.pp`, then you
will get a file named `MyPlotFile.pp` which
may be opened in [[gmsh]] for visualization purposes.

By default, all conductor surfaces will be held
at zero potential. If you wish to set one or more
conductor surfaces to non-zero potential, you may do this
by saying `--PotFile MyPotFile` (here "pot" is short
for "potential). `MyPotFile` should be simply a
list of (surface label, potential value) pairs,
like this:

  ````
 UpperSurface  1.2
 LowerSurface -3.4
  ````
{.toc}

where `UpperSurface` and `LowerSurface` are the
labels you assigned to the surfaces in question
in the `.scuffgeo` file. (The label is the string
following the `OBJECT` or `SURFACE` keyword
in the `.scuffgeo` file.) This would set the
conductor surface `UpperSurface` to a potential of
1.2 volts and the conductor surface `LowerSurface`
to a potential of -3.4 volts.

By default, the calculation will be performed with
no external electrostatic field. You can say e.g.
`--ConstField X` to request that the calculation
be performed in the presence of a constant unit-strength
electrostatic field pointing in the positive $x$
direction.

Alternatively, you may use the `--PhiExt` option
to specify the name of a file describing a more
complicated (non-constant) external field.

--------------------------------------------------
<a name="Examples"></a>
# 2. Examples of calculations using <span class="SC">scuff-static</span>

+ [Polarizability of platonic solids][PlatonicSolids]

+ [Self- and mutual-capacitance of pairs of conductors][TwoBodyCapacitors]

+ [Electrostatic fields in the vicinity of a complicated gate array][PlatonicSolids]

[PlatonicSolids]: ../../examples/PlatonicSolids/PlatonicSolids.md
[TwoBodyCapacitors]: ../../examples/TwoBodyCapacitors/TwoBodyCapacitors.md
[CommonOptions]: ../GeneralReference.md#CommonOptions
