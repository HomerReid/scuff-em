<h1>Computing the photonic local density of states with
    <span class="SC">scuff-ldos</span>
</h1>
--------------------------------------------------

[[scuff-ldos]]
is a tool for computing the electromagnetic local density 
of states (LDOS) at points inside or outside compact or 
extended material bodies.

The inputs you supply to [[scuff-ldos]] calculation are

+ A `.scuffgeo` file describing your geometry.

+ A list of evaluation points $\mathbf x$ at which you want to know
  the LDOS.

+ One or more angular frequencies $\omega$ at which to perform
  calculations.

+ Optionally, for periodic geometries, you may additionally specify
  a list of Bloch wavevectors at which to evaluate wavevector-resolved
  contributions to the LDOS. If you do not specify such a list,
  [[scuff-ldos]] will evaluate an integral over the Brillouin zone
  to compute the total LDOS at each $(\omega, \mathbf x)$ point.

The outputs you get back from a [[scuff-ldos]] calculation may
include

+ The LDOS at each $(\omega, \mathbf x)$ point.

+ For periodic geometries, the contributions of individual
  Bloch wavevectors $\mathbf k_{\text{Bloch}}$ to the LDOS.
  If you supplied a list of Bloch vectors as an input,
  you will get wavevector-resolved information for each point
  in your list; otherwise, you will get wavevector-resolved 
  information for each point chosen automatically by 
  [[scuff-ldos]] in its numerical-cubature evaluation of 
  the Brillouin-zone integral.

+ Optionally, the full Cartesian components of the scattering
  parts of the dyadic Green's functions (DGFs) [the 
  electric / magnetic LDOS is proportional to the imaginary part 
  of the trace of the electric / magnetic DGFs].

For testing purposes, [[scuff-ldos]] also includes built-in
functionality to compute the LDOS for some geometries that
may be handled analytically (specifically, infinite-area PEC
ground planes and infinite-area dielectric half-spaces).

[TOC]

--------------------------------------------------
# 1. What <span class="SC">scuff-ldos</span> actually computes

Some technical details on the calculation performed by
[[scuff-ldos]] may be found in [this memo](../../tex/scuff-ldos.pdf).
The long story short is as follows:

What [[scuff-ldos]] actually computes is the scattering
part of the dyadic Green's functions (DGFs) of the geometry 
in question. For a given (angular frequency, evaluation point)
pair $(\omega, \mathbf x)$, these are $3\times 3$ matrices
giving the scattered fields at $\mathbf x$ due to point sources
at $\mathbf x$, with all fields and sources having time
dependence $e^{-i\omega t}$. (The full definition may be
found in the memo above). The LDOS is obtained from the 
imaginary part of the traces of the DGFs.

For non-periodic geometries, [[scuff-ldos]] does 6 scattering
calculations for each $(\omega, \mathbf x)$ point---specifically,
scattering calculations in which the incident field is the
field of an electric or magnetic point source oriented in 
each of the 3 cartesian directions. (Because these calculations
involve the same scattering geometry at the same frequency,
just with different incident fields, they are fast in a 
BEM solver like [[scuff-em]] because the BEM matrix need
not be recomputed anew for each new incident field.) The
results for the LDOS at each $(\omega,\mathbf x)$ point are
reported in the `.LDOS` output file.

For periodic geometries, [[scuff-ldos]] does many scattering
calculations for each $(\omega, \mathbf x)$ point.
Indeed, the DGFs at $\mathbf x$ are defined as the response of 
the system to a single point source at $\mathbf x$; however, in
in [[scuff-em]] calculations for periodic geometries,
*all* currents and fields, including incident fields, 
must be Bloch-periodic, a condition which is not satisfied
by the fields of a single point source at $\mathbf x$. Instead,
what [[scuff-em]] can compute is the response of the system
to a phased *array* of point sources---that is, an infinite
collection of point sources located at points $\mathbf x+\mathbf L$
with phases $e^{ik_BL}$; here $\mathbf L$ ranges over all 
lattice vectors in a 1D or 2D lattice, and $k_B$ is 
a 1D or 2D Bloch wavevector that ranges over the
Brillouin zone (BZ) of the reciprocal lattice. By performing
these calculations at *all* possible Bloch vectors
$k_B$ and adding up the results---that is, by
performing an integration over the BZ, effectively an 
inverse Fourier transform---we recover 
the fields of just the single point source at $\mathbf x$.
For periodic geometries, [[scuff-ldos]] performs this 
BZ integration by numerical cubature for
each $(\omega,\mathbf x)$ point. This involves sampling
the integrand (that is, computing Bloch-periodic DGFs) 
at large numbers of $k_B$ points; these samples, 
which provide Bloch-vector-resolved information on the 
LDOS and DGFs of the system, are reported by [[scuff-ldos]] 
in the`.byOmegakBloch` output file, while results 
for the full (BZ-integrated) LDOS are written to the 
`.LDOS` output file.

--------------------------------------------------
<a name="CommandLineOptions"></a>
# 2. <span class="SC">scuff-ldos</span> command-line options

### Common options

[[scuff-ldos]] recognizes the following subset of the 
[list of commonly accepted options to <span class="SC">scuff-em</span> command-line codes][CommonOptions].

  ````
--geometry
--EPFile
--Omega
--OmegaFile
--OmegakBlochFile
--AbsTol
--RelTol
--FileBase
--Cache
--ReadCache
--WriteCache
  ````
{.toc}

Of these options, `--geometry` and `--EPFile` are
always mandatory, while one of 
`--Omega`, `--OmegaFile`, or `--OmegakBlochFile` 
must also be specified. All
other command-line arguments are optional.

If you specify `--Omega` or `--OmegaFile` in
a calculation involving a periodic geometry,
then [[scuff-ldos]] will perform a numerical
cubature over the Brillouin zone for each
$\omega$ value. (The options `--BZSymmetry`,
`--AbsTol`, and `--RelTol` control the 
parameters of this cubature.) Samples
of the integrand at the cubature points 
will be written to the `.byOmegakBloch` 
file, while the full integrated resuts
will be written to the `.LDOS` file.

Alternatively, if you use `--OmegakBlochFile`
to specify a list of ($\omega, \mathbf k_B$)
points, then [[scuff-ldos]] will skip the
numerical BZ cubature and instead perform
computations at just the points you 
specified. In this case you wil get back a
`.byOmegakBloch` file, but not an
`.LDOS` file.

### Options requesting analytical LDOS calculations

  ````
--GroundPlane

--HalfSpace PEC
--HalfSpace Aluminum

--SkipBZIntegration
  ````
{.toc}

As noted above, for testing purposes [[scuff-ldos]] incorporates

The first option here instructs [[scuff-ldos]] to 
bypass the usual LDOS calculation it would otherwise 
perform and instead to compute the LDOS of an
auxiliary geometry in which the half-space 
region lying below the $xy$ plane (the region
$z<0$) is filled with a homogeneous material
described by the given 
[<span class="SC">scuff-em</span> material designation][Materials].
If this material is `PEC`, then the calculation
is performed using the image-source method.
Otherwise, the calculation is performed using the 
analytical plane-wave decomposition 
outlined in the 
[<span class="SC">scuff-ldos</span> memo](scuff-ldos.pdf).

As illustrated by [this example][HalfSpaceLDOS],
the `--HalfSpace` option is intended to be
tacked on to an otherwise complete [[scuff-ldos]]
command-line containing options such as 
`--geometry` and `--EPFile.` With `--HalfSpace`,
[[scuff-ldos]] performs the same calculation that 
it would do without that option---using the same
evaluation points and the same frequency options---but
just does the calculation a different way. (Although
the surface meshes specified in the `.scuffgeo`
file are not referenced in this case, a 
valid `.scuffgeo` file must still be supplied;
the `LATTICE...ENDLATTICE` section of this 
file is used to define the lattice used 
for the analytical calculation.)

--------------------------------------------------
# 3. <span class="SC">scuff-ldos</span> output files

### The `.log` file

Like all command-line codes in the [[scuff-em]] suite,
[[scuff-ldos]] produces a text output file named
`FileBase.log` that you can follow to monitor the
status of your calculation.

### The `.LDOS` file

This file reports values of the electric and magnetic
LDOS for each angular frequency and each evaluation
point you requested. This file is always produced
for calculations on non-periodic geometries. For
calculation on periodic geometries, this file is 
produced only if you specified the 
`--Omega` and/or `--OmegaFile` command-line options.

### The `.byOmegakBloch` file

This file reports Bloch-vector resolved versions
of the information reported by the `.LDOS` file.
This file is only produced for calculations on
periodic geometries.

--------------------------------------------------
<a name="Examples"></a>
# 4. Examples of calculations using <span class="SC">scuff-ldos</span>

+ [LDOS above a dielectric half-space][HalfSpaceLDOS]

[HalfSpaceLDOS]:   ../../examples/HalfSpaceLDOS/HalfSpaceLDOS.md
[Materials]:       ../../reference/Materials.md
