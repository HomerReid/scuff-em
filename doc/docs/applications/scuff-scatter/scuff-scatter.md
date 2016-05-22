<h1> Solving electromagnetic scattering problems with 
     <span class="SC">scuff-scatter</span>
</h1>

[[scuff-scatter]] is a tool within the [[scuff-em]] code suite
for solving classical scattering problems involving
user-specified incident fields impinging on a material
geometry.

To run a scattering calculation using [[scuff-scatter]], you will

+ Create a [<span class="SC">scuff-em</span> geometry file][Geometries]
describing the shapes and material properties of the scattering objects in your geometry

+ Choose the [incident field][IncidentFields] that will scatter off your objects: a plane wave, a gaussian beam, a point dipole source, or some combination thereof

+ Run [[scuff-scatter]] with command-line options specifying the geometry, the frequencies, the incident field, and the type of output you wish to get back.

The various output quantities that you can ask [[scuff-scatter]] to generate include the following:

+ The components of the scattered and total electric and magnetic fields at arbitrary user-specified points away from scattering surfaces. (The points may lie inside or outside the scattering objects).

+ The components of the total electric and magnetic fields on the scattering surfaces. (These quantities may alternatively be interpreted as effective surface currents and charges that give rise to the scattered fields.)

+ The electric and magnetic dipole moments induced by the incident field on the scattering objects. (These are obtained from the interpretation of the tangential fields as effective sources that radiate the scattered fields.)

+ The total power scattered by, and the total power absorbed by, the scattering objects from the incident field.

+ The total force and/or torque exerted on the scattering objects by the incident fields (radiation pressure).

+ Visualization files plotting the electric and magnetic surface currents, and the associated charge densities, 
  induced by the incident fields on the scattering objects.

+ Visualization files plotting field components and Poynting fluxes on arbitrary user-specified surface meshes.

For more sophisticated users, [[scuff-scatter]] also offers
an *advanced mode* of operation that exposes---at the
command-line level---some of the key efficiencies of the
[surface-integral-equation formulation][Implementation]
implemented by the
[<span class="SC">scuff-em</sc> core library][libscuff].
This offers significant speedup for certain types of computations,
at the expense of a slightly more effort required to set
up your calculation.

[TOC]

<a name="Options"></a>
# 1. <span class="SC">scuff-scatter</span> command-line options

### Common options

[[scuff-scatter]] recognizes the following subset of the 
[list of commonly accepted options to <span class="SC">scuff-em</span> command-line codes][CommonOptions].

````bash
--geometry
--EPFile
--Omega
--OmegaFile
--Cache
--ReadCache
--WriteCache
````

## Options defining the scattering problem

````bash
--geometry MyGeometry.scuffgeo
````

Specifies the geometry input file.

````bash
--Omega      3.1415
--OmegaFile  MyOmegaFile
--Lambda     0.5
--LambdaFile MyLambdaFile
````

Specifies the angular frequencies at which to
run calculations. (Angular frequencies are interpreted
in units of $c/1\,\mu\text{m}=3\cdot 10^{14}$ rad/sec.)
The `--Omega` option may be used more than once 
to specify multiple frequencies. Alternatively,
the `--OmegaFile` option may be used to specify the
name of a file containing a list of frequencies (one per
line) at which to run calculations.

The options `--Lambda` and `--LambdaFile` may alternatively
be used to define the frequencies at which to run calculations
in terms of the corresponding free-space wavelength
$\lambda=\frac{2\pi}{\omega}$, interpreted in units
of microns. Thus `--Omega 3.1415` and `--Lambda 0.5`
are equivalent; both specify an angular frequency
$\omega=\pi \cdot 3\cdot 10^{14}$rad/sec, 
corresponding
to a free-space wavelength of $\lambda=0.5\,\mu$m.

Note: Even if you use `--Lambda / --LambdaFile`
instead of `--Omega / --OmegaFile` to specify the
computational frequencies, the output files
will still report data in terms of the frequency
(the equivalent `--Omega` values), not the 
wavelength. To plot your data versus wavelength
instead of frequency, just plot versus the
quantity $\frac{2\pi}{\omega}$.

## Options defining the incident field

The options for specifying incident fields in
[[scuff-em]] are described in detail on the page
[Incident fields in <span class="SC">scuff-em</span>][IncidentFields];
here we just list the
available options without commentary.

````bash
--pwDirection    nx ny nz
--pwPolarization Ex Ey Ez
````


````bash
--psStrength Px Py Pz
--psLocation xx yy zz
````


````bash
--gbDirection nx ny nz
--gbPolarization Ex Ey Ez
--gbCenter Cx Cy Cz
--gbWaist W
````

(As in [[scuff-scatter]], these options may occur multiple times 
to define superpositions of multiple types of incident field.)

## Options requesting scattered and total fields

````bash
 --EPFile MyEPFile
````

Specifies a list of evaluation points at which to
compute and report components of the scattered and total
fields. This option may be specified more than once to 
define multiple sets of field evaluation points. 

## Options requesting power, force, and torque data

````bash
 --PFTFile     MyGeometry.PFT
 --EMTPFTFile  MyGeometry.EMTPFT
 --OPFTFile    MyGeometry.OPFT
 --DSIPFTFile  MyGeometry.DSIPFT
````

Each of these options requests that power, force, and torque (PFT)
data be written to a file of the specified name. The resulting
files all have the same file format---reporting absorbed and
scattered power, force (radiation pressure), and torque
for all objects in the geometry at all frequencies you
requested (see the file header for details)---but differ 
in the algorithm used to compute the force:

+ The "energy-momentum-transfer" PFT (EMTPFT) method
computes powers, forces and torques by considering the
Joule heating of, and Lorentz force on, the
surface currents in the presence of the total 
fields. (This is the default, so the ``--PFTFile`` 
option is synonymous with ``--EMTPFTFile``.)

+ The "displaced-surface-integral" PFT (DSIPFT) method
computes PFTs by integrating the Poynting vector
and Maxwell stress tensor over a bounding surface 
surrounding the body.

+ The "overlap" PFT (OPFT) method computes PFTs
directly from the surface currents by exploiting
the relationship between the surface currents
and the total electric and magnetic fields at 
body surfaces.

<a name="AdvancedMode"></a>
# 2. <span class="SC">scuff-scatter</span> advanced mode

For some types of calculation it is possible to achieve
significant computational accelerations by taking advantage
of certain efficiencies inherent in the particular 
mathematical strategy used by the 
[<span class="SC">scuff-em</sc> core library][libscuff]
to solve Maxwell's equations---namely, the
discretized surface-integral-equation (SIE) formulation.

You can read about all the gory details of SIE solvers 
[here][Implementation], but for the purposes of 
this discussion all you really need to know is this: For a
given material geometry irradiated by a given incident field 
at a given frequency, [[scuff-em]] assembles and solves a 
linear system of the form
$$ \mathbf{M}(\omega) \mathbf{c} = \mathbf{f}^{\text{inc}}$$
where 

+ **c** represents the unknown surface currents for which we are solving,

+ the RHS vector **f** depends on the geometry, the frequency, and the incident field,

+ the matrix **M** depends on the scattering geometry and the frequency but not on the incident field.
More specifically, for a scattering geometry consisting of *N* objects
(or *N* surfaces in a 
[regions-and-surfaces geometry specification][Geometries]),
the matrix **M** has an *N&times;N* block structure in which the *(m,n)* block
describes the interactions of object *m* with object *n*.

Armed with just this much knowledge, we can understand the
two key efficiencies possible in SIE scattering calculations:

+ **(1)** First, suppose that, in a geometry consisting of 2 or more bodies,
we would like to perform calculations for various different
relative geometric configurations of the bodies---for example,
different separation distances or rotation angles between bodies---at
the same frequency.
The diagonal blocks of the **M** matrix, which represent 
the self-interactions of objects and are the most costly
blocks to compute, are *independent* of the relative 
configuration of the various objects in the geometry,
and thus need only be computed *once* for a given
geometry at a given frequency, after which they may be
reused for any number of calculations involving 
rearrangements of the relative positions and orientations
of the bodies.    


    Thus, if we are interested in running calculations
for a sphere-cube geometry at (say) 7 different
values of the surface-surface separation, it greatly
behooves us to assemble the diagonal (self-interaction)
blocks just *once* per frequency, then reuse 
these blocks for each of the 7 separation distances.
The sphere-cube interaction block of the matrix
must be recomputed at each separation distance, but 
this is relatively cheap compared to the cost of 
computing the sphere-sphere and cube-cube 
self-interaction blocks.   


&nbsp;


+ **(2)** In the equation above, the LHS is *independent* of
the incident field. This means that, once we have 
assembled and LU-factorized the **M** matrix 
for given geometry at a given frequency (a procedure
which scales asymptotically like $\sim T^3$ with $T$ 
the total number of triangles in our surface meshes) 
we can solve scattering problems for any number of
incident fields with cost $\sim T^2$ per incident
field---that is, essentially *for free* compared
to the cost of assembling and factorizing the matrix.


    Thus, if we are interested in observing the 
scattering properties of our geometry under irradiation
by 7 different types of incident field (say,
plane waves originating from 7 different angles)
it greatly behooves us to form and LU-factorize
the **M** matrix just *once* for this frequency,
then reuse the factorized matrix to solve the linear
system above for the 7 different types of incident field.

To take advantage of efficiency **(a)**, [[scuff-scatter]]
supports the command-line option

````
  --transfile MyTransFile
````

where ``MyTransFile`` is a 
[list of geometrical transformations][Transformations].

To take advantage of efficiency **(b)**, [[scuff-scatter]]
supports the command-line option

````
  --IFFile    MyIFFile
````

where ``MyIFFile`` is a 
[list of incident fields][IFList].

<a name="Examples"></a>
# 3. <span class="SC">scuff-scatter</span> examples

+ [Mie scattering][MieScattering]
+ [Electrostatics of a spherical dielectric shell][DielectricShell]
+ [Spatially-resolved study of plane-wave transmission through a infinite-area thin dielectric film][ThinFilm]
+ [Diffraction of a plane wave by a discs, disc arrays, and hole arrays][DiffractionPatterns]

[CommonOptions]:               ../GeneralReference.md#CommonOptions
[Geometries]:                  ../../reference/Geometries.md
[Transformations]:              ../../reference/Transformations.md
[IncidentFields]:              ../../reference/IncidentFields.md
[Implementation]:              ../../forDevelopers/Implementation.md
[libscuff]:                    ../../API/libscuff.md
[MieScattering]:               ../../examples/MieScattering/MieScattering.md
[DielectricShell]:             ../../examples/DielectricShell/DielectricShell.md
[ThinFilm]:                    ../../examples/ThinFilm/ThinFilm.md
[DiffractionPatterns]:         ../../examples/DiffractionPatterns/DiffractionPatterns.md
