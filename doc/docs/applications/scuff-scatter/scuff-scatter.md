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
--Omega      3.14
--OmegaFile  MyOmegaFile
````

Specifies the angular frequencies at which to
run calculations. (Angular frequencies are interpreted
in units of $c/1\,\mu\text{m}=3\cdot 10^{14}$ rad/sec.)
The `--Omega` option may be used more than once 
to specify multiple frequencies. Alternatively,
the `--OmegaFile` option may be used to specify the
name of a file containing a list of frequencies (one per
line) at which to run calculations.

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

<a name="Examples"></a>
# 2. <span class="SC">scuff-scatter</span> examples

+ [Mie scattering][MieScattering]
+ [Electrostatics of a spherical dielectric shell][DielectricShell]
+ [Spatially-resolved study of plane-wave transmission through a infinite-area thin dielectric film][ThinFilm]
+ [Diffraction of a plane wave by a discs, disc arrays, and hole arrays][DiffractionPatterns]

[CommonOptions]:               ../GeneralReference.md#CommonOptions
[Geometries]:                  ../../reference/Geometries.md
[IncidentFields]:              ../../reference/IncidentFields.md
[MieScattering]:               ../../examples/MieScattering/MieScattering.md
[DielectricShell]:             ../../examples/DielectricShell/DielectricShell.md
[ThinFilm]:                    ../../examples/ThinFilm/ThinFilm.md
[DiffractionPatterns]:         ../../examples/DiffractionPatterns/DiffractionPatterns.md
