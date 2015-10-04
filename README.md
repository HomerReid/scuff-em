[![Build Status](https://travis-ci.org/HomerReid/scuff-em.svg?branch=master)](https://travis-ci.org/HomerReid/scuff-em)

SCUFF-EM
========

A comprehensive and full-featured computational physics suite for 
boundary-element analysis of electromagnetic scattering, 
fluctuation-induced phenomena (Casimir forces and radiative 
heat transfer), nanophotonics, RF device engineering, 
electrostatics, and more. Includes a core library with C++ and 
python APIs as well as many command-line applications.

SCUFF-EM stands for

 Surface-CUrrent-Field Formulation of ElectroMagnetism

For documentation and further information on SCUFF-EM visit the 
[SCUFF-EM documentation homepage:](http://homerreid.github.io/scuff-em-documentation)

http://homerreid.github.io/scuff-em-documentation

An older version of the SCUFF-EM documentation is available here:

http://www.homerreid.com/scuff-em

# General reference

* [Top-level overview][TopLevel]
* [Installation][Installing]
* [Geometry files][Geometries]
* [Material descriptions][Materials]
* [Geometrical transformations][Transformations]
* [Incident Fields][IncidentFields]

# Tutorial Examples

* [Mie scattering][MieScattering]
* [Electrostatics of a spherical dielectric shell][DielectricShell]
* [Spatially-resolved study of plane-wave transmission through an infinite-area thin dielectric film][ThinFilm]
* [Imaging diffraction patterns of discs, disc arrays, and hole arrays in metal screens][DiffractionPatterns]
* [Casimir forces in a compact geometry][CubeTorus]
* [Casimir forces in a 1D extended geometry][SiliconBeams]
* [Casimir forces in a 2D extended geometry][SiliconSlabs]
* [Thermal radiation, heat transfer, and non-equilibrium Casimir forces between silicon dioxide spheres][Spheres]
* [Spatial distribution of heat flux from a warm tip over a cold substrate][TipSubstrate]
* [LDOS and dyadic Green functions above an infinite aluminum half-space][HalfSpaceLDOS]
* [Electrostatic polarizability of platonic solids][PlatonicSolids]

# Command-line application reference

+ [General reference for SCUFF-EM command-line applications][GeneralReference]

+ *Nanophotonics codes*

    + [SCUFF-SCATTER][scuff-scatter]: general-purpose electromagnetic scattering
    + [SCUFF-LDOS][scuff-ldos]                  - photonic **l**ocal **d**ensity **o**f **s**tates
    + [SCUFF-TRANSMISSION][scuff-transmission]  - plane-wave transmission through 2D extended structures
    + [SCUFF-TMATRIX][scuff-tmatrix]            - T-matrices of arbitrary compact objects in the spherical-wave basis  

+ *Equilibrium Casimir codes*
    + [SCUFF-CAS3D][scuff-cas3D]   - Casimir energies, forces, and torques
    + [SCUFF-CASPOL][scuff-caspol] - Casimir-Polder potentials
       
+ *Non-equilibrium Casimir/ heat-transfer code*
    + [SCUFF-NEQ][scuff-neq]       - radiative heat transfer and non-equilibrium Casimir forces/torques
  
+ *RF / microwave engineering code*
    + [SCUFF-RF][scuff-RF]         - multiport network parameters
                                     (S- and impedance parameters),
                                     and radiated fields, for passive RF
                                     and microwave structures.

+ *Electrostatics code*

    + [SCUFF-STATIC][scuff-static] - pure electrostatics problems:
                                     capacitance matrices, DC polarizabilities,
                                     electrostatic potentials and fields

+ *Utility code*
    + [SCUFF-ANALYZE][scuff-analyze] - diagnostic tool to print info on SCUFF-EM geometries

# API reference

* [LIBSCUFF][libscuff] - Accessing SCUFF-EM from C++ and python programs

# Developer reference

* [Implementation][Implementation] - how SCUFF-EM works
* [DataStructures][DataStructures] - data structures and methods inside the SCUFF-EM core library
* [Documentation][Documentation]   - about the SCUFF-EM documentation

# Technical memos

* [LIBSCUFF Implementation and Technical Details](tex/lsInnards.pdf) - a technical memo describing many details of the core library implementation
* [Computation of power, force, and torque in SCUFF-EM](tex/PFT.pdf) - a technical memo describing methods for computing power, force and torque, with applications to both classical scattering and non-equilibrium fluctuational electrodynamics
* [Computation of fields near panels in SCUFF-EM](tex/NearFields.pdf) - a technical memo describing the computation of fields near triangular panels in discretized surface meshes
* [Computation of Green's Functions and LDOS in SCUFF-LDOS](applications/scuff-ldos) code for computing dyadic Green's functions and local photonic densities of states 


[TopLevel]:                http://homerreid.github.io/scuff-em-documentation/reference/TopLevel
[Installing]:              http://homerreid.github.io/scuff-em-documentation/reference/Installing
[Geometries]:              http://homerreid.github.io/scuff-em-documentation/reference/Geometries
[Materials]:               http://homerreid.github.io/scuff-em-documentation/reference/Materials
[Transformations]:         http://homerreid.github.io/scuff-em-documentation/reference/Transformations
[IncidentFields]:          http://homerreid.github.io/scuff-em-documentation/reference/IncidentFields
[MieScattering]:           http://homerreid.github.io/scuff-em-documentation/examples/MieScattering/MieScattering
[DielectricShell]:         http://homerreid.github.io/scuff-em-documentation/examples/DielectricShell/DielectricShell
[ThinFilm]:                http://homerreid.github.io/scuff-em-documentation/examples/ThinFilm/ThinFilm
[DiffractionPatterns]:     http://homerreid.github.io/scuff-em-documentation/examples/DiffractionPatterns/DiffractionPatterns
[CubeTorus]:               http://homerreid.github.io/scuff-em-documentation/examples/CubeTorus
[SiliconBeams]:            http://homerreid.github.io/scuff-em-documentation/examples/SiliconBeams/SiliconBeams
[SiliconSlabs]:            http://homerreid.github.io/scuff-em-documentation/examples/SiliconSlabs/SiliconSlabs
[Spheres]:                 http://homerreid.github.io/scuff-em-documentation/examples/SiO2Spheres/SiO/Spheres
[TipSubstrate]:            http://homerreid.github.io/scuff-em-documentation/examples/TipSubstrate/TipSubstrate
[HalfSpaceLDOS]:           http://homerreid.github.io/scuff-em-documentation/examples/HalfSpaceLDOS/HalfSpaceLDOS
[PlatonicSolids]:          http://homerreid.github.io/scuff-em-documentation/examples/PlatonicSolids/PlatonicSolids
[scuffEMLogo]:             http://homerreid.github.io/scuff-em-documentation/img/scuffEMLogo.png
[GeneralReference]:        http://homerreid.github.io/scuff-em-documentation/applications/GeneralReference
[scuff-scatter]:           http://homerreid.github.io/scuff-em-documentation/applications/scuff-scatter/scuff-scatter
[scuff-ldos]:              http://homerreid.github.io/scuff-em-documentation/applications/scuff-ldos/scuff-ldos
[scuff-transmission]:      http://homerreid.github.io/scuff-em-documentation/applications/scuff-transmission/scuff-transmission
[scuff-tmatrix]:           http://homerreid.github.io/scuff-em-documentation/applications/scuff-tmatrix/scuff-tmatrix
[scuff-cas3D]:             http://homerreid.github.io/scuff-em-documentation/applications/scuff-cas3D/scuff-cas3D
[scuff-caspol]:            http://homerreid.github.io/scuff-em-documentation/applications/scuff-caspol/scuff-caspol
[scuff-neq]:               http://homerreid.github.io/scuff-em-documentation/applications/scuff-neq/scuff-neq
[scuff-RF]:                http://homerreid.github.io/scuff-em-documentation/applications/scuff-RF/scuff-RF
[scuff-static]:            http://homerreid.github.io/scuff-em-documentation/applications/scuff-static/scuff-static
[scuff-analyze]:           http://homerreid.github.io/scuff-em-documentation/applications/scuff-analyze/scuff-analyze
[libscuff]:                http://homerreid.github.io/scuff-em-documentation/API/libscuff
[Implementation]:          http://homerreid.github.io/scuff-em-documentation/forDevelopers/Implementation
[DataStructures]:          http://homerreid.github.io/scuff-em-documentation/forDevelopers/DataStructures
[Documentation]:           http://homerreid.github.io/scuff-em-documentation/forDevelopers/Documentation
