<p align="center"><img align="center" src="img/scuffEMLogo.png"></p>

<p align="center"><h1 align="center">
 <span class="SmallCaps">scuff-em</span> documentation: Table of contents
</h1> 
</p>

**If you find any inconsistencies or missing bits in the documentation, please
  [file an issue on the <span class="SC">scuff-em</span> GitHub page][GitHub].**

## General reference

* [Top-level overview](reference/TopLevel.md)
* [Installation](reference/Installing.md)
* [Geometry files](reference/Geometries.md)
* [Material descriptions](reference/Materials.md)
* [Geometrical transformations](reference/Transformations.md)
* [Incident Fields](reference/IncidentFields.md)
* [Brillouin-zone integration](reference/BrillouinZoneIntegration.md)
* [FAQ](reference/FAQ.md)

# Tutorial Examples

### Photonics, imaging 

* [Mie scattering](examples/MieScattering/MieScattering.md)
* [Electrostatics of a spherical dielectric shell](examples/DielectricShell/DielectricShell.md)
* [LDOS and polarization-sensitive response in plasmonic bowtie antennas](http://homerreid.github.io/SCUFFEMTutorialSymposium/BowtieAntennas/)
* [LDOS and dyadic Green functions above an infinite aluminum half-space](examples/HalfSpaceLDOS/HalfSpaceLDOS.md)
* [Spatially-resolved study of plane-wave transmission through an infinite-area thin dielectric film](examples/ThinFilm/ThinFilm.md)
* [Plane-wave transmission with rotation of polarization: Chiral metasurface rotator](http://homerreid.github.io/SCUFFEMTutorialSymposium/ChiralMetasurfaceRotator)
* [Imaging diffraction patterns of discs, disc arrays, and hole arrays in metal screens](examples/DiffractionPatterns/DiffractionPatterns.md)
* [Imaging diffraction patterns of aperiodic arrays: Vogel spirals](http://homerreid.github.io/SCUFFEMTutorialSymposium/SpiralArrays/)

### Fluctuation-induced phenomena

* [Casimir forces in a compact geometry](examples/CubeTorus.md)
* [Casimir forces in a 1D extended geometry](examples/SiliconBeams/SiliconBeams.md)
* [Casimir forces in a 2D extended geometry](examples/SiliconSlabs/SiliconSlabs.md)
* [Thermal radiation, heat transfer, and non-equilibrium Casimir forces between silicon dioxide spheres](examples/SiO2Spheres/SiO2Spheres.md)
* [Spatial distribution of heat flux from a warm tip over a cold substrate](examples/TipSubstrate/TipSubstrate.md)

### Pure electrostatics

* [Electrostatic polarizability of platonic solids](examples/PlatonicSolids/PlatonicSolids.md)
* [Capacitance of two-body capacitors](examples/TwoBodyCapacitors/TwoBodyCapacitors.md)
* [Capacitance of PCB stripline trace](examples/StriplineCapacitor/StriplineCapacitor.md)
* [Implicit handling of multilayer dielectric substrates](examples/ImplicitSubstrate/ImplicitSubstrate.md)

### BU Symposium on Open-Source CAD Tools for Photonic Design Modeling

* [BU Symposium on Open-Source CAD tools: <span class=SC>scuff-em</span> for photonics](http://homerreid.github.io/SCUFFEMTutorialSymposium/)

# Command-line application reference

* [General reference for <span class="SC">scuff-em</span> command-line applications][GeneralReference]

+ *Nanophotonics codes*
    + [<span class="SC">scuff-scatter</span>][scuff-scatter]: general-purpose electromagnetic scattering
    + [<span class="SC">scuff-ldos</span>][scuff-ldos]                  - photonic **l**ocal **d**ensity **o**f **s**tates
    + [<span class="SC">scuff-transmission</span>][scuff-transmission]  - plane-wave transmission through 2D extended structures
    + [<span class="SC">scuff-tmatrix</span>][scuff-tmatrix]            - T-matrices of arbitrary compact objects in the spherical-wave basis  

* *Equilibrium Casimir codes*
    + [<span class="SC">scuff-cas3d</span>][scuff-cas3D]   - Casimir energies, forces, and torques
    + [<span class="SC">scuff-caspol</span>][scuff-caspol] - Casimir-Polder potentials
       
* *Non-equilibrium Casimir/ heat-transfer code*
    + [<span class="SC">scuff-neq</span>][scuff-neq]       - radiative heat transfer and non-equilibrium Casimir forces/torques
  
* *RF / microwave engineering code*
    + [<span class="SC">scuff-rf</span>][scuff-RF]         - multiport network parameters
                                     (S- and impedance parameters),
                                     and radiated fields, for passive RF
                                     and microwave structures.

* *Electrostatics code*

    + [<span class="SC">scuff-static</span>][scuff-static] - pure electrostatics problems:
                                     capacitance matrices, DC polarizabilities,
                                     electrostatic potentials and fields

* *Utility codes*
    + [<span class="SC">scuff-analyze</span>][scuff-analyze] - diagnostic tool to print info on [[scuff-em]] geometries
    + [<span class="SC">scuff-integrate</span>][scuff-integrate] - utility tool to integrate functions using samples tabulated in data files

## Validation Test Suite

* [Overview of the <span class="CodeName">scuff-em</span> test suite](tests/Overview.md)
* [Mie scattering](tests/MieScattering/MieScattering.md)
* [Fresnel scattering](tests/FresnelScattering/FresnelScattering.md)
* [Equilibrium Casimir forces between spheres](tests/CasimirSpheres/CasimirSpheres.md)
* [Equilibrium Casimir forces between plates](tests/CasimirPlates/CasimirPlates.md)
* [Equilibrium Casimir-Polder potential near a sphere](tests/CPSphere/CPSphere.md)
* [Equilibrium Casimir-Polder potential near a plate](tests/CPPlate/CPPlate.md)
* [Heat transfer and non-equilibrium Casimir forces between spheres](tests/NEQSpheres/NEQSpheres.md)
* [Low-level tests of the <span class="CodeName">scuff-em</span> core library](tests/libscuff/libscuff.md)

## API reference

* [<span class="SC">libscuff</span>][libscuff] - Accessing [[scuff-em]] from C++ and python programs

## Developer reference

* [Implementation][Implementation]        - how [[scuff-em]] works
* [DataStructures][DataStructures]        - data structures and methods inside the [[scuff-em]] core library
* [Documentation][Documentation]          - about the [[scuff-em]] documentation
* [Singular integrals][SingularIntegrals] - [[scuff-em]]'s code for computing singular triangle-product integrals

## FAQ

* [FAQ][FAQ] - frequently asked questions about [[scuff-em]]

## Technical memos

* [<span class="SC">libscuff</span> Implementation and Technical Details](tex/lsInnards.pdf) - a technical memo describing many details of the core library implementation
* [Computation of power, force, and torque in <span class="SC">scuff-em</span>](tex/PFT.pdf) - a technical memo describing methods for computing power, force and torque, with applications to both classical scattering and non-equilibrium fluctuational electrodynamics
* [Computation of fields near panels in <span class="SC">scuff-em</span>](tex/NearFields.pdf) - a technical memo describing the computation of fields near triangular panels in discretized surface meshes
* [Computation of Green's Functions and LDOS in <span class="SC">scuff-em</span>](tex/scuff-ldos.pdf) - a technical memo describing the implementation of the [<span class="SC">scuff-ldos</span>](applications/scuff-ldos) code for computing dyadic Green's functions and local photonic densities of states
* [Computation of reflection and transmission coefficients in <span class="SC">scuff-em</span>](tex/scuff-transmission.pdf) - a technical memo describing the implementation of the [<span class="SC">scuff-transmission</span>](applications/scuff-transmission) code for computing plane-wave transmission and reflection coefficients
* [<span class=SC>scuff-static:</span> Pure electrostatics in <span class=SC>scuff-em</span>](tex/scuff-static.pdf) -- a technical memo describing the implementation of the [<span class=SC>scuff-static</span>](applications/scuff-static) code for electrostatics
* [Electromagnetism in the vector-spherical-wave (VSW) basis](tex/scuffSpherical.pdf) -- a technical memo collecting results of various classical electromagnetism calculations in the vector-spherical-wave basis

[scuffEMLogo]:        img/scuffEMLogo.png
[GeneralReference]:   applications/GeneralReference.md
[scuff-scatter]:      applications/scuff-scatter/scuff-scatter.md
[scuff-ldos]:         applications/scuff-ldos/scuff-ldos.md
[scuff-transmission]: applications/scuff-transmission/scuff-transmission.md
[scuff-tmatrix]:      applications/scuff-tmatrix/scuff-tmatrix.md
[scuff-cas3D]:        applications/scuff-cas3D/scuff-cas3D.md
[scuff-caspol]:       applications/scuff-caspol/scuff-caspol.md
[scuff-neq]:          applications/scuff-neq/scuff-neq.md
[scuff-RF]:           applications/scuff-RF/scuff-RF.md
[scuff-static]:       applications/scuff-static/scuff-static.md
[scuff-analyze]:      applications/scuff-analyze/scuff-analyze.md
[scuff-integrate]:    applications/scuff-integrate/scuff-integrate.md
[libscuff]:           API/libscuff.md
[Implementation]:     forDevelopers/Implementation.md
[DataStructures]:     forDevelopers/DataStructures.md
[Documentation]:      forDevelopers/Documentation.md
[SingularIntegrals]:  forDevelopers/SingularIntegrals.md
[GitHub]:             https://github.com/HomerReid/scuff-em/
[FAQ]:                reference/FAQ.md
