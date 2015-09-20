<p align="center"><img align="center" src="img/scuffEMLogo.png"></p>

<p align="center"><h1 align="center">
 <span class="SmallCaps">scuff-em</span> documentation: Table of contents
</h1> 
</p>

** Note:** As of June 2015, the [[scuff-em]] documentation is the process
   of being ported from an older format (pure HTML, maintained separately
   from the code repository) to a newer format (plain-text markdown files,
   maintained together with the repository, and built by
   [mkdocs](http://www.mkdocs.org) into the documentation website).
   See [here](forDevelopers/Documentation.md)
   for more on the new documentation system and information
   on how to contribute.

   While the porting process is underway, some of the links on these  
   pages will take you to earlier versions of the documentation. Thanks
   for bearing with us as we complete the update!

## General reference

* [Top-level overview](reference/TopLevel.md)
* [Installation](reference/Installing.md)
* [Geometry files](reference/Geometries.md)
* [Material descriptions](reference/Materials.md)
* [Geometrical transformations](reference/Transformations.md)
* [Incident Fields](reference/IncidentFields.md)

## Command-line application reference

+ [General reference for <span class="SC">scuff-em</span> command-line applications][GeneralReference]

<br>

+ *Nanophotonics codes*
    + [<span class="SC">scuff-scatter</span>][scuff-scatter]: general-purpose electromagnetic scattering
    + [<span class="SC">scuff-ldos</span>][scuff-ldos]                  - photonic **l**ocal **d**ensity **o**f **s**tates
    + [<span class="SC">scuff-transmission</span>][scuff-transmission]  - plane-wave transmission through 2D extended structures
    + [<span class="SC">scuff-tmatrix</span>][scuff-tmatrix]            - T-matrices of arbitrary compact objects in the spherical-wave basis  

<br>

+ *Equilibrium Casimir codes*
    + [<span class="SC">scuff-cas3d</span>][scuff-cas3D]   - Casimir energies, forces, and torques
    + [<span class="SC">scuff-caspol</span>][scuff-caspol] - Casimir-Polder potentials
       
<br>
      
+ *Non-equilibrium Casimir/ heat-transfer code*
    + [<span class="SC">scuff-neq</span>][scuff-neq]       - radiative heat transfer and non-equilibrium Casimir forces/torques
  
<br>
  
+ *RF / microwave engineering code*
    + [<span class="SC">scuff-rf</span>][scuff-RF]         - multiport network parameters
                                     (S- and impedance parameters),
                                     and radiated fields, for passive RF
                                     and microwave structures.

<br>

+ *Electrostatics code*

    + [<span class="SC">scuff-static</span>][scuff-static] - pure electrostatics problems:
                                     capacitance matrices, DC polarizabilities,
                                     electrostatic potentials and fields

<br>

+ *Utility code*
    + [<span class="SC">scuff-analyze</span>][scuff-analyze] - diagnostic tool to print info on [[scuff-em]] geometries

<br>

## API reference

* [<span class="SC">libscuff</span>][libscuff] - Accessing [[scuff-em]] from C++ and python programs

## Developer reference

* [Implementation][Implementation] - how [[scuff-em]] works
* [DataStructures][DataStructures] - data structures and methods inside the [[scuff-em]] core library
* [Documentation][Documentation]   - about the [[scuff-em]] documentation

## FAQ

## Technical memos

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
[libscuff]:           API/libscuff.md
[Implementation]:     forDevelopers/Implementation.md
[DataStructures]:     forDevelopers/DataStructures.md
[Documentation]:      forDevelopers/Documentation.md
