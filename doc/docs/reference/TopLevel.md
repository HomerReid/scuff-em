# Top-level overview of [[scuff-em]]

[[scuff-em]] is a free, open-source software 
implementation of the boundary-element method (BEM)
(or the "method of moments") of electromagnetic
scattering. (More specifically, [[scuff-em]]
implements the EFIE and PMCHWT formulations
of the BEM using RWG basis functions.)

## Calculating with [[scuff-em]]

Access to the [[scuff-em]] computational engine is available
via multiple interfaces.

The *command-line interface* consists of a large number
of [command-line applications](#AvailableApplications) for
running various types of standard calculations in computational
physics. Using [[scuff-em]] in this way requires only 
learning how to specify command-line options to 

The *application programming interface* consists of 
[C++ and python APIs](../API/libscuff.md)
that allow access to internal [[scuff-em]] data structures
and methods for maximal flexibility in implementing your
own custom-designed physics codes.

## Inputs to [[scuff-em]] calculations

Typical inputs to [[scuff-em]] calculations include

+ A [geometry file](Geometries.md) describing the scattering geometry

+ An optional [list of geometric transformations](Transformations.md) 
  to be applied to the geometry, with calculations generally repeated
  at each transformation

+ Specifications of the frequencies (and, for extended geometries,
  the Bloch vectors) at which you want to perform calculations

+ For scattering codes: a specification of the incident fields

+ Specifications of the output quantities you wish to get back: 
  field components at individual points in space, power/force/torque
  information, Casimir quantities, heat-transfer rates, impedance 
  parameters, capacitances, polarizabilities, etc.

## Outputs from [[scuff-em]] calculations

Typical outputs from [[scuff-em]] calculations include

+ text-based data files reporting output quantities

+ Visualization files written in 
  [<span class="SC">GMSH</span>][GMSH] post-processing
  format.

<a name="AvailableApplications"></a>
## Command-line Applications

### Nanophotonics / electromagnetic scattering 

 + [scuff-scatter][scuff-scatter]
> A general-purpose solver forproblems involving
> Available outputs include: scattered and total fields
> at arbitrary points in space; visualization of fields 
> and surface currents; absorbed and scattered power;
> force and torque (radiation pressure); induced dipole
> or spherical multipole moments; and more.
> 

 + scuff-transmission: 
> A specialized solver for computing plane-wave transmission
> in 2D extended geometries: thin films, perforated screens,
> nanoparticle arrays, etc. 

 + scuff-tmatrix:
> A specialized code for computing the elements
> in 2D extended geometries: thin films, perforated screens,
> nanoparticle arrays, etc. 

### Fluctuation-induced interactions

### RF / microwave engineering

### Electrostatics

##Citing [[scuff-em]]

If you find [[scuff-em]] useful for generating
results included in publications, please consider citing both 
**(a)** one of the papers discussing the implementation of
[[scuff-em]], and 
**(b)** the URL for the code. For example, if you are writing
in LaTeX, you might write something like this:

````tex
Numerical computations were performed using {\sc scuff-em}, a free,
open-source software implementation of the boundary-element 
method~\cite{SCUFF1, SCUFF2}.
````

Here the ``SCUFF1`` and ``SCUFF2``
references refer to the following ``.bibtex`` entries:

````tex
@ARTICLE{SCUFF1,
author = {{Homer Reid}, M.~T. and {Johnson}, S.~G.},
title = "{Efficient Computation of Power, Force, and Torque in 
BEM Scattering Calculations}",
journal = {ArXiv e-prints},
archivePrefix = "arXiv",
eprint = {1307.2966},
primaryClass = "physics.comp-ph",
keywords = {Physics - Computational Physics, Physics - Classical Physics},
year = 2013,
month = jul,
}

@ARTICLE{SCUFF2,
note="\texttt{http://homerreid.com/scuff-EM}"
}
````
