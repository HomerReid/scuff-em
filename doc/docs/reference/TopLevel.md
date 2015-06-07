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
of [command-line applications][#AvailableApplications] for
running various types of standard calculations in computational
physics. Using [[scuff-em]] in this way requires only 
learning how to specify command-line options to 

The code-level interface consists of 
[C++ and python APIs][../API/libscuff.md]
that allow access to internal [[scuff-em]] data structures
and methods for maximal flexibility in implementing your
own custom-designed physics codes.

## Inputs to [[scuff-em]] calculations

Typical inputs to [[scuff-em]] calculations include

+ A [geometry file](Geometries.md) describing the scattering geometry

+ Specifications of the output quantities you wish to get back: 
  field components at individual points in space, power/force/torque
  information, Casimir quantities, heat-transfer rates, impedance 
  parameters, capacitances, polarizabilities, etc.

+ For scattering codes: a specification of the incident fields

+ An optional [list of geometric transformations](Transformations.md) 
  to be applied to the geometry, with calculations generally repeated

+ Specifications of the frequencies and/or Bloch vectors at which 
  you want to perform calculations

## Outputs from [[scuff-em]] calculations

Typical outputs from [[scuff-em]] calculations include

+ text-based data files reporting output quantities

+ Visualization files written in 
  [<span class="SC">GMSH</span>][GMSH] post-processing

<a name="AvailableApplications"></a>

## Command-line Applications

## Installing [[scuff-em]]

[GMSH]: http://www.geuz.org/gmsh
