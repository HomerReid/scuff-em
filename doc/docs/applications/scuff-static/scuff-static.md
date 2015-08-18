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

# 1. What <span class="SC">scuff-neq</span> actually computes

--------------------------------------------------

<a name="CommandLineOptions"></a>
# 2. <span class="SC">scuff-static</span> command-line options

--------------------------------------------------

<a name="OutputFiles"></a>
# 3. <span class="SC">scuff-static</span> output files

<a name="Examples"></a>
# 4. Examples of calculations using <span class="SC">scuff-static</span>

+ [Polarizability of the platonic solids](../../examples/PlatonicSolids.md)

+ [Self- and mutual-capacitance of irregularly shaped conductors](../../examples/PlatonicSolids.md)

+ [Electrostatic fields in the vicinity of a complicated gate array](../../examples/
