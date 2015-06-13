# Computing equilibrium Casimir energies, forces, and torques with [[scuff-cas3d]]

[[scuff-cas3d]] is a command-line application within the [[scuff-em]] suite
for modeling Casimir interactions between compact or extended homogeneous
bodies of arbitary shape and arbitrary (linear, isotropic, piecewise homogeneous)
frequency-dependent permittivity and permeability. [[scuff-cas3d]] implements the 
[*fluctuating-surface current (FSC)* approach][FSCPaper]
to numerical Casimir modeling.]

The basic flow of a typical [[scuff-cas3d]] run goes something like this:

+ You create a [<span class="SC">scuff-em</span> geometry file][Geometries]
describing the interacting objects or surfaces in your geometry.

+ Optionally, you define a 
[list of geometric transformations][Transformations] 
to be applied to the geometry for Casimir computations.
For example, if your geometry consists of two nanoparticles, you
might ask for the Casimir force between the particles at 10 different
values of the surface--surface separation.

+ Then you 
+ 
You create a text file with extension .scuffgeo that lists all the objects in your geometry and assigns to each a "material property" specification (something like GOLD or SILICON).
You create a second text file with extension .trans containing a list of transformations to be applied to the geometry, where each transformation is a series of displacements and rotations applied to one or more of the objects in your geometry.
You run scuff-cas3d with appropriate command-line settings to compute Casimir energies, forces, and/or torques for your geometry under each of the transformations you described.
Finally, you interpret the variety of output files that scuff-cas3d emits. In general the one you will care most about is the .out file, which simply lists the computed energy, force, and/or torque at each of the transformations you specified.

[TOC]

# 1. What [[scuff-cas3d]] actually computes

## 1a. Compact objects

The Casimir energy $\mc E$ of a collection of compact bodies,
and the (say) $z$-directed Casimir force on one of those 
bodies, are computed in the FSC approach as integrals 
over the 

$$ \begin{array}{ccc}
   \mc E   &=& \int_0^\infty E(\xi) \, d\xi \\
   \mc F_z &=& \int_0^\infty F(\xi) \, d\xi \\
   \end{array}{ccc}
$$

## 1b. Extended objects

For an extended material configuration described by a
[Bloch-periodic <]

the energy and force \textit{density}
(that is, the energy per unit length for a 1D extended geometry
or per unit area for a 2D extended geometry)

$$ \begin{array}{ccc}
   \mc E 
     &=& 
   \int_0^\infty d\xi \, \int_{\int_{\hbox{BZ}} e(\xi, \mathbf{k}) d\mathbf{k}}_{E(\xi)}
\\
   \end{array}{ccc}
$$

# 2. [[scuff-cas3d]] command-line options

# 3. [[scuff-cas3d]] output files

# 4. [[scuff-cas3d]] 

[EarlierVersion]: http://homerreid.com/scuff-em/scuff-cas3d
[FSCPaper]: http://dx.doi.org/10.1103/PhysRevA.88.022514
[Geometries]: ../reference/Geometries.md
[Transformations]: ../reference/Transformations.md
