# Modeling non-equilibrium electromagnetic fluctuations with [[scuff-neq]]

[[scuff-neq]] is an application code in the [[scuff-em]] suite for 
studying non-equilibrium (NEQ) electromagnetic fluctuation 
phenomena--specifically, for computing *radiative heat-transfer rates* 
and *non-equilibrium Casimir forces* between bodies of arbitrary shapes 
and arbitrary (linear, isotropic, piecewise homogeneous)
frequency-dependent permittivity and permeability.

Mechanically, working with [[scuff-neq]] is similar in many ways to 
working with the equilibrium Casimir code [scuff-cas3d][scuff-cas3d]{.SC}. 
In particular,

+ As in [[scuff-cas3d]], you can request either *(a)* frequency-resolved 
information on heat-transfer rates and NEQ Casimir forces (in which case 
you will specify a list of frequencies and will get back a list of 
frequency-specific values of energy and momentum fluxes) or 
*(b)* frequency-integrated information, in which case you will assign 
temperatures to each body in your geometry and [[scuff-neq]] will 
numerically integrate the fluxes, weighted by appropriate Bose-Einstein 
factors, to obtain the total heat-transfer rate or NEQ Casimir force. 
(For more details, see 
[What <span class="SC">scuff-neq</span> actually computes](#WhatItComputes).)

+ As in [[scuff-cas3d]], you can specify an optional list of 
[geometrical transformations](../../reference/Transformations.md) 
describing various displacements and rotations of the bodies 
in your geometry; in this case you will get back values of the 
frequency-resolved or frequency-integrated quantities for each 
transformation you specify.

A bonus feature of [[scuff-neq]] that is **not** present in 
[[scuff-cas3d]] is the ability to obtain spatially-resolved 
information on energy and momentum fluxes. More specifically, 
you can specify to [[scuff-neq]] a list of evaluation points, 
and you will get back values of the (thermally and 
temporally averaged) Poynting flux and Maxwell stress tensor 
at each point you requested.

In addition to numerical output on heat-transfer rates and
Casimir quantities, you can also request visualization outputs
that plot the spatial distribution of the heat or momentum flux.

For Casimir forces, the quantities computed by [[scuff-neq]] are 
only the *non-equilibrium* contributions to the total force---that 
is, the contributions to the force arising from the temperature 
differences between individual bodies and the surrounding medium.
To get the total force, these must be added to the equilibrium 
contributions, which are the Casimir forces for the case in which 
all bodies are at the temperature of the surrounding medium. 
These contributions must be computed by doing a separate 
[[scuff-cas3d]] calculation, on the same geometry, at the 
temperature of the external medium. 
(For heat-transfer rates there is of course no equilibrium 
contribution, as there is no net power transfer between 
bodies at thermal equilibrium.)

One difference between [[scuff-cas3d]] and [[scuff-neq]] is that,
whereas [[scuff-cas3d]] reports only the Casimir force on one body in 
a geometry (namely, the first body listed in the `.scuffgeo' file), 
[[scuff-neq]] reports forces and heat-transfer rates for *all* bodies
in the geometry. The extra information would typically be redundant 
in an equilibrium Casimir calculation, since the equilibrium Casimir 
force on the second body (in a two-body geometry) is just equal and 
opposite to the force on the first body; but in general no such 
relation holds in the non-equilibrium case.

[TOC]

# What <span class="SC">scuff-neq</span> actually computes

# <span class="SC">scuff-neq</span> command-line options

# <span class="SC">scuff-neq</span> output files

[scuff-cas3d]: ../scuff-cas3d/scuff-cas3d.md
