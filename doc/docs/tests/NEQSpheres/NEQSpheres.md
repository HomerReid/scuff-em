# Non-equilibrium fluctuation-induced interactions between spheres: heat radiation and non-equilibrium Casimir forces

This test validates the [[scuff-neq]] application module
of the [[scuff-em]] code suite by using it to compute
**(a)** the temperature-dependent rate of heat radiation
from individual spheres,
**(b)** radiative heat-transfer rates and non-equilibrium
Casimir forces between spheres at various separation 
distances.

## Analytical solution

Analytical formulas for single-sphere heat radiation
and sphere-sphere heat-transfer rates and non-equilibrium 
Casimir forces were obtained by M. Krueger and are
discussed in this paper:

* [M. Krueger, G. Bimonte, T. Emig, and M. Kardar, "Trace formulas for nonequilibrium {C}asimir interactions, heat radiation, and heat transfer for arbitrary objects", Physical Review B \textbf{86} 115423 (2012)]({http://link.aps.org/doi/10.1103/PhysRevB.86.115423)

The formulas are a bit complicated to reproduce here, but here
is a simple [<span class="SC">julia</span>](http://julialang.org)
code that implements them:

* [`KruegerFormulas.jl`](KruegerFormulas.jl)

This code may be used to

## <span class="SC">scuff-em</span> solution

Values of the flux quantity $\phi^{\text{\small 
heat radiation
and between two dielectric spheres may be computed
using [[scuff-cas3d]] as follows:

````bash
 % scuff-neq --geometry SiO2Sphere_501.scuffgeo --OmegaFile OmegaFile
````

````bash
 % scuff-neq --geometry SiO2Spheres_501.scuffgeo --OmegaFile OmegaFile
````

Here the two `.scuffgeo` files 
([`PECSpheres_501.scuffgeo`](PECSpheres_501.scuffgeo) and [`E10Spheres_501.scuffgeo`](E10Spheres_501.scuffgeo)]
describe the two geometric configurations
(two PEC spheres and two dielectric spheres of radius $R=1\, \mu$m 
separated by an initial center-center distance of $d$=3 $\mu$m) while
[`Spheres.trans`](Spheres.trans) specifies the list of center-center 
separation distances $d$ at which we compute the energy and force.
(Both geometries refer to the same surface mesh file for the 
sphere, [`Sphere_327.msh`](Sphere_327.msh).

The above calculations produce output files named 
[`PECSpheres_327.out`](PECSpheres_327.out) and 
[`E10Spheres_327.out`](E10Spheres_327.out). Plotting against
the theoretical predictions of Emig et. al (referenced above)
yields good agreement:

![CasimirSphereData.png](CasimirSphereData.png).

Here's the [[gnuplot]] script I used to produce this 
plot: [`Plotter.gp`](Plotter.gp).

[EmigPaper]:                          http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.170403

{!Links.md!}
