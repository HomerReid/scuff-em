<h1> Modeling non-equilibrium electromagnetic fluctuations with
     <span class="SC">scuff-neq</span>
</h1>

[[scuff-neq]] is an application code in the [[scuff-em]] suite for 
studying non-equilibrium (NEQ) electromagnetic-fluctuation-induced 
phenomena--specifically, for computing *radiative heat-transfer rates* 
and *non-equilibrium Casimir forces and torques* for bodies of 
arbitrary shapes and arbitrary (linear, isotropic, piecewise 
homogeneous) frequency-dependent permittivity and permeability.
[[scuff-cas3d]] implements the 
[*fluctuating-surface current (FSC)* approach][FSCPaper] 
to numerical modeling of non-equilibrium fluctuation phenomena.

Mechanically, working with [[scuff-neq]] is similar in many ways to 
working with the equilibrium Casimir code [scuff-cas3d][scuff-cas3D]{.SC}. 
In particular,

+ As in [[scuff-cas3d]], you can request either **(a)** frequency-resolved 
information on heat-transfer rates and NEQ Casimir forces (in which case 
you will specify a list of frequencies and will get back a list of 
frequency-specific values of energy and momentum fluxes) or 
**(b)** frequency-integrated information, in which case you will assign 
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
you can specify to [[scuff-neq]] a 
[list of evaluation points][EPFile]
and you will get back values of the (thermally and 
temporally averaged) Poynting flux and Maxwell stress tensor 
at each point you requested.

In addition to numerical output on heat-transfer rates and
Casimir quantities, you can also request visualization outputs
that plot the spatial distribution of the heat or momentum flux.

For Casimir forces and torques, the quantities computed by 
[[scuff-neq]] are only the *non-equilibrium* contributions to 
the total force and torque---that
is, the contributions arising from the temperature 
differences between individual bodies and the surrounding environment.
To get the total force, these must be added to the *equilibrium*
contributions, which are the Casimir forces and torques for the 
case in which all bodies are at the temperature of the environment.
These contributions must be computed by doing a separate 
[[scuff-cas3d]] calculation, on the same geometry, at the 
temperature of the external medium. 
(For heat-transfer rates there is of course no equilibrium 
contribution, as there is no net power transfer between 
bodies at thermal equilibrium.)

One difference between [[scuff-cas3d]] and [[scuff-neq]] is that,
whereas [[scuff-cas3d]] reports only the Casimir force on one body in 
a geometry (namely, the first body listed in the `.scuffgeo` file), 
[[scuff-neq]] reports forces and heat-transfer rates for *all* bodies
in the geometry. <span>[</span>The extra information would typically 
be redundant in an equilibrium Casimir calculation, since the equilibrium 
Casimir force on the second body (in a two-body geometry) is just equal
and opposite to the force on the first body; but in general no such 
relation holds in the non-equilibrium case.<span>]</span>

In fact, the output from [[scuff-neq]] is even more detailed
than that: in addition to the total power/force/torque on
each body, you also get the contributions of each individual
source body to those quantities. All of this means that the
output from [[scuff-neq]] requires some effort to interpret,
as discussed in more detail below.

[TOC]

# 1. What <span class="SC">scuff-neq</span> actually computes

[[scuff-neq]] implements the
[FSC approach to non-equilibrium fluctuation phenomena][FSCPaper],
an algorithm for computing the thermal averages of power,
force, and torque (PFT) quantities in geometries consisting
of homogeneous material bodies at various temperatures embedded
in an finite-temperature or zero-temperature environment.
[[scuff-neq]] can compute both
spatially-*resolved* and spatially-*integrated* PFT quantities.
(Examples of spatially-resolved quantities include components
of the average Poynting flux or Maxwell stress tensor at individual
points in space. Examples of spatially-integrated quantities 
include the total power absorbed by, or the total force or torque 
on, a compact homogeneous body. Spatially-integrated quantities 
are generally obtained by integrating spatially-resolved quantities
over closed bounding surfaces, although this is not necessarily
the way they are computed by [[scuff-neq]].)

In general, for a geometry consisting of multiple homogeneous
bodies, PFT quantities receive contributions from source
fluctuations in all bodies and at all frequencies, and the
the thermal average of a PFT quantity *Q* may be
written in the form

$$ \big\langle Q\big\rangle
   =\int_0^\infty 
      \left\{ \sum_s 
       \Big[ \Theta_s \left( T_s, \omega \right)
             -
             \Theta_s \left( T_{\hbox{\scriptsize{env}}}, \omega
                      \right)
       \Big] \Phi_{s}(\omega\right) \right\} d\omega
$$

<a name="CommandLineOptions"></a>
# 2. <span class="SC">scuff-neq</span> command-line options

## Common options

[[scuff-neq]] recognizes the following subset of the 
[list of commonly accepted options to <span class="SC">scuff-em</span> command-line codes][CommonOptions].

````
`--geometry`
`--TransFile`
` `
`--Omega`
`--OmegaFile`
`--OmegaQuadrature`
`--OmegaMin`
` `
`--AbsTol`
`--RelTol`
` ` 
`--FileBase`
` `    
`--Cache`
`--ReadCache`
`--WriteCache`
````
{.toc}

## Options requesting output quantities

````
`--PAbs`  
`--PRad`  
`--XForce`  
`--YForce`  
`--ZForce`  
`--XTorque`  
`--YTorque`  
`--ZTorque`  
````{.toc}

Specifies the quantities in which you are interested:
absorbed power (`--PAbs`), radiated power (`--PRad`),
Cartesian force components, or Cartesian torque components.
You may specify none, all, or any subset of these options,
but each option you specify will generally increase
the computation time (you can scrutinize the
[`.log` file](#LogFile) to see how *much* additional time each
extra output quantity takes to compute).

## Option requesting visualization output

+ `--PlotFlux`
{.toc}

## Options specifying object temperatures

````
--Temperature UpperSphere 300
--Temperature LowerSphere 100

--Temperature ENVIRONMENT 100
````

> The first two options here set the temperatures
> of the objects labeled `UpperSphere` and
> `LowerSphere` in the `.scuffgeo` file.
> **Temperature specifications are interpreted in 
> units of Kelvin**, so `300` corresponds to 
> room temperature.
>
> The third option here sets the temperature of
> the environment in which the objects are embedded.
> (The keywords `MEDIUM` and `EXTERIOR` may be used
> here interchangeably with `ENVIRONMENT`).
>
> Note that the temperatures of all objects, and of
> the environment, are 0 by default. This means that,
> if you request a full frequency-integrated calculation
> (which you do by omitting the `--omega` or `--omegaFile`
> option) and you do not specify any `--temperature` 
> options, the code will chug for a while (computing 
> temperature-independent fluxes at various frequencies)
> before reporting strictly zero values for all  
> quantities! This is probably not what you want.

## Options controlling the computation of power, force, and torque

`--DSIPoints 302`
`--DSIRadius 5.0`
`--DSIMesh BoundingMesh.msh`
`--DSIFarField`
`--ForceDSI`
{.toc}

## Other options

`--OmitSelfTerms`
{.toc}

> This h

--------------------------------------------------

# <span class="SC">scuff-neq</span> output files

<a name="#LogFile"></a>
## The `.log` file 

Like all command-line codes in the [[scuff-em]] suite,
[[scuff-cas3d]] writes a [`.log` file][LogFiles] that you
can monitor to keep track of your calculation's progress.

## Output files for spatially-integrated PFTs: The `.SIFlux`, `.SIIntegrand`, and `.NEQPFT` files

If you requested the computation of any spatially-integrated
PFTs (by setting command-line options such as `--PAbs` or `--YForce`),
you will get back

To understand what is written to these files, let $Q_d$ be
the spatially-integrated PFT on a destination body $d$,
and write the FSC decomposition of the thermal average
of $Q_d$ in the form

$$ \big\langle Q_d\big\rangle
   = 
   \underbrace{ 
    \Bigg[ \int_0^\infty
     \underbrace{ 
      \bigg\{ \hbar\omega_0^2 \sum_s \, \widehat \Theta_s(u)
       \underbrace{ \Phi_{s\to d}(u)}_{\texttt{.SIFlux}}
      \bigg\}
                }_{\texttt{.SIIntegrand}}
    \,\,du \Bigg]
              }_{\texttt{.NEQPFT}}
$$

## Output files for spatially-resolved PFTs: The `.SRFlux`, `.SRIntegrand`, and `.PVMST` files

<a name="Examples"></a>
# 4. Examples of calculations using <span class="SC">scuff-neq</span>

+ [Heat radiation from a warm sphere in a cold environment][SiO2Spheres]

+ [Heat transfer and non-equilibrium Casimir forces between warm and cold spheres][SiO2Spheres]

+ [Spatial distribution of poynting flux from a warm tip above a cold substrate][SiO2Spheres]

[scuff-cas3D]:  ../scuff-cas3D/scuff-cas3D.md
[EPFile]:       ../../applications/GeneralReference.md#EvaluationPoints
[FSCPaper]:     http://link.aps.org/doi/10.1103/PhysRevB.88.054305
[LogFiles]:     ../GeneralReference.md#LogFiles
[SiO2Sphere]:   ../../examples/SiO2Spheres/SiO2Spheres.md#SingleSphere
[SiO2Spheres]:   ../../examples/SiO2Spheres/SiO2Spheres.md#TwoSpheres
[TipSubstrate]: ../../examples/TipSubstrate/TipSubstrate.md
