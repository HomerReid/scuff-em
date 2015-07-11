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
material regions, PFT quantities receive contributions from source
fluctuations in all regions and at all frequencies, and the
the thermal average of a PFT quantity *Q* may be
written in the form

$$ \big\langle Q\big\rangle
    = \int_0^\infty \, \sum_r \, \Theta(T_r,\omega) \Phi_r(\omega)\,d\omega
$$
where $T_r$ is the temperature of region $r$, $\Phi_r(\omega)$ is
a temperature-independent *generalized flux* describing the 
contribution of frequency-$\omega$ source fluctuations in region $r$,
and 
$$\Theta(T_r,\omega) = \frac{\hbar\omega}{e^{\hbar \omega/kT_r} - 1}$$
is the Bose-Einstein factor.

The sum over regions $r$ in this equation includes the
contributions of the external environment. To isolate these
contributions it is convenient to decompose $\langle Q \rangle$
into a sum of two terms:
$$ \begin{array}{rcl}
 \big\langle Q\big\rangle
&=&
 \big\langle Q\big\rangle^{\small EQ} + 
   \big\langle Q\big\rangle^{\small NEQ}
\\[8pt]
 \big\langle Q\big\rangle^{\small EQ}
&\equiv&
  \displaystyle{\int_0^\infty} \Theta(T_{\small env},\omega) 
                \sum_r \, \Phi_r(\omega)\,d\omega 
\\[8pt]
 \qquad \big\langle Q\big\rangle^{\small NEQ} 
&\equiv& 
  \displaystyle{\int_0^\infty \sum_s }
  \Big[ \Theta(T_s, \omega) - \Theta(T_{\small env},\omega)\Big]
        \Phi_s(\omega)\,d\omega 
\\[4pt]
&=&
  \displaystyle{\int_0^\infty \sum_s}
  \Delta \Theta(T_s, \omega) \Phi_s(\omega)\,d\omega 
\end{array}
$$

where 

$$ \Delta \Theta(T_s, \omega) \equiv 
   \Theta(T_s, \omega) - \Theta(T_{\small env},\omega).
$$

The quantity $\langle Q\rangle^{\small EQ}$
is the average value of $Q$ that would obtain if 
the temperature in all material regions were equal
to the environment temperature $T_{\small env}$---that is,
it is the *equilibrium* value of $\langle Q\rangle$
at temperature $T_{\small env}$. The equilibrium value
of PFT quantities may be computed by methods that are 
less costly than [[scuff-neq]]. (For example,
if $Q$ is a spatially-integrated force or torque, then 
$\langle Q\rangle^{\small EQ}$ is just the equilibrium
Casimir force, which is computed efficiently by 
[scuff-cas3d][scuff-cas3D]{.SC}.
On the other hand, if $Q$ is a spatially-integrated
power transfer quantity, then 
$\langle Q\rangle^{\small EQ}=0$ identically.)
Thus this contribution is not computed by [[scuff-neq]].

The quantity $\langle Q\rangle^{\small NEQ}$
is the extent to which $\langle Q\rangle$
*deviates* from its equilibrium value, and
the sum in its definition ranges only over 
the source bodies in the geometry, not including
the environment contribution.
$\langle Q \rangle^{\small NEQ}$ is the
quantity that is computed by [[scuff-neq]].

<a name="CommandLineOptions"></a>
# 2. <span class="SC">scuff-neq</span> command-line options

### Common options

[[scuff-neq]] recognizes the following subset of the 
[list of commonly accepted options to <span class="SC">scuff-em</span> command-line codes][CommonOptions].

 
  ````
--geometry
--TransFile
--Omega
--OmegaFile
--OmegaQuadrature
--OmegaMin
--AbsTol
--RelTol
--FileBase
--Cache
--ReadCache
--WriteCache
  ````
{.toc}

### Options requesting output quantities

  ````
--PAbs  
--PRad  
--XForce  
--YForce  
--ZForce
--XTorque  
--YTorque
--ZTorque
  ````
{.toc}

Specifies the quantities in which you are interested:
absorbed power (`--PAbs`), radiated power (`--PRad`),
Cartesian force components, or Cartesian torque components.
You may specify none, all, or any subset of these options,
but each option you specify will generally increase
the computation time (you can scrutinize the
[`.log` file](#LogFile) to see how *much* additional time each
extra output quantity takes to compute).

### Option requesting visualization output

  ````
--PlotFlux
  ````
{.toc}

This option directs [[scuff-neq]] to produce visualization
files (in addition to its usual output files) which may be 
opened in [[gmsh]] to visualize,
for each spatially-integrated PFT quantity you requested,
the spatial distribution of the Poynting flux or 
Maxwell stress on the surfaces of objects or 
the displaced bounding surfaces over which those quantities
are integrated to compute the total PFT quantity.

### Options specifying object temperatures

  ````
--Temperature UpperSphere 300
--Temperature LowerSphere 100
--Temperature ENVIRONMENT 100
  ````
{.toc}

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
> the environment, are zero by default. This means that,
> if you request a full frequency-integrated calculation
> (which you do by omitting the `--omega` or `--omegaFile`
> option) and you do not specify any `--temperature` 
> options, the code will chug for a while (computing 
> temperature-independent fluxes at various frequencies)
> before reporting strictly zero values for all
> quantities! This is probably not what you want.

### Options controlling the computation of power, force, and torque

  ````
--ForceDSI
  ````
{.toc}


  ````
--DSIPoints 302
--DSIRadius 5.0
--DSIMesh BoundingMesh.msh
--DSIFarField
  ````
{.toc}

As detailed in [this paper][SIEPFTPaper], there are several
ways to compute PFTs in surface-integral formulations,
including the "displaced-surface-integral" (DSIPFT),
"equivalence principle" (EPPFT), and "overlap" (OPFT) 
methods.

By default, [[scuff-neq]] uses different algorithms for 
different cases of the PFT computation:

* Power computation (self term): EPPFT
* Power computation (non-self terms): EPPFT
* Force/Torque computation (self term): DSIPFT
* Force/Torque computation (non-self term): OPFT

However, you can override this default behavior 
by specifying `--ForceDSI`, in which case DSIPFT
will be used in all cases.

The other options here set specific parameters 
that are only used for DSIPFT calculations.

`--DSIMesh` specifies the name of a mesh file
(such as a [[GMSH]]-produced `.msh` file) that 
defines the bounding surface over which the 
surface integrals are computed. (This surface is
automatically displaced and rotated appropriately
for each object in accordance with any
geometrical transformations that may be 
specified in the `.scuffgeo` file and/or the
`.trans` file.) The surface integral is evaluated
via a one-point cubature over the surface of 
each panel in the bounding mesh; thus, the
finer the bounding mesh, the more accurate and
the more expensive the computation).

If you do not specify a `--DSIMesh`, then the
surface integral for each object is computed 
using Lebedev cubature over a bounding sphere 
centered at the origin of coordinates of the 
object mesh (appropriately displaced if the 
object has been displaced via statements in 
the `.scuffgeo` or `.trans` files.) In this case,
you may use `--DSIRadius` and `--DSIPoints` to
set the radius of this sphere (in microns) and 
the number of Lebedev cubature points (the more 
points, the more accurate and expensive the calculation).
To see the possible values that may be specified 
for `--DSIPoints,` type `scuff-neq --help.`

Finally, you may use `--DSIFarField` to request
that the Poynting vector and Maxwell tensor
on the bounding surface be computed using only
the far-field (radiation-zone) contributions 
of the surface currents to the fields.

### Other options

  ````
--OmitSelfTerms
  ````
{.toc}

Omit the contributions of sources in individual bodies
to the total PFTs on those bodies themselves.

--------------------------------------------------

# 3. <span class="SC">scuff-neq</span> output files

<a name="#LogFile"></a>
### The `.log` file 

Like all command-line codes in the [[scuff-em]] suite,
[[scuff-cas3d]] writes a [`.log` file][LogFiles] that you
can monitor to keep track of your calculation's progress.

### Output files for spatially-integrated PFTs: The `.SIFlux`, `.SIIntegrand`, and `.NEQPFT` files

If you requested the computation of any spatially-integrated
PFTs (by setting command-line options such as `--PAbs` or `--YForce`),
you will get back files reporting various contributions to 
these quantities.
To understand what is written to these files, let $Q_d$ be
the spatially-integrated PFT on a destination body $d$,
and write the FSC decomposition of the thermal average
of $Q_d$ in the form

$$ \big\langle Q_d\big\rangle
   = 
   \underbrace{ 
    \Bigg[ \int_0^\infty
     \underbrace{ 
      \bigg\{ \hbar\omega_0^2 \sum_s \, \Delta \widehat \Theta_s(u)
       \underbrace{ \Phi_{s\to d}(u)}_{\texttt{.SIFlux}}
      \bigg\}
                }_{\texttt{.SIIntegrand}}
    \,\,du \Bigg]
              }_{\texttt{.NEQPFT}}
$$

In this equation,

+ $u$ is a dimensionless frequency variable: $u=\omega/\omega_0$, 
where $\omega_0=3\cdot 10^{14}$ rad/sec. (Thus $u$ agrees numerically
with the arguments to the `--omega` option.)

+ $\widehat\Theta(u)$ is a dimensionless version of the usual
Bose-Einstein factor, defined by $\Theta(\omega)=\hbar \omega_0 \widehat\Theta(u)$.

+ $\Delta \widehat \Theta_s(u)=\widehat \Theta_s(u) - \widehat \Theta_{\small env}(u)$ 
is the difference between the dimensionless Bose-Einstein factors of source body $s$ 
and the environment.

+ $\Phi_{s\to d}$ is the temperature-independent generalized flux
describing the contributions of fluctuations in source body $s$
to the power, force, or torque on destination body $d$.
Values of this quantity are written to the `.SIFlux` output
file.

+ $\hbar\omega_0^2 \Delta \widehat \Theta_s(u) \Phi_{s\to d}(u)$ is
the spectral density of temperature-weighted contributions
from fluctuations in source body $s$ to the PFT on destination
body $d$. 
Values of this quantity are written to the `.SIIntegrand` output
file.

+ Finally, $\langle Q_d \rangle$ is the total thermally-averaged
PFT on body $d$. Values of this quantity are written to the 
`.NEQPFT` output file.

In all of these files, each single line corresponds to
a single frequency, a single geometric transformation,
and a single pair of (source,destination) objects.

At the top of each output file you will find a file header
explaining the significance of each of the various
columns in the file. One of the columns will be
described in the header as `# (source object, dest object),`
and will take values like `12`, `22`, or `02.`
The first case (`12`) indicates that the data on that
line correspond to the contributions of object 1 to the 
PFT on object 2. (The ordering of objects corresponds
with the order of their appearance in the `.scuffgeo`
file).  The second case (`22`) indicate that the data
on that line correspond to the self-contributions of 
object 2 to its own PFT. The third case (`02`) 
indicates that the data on that line correspond to 
the *total* PFT on object 2---that is, the sum of 
contributions from all source objects.

### Output files for spatially-resolved PFTs: The `.SRFlux`, `.SRIntegrand`, and `.PVMST` files

If you requested the computation of spatially-resolved
power and momentum flux (by specifying the `--EPFile` 
command-line option), you will get back files reporting 
various contributions to these quantities. The breakdown
here is similar to that described above for spatially-integrated
quantities. To understand this, let $Q(\mathbf{x})$ be
a spatially-resolved PFT quantity (a component 
of the Poynting vector or Maxwell stress tensor)
at a point $\mathbf{x}$. Then the thermal average of $Q$ 
may be written in the form

$$ \big\langle Q(\mathbf{x})\big\rangle
   = 
   \underbrace{ 
    \Bigg[ \int_0^\infty
     \underbrace{ 
      \bigg\{ \hbar\omega_0^2 \sum_s \, \Delta \widehat \Theta_s(u)
       \underbrace{ \Phi_{s\to\mathbf x}(u)}_{\texttt{.SRFlux}}
      \bigg\}
                }_{\texttt{.SRIntegrand}}
    \,\,du \Bigg]
              }_{\texttt{.PVMST}}
$$

In this equation,

+ $\Phi_{s\to \mathbf{x}}$ is the temperature-independent generalized flux
describing the contributions of fluctuations in source body $s$
to the Poynting flux or Maxwell stress at $\mathbf{x}$.
Values of this quantity are written to the `.SRFlux` output
file.

+ $\hbar\omega_0^2 \Delta \widehat \Theta_s(u) \Phi_{s\to \mathbf{x}}(u)$ is
the spectral density of temperature-weighted contributions
from fluctuations in source body $s$ to the Poynting flux or 
Maxwell stress at $\mathbf{x}$. Values of this quantity are 
written to the `.SRIntegrand` output
file.

+ Finally, $\langle Q(\mathbf{x})\rangle$ is the total thermally-averaged
Poynting vector or Maxwell stress tensor at $\mathbf{x}$.
Values of this quantity are written to the `.PVMST` output file.

In all of these files, each line corresponds to
a single frequency, a single geometric transformation,
and a single source object.
At the top of each output file you will find a file header
explaining how to interpret the various data columns
on each line.

### Units of output quantities

* The units of the total (frequency-integrated)
spatially-integrated output quantities reported in
the `.NEQPFT` file are *watts* for power, *nanoNewtons*
for force, and *nanoNewtons $\times$ microns* for torque.

* The quantities in the `.SIIntegrand` output file
are the PFT quantities per unit *dimensionless* frequency,
so have the same units as the corresponding quantities
in the `.NEQPFT` file.

* The quantities in the `.SIFlux` output file are
the quantities per unit dimensionless frequency
per watt of thermal energy, so these quantities
have the same units as the quantities in 
the `.NEQPFT` and `.SIIntegrand` file, but divided
by watts: thus the power flux is *dimensionless*,
the force flux has units of *nanoNewtons / watts*,
and the torque flux has units of *nanoNewtons microns/watts.*

<a name="Examples"></a>
# 4. Examples of calculations using <span class="SC">scuff-neq</span>

+ [Heat radiation from a warm sphere in a cold environment][SiO2Spheres]

+ [Heat transfer and non-equilibrium Casimir forces between warm and cold spheres][SiO2Spheres]

+ [Spatial distribution of poynting flux from a warm tip above a cold substrate][TipSubstrate]

[scuff-cas3D]:   ../scuff-cas3D/scuff-cas3D.md
[EPFile]:        ../../applications/GeneralReference.md#EvaluationPoints
[FSCPaper]:      http://link.aps.org/doi/10.1103/PhysRevB.88.054305
[SIEPFTPaper]:	 http://dx.doi.org/10.1109/TAP.2015.2438393
[LogFiles]:      ../GeneralReference.md#LogFiles
[SiO2Sphere]:    ../../examples/SiO2Spheres/SiO2Spheres.md
[SiO2Spheres]:   ../../examples/SiO2Spheres/SiO2Spheres.md
[TipSubstrate]:  ../../examples/TipSubstrate/TipSubstrate.md
[CommonOptions]: ../GeneralReference.md#CommonOptions
