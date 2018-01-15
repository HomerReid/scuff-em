<h1> Modeling non-equilibrium electromagnetic fluctuations with
     <span class="SC">scuff-neq</span>
</h1>

[[scuff-neq]] is an application code in the [[scuff-em]] suite for 
studying non-equilibrium (NEQ) electromagnetic-fluctuation-induced 
phenomena--specifically, for computing *radiative heat-transfer rates* 
and *non-equilibrium Casimir forces and torques* for bodies of 
arbitrary shapes and arbitrary (linear, isotropic, piecewise 
homogeneous) frequency-dependent permittivity and permeability.

[[scuff-neq]] implements the
[*fluctuating-surface current (FSC)* approach][FSCPaper] approach
to numerical modeling of non-equilibrium fluctuation phenomena.
The FSC method was first developed to compute *equilibrium*
Casimir forces (as implemented by [<span class=SC>scuff-cas3d</span>][scuff-cas3D])
[[1]](scuff-neq.md#bibliography)
and was subsequently extended to handle non-equilibrium
radiative heat-transfer problems 
[[2]](scuff-neq.md#bibliography)
and later non-equilibrium
forces and torques
[[3]](scuff-neq.md#bibliography).

As discussed in more detail below,
[what <span class="SC">scuff-neq</span> actually computes](#WhatItComputes)
are temperature-independent quantities known as *generalized fluxes*
that describe frequency-resolved contributions to various quantities of
interest in fluctuational electromagnetism, 
including both **(a)** *spatially-integrated* (SI) quantities such as
the total rate of heat absorbed or radiated by, or the total force
or torque on, a material body, and **(b)** *spatially-resolved* (SR)
quantities such as the Poynting flux or Maxwell stress tensor at
individual points in space. [[scuff-neq]] outputs data files
giving values of generalized fluxes at various frequencies you specify;
these data may then be integrated by the separate utility 
code [<span class=SC>scuff-integrate</span>][scuffIntegrate]
to compute full thermally-averaged values of SI and SR quantities
for various body temperatures.

When requesting spatially-integrated quantities, you will
specify which of the various computational techniques you
want the code to use to compute total energies and forces/torques
on bodies; the main options are the *displaced surface-integral* (DSI)
method, which integrates the Poynting vector or Maxwell stress
tensor over a bounding surface surrounding the body,
or the *energy-momentum transfer* (EMT) method, which
integrates the local power absorption and force over the 
volume of the body. (Algorithms for spatially-integrated
PFT calculation are discussed in 
Refs. [ [3] and [4] ](scuff-neq.md#bibliography)
below.


When requesting spatially-resolved quantities, you will
specify a [list of evaluation points][EPFile], and you
will get back a data file giving (generalized fluxes for)
the Poynting vector and Maxwell stress tensor at each point.

[[scuff-neq]] also offers options for generating *visualization*
files showing the spatial distribution of energy and momentum
absorption by bodies.

As in [[scuff-cas3d]], you can specify an optional list of 
[geometrical transformations](../../reference/Transformations.md) 
describing various displacements and rotations of the bodies 
in your geometry; in this case you will get back multiple 
copies of all output quantities, one for each transformation
in your list.

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
whereas in [[scuff-cas3d]] you use command-line options like 
`--zForce` or `--xTorque` to request computation of specific 
SI quantities, in [[scuff-neq]] you automatically get all 8
SI quantities (heat radiated, heat absorbed, plus all 3
components of the force and torque). 

Another difference between [[scuff-cas3d]] and [[scuff-neq]] is that,
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

Consider a collection of one or more homogeneous material
bodies, each maintained at a given temperature, embedded
in an finite- (or zero-) temperature environment.
The electromagnetic fields radiated by thermally-induced
fluctuating sources in each body carry energy and momentum,
which we may characterize either in terms 
of **(a)** spatially-resolved (SR) quantities like the 
time-average Poynting vector or Maxwell stress tensor
at given points in space, or **(b)** spatially-integrated (SI)
quantities like the total time-average heat absorbed or
radiated by, or force or torque on, entire bodies.
In the 
[FSC approach to non-equilibrium fluctuation phenomena][FSCPaper]
implemented by [[scuff-neq]],
the total thermally-averaged value of any SR or SI quantity $Q$
is computed by summing the contributions of fluctuations from
all bodies at all frequencies:

$$ \big\langle Q\big\rangle
    = \int_0^\infty \, \sum_b \, \Theta(T_b,\omega) \Phi_b(\omega)\,d\omega
$$
where $T_b$ is the temperature of body$b$, $\Phi_b(\omega)$ is
a temperature-independent *generalized flux* describing the 
contribution of frequency-$\omega$ source fluctuations in body $b$,
and 
$$\Theta(T_b,\omega) = \frac{\hbar\omega}{e^{\hbar \omega/kT_b} - 1}$$
is the Bose-Einstein factor.

The sum over bodies $b$ in this equation includes the
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
                \sum_b \, \Phi_b(\omega)\,d\omega 
\\[8pt]
 \qquad \big\langle Q\big\rangle^{\small NEQ} 
&\equiv& 
  \displaystyle{\int_0^\infty \sum_b }
  \Big[ \Theta(T_b, \omega) - \Theta(T_{\small env},\omega)\Big]
        \Phi_b(\omega)\,d\omega 
\\[4pt]
&=&
  \displaystyle{\int_0^\infty \sum_b}
  \Delta \Theta(T_b, \omega) \Phi_b(\omega)\,d\omega 
\end{array}
$$

where 

$$ \Delta \Theta(T_b, \omega) \equiv 
   \Theta(T_b, \omega) - \Theta(T_{\small env},\omega).
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

The job of [[scuff-neq]] is to compute
the quantity $\Phi_b(\omega)$---known as a *generalized flux*---
that enters the integral defining $\langle Q\rangle^{\small NEQ}$;
[[scuff-neq]] will do this computation at each of multiple
frequencies that you specify in advance. These data are then used by 
the separate [[scuff-integrate]] utility to evaluate
the actual $\omega$ integrals and compute thermally-averaged 
power, force, and torque quantities at various temperatures.

<a name="CommandLineOptions"></a>
# 2. <span class="SC">scuff-neq</span> command-line options

### Options specifying the geometry

````--geometry  MyGeometry.scuffgeo
--TransFile MyTransFile
````

{.toc}

The mandatory `--geometry` option specifies the
[<span class=SC>scuff-em</span> geometry file][Geometries]
describing your geometry.

The optional `--TransFile` option specifies a
[list of geometrical transformations][Transformations]
to be applied to your geometry. Each output quantity you
request will be separately computed and reported for 
each transform you request.

### Options specifying frequencies

````--omega 1.23
--OmegaFile MyOmegaFile
````

The first option here requests computation at just a
single angular frequency (interpreted in 
[standard <span class=SC>scuff-em</span> units][Units]
of $\omega_0=3\cdot 10^{14}$ rad/sec).

The second option specifies a file containing a list
of frequencies $\omega$ (one per line) at which to compute
generalized fluxes $\Phi(\omega)$.

### Options requesting spatially-integrated output quantities

````--EMTPFT
--DSIPFT
--OPFT   
````

Request calculation of spatially-integrated quantities
via one of the available methods: 
 the energy/momentum transfer (EMT) method,
 the displaced-surface-integral (DSI) method,
 or the overlap (O) method.
You may specify none, any combination, or all of these options.

Note that all of these methods are computing the same quantities,
so in principle their results should be equivalent, and specifying
more than one should be redundant. However, in practice it is
often helpful to compute SI quantities via multiple methods
as a confirmation and accuracy check.

For the particular case of `--DSIPFT,` there are several
options controlling the particular implementation of the 
algorithm; see below.

### Option requesting spatially-resolved output quantities

````
--EPFile MyEPFile
````

Specifies a file containing a [list of evaluation points][EPFile],
one per line (three Cartesian coordinates),
at which values of the Poynting vector and Maxwell stress tensor
are to be computed.

#### Options controlling DSI calculations

The displaced-surface-integral (DSI) method
total absorbed/radiated power,
force, and torque on each body in a geometry by integrating
the Poynting flux and Maxwell stress over a bounding surface
surrounding the body. You have two options for *what* bounding
surface is used: 
**(a)** a [<span class=SC>gmsh</span>](gmsh.info) mesh file you provide
giving a discretized version of your bounding surface, or 
**(b)** a sphere centered at the body. (In both cases,
the surface is displaced and re-centered at the origin of 
each body to compute DSIPFT for that body.)

To use option **(a)**, say

```
 --DSIMesh    MyDSIMesh.msh
```

If the `--DSIMesh` option is not specified, [[scuff-neq]]
defaults to a spherical bounding surface. In this case you may 
use the command-line options

```
 --DSIRadius  2.0
 --DSIPoints  302
```

to set the radius of this sphere and the number of cubature
points used to approximate the surface integral over it.
(Use `scuff-neq --help` to see a list of the available values
of `--DSIPoints`). If you don't specify values for these
options, they will be set automatically.

```
 --DSIPoints2 770
```

For the case of a spherical bounding surface, you may specify
`-DSIPoints2` to request a "second opinion" DSI calculation
with a different number of cubature points to check the 
accuracy of the numerical cubature.  The second-opinion DSI
data will be written to a file with file extension `.SIFlux.DSI2`
(This is optional; if you
don't specify `--DSIPoints2` no second-opinion calculation
is done.)

```
 --DSIOmegaFile DSIOmegaFile
```

Because the DSI calculation can be slow, in some cases you
may wish to request that it only be performed at a subset
of the frequencies in your `--OmegaFile`. You specify
that subset using the `-DSIOmegaFile` option. For example,
you might say `--EMTPFT` to request EMT calculations at 
all frequencies, but use `--DSIPFT` with a `--DSIOmegaFile`
containing only every 5th or 10th frequency in your full 
`OmegaFile` to obtain a sanity check on the EMTPFT results.

### Option requesting visualization output

````
--PlotFlux
````

This option directs [[scuff-neq]] to produce visualization
files (in addition to its usual output files) which may be 
opened in [[gmsh]] to visualize,
for each spatially-integrated PFT quantity you requested,
the spatial distribution of the Poynting flux or 
Maxwell stress on the surfaces of objects or 
the displaced bounding surfaces over which those quantities
are integrated to compute the total PFT quantity.

### Other options

````
--FileBase MyFileBaSe
````

Specifies the base filename for all output files.
If this option is not specified, `--FileBase` is set
to the base filename of the `.scuffgeo` file.

````
--SourceObject MySourceObject
--DestObject   MyDestObject 
````

For a geometry containing $N$ bodies,
by default <span class=SC>scuff-neq</span> computes
the full $N\times N$ matrix of generalized flux
quantities $Q_{s\to d}$ describing the
energy/momentum transfer to destination body $d$ due
to fluctuations in source body $s$.
You may use either or both of `--SourceObject`
and `--DestObject` to compute
only quantities sourced by a given source body
and/or only energy/momentum transfer to a given 
destination body; the arguments to these options are the
labels you specified for the objects in question in the `.scuffgeo` file.

Note: Because of the way the EMTPFT method works, this
algorithm always computes results for all destination
bodies, so `--DestObject` has no effect on EMTPFT 
calculations. (But `--SourceObject` can still be used
to reduce the number of EMTPFT calculations performed.)

````
--OmitSelfTerms
````

This flag may be set to bypass the calculations
of self-contributions
(i.e. fluxes of the form $\Phi_{s\to s}$).

--------------------------------------------------
# 3. <span class="SC">scuff-neq</span> output files

<a name="#LogFile"></a>
### The `.log` file 

Like all command-line codes in the [[scuff-em]] suite,
[[scuff-cas3d]] writes a [`.log` file][LogFiles] that you
can monitor to keep track of your calculation's progress.

<a name="#SIFluxFile"></a>
### Output files for spatially-integrated PFTs: The `.SIFlux.METHOD` file

If you used an option like `--EMTPFT` or `--DSIPFT` to request
computation of spatially-integrated PFT data,
you will get back a file with extension `.SIFlux.METHOD`
(where `SIFlux` stands for "spatially-integrated flux" and 
where `METHOD` is replaced by the method used, i.e., EMTPFT, DSIPFT,
etc.)
reporting values
of the *generalized fluxes *$\Phi(\omega)$ at each frequency you requested. 

Generalized fluxes represent temperature-independent
frequency-resolved contributions to
thermally-averaged rates of energy and momentum transfer; 
as noted, above, they are multiplied by temperature-dependent
Bose-Einstein factors and integrated over frequency to yield
full thermally and temporally-averaged data. This integration
operation is carried out by the utility
code [<span class=SC>scuff-integrate</span>][scuffIntegrate],
as described [below](scuff-neq.md/#scuffIntegrate).

### Output files for spatially-resolved PFTs: The `.SRFlux` file

If you requested the computation of spatially-resolved
power and momentum flux (by specifying the `--EPFile` 
command-line option), you will get back a file with
extension `.SRFlux` reporting temperature-dependent
frequency-resolving contributions to the Poynting flux
and Maxwell stress at each point in your `EPFile`.
(See below for more information on what these values mean.)

--------------------------------------------------
# 4. Using <span class="SC">scuff-integrate</span> to perform frequency integrals

After running [[scuff-neq]] to get temperature-independent 
generalized flux data at various frequencies, we
run [[scuff-integrate]] to evaluate the frequency integrals
that compute total thermally-averaged rates of energy
and momentum transfer.

## Specifying temperatures

To convert frequency-resolved data into thermally-averaged
data, we need to know the temperature of all bodies in
the geometry, as well as that of the environment; in 
general, we will want results for multiple 
 *different* combinations of body and environment temperatures
This information is specified to [[scuff-integrate]] in the
form of a simple text file (I usually call it simply `TemperatureFile`)
in which each line describes
a separate assignment of temperatures to the environment
and to the various bodies in the geometry, with format
```
 TEnv   T1   T2   ...   TN
```
where `TEnv` is the environment temperature 
and `T1`, `T2`, ..., `TN` are the temperatures of the
first, second, ..., `N`th bodies in your geometry 
(indexed in order of their appearance in the `.scuffgeo`
file).

For example, here's a `TemperatureFile` for a two-body
geometry in which we request data for 
the first body held at temperatures of $T_1=\{50,100,150\}$ K
with the second body held fixed at $T_2=100 $ K and the 
environment kept at $T_{\text{env}}=0$ K
(temperatures are always specified in units of Kelvin):

```bash
0 50  100
0 100 100 
0 150 100 
```

## Spatially-integrated quantities

For spatially-integrated quantities (total heat radiation/absorption,
force, and torque on bodies),
[[scuff-integrate]] inputs a `.SIFlux` produced by [[scuff-neq]]
and produces two further output files with extensions
`.NEQPFT` and `.SIIntegrand,` which report respectively 
the total (frequency-integrated) thermally-averaged quantities 
and the temperature-weighted integrand of the frequency
integral that produced them. 

To understand the various quantities written to these 3 files, 
let $Q_d$ be
the spatially-integrated PFT on a destination body $d$,
and write the FSC decomposition of the thermal average
of $Q_d$ in the form

$$ \big\langle Q_d\big\rangle
   = 
   \underbrace{ 
    \Bigg[ \int_0^\infty
     \underbrace{ 
      \bigg\{ \hbar\omega_0^2 \sum_b \, \Delta \widehat \Theta_b(u)
       \underbrace{ \Phi_{s\to d}(u)}_{\texttt{.SIFlux.XXXPFT}}
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

+ $\Delta \widehat \Theta_b(u)=\widehat \Theta_b(u) - \widehat \Theta_{\small env}(u)$ 
is the difference between the dimensionless Bose-Einstein factors of source body $s$ 
and the environment.

+ $\Phi_{s\to d}$ is the temperature-independent generalized flux
describing the contributions of fluctuations in source body $s$
to the power, force, or torque on destination body $d$.
Values of this quantity are written to the `.SIFlux.METHOD` output
file (where `METHOD` is a string like `EMTPFT` identifying the 
spatial-integration method used).

+ $\hbar\omega_0^2 \Delta \widehat \Theta_b(u) \Phi_{s\to d}(u)$ is
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

## Spatially-resolved quantities

For spatially-resolved data,
[[scuff-integrate]] inputs a `.SRFlux` produced by [[scuff-neq]]
and produces two further output files with extensions
`.PVMST` and `.SRIntegrand,` which report respectively 
values of the total (frequency-integrated) thermally-averaged
 **P**oynting **V**ector and **M**axwell **S**tress **T**ensor
and the temperature-weighted integrand of the frequency
integral that produced them. 

The breakdown
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
      \bigg\{ \hbar\omega_0^2 \sum_b \, \Delta \widehat \Theta_b(u)
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

+ $\hbar\omega_0^2 \Delta \widehat \Theta_b(u) \Phi_{s\to \mathbf{x}}(u)$ is
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

<a name="bibliography"></a>

## Bibliography

Here are the original papers cited above describing the FSC approach to 
fluctuational electromagnetism:

1. [M. T. H. Reid et al, "Efficient Computation of Casimir Interactions between Arbitrary 3D Objects." *Physical Review Letters* **103** 040401 (2009). DOI: 10.1103/PhysRevLett.103.040401](http://link.aps.org/doi/10.1103/PhysRevLett.103.040401)

2. [A. Rodriguez, M. T. H. Reid, S.G. Johnson, "Fluctuating-surface-current formulation of radiative heat transfer: Theory and applications." *Physical Review B* **88** 054305 (2013). DOI: doi/10.1103/PhysRevB.88.054305](http://link.aps.org/doi/10.1103/PhysRevB.88.054305)

3. [M. T. H. Reid et al, "Photon Torpedoes and Rytov Pinwheels: Integral-Equation Modeling of Non-Equilibrium Fluctuation-Induced Forces and Torques on Nanoparticles." arXiv:1708.01985](https://arxiv.org/abs/1708.01985)

The various algorithms for computing spatially-integrated data
(DSI, EMT, etc.) are described in Ref. 3 here; see also [this paper][SIEPFTPaper].


[Geometries]:       ../../reference/Geometries.md
[Transformations]:  ../../reference/Transformations.md
[Units]:            ../../reference/FAQ.md/#Units
[scuff-cas3D]:      ../scuff-cas3D/scuff-cas3D.md
[scuffIntegrate]:   ../scuff-integrate/scuff-integrate.md
[EPFile]:           ../../applications/GeneralReference.md#EvaluationPoints
[FSCPaper]:         scuff-neq.md#bibliography
[SIEPFTPaper]:	    http://dx.doi.org/10.1109/TAP.2015.2438393
[LogFiles]:        ../GeneralReference.md#LogFiles
[SiO2Sphere]:      ../../examples/SiO2Spheres/SiO2Spheres.md
[SiO2Spheres]:     ../../examples/SiO2Spheres/SiO2Spheres.md
[TipSubstrate]:    ../../examples/TipSubstrate/TipSubstrate.md
[CommonOptions]:   ../GeneralReference.md#CommonOptions
