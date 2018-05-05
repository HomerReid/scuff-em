# Input impedance, mutual coupling, and radiated fields of microstrip RF devices

In this example, we use the
[<span class=SC>scuff-em</span> RF module][scuff-rf]
together with
[<span class=SC>scuff-em</span>'s built-in support for
implicit multilayer dielectric substrates in full-wave calculations][libSubstrate]
to compute multiport $Z$-parameters, $S$-parameters, and radiated
field pattern of antennas, coplanar waveguides, and other *microstrip RF
devices*---arbitrarily-shaped metal traces lying on the surface of a
dielectric substrate, optionally terminated below by a perfectly-conducting
ground plane. More specifically, we will consider three particular examples
(click the pictures for larger images):

|                                                  |                                                                                                                                                                             |
|:------------------------------------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|[![](EFAntenna.png)](EFAntenna.png)               | an *edge-fed patch antenna*, for which we will compute input impedance over a wide bandwidth and radiation patterns at specific resonance frequencies                       |
|[![](CoupledAntennas.png)](CoupledAntennas.png)   | two coupled *center-fed patch antennas*, whose mutual coupling we will study as a function of their separation, and                                                         |
|[![](CPW.png)](CPW.png)                           | a section of a *coplanar waveguide*, which we will treat as a two-port device and compare computed multiport S-parameters to the predictions of transmission-line theory.   |

In all of these calculations, the metal traces (red regions) are described by
the usual [discretized triangle mesh files common to in all <span class=SC>scuff-em</span>
geometries][scuffGeometries], but the dielectric substrate (green) and ground plane
(black bottom layer) are *not meshed*; instead, the presence of the substrate layers
is taken into account 
[*implicitly* through the use of an appropriate Green's function.][libSubstrate]

The input files for all of these examples may be found in
the `share/scuff-em/examples/MicrostripDevices`
subdirectory of your `scuff-em` installation.

# The <span class=SC>scuff-em</span> RF module

The [<span class=SC>scuff-em</span> RF module][scuff-rf]
that we will use for these calculations is an extension of the
[<span class=SC>scuff-em</span> core library][libscuff] to allow
RF device modeling. More specifically, the RF module extends
the surface-integral-equation (SIE) EM solver implemented by
<span class=SC>scuff-em</span> in two key ways:

+ As a new type of incident-field source, <span class=SC>scuff-rf</span>
introduces the notion of *RF port.* An RF port is simply
a localized region of a material body into which an RF current
may be injected. The fields radiated by currents forced into RF ports 
act like the incident fields in an SIE scattering problem, which 
<span class=SC>scuff-em</span> solves like any other scattering problem
 to determine the induced currents everywhere
*else* (outside the ports) on the surfaces in a geometry.

+ As a new type of post-processing calculation, <span class=SC>scuff-rf</span>
implements an algorithm for computing the matrix of *impedance parameters*
for a multiport RF structure.

Technical details regarding the implementation of both of these features
are discussed in the [<span class=SC>scuff-rf</span> technical memo][scuffRFMemo].

# Describing microstrip geometries in <span class=SC>scuff-em</span>

Describing

{!Links.md!}
