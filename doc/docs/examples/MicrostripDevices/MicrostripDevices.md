<h1> Input impedance, mutual coupling, and radiated fields of microstrip RF devices </h1>

In this example, we use the
[<span class=SC>scuff-em</span> RF module][scuff-rf]
together with
[<span class=SC>scuff-em</span>'s built-in support for
implicit multilayer dielectric substrates in full-wave calculations][libSubstrate]
to compute multiport $Z$-parameters, $S$-parameters, and radiated
field pattern of antennas, coplanar waveguides, and other *microstrip RF
devices*---arbitrarily-shaped finite-sized metal traces lying on the surface of
an (infinite-area) dielectric substrate, optionally terminated below by a
perfectly-conducting ground plane. More specifically, we will consider 
three particular examples (click the pictures for larger images):

|                                                  |                                                                                                                                                                             |
|:------------------------------------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|[![](EFAntenna.png)](EFAntenna.png)               | an *edge-fed patch antenna*, for which we will compute input impedance over a wide bandwidth and radiation patterns at specific resonance frequencies                       |
|[![](CoupledAntennas.png)](CoupledAntennas.png)   | two coupled *center-fed patch antennas*, whose mutual coupling we will study as a function of their separation, and                                                         |
|[![](CPW.png)](CPW.png)                           | a section of a *coplanar waveguide*, which we will treat as a two-port device and compare computed multiport S-parameters to the predictions of transmission-line theory.   |

In all of these calculations, the metal traces (red regions) are described by
the usual [discretized triangle mesh files common to in all <span class=SC>scuff-em</span>
geometries][scuffEMGeometries], but the dielectric substrate (green) and ground plane
(black bottom layer) are *not meshed*; instead, the presence of the substrate layers
is taken into account 
[*implicitly* through the use of an appropriate Green's function.][libSubstrate]
(Note that this implicitly-described substrate is always infinitely extended in the
$x$ and $y$ directions; do not be misled by the apparently finite transverse size of
the substrates in the pictures above.)

The input files for all of these examples may be found in
the `share/scuff-em/examples/MicrostripDevices`
subdirectory of your `scuff-em` installation.

[TOC]

## The <span class=SC>scuff-em</span> RF module: A brief overview

The <span class=SC>scuff-em</span> RF module is discussed in detail
at [its main documentation page][scuff-rf], but for busy readers
we summarize here some salient points.

+ The module extends the [<span class=SC>scuff-em</span> core library][libscuff] to
  allow RF modeling problems to be studied within the framework of the
  [surface-integral-equation (SIE) approach to electromagnetic scattering implemented by <span class=SC>scuff-em</span>][libscuff].

+ The functionality of the module may be accessed either from C++/python
  API codes or via the `scuff-rf` command-line module distributed with
  <span class=SC>scuff-em</span>; in particular, all of the microstrip-device
  calculations covered in this example may be done equally well
  in python or from the command line, and we will show how to do everything
  both ways.

+ The primary way in which the RF module extends the core library is by adding
  support for [*RF ports*][RFPorts], regions of meshed surfaces into which RF currents are
  forced. This means that, in addition to specifying the usual
  [collection of meshed surfaces constituting a <span class=SC>scuff-em</span> geometry][scuffEMGeometries],
  you must also define one or more ports for your structure.
  Port terminals in <span class=SC>scuff-em</span> are defined by geometric objects:
  lines or polygons (identifying regions on the boundaries of meshed surfaces)
  or points (identifying source/sink points inside meshed surfaces).
  Port definitions may be communicated to <span class=SC>scuff-em</span> by
  **(a)** writing a simple text file (`.ports` file),
  **(b)** making API calls from C++/python,
  **(c)** creating a GDSII file containing points, lines, or polygons tagged by special text strings.
  (Examples of all of these methods are discussed below).

## Example 1: Edge-fed patch antenna

Our first example is the *edge-fed patch antenna* originally studied in this paper:

+ [Wu et al., "Feeding Structure Contribution to Radiation by Patch Antennas with Rectangular Boundaries," IEEE Transactions on Antennas and Propagation **40** 1245 (1992). DOI: 10.1109/8.18245](https://doi.org/10.1109/8.182458)

This antenna is also used as a demonstration example for the commercial solver FEKO:

+ [http://feko.info/fekoweb/applications/white-papers/edge_fed_rectangular_microstrip_patch_antenna/document](http://feko.info/fekoweb/applications/white-papers/edge_fed_rectangular_microstrip_patch_antenna/document)

A

### Python calculation

### Command-line solution

Here's how the same calculation could be run from the command line using `scuff-rf`:

## Example 2: Coupled center-fed patch antennas

<a name="CoupledAntennaSCUFFGEOFile">

The <span class=SC>scuff-em</span> geometry files

## Example 3: Coplanar waveguide section

{!Links.md!}
