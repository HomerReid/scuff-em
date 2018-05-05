# Input impedance, mutual coupling, and radiated fields of microstrip RF devices

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

# Describing microstrip geometries in <span class=SC>scuff-em</span>

The description of microstrip geometries in <span class=SC>scuff-em</span>
consists of three ingredients:

1. one or more meshed polygons defining metal traces (red regions in the
pictures above)

2. a description of the substrate (its relative permittivity $\epsilon_r$ and
thickness $h$).

3. a definition of where the *ports* are located on the traces (indicated
schematically by white lines and arrows in the pictures above.)

We'll consider each of these items in turn.

## 1. Specifying meshed metal traces

### 1A. Write a `.scuffgeo` file

One way to define metal traces is simply to write a standard [<span class=SC>scuff-em</span> geometry file][scuffEMGeometries]
specifying one or more mesh files (such as [<span class=SC>gmsh</span>][GMSH] `.msh` files or [<span class=SC>comsol</span>][COMSOL] `.mphtxt` files),
optionally subjected to rigid geometrical transformations (translations and rotations). For example, [<span class=SC>scuff-em</span>
geometry files for the coupled-antenna geometry](#CoupledAntennaSCUFFGEOFile) shown above describe the geometry as containing two copies
of a mesh file for the single antenna, the second of which is displaced relative to the first.

### 1B. Use API routines to add polygon meshes

Alternatively, you can skip writing a `.scuffgeo` file altogether and instead just build up your
geometry from a C++ or python code by making one or more calls to the `AddMetalTraceMesh()`
API routine. Each call to this routine adds one copy of the meshed polygon to your geometry,
optionally displaced or rotated. Thus, the following two lines in a python code
yield the same effect as the `.scuffgeo` file above:

```python
  AddMetalTraceMesh('

### 1C. Use API routines to import layers from GDSII files

As an alternative to calling `AddMetalTraceMesh` to specify pre-meshed polygons,
you can call the `ImportGDSIILayer()` API routine to request that all polygons
be instantiated as metal traces in your geometry.

## 2. Define the substrate

### 2A. Write a `.substrate` file

One way to specify the substrate is to write a simple text file (conventionally 
assigned file extension `.substrate`) describing the substrate in your geometry.
As discussed in more detail in the technical documentation for the 
[<span class=SC>scuff-em</span> full-wave substrate module][FullWaveSubstrate],
each line of this file adds either **(a)** a single material layer to the substrate,
specified by the $z$-coordinate of its upper surface and a
[<span class=SC>scuff-em</span> material descriptor][scuffEMMaterials] describing
its bulk permittivity and permeability, or **(b)** a perfectly conducting ground plane,
specified by its height ($z$-coordinate) and the keyword `GROUNDPLANE`.
Here's a `.substrate` file defining a dielectric layer of relative
permittivity $\epsilon_r = 2.2$ and thickness 0.794 mm terminated 
below by a ground plane.

```python
0.0    CONST_EPS_2.2
-0.794 GROUNDPLANE
```

Note that, if the `GROUNDPLANE` line were omitted, this file would
describe an infinitely thick substrate (half-space) filling the region $z<0$.

### 2B. Include a `SUBSTRATE....ENDSUBSTRATE` section in your `.scuffgeo` file

If you chose to write a `.scuffgeo` file to define your metal traces, you can
include the substrate definition in this file by adding a `SUBSTRATE...ENDSUBSTRATE`
clause. The text between the opening and closing keywords is processed
as if it had been read in from the `.substrate` file; thus, to specify the 
substrate described above we would say

```python
SUBSTRATE
     0.0   CONST_EPS_2.2
  -0.794   GROUNDPLANE
ENDSUBSTRATE
```

### 2C. Use API routines to define substrate properties

Finally, you can specify substrate properties directly from C++ or
python codes by calling API routines:

```python
 SetSubstratePermittivity(2.2);

 # the following are equivalent
 SetSubstrateThickness(+0.794);
 AddGroundPlane(-0.794);

 AddSubstrateLayer(0.0, 2.2);

 SetSubstrateFile("MyFile.substrate");
```

## 3. Specifying ports

### 3A. Write a `.Ports` file

### 3B. Import port definitions from a `.GDSII` file


# Double-checking your geometry specification: The `PlotGeometry` API routine

# Example 1: Edge-fed patch antenna

# Example 2: Coupled center-fed patch antennas

<a name="CoupledAntennaSCUFFGEOFile">

The <span class=SC>scuff-em</span> geometry files

# Example 3: Coplanar waveguide section

{!Links.md!}
