# Diffraction patterns for discs, disc arrays, and hole arrays in metal screens

In this example, we shine a laser beam (or a plane wave) on an 
(infinite-area) metal screen perforated by a square-lattice
array of circular holes, and produce images of the diffraction 
patterns as observed on a visualization surface located behind the 
perforated screen. Here's a schematic depiction of the configuration:

![Diffraction experiment image](DiffractionSchematic.png)

The files for this example may be found in the
`share/scuff-em/examples/DiffractionPatterns` subdirectory
of your [[scuff-em]] installation.

--------------------------------------------------
## [[gmsh]] geometry file and surface mesh for the screen unit cell

The [[gmsh]] geometry file [`HoleyScreenUnitCell.geo`](HoleyScreenUnitCell.geo)
describes an (infinitely thin) square metallic screen, 
of dimensions 1&mu;m &times; 1&mu;m, with a hole of radius 0.25 &mu;m
centered at the center of the square. I produce a surface mesh for
this geometry by saying

````bash
% gmsh -2 -clscale 1 HoleyScreenUnitCell.geo
% RenameMesh HoleyScreenUnitCell.msh
````
(where [`RenameMesh`](../SiO2Spheres/RenameMesh) is a simple 
`bash` script that uses [[scuff-analyze]] to count the number 
of interior edges in a surface mesh and rename the mesh file 
accordingly.)
This produces the file `HoleyScrenUnitCell_1228.msh,`
which you can visualize by opening in [[gmsh]]::

````bash
% gmsh HoleyScreenUnitCell_1228.msh
````
![Holey screen mesh image](HoleyScreenUnitCell_320.png)

If you need finer meshing resolution, you can tweak
the `-clscale ` argument to [[gmsh]] (it stands
for "characteristic length scale"; for example,
running with `-clscale 0.5` will do the meshing
with 2x finer lengthscale (corresponding to a 
roughly 4x increase in the number of triangles).

--------------------------------------------------
## [[scuff-em]] geometry file

The [[scuff-em]] geometry file
[`HoleyScreen_1228.scuffgeo`](HoleyScreen_1228.scuffgeo)
describes an infinite square lattice whose unit
cell is defined by the unit-cell mesh we created
above. It looks like this:

````bash
LATTICE
        VECTOR 1 0
        VECTOR 0 1
ENDLATTICE    

OBJECT HoleyScreen
        MESHFILE HoleyScreenUnitCell_1228.msh
ENDOBJECT
````

Note that we don't have to specify a `MATERIAL`
for the screen, since PEC is the default.

We can use [[scuff-analyze]] to produce an image
of what the full geometry looks like, including
the lattice repetitions:

````bash
% scuff-analyze --geometry HoleyScreen_1228.scuffgeo --WriteGMSHFiles --Neighbors 2
````

This produces the file `HoleyScreen_1228.pp`, which you 
can view by opening it in [[gmsh]]:

````bash
% gmsh HoleyScreen_1228.pp
````

![HoleyScreen geometry image](HoleyScreenGeometry.png)

--------------------------------------------------
## Field visualization mesh

The next step is to create a meshed representation of the
surface on which we will visualize the diffraction patterns.
Here's a [[gmsh]] file called
[`FVMesh.geo`](FVMesh.geo) that describes a square of
side length 1 micron, parallel to the *xy* plane and
located at a height of *z*=1 micron, thus corresponding
to the region enclosed by the dotted line in the schematic
figure above. ("FVMesh" stands for `field-visualization
mesh.) This `.geo` file contains a user-specifiable
parameter `N` that sets the number of triangle edges per
unit length in the mesh representation; I would
like to set this number to 50, so I say

````bash
% gmsh -2 -setnumber N 50 FVMesh.geo -o FVMesh.msh
% RenameMesh FVMesh.msh
````

This produces the file `FVMesh_7400.msh`:

![FVMesh image](FVMesh_7400.png)

--------------------------------------------------
## Running [[scuff-scatter]]

Now all that's left is to run the calculation.
Put the following content into a little text
file called `scuff-scatter.args`:


--------------------------------------------------

[scuff-neq]:              ../../applications/scuff-neq/scuff-neq.md
[Transformations]:        ../../reference/Transformations
[KruegerPaper]:           http://dx.doi.org/10.1103/PhysRevB.86.115423
