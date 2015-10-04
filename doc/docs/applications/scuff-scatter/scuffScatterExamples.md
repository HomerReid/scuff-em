<table>
<col width="33%" />
<col width="33%" />
<col width="33%" />
<tbody>
<tr class="odd">
<td align="left"><img src="scuff-EM/scuff-scatter/scuff-scatter.png" /></td>
<td align="left"></td>
<td align="left"><h1>Example calculations in <a href="scuff-EM/scuff-scatter">scuff-scatter</a></h1></td>
</tr>
</tbody>
</table>

Here are some examples of calculations you can do with scuff-scatter. Input files and command-line runscripts for all these examples are included in the `share/scuff-em/examples` subdirectory of the scuff-em installation.

* * * * *

### 4a. Mie scattering

We first demonstrate how to use scuff-scatter to solve the canonical textbook problem of *Mie scattering* -- the scattering of a plane wave from a dielectric sphere. The files for this example are in the `share/scuff-em/examples/SolidSphere` subdirectory of the `scuff-em` source distribution.

We begin by creating a [gmsh](http://geuz.org/gmsh) geometry file for the sphere ([`Sphere.geo`](scuff-em/scuff-scatter/Sphere.geo)). We turn this geometry file into a mesh file by running the following command:

~~~~ {.listing}
      % gmsh -2 -clscale 1.0 Sphere.geo
    
~~~~

This produces a file named [`Sphere.msh,`](scuff-em/scuff-scatter/Sphere.msh) which looks like this:

![](scuff-em/scuff-scatter/SphereMesh.png)

You can adjust the fineness of the surface mesh by varying the `-clscale` parameter (which stands for "characteristic length scale"); finer meshes will be more accurate but will take longer to simulate.

Next we create a [scuff-em geometry file](scuff-em/reference/scuffEMGeometries.shtml) that will tell scuff-scatter about our geometry, including both the surface mesh and the material properties (dielectric function) of the sphere. As a first example, we'll use a dielectric model for silicon carbide that expresses the relative permittivity as a rational function of ω; in this case we'll call the geometry file `SiCSphere.scuffgeo.`

~~~~ {.listing}
MATERIAL SiliconCarbide
   
   EpsInf = 6.7;
   a0     = -3.32377e28;
   a1     = +8.93329e11;
   b0     = -2.21677e28;
   b1     = 8.93329e11;
   Eps(w) = EpsInf * (a0 + i*a1*w + w*w) / ( b0 + i*b1*w + w*w);

ENDMATERIAL 

OBJECT TheSphere
        MESHFILE Sphere.msh
        MATERIAL SiliconCarbide
ENDOBJECT
 
~~~~

We create a simple file called [`OmegaValues.dat`](scuff-em/scuff-scatter/OmegaValues.dat) containing a list of angular frequencies at which to run the scattering problem:

~~~~ {.listing}
    0.010
    0.013
    ...
    10.0
    
~~~~

(We pause to note one subtlety here: Angular frequencies specified using the `--Omega` or `--OmegaFile` arguments are interpreted in units of **c / 1 μm = 3•10<sup>14</sup> rad/sec**. These are natural frequency units to use for problems involving micron-sized objects; in particular, for Mie scattering from a sphere of radius 1 μm, as we are considering here, the numerical value of `Omega` is just the quantity `kR` (wavenumber times radius) known as the \`\`size parameter'' in the Mie scattering literature. In contrast, when specifying functions of angular frequency like `Eps(w)` in `MATERIAL...ENDMATERIAL` sections of geometry files or in any other [scuff-em material description,](scuff-em/reference/scuffEMMaterials.shtml) the `w` variable is always interpreted in units of **1 `rad/sec`**, because these are the units in which tabulated material properties and functional forms for model dielectric functions are typically expressed.)

Finally, we'll create a little text file called `Args` that will contain a list of command-line options for scuff-scatter; these will include **(1)** a specification of the geometry, **(2)** the frequency list, **(3)** the name of an output file for the power, force, and torque **(4)** the name of a cache file for geometric data (this file doesn't exist yet, but will be created by our first run of scuff-scatter), and **(5)** a specification of the incident field.

~~~~ {.listing}
    geometry SiCSphere.scuffgeo
    OmegaFile OmegaValues.dat
    PFTFile SiCSphere.PFT
    Cache Sphere.cache
    pwDirection 0 0 1
    pwPolarization 1 0 0
    
~~~~

And now we just pipe this little file into the standard input of scuff-scatter:

~~~~ {.listing}
    % scuff-scatter < Args 
    
~~~~

This produces the file `SiCSphere.PFT`, which contains one line per simulated frequency; each line contains data on the scattered and total power, the force, and the torque on the particle at that frequency.

(On my fairly standard workstation (8 Xeon E5420 cores, 2.5 GHz), this calculation takes a few minutes to run. You can monitor its progress by following the `scuff-scatter.log` file. Note that, during computationally-intensive operations such as the BEM matrix assembly, the code should be using all available CPU cores on your workstation; if you find that this is not the case (for example, by monitoring CPU usage using [htop](http://htop.sourceforge.net)) you may need to [reconfigure and recompile with different openmp/pthreads configuration options.](scuff-em/reference/scuffEMInstallation.shtml))

Here's a comparison of the scuff-scatter results with the analytical Mie series, as computed using [this Mathematica script.](scuff-em/scuff-scatter/Mie.math) [Like most Mie codes, this script computes the absorption and scattering *cross-sections*, which we multiply by the incoming beam flux (1/2*Z<sub>0</sub>* for a unit-strength plane wave) to get values for the absorbed and scattered *power*.]

![](scuff-em/scuff-scatter/SiCData.png)

Now let's redo the calculation for a sphere made of gold instead of silicon carbide. In this case we will name our scuff-em geometry file `GoldSphere.scuffgeo`:

~~~~ {.listing}
    MATERIAL Gold
      wp = 1.37e16;
      gamma = 5.32e13;
      Eps(w) = 1 - wp^2 / (w * (w + i*gamma));
    ENDMATERIAL

    OBJECT TheSphere
        MESHFILE Sphere.msh
        MATERIAL Gold
    ENDOBJECT
    
~~~~

Since most of the command-line arguments to scuff-scatter will be the same as before, we can reuse the same `Args` file, with the options that need to be given new values specified on the command line:

~~~~ {.listing}
    % scuff-scatter --geometry GoldSphere.scuffgeo --PFTFile GoldSphere.PFT < Args
    
~~~~

(Note that we don't have to change the name of the cache file specified with the `Cache` option; because we are using the same surface mesh as before, and because [cached geometric data in scuff-em are independent of material properties](scuff-em/reference/scuffEMMisc.shtml#Caching), we can take advantage of geometric data computed during the earlier run for the silicon carbide sphere.)

Now our data look like this:

![](scuff-em/scuff-scatter/GoldData.png)

In some cases it's useful to look at how the induced surface currents vary over the surface of the object. Let's re-run the SiC example, now at just the single angular frequency of ω=0.1, and ask for a surface current plot.

~~~~ {.listing}
     scuff-scatter --geometry SiCSphere.scuffgeo --Omega 0.1 --Cache Sphere.cache --pwDirection 0 0 1 --pwPolarization 1 0 0 --PlotSurfaceCurrents
    
~~~~

This produces a file named `SiCSphere.0.1.pp`, which we can open in gmsh like this:

~~~~ {.listing}
     gmsh SiCSphere.0.1.pp
    
~~~~

![](scuff-em/scuff-scatter/SphereSurfaceCurrents.png)

### 

* * * * *

### 4b. Electrostatics of a spherical dielectric shell

For our next trick, we'll consider a spherical shell of dielectric material illuminated by a plane wave of such low frequency that we may think of the incident field as a spatially constant DC electric field. In this case it is easy to obtain an [exact analytical solution of the scattering problem,](scuff-em/scuff-scatter/SphericalShellElectrostatics.pdf) which we will reproduce numerically using scuff-scatter. We will take the outer and inner radii to be *R<sub>out</sub>=1* and *R<sub>in</sub>=0.5.* The files for this example are in the `share/scuff-em/examples/SphericalShell` subdirectory of the `scuff-em` source distribution.

To represent a spherical shell in scuff-em, we need two surface meshes, one each for the inner and outer spherical surfaces. These are described by GMSH mesh files `Sphere_R1P0.msh` and `Sphere_R0P5.msh` We describe the shell as an inner vacuum sphere embedded in the outer sphere; the geometry file for this situation is `SphericalShell.scuffgeo:`

~~~~ {.Listing}
 OBJECT OuterSphere 
     MESHFILE Sphere_R1P0.msh
     MATERIAL CONST_EPS_10
 ENDOBJECT 
  
 OBJECT InnerSphere 
     MESHFILE Sphere_R0P5.msh
     MATERIAL Vacuum
 ENDOBJECT 
      
~~~~

[![](scuff-em/scuff-scatter/SphericalShell.png)](scuff-em/scuff-scatter/SphericalShell.png)

We will run two separate calculations. First, we will fix the relative permittivity of the shell at ε<sup>r</sup>=10 and look at the *z* component of the electric field at points on the *z* axis ranging from the origin (the center of the concentric spheres) to the exterior medium. We create a file called [`LineOfPoints`](scuff-em/scuff-scatter/LineOfPoints) which lists the Cartesian coordinates of each evaluation point:

~~~~ {.Listing}
      0.0 0.0 0.000
      0.0 0.0 0.025
      0.0 0.0 0.050
      ...
      0.0 0.0 2.000
    
~~~~

We will pass this file to scuff-scatter using the `--EPFile` option:

~~~~ {.Listing}
     % scuff-scatter --EPFile LineOfPoints < Args
    
~~~~

where the `Args` file looks like this:

~~~~ {.Listing}
      geometry       SphericalShell.scuffgeo
      cache          SphericalShell.cache
      omega          0.001
      pwDirection    1.0 0.0 0.0
      pwPolarization 0.0 0.0 1.0
    
~~~~

Note that we choose a frequency low enough to ensure we are well within the electrostatic limit, and that the constant *z*-directed electrostatic field described in the memo above becomes a linearly polarized plane wave traveling in the *x-* direction.

This run of the code produces files `SphericalShell.scattered` and `SphericalShell.total`. Plotting the 8th vs. the 3rd column of the latter (real part of *E<sub>z</sub>* vs. *z*, as noted [here)](scuff-em/scuff-scatter/scuffScatterFiles.shtml#EPFileTable) yields good agreement with the analytical calculation, modulo some funkiness at points on or near the boundary surfaces which is to be expected in an SIE/BEM calculation:

![](scuff-em/scuff-scatter/EzVsZ.png)

Next, we will vary the shell permittivity and look at the electric field at the center of the shell. In this case the analytical solution makes the interesting prediction

![](scuff-em/scuff-scatter/EzVsEpsEq.png)

which we will try to verify numerically.

This calculation is slightly trickier than the last one, because scuff-scatter doesn't offer command-line options for varying the dielectric constant of an object. One way around this is to use the python interface to scuff-em, as discussed [on this page.](scuff-em/libscuff/python.shtml) Here we will pursue a different solution involving a shell script that modifies the `.scuffgeo` file for each different value of ε we want to simulate. That script looks like this:

~~~~ {.listing}
     #!/bin/bash

     cat EpsValues | while read EPS
     do

       # copy the .scuffgeo file with EPS_10 replaced by EPS_xx
       sed "s/EPS_10/EPS_${EPS}/" SphericalShell.scuffgeo > temp.scuffgeo

       # run scuff-scatter to get E-field at origin
       /bin/rm -f CenterPoint.total
       /bin/rm -f CenterPoint.scattered
       scuff-scatter --geometry temp.scuffgeo --EPFile CenterPoint < Args

       # extract the z-component of the field from the output file
       EZ=`awk '{print $8}' CenterPoint.total`
       echo "${EPS} ${EZ}" >> EzVsEps.dat

     done
    
~~~~

(Here `EpsValues` is a file containing the values of epsilon that we want to simulate, and `CenterPoint` is a file containing just the first line of the file `LineOfPoints` for the cartesian coordinates of the origin.)

The result of the calculation looks like this:

![](scuff-em/scuff-scatter/EzVsEps.png)

* * * * *

### 4c. Spatially-resolved study of plane-wave transmission through a infinite-area thin dielectric film

The previous examples dealt with *compact* scatterers. We'll next consider an [*extended*](scuff-em/reference/scuffEMGeometries.shtml#ExtendedGeometries) geometry -- namely, a thin dielectric film of finite thickness in the *z* direction but infinitely extended in the *x* and *y* directions. This is the same geometry for which we used scuff-transmission to look at the plane-wave transmission and reflection coefficients as a function of frequency in [this example,](scuff-em/scuff-transmission/index.shtml#ThinFilm) but here we'll do a different calculation -- namely, we'll pick a single frequency and look at how the electric and magnetic fields vary in space, both inside and outside the thin film. (The files for this example may be found in the `share/scuff-em/examples/ThinFilm` directory of the scuff-em installation.)

The mesh file and `.scuffgeo` file for this geometry are discussed in the [documentation for scuff-transmission.](scuff-em/scuff-transmission/index.shtml#ThinFilm) The geometry consists of a film of thickness *T*=1μm, with relative dielectric constant ε<sup>*r*</sup>=100, illuminated from below by a plane wave at normal incidence. (We'll take the incident field to be linearly polarized with **E** field pointing in the *x* direction.) The lower and upper surfaces of the film are at *z=0* and *z=T.* For this geometry it is easy to solve Maxwell's equation directly to obtain the **E** and **H** fields directly at points below, within, and above the film:

![](scuff-em/scuff-scatter/EHvsz.png)

We will try to reproduce this behavior using scuff-scatter. First create a little text file ([ThinFilm.EvalPoints](scuff-em/scuff-scatter/ThinFilm.EvalPoints)) containing the coordinates of a bunch of points on a straight line passing from below the film to above the film. Then put the following command-line arguments into a file called `ThinFilm_58.args:`

~~~~ {.Listing}
 geometry ThinFilm_58.scuffgeo
 cache ThinFilm_58.scuffcache
 omega 1.0
 EPFile ThinFilm.EvalPoints
 pwDirection 0 0 1
 pwPolarization 1 0 0
    
~~~~

and pipe it into scuff-scatter:

~~~~ {.Listing}
 scuff-scatter < scuff-scatter.args
    
~~~~

This produces files named `ThinFilm.scattered` and `ThinFilm.total`. Plotting the 4th vs the 3rd column [(Re *E<sub>x</sub>* vs. *z*)](scuff-em/scuff-scatter/scuffScatterFiles.shtml#FieldOutput) as well as the the 12th vs the 3rd column [(Re *H<sub>y</sub>* vs. *z*)](scuff-em/scuff-scatter/scuffScatterFiles.shtml#FieldOutput) of the `.total` file yields excellent agreement with theory:

![](scuff-em/scuff-scatter/ThinFilmEField.png)

![](scuff-em/scuff-scatter/ThinFilmHField.png)

* * * * *

### 4d. Diffraction of a gaussian beam by a finite disc

As a final example, we shine a laser beam on a small metallic disc and observe the resulting diffraction pattern on an observation screen located behind the disc. The files for this example are in the `share/scuff-em/examples/DiffractionPatterns` subdirectory of the `scuff-em` source distribution.

|[scuff-scatter](scuff-em/scuff-scatter) Documentation|
|:----------------------------------------------------|
|[1. scuff-scatter Command-Line Options](scuff-EM/scuff-scatter/scuffScatterOptions.shtml)|
|[2. scuff-scatter Output Files](scuff-EM/scuff-scatter/scuffScatterFiles.shtml)|
|3. scuff-scatter Examples|


