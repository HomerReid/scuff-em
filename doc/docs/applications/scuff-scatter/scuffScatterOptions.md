||
|![](scuff-EM/scuff-scatter/scuff-scatter.png)||[scuff-scatter](scuff-EM/scuff-scatter) Command-Line Reference|

The following table summarizes all command-line options currently available in [scuff-scatter.](scuff-EM/scuff-scatter)

As is true for all programs in the scuff-em suite, command-line options may be specified in a text file catted to standard input; see [here](scuff-EM/reference/scuffEMMisc.shtml#Options) for an example of how this works.

Note that, if you find yourself needing more flexibility than can be achieved with these command-line options, then you have probably graduated from using prefab command-line programs to writing your own applications using the [c++,](scuff-EM/libscuff/cpp.shtml) [python,](scuff-EM/libscuff/python.shtml) or [matlab](scuff-EM/libscuff/matlab.shtml) interfaces to libscuff.

Option

Description

*Options controlling the scattering geometry*

    --geometry MyGeometry.scuffgeo

Specifies the [geometry file](scuff-EM/reference/scuffEMGeometries.shtml) describing the scattering geometry. This option is always mandatory.

*Options specifying the incident field*

<table>
<col width="100%" />
<tbody>
<tr class="odd">
<td align="left"><pre><code>--pwDirection nx ny nz </code></pre></td>
</tr>
<tr class="even">
<td align="left"><pre><code>--pwPolarization Ex Ey Ez </code></pre></td>
</tr>
</tbody>
</table>

Selects the incident field to be a plane wave, with propagation vector **n**=*(nx,ny,nz)* and **E**-field polarization vector **E**=*(Ex,Ey,Ez).*

More specifically, the fields of a plane wave are

![](scuff-em/libscuff/PlaneWave.png)

where the components of the vectors **n** and **E<sub>0</sub>** are what you specify with the `--pwDirection` and `--pwPolarization` options. (The frequency ω is specified using other options described below; the quantities ε and μ are the material properties of the exterior medium at this frequency, which are determined by the material property designation you give the external medium in the `.scuffgeo` file; the wave impedance of the medium is *Z*=√μ/ε≈377 Ω in vacuum.)

The values specified for `--pwPolarization` may be [complex numbers.](scuff-em/reference/scuffEMMisc.shtml#Complex)

As an example, the command-line options

~~~~ {.listing}
             --pwDirection 0 0 1 --pwPolarization 0.7071 0.7071i 0.0
           
~~~~

will specify an incident field consisting of a circularly polarized plane wave traveling in the positive *z* direction.

<table>
<col width="100%" />
<tbody>
<tr class="odd">
<td align="left"><pre><code>--gbDirection nx ny nz </code></pre></td>
</tr>
<tr class="even">
<td align="left"><pre><code>--gbPolarization Ex Ey Ez </code></pre></td>
</tr>
<tr class="odd">
<td align="left"><pre><code>--gbCenter Cx Cy Cz </code></pre></td>
</tr>
<tr class="even">
<td align="left"><pre><code>--gbWaist W</code></pre></td>
</tr>
</tbody>
</table>

Selects the incident field to be a focused Gaussian beam, traveling in the direction defined by the unit vector **n**=`(nx,ny,nz)`, with **E**-field polarization vector **E**=`(Ex,Ey,Ez)`, beam center point with cartesian coordinates **C**=`(Cx,Cy,Cz)`, and beam waist `W`.

The values specified for `--gbPolarization` may be [complex numbers.](scuff-em/reference/scuffEMMisc.shtml#Complex)

For more information on gaussian beams in scuff-EM, see the [`IncField` documentation.](scuff-em/libscuff/IncField.shtml#GaussianBeam) For an example of how to use a gaussian beam in scuff-scatter, see [this example.](scuff-em/scuff-scatter/scuffScatterExamples.shtml#Diffraction)

<table>
<col width="100%" />
<tbody>
<tr class="odd">
<td align="left"><pre><code>--psStrength Px Py Pz</code></pre></td>
</tr>
<tr class="even">
<td align="left"><pre><code>--psLocation xx yy zz</code></pre></td>
</tr>
</tbody>
</table>

Selects the incident field to be the field of a pointlike electric dipole radiator with dipole moment **P**=`(Px,Py,Pz)` and located at cartesian coordinates (`xx,yy,zz`).

The values specified for `--psStrength` may be [complex numbers.](scuff-em/reference/scuffEMMisc.shtml#Complex)

You may define the incident field to be a superposition of the fields of multiple point radiators by specifying these options more than once. (The *n*th occurrence of `--psStrength` will be paired with the *n*th occurrence of `--psLocation.`)

*Options specifying the frequency*

     --Omega xx 

Specifies the angular frequency of the scattering problem, in units of 3•10<sup>14</sup> rad/s (=c / 1 μm).

The value specified for `--Omega` may be a [complex number.](scuff-em/reference/scuffEMMisc.shtml#Complex)

You may request computations at more than one frequency by using the `--Omega` option more than once, i.e. you may say

~~~~ {.listing}
            --Omega 1e-2 --Omega 1e-1 --Omega 1 --Omega 10
           
~~~~

(However, for more than a few frequencies it is more convenient to use the `--OmegaFile` option discussed below.)

     --OmegaFile xx 

Specifies the name of a file containing a list of frequencies at which to do computations.

The file should contain one frequency per line; blank lines and comments (lines beginning with `#`) are ignored.

*Options requesting output data on fields away from scattering surfaces*

<table>
<col width="100%" />
<tbody>
<tr class="odd">
<td align="left"><pre><code> --EPFile MyEPFile </code></pre></td>
</tr>
</tbody>
</table>

Specifies a list of points at which the scattered and total fields are to be calculated.

The file `MyEPFile` should contain 3 numbers on each line (the cartesian components of the evaluation points); blank lines and comments (lines beginning with `#`) are ignored.

You may use the `--EPFile` option more than once to specify multiple lists of field evaluation points. Each input file specified using `--EPFile` will generate separate output files tabulating the fields at the requested points (see below).

Note that the evaluation points in the `--EPFile` should not lie close to the surfaces of scattering objects (where \`\`close'' means roughly \`\`within a few discretization lengths,'' where the discretization length is the linear dimension of the panels in the surface mesh). If you need information on the fields at the surfaces of scatterers, use the alternative output options described below.

<table>
<col width="100%" />
<tbody>
<tr class="odd">
<td align="left"><pre><code> --FluxMesh FluxMesh.msh</code></pre></td>
</tr>
</tbody>
</table>

Requests generation of graphical data files plotting field strengths and Poynting fluxes on user-supplied surface meshes. (See examples below).

You may use the `--FluxMesh` option more than once to specify multiple flux mesh files.

Note that the visualization surfaces specified with `--FluxMesh` options must not intersect the surfaces of any objects in the scattering geometry, and in fact no point on a `--FluxMesh` surface should lie within a few discretization lengths of any scattering object. If you would like graphical information on the fields at scattering surfaces, use `--PlotSurfaceCurrents` instead of `--FluxMesh.`

*Options requesting output data on power, force, and torque on scattering objects*

<table>
<col width="100%" />
<tbody>
<tr class="odd">
<td align="left"><pre><code> --PFTFile MyPFTFile.dat</code></pre></td>
</tr>
</tbody>
</table>

Requests that values of scattered and absorbed power, force (radiation pressure), and torque for all scattering objects be written to file `MyPFTFile.dat.` (`PFT` stands for **p**ower, **f**orce, and **t**orque.)

The absorbed power, force, and torque represent respectively transfers of energy, linear momentum, and angular momentum from the incident field sources to material objects. In principle, these quantities could be evaluated by integrating the Poynting vector and the Maxwell stress tensor over a sphere bounding the object, with the field values needed obtained by specifying an `EPFile` as above. However, this procedure would be both unwieldy and numerically ill-behaved. Instead, the computation requested by the `--PFTFile` option makes use of special formulas specific to the SIE/BEM formulation that are both faster and more accurate than such an approach.

*Options requesting output data on fields/sources at scattering surfaces*

<table>
<col width="100%" />
<tbody>
<tr class="odd">
<td align="left"><pre><code> --MomentFile MyMomentFile.dat</code></pre></td>
</tr>
</tbody>
</table>

Requests that values of the electric and magnetic dipole moments induced on each scattering object by the incident field be written to the file `MyMomentFile.dat.`

<table>
<col width="100%" />
<tbody>
<tr class="odd">
<td align="left"><pre><code> --PSDFile MyPSDFile.dat</code></pre></td>
</tr>
</tbody>
</table>

Requests that values of *panel source densities* (PSDs) be written to the file `MyPSDFile.dat`. The panel source densities may be interpreted either as **(1)**effective electric and magnetic surface currents and charges induced by the incident field on the surfaces of scattering objects, or as **(2)** the tangential and normal components of the total **E-** and **H-** fields at the surfaces of the scattering objects.

You can think of the information contained in the `--PSDFile` output as a spatially-resolved version of the information contained in the `--MomentFile` output.

The `--PSDFile` file produces output in a text file format that is intended for use as input to various post-processing scripts. If instead you would prefer to get the output in a format suitable for graphical visualization, use `--PlotSurfaceCurrents` instead.

<table>
<col width="100%" />
<tbody>
<tr class="odd">
<td align="left"><pre><code> --PlotSurfaceCurrents </code></pre></td>
</tr>
</tbody>
</table>

Requests the creation of graphical data files which you can open in gmsh to visualize the electric and magnetic current and charge distributions induced by the incident field on the surfaces of the scattering objects. See examples below.

*Cache options*

     --Cache MyCache.cache 

     --ReadCache InCache1.cache 

     --ReadCache InCache2.cache 

     --WriteCache OutCache.cache 

Specifies the names of [cache files for geometric data.](scuff-em/reference/scuffEMMisc.shtml#Caching)

The `--ReadCache` option allows you to specify a file from which geometric data will be preloaded before the calculation begins. This option may be specified any number of times. If the specified file does not exist, the option is silently ignored.

The `--WriteCache ` option allows you to specify a file to which geometric data will be written after the calculation has completed. The resulting cache file will contain any data that were preloaded from `--ReadCache` files, plus any data that were newly generated during the course of the scuff-scatter run.

The `--Cache XX` option is equivalent to saying `--ReadCache XX --WriteCache XX.`

For more information on geometric data caching in scuff-em, see [here.](scuff-em/reference/scuffEMMisc.shtml#Caching)

For examples of how caching is used in practical scuff-scatter runs, see [this example](scuff-em/reference/scuffEMMisc.shtml#Mie) or [this example.](scuff-em/reference/scuffEMMisc.shtml#SphericalShell)

*Other options*

     --HDF5File MyFile.hdf5 

Export numerical matrices and vectors [including the BEM matrix, the vector of incident-field projections (RHS vector), and the vector of surface-current coefficients (solution vector)] to an HDF5 data file named `MyFile.hdf5.` See [here](scuff-EM/scuff-scatter/scuffScatterFiles.shtml) for more information.

<table>
<col width="100%" />
<thead>
<tr class="header">
<th align="left"><a href="scuff-em/scuff-scatter">scuff-scatter</a> Documentation</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">1. scuff-scatter Command-Line Options</td>
</tr>
<tr class="even">
<td align="left"><a href="scuff-EM/scuff-scatter/scuffScatterFiles.shtml">2. scuff-scatter Output Files</a></td>
</tr>
<tr class="odd">
<td align="left"><a href="scuff-EM/scuff-scatter/scuffScatterExamples.shtml">3. scuff-scatter Examples</a>
<blockquote>
<table>
<tbody>
<tr class="odd">
<td align="left"><a href="scuff-EM/scuff-scatter/scuffScatterExamples.shtml#Mie">4a. Mie scattering</a></td>
</tr>
<tr class="even">
<td align="left"><a href="scuff-EM/scuff-scatter/scuffScatterExamples.shtml#SphericalShell">4b. Electrostatics of a spherical dielectric shell</a></td>
</tr>
<tr class="odd">
<td align="left">4c. Spatially-resolved study of plane-wave transmission through a infinite-area thin dielectric film</td>
</tr>
<tr class="even">
<td align="left">4d. Diffraction of a Gaussian laser beam by a finite disc</td>
</tr>
</tbody>
</table>
</blockquote></td>
</tr>
</tbody>
</table>


