||
|![](scuff-EM/scuff-scatter/scuff-scatter.png)||[scuff-scatter Output File Reference](scuff-EM/scuff-scatter)|

### 1. Field components at evaluation points

If you use the `--EPFile` command-line argument to specify a file named `MyPoints.dat`, scuff-scatter will produce output files named `MyPoints.scattered` and `MyPoints.total`. The former will contain information on the *scattered* electric and magnetic fields, while the latter will contain information on the *total* fields.

(Note that the scattered and total fields differ only in regions containing incident field sources. For example, if we have a dielectric sphere illuminated by a plane wave, the scattered and total fields agree for evaluation points inside the sphere; outside the sphere, the two fields differ, with the difference being exactly the field of the incident plane wave.)

Each line of the `.scattered` and `.total` files will correspond to a single evaluation point (that is, a single line of the file `MyPoints.dat`). There will be 15 numbers on each line, separated by spaces, as follows:

||
|*x y z re(Ex) im(Ex) re(Ey) im(Ey) re(Ez) im(Ez) re(Hx) im(Hx) re(Hy) im(Hy) re(Hz) im(Hz)*|

More specifically, the 15 columns of the output files are as follows:

|Column|Significance|Column|Significance|Column|Significance|
|:-----|:-----------|:-----|:-----------|:-----|:-----------|
|1|*x*|2|*y*|3|*z*|
|4|*Re(E<sub>x</sub>)*|5|*Im(E<sub>x</sub>)*|6|*Re(E<sub>y</sub>)*|
|7|*Im(E<sub>y</sub>)*|8|*Re(E<sub>z</sub>)*|9|*Im(E<sub>z</sub>)*|
|10|*Re(H<sub>x</sub>)*|11|*Im(H<sub>x</sub>)*|12|*Re(H<sub>y</sub>)*|
|13|*Im(H<sub>y</sub>)*|14|*Re(H<sub>z</sub>)*|15|*Im(H<sub>z</sub>)*|

### 2. Power, force, and torque

If you use the `--PFTFile` command-line argument to specify a file named `MyPFTFile.out`, then scuff-scatter will *append* to this file a single line for each frequency at which you requested scattering calculations.

This line begins with the angular frequency, followed by the label of the first object in the geometry, and then 8 numbers for that object: the absorbed power, the scattered power, the three cartesian force components, and the three cartesian torque components. Next comes the label of the second object in the geometry, and then 8 numbers for that object, etc. Thus the format of each line is:

||
|*Omega Object-Label-1 Pabs Pscat Fx Fy Fz Tx Ty Tz Object-Label-2 Pabs Pscat Fx ...*|

Here *Pabs* is the power absorbed by the object from the incident fields, while *Pscat* is the total power carried away by the scattered fields. (The sum of these two quantities is sometimes called the the "extinction"; it is the total power taken out of the incident fields by the scattering process.) *Pabs* and *Pscat* are reported in units of *watts.*

Note that these are powers, not cross-sections. To convert these numbers into cross-sections, you would divide by the incident flux (power-per-unit-area) of your field source. For a plane wave in vacuum, the incident flux is *|**E<sub>0</sub>**|<sup>2</sup>/2Z<sub>0</sub>* ≈ 1/(2\*377) watts / micron<sup>2</sup> for *E<sub>0</sub>*=1 V/micron.

`Fx,Fy,Fz` are the cartesian components of the force (radiation pressure) on the object; these are reported in units of *nanoNewtons.*

`Tx,Ty,Tz` are the cartesian components of the torque on the object; these are reported in units of *nanoNewtons • microns.*

Actually, the units quoted here are based on the assumption that you were using units of volts/micron when you specified the numerical values for the incident-field magnitudes. For example, if you declared a plane wave with **E-**field polarization vector `--pwPolarization 1 0 0`, and you meant by this that the *x*-component of the field had units of `1` volt/micron, then the units of the `PFT` quantities will be as indicated above. If, on the other hand, you actually intended units of volts/meter, then the numbers reported in the `PFTFile` for power, force, and torque should be interpreted in gigawatts, milliNewtons, and milliNewtons • meters.

### 3. Induced dipole moments

If you use the `--MomentFile` command-line argument to specify a file named `MyMoments.out`, then scuff-scatter will *append* to this file a single line for each frequency at which you requested scattering calculations.

This line begins with the angular frequency, followed by the label of the first object in the geometry, and then 12 numbers for that object: the real and imaginary parts of the cartesian components of the induced electric and magnetic dipole moments for that object. Then comes the label of the second object, followed by another 12 numbers, etc. The format of each line is

||
|*Omega Object-Label-1 Re(Px) Im(Px) Re(Py) Im(Py) Re(Pz) Im(Pz) Re(Mx) Im(Mx) Re(My) Im(My) Re(Mz) Im(Mz) Object-Label-2 Re(Px) ...*|

Here `Px,Py,Pz` and `Mx,My,Mz` are the cartesian components of the electric and magnetic dipole moments induced on the object by the incident field. If the incident field units were volts/micron, then the units of the electric dipole moment are 8.85•10<sup>-6</sup> coulomb•meter. (The numerical constant here is the value of ε<sub>0</sub> in units of farads/micron.)

* * * * *

### 4. Panel source densities in text format

If you use the `--PSDFile` command-line argument to specify a file named `MyPSDFile.out`, then scuff-scatter will append data to this file describing the densities of electric and magnetic current and charge induced by the incident fields on the surfaces of your scattering objects. These data essentially constitute a spatially-resolved version of the same information contained in the file written by the `--MomentFile` option.

The `--PSDFile` option is designed to output data in an unadorned text file format for use as input to subsequent post-processing scripts you might cook up. If instead you want to visualize the information graphically, you probably want to use `--PlotSurfaceCurrents` instead.

Each line in the `MyPSDFile.out` will contain data on the source densities at the centroid of a single triangular panel in your surface mesh. The total number of lines in the file will be the total number of panels in the geometry (times the number of frequencies you simulate, if you do multiple frequencies).

The first number on the line is the angular frequency. The next four numbers on the line are the coordinates of the panel centroid and the panel area.

The next 2 numbers are the real and imaginary components of the electric surface charge density σ=∇•**K**/iω at the panel centroid.

The next 6 numbers are the real and imaginary components of the *x,y,z* components of the electric surface current density **K** at the panel centroid.

The next 2 numbers are the real and imaginary components of the magnetic surface charge density η=∇•**N**/iω at the panel centroid.

The next 6 numbers are the real and imaginary components of the *x,y,z* components of the magnetic surface current density **N** at the panel centroid.

Finally, the last number on the line is the inward-directed normal Poynting flux (power per unit area flowing into your scattering object) at the panel centroid. This is a spatially-resolved measure of power absorption by your object; if you multiply this number by the panel area (column 4) and sum over all panels on the surface of a scattering object, you will get the net power absorbed by that object from the incident field. (Of course, you don't need to do this calculation yourself -- it will be done automatically for you if you specify the `--PFTFile` option.)

To summarize, the columns of each line of the output file produced by `--PSDFile` are as follows:

Column

Significance

1

angular frequency ω

2

*x*-coordinate of panel centroid

3

*y*-coordinate of panel centroid

4

*z*-coordinate of panel centroid

5

panel area <sup>2</sup>

6,7

real,imag parts of *σ*

8,9

real,imag parts of *K\_x*

10,11

real,imag parts of *K\_y*

12,13

real,imag parts of *K\_z*

14,15

real,imag parts of *η*

16,17

real,imag parts of *N\_x*

18,19

real,imag parts of *N\_y*

20,21

real,imag parts of *N\_z*

22

inward-directed normal Poynting flux

**Units.** The easiest way to think about the units of the PSD quantities is that the panel area is in μm<sup>2</sup>, electric surface currents have units of amperes/μm, magnetic surface currents have units of volts/μm, and the Poynting flux has units of watts / μm<sup>2</sup>.

* * * * *

### 5. Panel source densities in graphical format

If you specify the `--PlotSurfaceCurrents` flag, then for each computational frequency you will get a file called something like `MyGeometry.1.234.pp` (where `MyGeometry.scuffgeo` was your `--geometry` file and 1.234 is the angular frequency) containing graphical data in gmsh for plotting the electric and magnetic surface currents and charges induced on your scattering geometry by the incident field. (The quantities plotted are the *real parts* of the currents and charges; you can think of the resulting plot as giving you a snapshot of the physical densities at a single instant in time.)

The visualization file produced by `--PlotSurfaceCurrents` also contains spatially-resolved data on the time-averaged inward-directed normal Poynting flux. This is the quantity whose integral over the surface of a closed object gives the total power absorbed by that object, i.e. the rate at which the object absorbs energy from the incident fields. (The spatial integral of this quantity, i.e. the total absorbed power, is reported for you in the output produced by the `--PFTFile` above; the gmsh visualization output here is useful for understanding the spatial distribution of power flow in and out of the object.)

* * * * *

### 5. Flux meshes

<table>
<col width="100%" />
<thead>
<tr class="header">
<th align="left"><a href="scuff-em/scuff-scatter">scuff-scatter</a> Documentation</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><a href="scuff-EM/scuff-scatter/scuffScatterOptions.shtml">1. scuff-scatter Command-Line Options</a></td>
</tr>
<tr class="even">
<td align="left">2. scuff-scatter Output Files</td>
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


