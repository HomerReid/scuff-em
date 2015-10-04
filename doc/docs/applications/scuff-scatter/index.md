||
|![](scuff-EM/scuff-scatter/scuff-scatter.png)||Solving electromagnetic scattering problems with scuff-scatter|

<table>
<col width="100%" />
<tbody>
<tr class="odd">
<td align="left"><p><br /></p>
<p>scuff-scatter is a tool within the scuff-em code suite for solving a broad class of electromagnetic scattering problems.</p>
<p>To run a scattering calculation using scuff-scatter, you will</p>
<ol>
<li>create a <a href="scuff-EM/reference/scuffEMGeometries.shtml">geometry file</a> describing the shapes and material properties of the scattering objects in your geometry;</li>
<li>choose the incident field that will scatter off your objects -- a plane wave, a gaussian beam, a point dipole source, or some combination thereof</li>
<li>run scuff-scatter with command-line options specifying the geometry, the frequencies, the incident field, and the type of output you wish to get back.</li>
</ol>
<p>The various output quantities that you can ask scuff-scatter to generate include the following:</p>
<ul>
<li>The components of the scattered and total electric and magnetic fields at arbitrary user-specified points away from scattering surfaces. (The points may lie inside or outside the scattering objects).</li>
<li>The components of the total electric and magnetic fields on the scattering surfaces. (These quantities may alternatively be interpreted as effective surface currents and charges that give rise to the scattered fields.)</li>
<li>The electric and magnetic dipole moments induced by the incident field on the scattering objects. (These are obtained from the interpretation of the tangential fields as effective sources that radiate the scattered fields.)</li>
<li>The total power scattered by, and the total power absorbed by, the scattering objects from the incident field.</li>
<li>The total force and/or torque exerted on the scattering objects by the incident fields (radiation pressure).</li>
<li>Visualization files plotting the electric and magnetic surface currents, and the associated charge densities, induced by the incident fields on the scattering objects.</li>
<li>Visualization files plotting field components and Poynting fluxes on arbitrary user-specified surface meshes.</li>
</ul>
<p></p>
<table>
<col width="100%" />
<thead>
<tr class="header">
<th align="left">Table Of Contents</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><a href="scuff-EM/scuff-scatter/scuffScatterOptions.shtml">1. scuff-scatter Command-Line Options</a></td>
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
</table></td>
</tr>
</tbody>
</table>


