# Power, force, and torque in Mie scattering

This test validates the algorithms implemented in
[[scuff-em]] for computing the power, force, and torque (PFT)
on bodies irradiated by external fields. The test 
uses the 
[<span class="SC">scuff-scatter</span>][scuff-scatter] application module
to compute the PFT for a lossy dielectric sphere
irradiated by a circularly-polarized plane wave
and compares the results to the predictions of 
Mie scattering theory.

## Exact solution

We consider a sphere centered at the origin irradiated by
a left-circularly-polarized plane wave at frequency $\omega$,
$$ \mathbf{E\sups{inc}}(\vb x)=
    \frac{1}{\sqrt{2}E_0\Big( \mathbf{\hat{x}} + i\mathbf{\hat{y}\Big)
    e^{ikz}
$$
with $k=\omega/c$ the vacuum wavenumber.

For a dielectric sphere of radius $R$, Mie theory
(see, for example,
[Bohren and Huffman][BohrenHuffman] and
[Marston and Crichton][MarstonCrichton]) predicts that
the cross sections for extinction, scattering, absorption,
force, and torque are
$$ \begin{array}{ccc}
  C_{\text{\tiny{ext}}}
  &=&
 \frac{2\pi}{k^2} \sum (2n+1)\text{Re }(a_n + b_n)
 \\[5pt]
 \end{array}
$$
where $a_n, b_n$ are the Mie scattering coefficients:

## <span class="SC">scuff-em</span> solution

The power, force, and torque on a gold sphere of radius
1 $\mu$ m may be computed with [[scuff-scatter]] by placing the 
following content into a little file called `scuff-scatter.args`:

````bash
 % geometry GoldSphere_501.scuffgeo
 % PFTFile  GoldSphere_501.PFT
 % OmegaFile MyOmegaFile
 % pwDirection    0 0 1
 % pwPolarization 0.7071 0.7071i 0
````

Here the file
[`E10HalfSpace_40.scuffgeo`](E10HalfSpace_40.scuffgeo)
describes the [<span class="SC">scuff-em"</span> geometry][scuffEMGeometries] 
(it refers to a mesh file named
[`Square_40.msh`](Square_40.msh)) 
and the command-line arguments ask for a calculation at 
angular frequency $\omega=1\cdot 3\times 10^{14}$ rad/sec
and at 20 incident angles in the range $0\le \theta\le 88$ degrees.

## Comparison

Running the above command yields the file
[`E10HalfSpace_40.transmission`](E10HalfSpace_40.transmission).
Plotting in [<span class="SC">gnuplot</sc>][gnuplot] yields 
a comparison of [[scuff-transmission]] data (point) to 
theoretical predictions (curves):

![FresnelData.png](FresnelData.png)

Here is the [[gnuplot]] script that I use to produce this 
plot: [PlotFresnelData.gp](PlotFresnelData.gp).

[scuffScatter]:               ../../applications/scuff-scatter/index.md
[scuffEMGeometries]:          ../../reference/Geometries.md
[scuffEMTransformations]:     ../../reference/Transformations.md
[scuffEMMaterials]:           ../../reference/Materials.md
[scuffEMInstallation]:        ../../reference/Installation.md
[EmigPaper]:                  http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.170403
[gnuplot]:                    https://www.gnuplot.info
[BohrenHuffman]:              http://onlinelibrary.wiley.com/book/10.1002/9783527618156
[MarstonCrichton]:            http://journals.aps.org/pra/abstract/10.1103/PhysRevA.30.2508
