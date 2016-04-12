# Fresnel Scattering

This test validates the [[scuff-transmission]] application module
of the [[scuff-em]] code suite by using it to study the textbook
case of *Fresnel scattering*: the transmission and reflection of
plane waves at a dielectric interface.

## Exact solution

The situation considered here is that of a plane wave impinging 
from below on a dielectric half-space (relative permittivity $\epsilon=10$)
filling the region $z>0$. (More details on the setup for
[[scuff-transmission]] calculations may be found in the document
[Computation of reflection and transmission coefficients in <span class="SC">scuff-em</span>](../../tex/scuff-transmission.pdf).)

For this case, the transmission and reflection coefficients
for the TE and TM polarizations read 

$$ t^{\scriptsize{\text{TE}}}
   =\frac{2\cos\theta}{\cos\theta + \cos\theta^\prime},
   \qquad
   r^{\scriptsize{\text{TE}}}
   =\frac{\cos\theta-\cos\theta^\prime}{\cos\theta + \cos\theta^\prime},
$$
$$
   t^{\scriptsize{\text{TM}}}
  =\frac{2n\cos\theta}{n^2 \cos\theta + \cos\theta^\prime},
   \qquad
   r^{\scriptsize{\text{TM}}}
  =\frac{n^2\cos\theta-\cos\theta^\prime}{n^2\cos\theta + \cos\theta^\prime}
$$

where $\theta$ is the incident angle ($\theta=0$ for normal incidence), $n=\sqrt{\epsilon}$ is the index of refraction, and 

$$ \cos\theta^\prime = \sqrt{ \epsilon - \sin^2\theta}.$$

## <span class="SC">scuff-em</span> solution

The transmission and reflection coefficients for the 
$\epsilon=10$ dielectric
half-space problem may be computed using [[scuff-transmission]]
as follows:

````bash
 % scuff-transmission --geometry E10HalfSpace_40.scuffgeo --Omega 1.0 --ThetaMin 0.0 --ThetaMax 88.0 --ThetaPoints 20"
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

[scuffEMGeometries]:                  ../../reference/Geometries.md
[scuffEMTransformations]:             ../../reference/Transformations.md
[scuffEMMaterials]:                   ../../reference/Materials.md
[scuffEMInstallation]:                ../../reference/Installation.md
[EmigPaper]:                          http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.170403
[gnuplot]:             ../https://www.gnuplot.info
