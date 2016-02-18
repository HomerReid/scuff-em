# Brillouin-zone integration in [[scuff-em]]

Many codes in the [[scuff-em]] suite require evaluating
integrals over the Brillouin zone (BZ) of a 1D or 2D
reciprocal lattice, i.e.

$$ Q(\omega) = \int_{\hbox{BZ}} \mathcal{Q}(\omega, \bf{k}) \,d\bf{k}.$$

Examples of calculations that require
Brillouin-zone integrations include

+ the Casimir force per unit imaginary frequency
$\omega=i\xi$ 
on an extended object in
[<span class="SC">scuff-cas3D</span>][scuffCas3D]

+ the Casimir-Polder potential per unit imaginary frequency
$\omega=i\xi$ 
experienced by a polarizable particle near an extended surface
in 
[<span class="SC">scuff-caspol</span>][scuffCaspol]

+ the local density of states at a given angular frequency
$\omega$
at user-specified evaluation points in 
[<span class="SC">scuff-ldos</span>][scuffldos].

In general, Brillouin-zone integrations are evaluated
by numerical cubature---that is, as weighted sums of
integrand samples:
$$ Q(\omega) \approx \sum W_n \mathcal Q(\omega, \bf k_n)$$
where $\{W_n, \bf k_n\}$ are the weights and points in
a cubature rule for the Brillouin zone of your reciprocal
lattice, and where each integrand sample
$\mathcal Q(\omega, \bf k_n)$ is computed by performing a
single [[scuff-em]] calculation at a fixed Bloch
wavevector. The [[scuff-em]] workflow offers two 
options for evaluating such cubatures:   

+ You can design and implement your own cubature scheme
involving your own custom-chosen weights and points
$\{W_n, \bf k_n\}$. In this case, you will typically
use the ``--byOmegakBloch`` command-line option to instruct
a [[scuff-em]] application code to report values of the
quantity $\mathcal Q(\omega, \bf k_n)$ at each of your points
(this output will typically be written to file with
extension ``.byOmegakBloch`` or ``.byXikBloch``),
then compute the weighted sums yourself in e.g.
[<span class="SC">julia</span>](http://julialang.org).    

+ Alternatively, you can ask [[scuff-em]] to perform the
BZ integration internally, using one of several
built-in cubature schemes. In this case the BZ-integrated
quantities $Q(\omega)$ will typically be written to 
an output file with extension `.byOmega` or `.byXi`.
You will *also* get an output file named `.byOmegakBloch`
or `.byXikBloch` that reports the Bloch-vector-resolved
integrand samples $\mathcal Q(\omega, \bf k_n)$
chosen internally by the [[scuff-em]] BZ integrator.   

### Command-line options for customizing built-in BZ integration

If you choose the second option above, you may use the following
set of command-line options, recognized by all application
codes in the [[scuff-em]] suite, to customize the cubature
scheme used by [[scuff-em]] to evaluate the Brillouin-zone integral.

  ````
--BZIMethod TC
--BZIMethod CC
--BZIMethod Adaptive
  ````
{.toc}

This option specifies the type of cubature rule used to
evaluate the Brillouin-zone integral:

+ For geometries with two-dimensional lattice periodicity,
  `--BZIMethod TC` selects a "triangle cubature" rule in which
  the Brillouin zone is divided into two triangles, over which
  the integral is evaluated using a fixed-order cubature scheme.
  This method is the default for geometries with 2D periodicity.

+ `--BZIMethod CC` selects fixed-order Clenshaw-Curtis quadrature.
  This option is available for geometries of both 1D and 2D
  periodicity.
  For 2D lattices, the BZ integral is evaluated
  using two nested 1D Clenshaw-Curtis quadrature rules.
  This method is the default for geometries with 1D periodicity.

+ `--BZIMethod Adaptive` selects an adaptive cubature scheme
  that automatically choose sample points and keeps going
  until a specified convergence tolerance is achieved.
  This option is available for geometries of both 1D and 2D
  periodicity.    

  ````
--BZIOrder NN
  ````
{.toc}
 
  For integration methods that use fixed-order cubature rules
  (`--BZIMethod TC` or `--BZIMethod CC`), this option may be 
  used to specify the order of the cubature rule.  

  For Clenshaw-Curtis cubature (`--BZIMethod CC`), the order `NN`
  may be any odd integer between 9 and 99. (In this case 
  the "order" actually refers to the number of sample points in the
  cubature rule.)    

  For triangle cubature (`--BZIMethod TC`), the order `NN`
  may be 1, 2, 4, 5, 7, 9, 13, 14, 16, 20, or 25. In this
  case the "order" is the degree of polynomial that is integrated
  exactly by the cubature rule, not the number of sample points.

  ````
--BZIRelTol     1.0e-2
--BZIMaxEvals   1000
  ````
{.toc}

  These options may be used to specify stopping criteria for
  adaptive integration methods (`--BZIMethod adaptive`).
  `--BZIRelTol` sets the relative convergence tolerance,
  while `--BZIMaxEvals` sets the upper limit on the number of 
  integrand samples the integrator may use.

  ````
--BZSymmetric
--FullBZ
  ````
{.toc}

  These options may be used to describe symmetries of the integrand function $\mathcal{Q}(\bf k)$.

  For geometries with 2D square-lattice periodicity, `--BZSymmetric` may be used to specify that the Brillouin-zone integrand is symmetric under interchange of $k_x$ and $k_y$ (i.e. $\mathcal{Q}(k_x, k_y)=\mathcal{Q}(k_y, k_x)$, saving a factor of approximately 2 in evaluating the BZ cubature. 

[Note that symmetry requires both that **(a)** the surface meshes in the geometry are symmetric under inversion about the line $x=y$, and **(b)** the integrand $\mathcal{Q}(\bf k)$ is symmetric under inversion about the line $k_x=k_y.$ Condition **(b)** is not always satisfied even for geometries that satisfy condition **(a)**; an example would be the calculation of LDOS at a point $\bf x$ that does not lie on the line $x=y.$]

  The default behavior of [[scuff-em]] is to integrate over just the *irreducible* portion of the Brillouin zone and double or quadruple the result to get the full BZ integral. You can override this behavior by specifying `--FullBZ,` in which case the integral over the full Brillouin zone will be computed at 2x or 4x greater cost.

For a 1D lattice with lattice constant $L$, 
the full Brillouin zone is the region $0\le k \le \frac{2\pi}{L}$, and the irreducible Brillouin zone is the region $0\le k \le \frac{\pi}{L}.$
For a 2D rectangular lattice with direct lattice constants $L_x, L_y$,
the full Brillouin zone is the rectangle in the $(k_x,k_y)$ plane with corners at the origin and at 
$(k_x,k_y)=(\frac{2\pi}{L_x},\frac{2\pi}{L_y})$,
while the irreducible Brillouin zone is the rectangle with corners at the origin and at 
$(k_x,k_y)=(\frac{\pi}{L_x},\frac{\pi}{L_y})$.

[scuff-analyze]:               ../applications/scuff-analyze/scuff-analyze.md
[scuffCas3D]:                  ../applications/scuff-cas3D/scuff-cas3D.md
[scuffCaspol]:                 ../applications/scuff-caspol/scuff-caspol.md
[scuffldos]:                   ../applications/scuff-ldos/scuff-ldos.md
