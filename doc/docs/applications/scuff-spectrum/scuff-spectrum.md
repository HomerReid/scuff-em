<h1> Computing electromagnetic modes with <span class="SC">scuff-spectrum</span>
</h1>

[[scuff-spectrum]] is a tool within the [[scuff-em]] code suite
for computing frequencies and field patterns of resonant modes
in arbitrary material geometries. More specifically, [[scuff-spectrum]]
implements [Beyn's contour-integral algorithm for nonlinear
eigenproblems][BeynMethod], described in more detail below,
to pinpoint frequencies $\omega$ at which the
[BEM system matrix $\mathbf{M}(\omega)$][PMCHWTSystem] is singular,
allowing nonzero surface currents to flow in the absence of
external fields. [[scuff-spectrum]] reports mode frequencies
and, optionally, visualization diagrams and various
types of information on the current and field patterns
of the corresponding eigenvectors.

Of course, the information reported by [[scuff-spectrum]] could
also be obtained *indirectly* from other codes in the
[<span class=SC>scuff-em</span> application suite][ApplicationSuite]---for
example, the resonant modes of the
[spherical dielectric cavity][SphericalDielectricCavity]
shown below would show up as blips in (say) the
scattering cross section as the frequency is scanned past
the vicinity of the resonance.
<!--# while the [circuit-QED resonator][CircuitQEDResonator]-->
However, this would be an inefficient and imprecise way to
determine resonance frequencies---we would probably wind up
computing scattering cross sections at hundreds of frequencies
on a dense grid just to approximate the resonance with 
one- or two-digit accuracy---and it would not yield
direct information on the current and field patterns
of the corresponding modes. In contrast, the method
used by [[scuff-spectrum]] is optimized to pinpoint
mode frequencies with many-digit accuracy---using
only a small number of samples---and, moreover, directly
furnishes all information we may need about the current a
and field patterns of the modes themselves.

[TOC]

<a name="BeynMethod"></a>
# 1. Beyn's algorithm: the computational engine behind <span class="SC">scuff-spectrum</span>

### Mode frequencies as nonlinear eigenvalues

The basic problem solved by <span class=SC>scuff-spectrum</span>
is to find values of the angular frequency $\omega$ at which the
[BEM system matrix $\mathbf{M}(\omega)$][PMCHWTSystem] is singular
(has non-empty nullspace), in which case we can find a 
nonzero vector $\mathbf{v}$ that solves the equation

$$ \mathbf{M}(\omega) \mathbf{v} = 0. \qquad (1) $$

This is just the 
[usual linear system of equations solved in <span class=SC>scuff-em</span> scattering problems][PMCHWTSystem],
but now with zero right-hand side, i.e. there is no incident field or
other external stimulus exciting the structure; instead,
the eigenvector $\mathbf{v}={\mathbf{k} \choose \mathbf{n}}$
describes a configuration of surface currents that can exist
in the absence of external fields---a resonance mode.

Equation (1) resembles an eigenvalue problem, but with the
complication that the matrix $\mathbf{M}(\omega)$ depends
in a complicated nonlinear way on the frequency $\omega$; 
this means we can't simply use standard solvers
for linear eigenproblems such as those implemented
in [<span class=SC>lapack</span>](http://www.netlib.org/lapack/).

### Beyn's contour-integral algorithm: A lightning overview

Instead, [[scuff-spectrum]] implements the contour-integral approach
to nonlinear eigenproblems proposed by W. Beyn in this 2012 paper:

+ Wolf-J&uuml;rgen Beyn, "An integral method for solving nonlinear eigenvalue problems." *Linear Algebra and its Applications* **436** 3839 (May 2012).

+ DOI: [https://doi.org/10.1016/j.laa.2011.03.030](https://doi.org/10.1016/j.laa.2011.03.030)

+ ArXiV: [https://arxiv.org/abs/1003.1580](https://arxiv.org/abs/1003.1580)

In a nutshell, Beyn's method locates eigenvalues of (1) by
evaluating a certain contour integral in the complex $\omega$ plane;
the evaluation proceeds via numerical quadrature, sampling
the integrand at a total of $N$ quadrature points
$\{\omega_n\},n=1,2,\cdots,N$ (see figure below).
The number of quadrature points, and the shape and location
of the contour, are user-tweakable parameters that you 
will specify as inputs to [[scuff-spectrum]]; the basic idea
is that Beyn's method will identify all eigenvalues located 
*inside* the contour, and will ignore any eigenvalues lying
*outside* the contour, so the challenge from the user's 
perspective is to have a rough idea of where your eigenvalues
will be---and to design contours that enclose the ones you
want while omitting the others.

We are being deliberately vague here about precisely
*what* contour integral we are evaluating, as this is not
crucial knowledge for using [[scuff-spectrum]] (for details,
see the original paper linked above). Pretty much all you need
to know is this: Evaluating the integrand at each
quadrature point $\omega_n$ involves **(a)** forming and factorizing 
the matrix $\mathbf{M}(\omega_n)$, and then **(b)** doing some simple
linear-algebra calculations involving the solution
of a small number (see below) of linear systems of equations
of the form $\mathbf{M}(\omega)\mathbf{x}=\mathbf{b}$
for certain given RHS vectors $\mathbf{b}.$ Of course, 
assembling the system matrix $\mathbf{M}(\omega)$ and solving
linear systems is exactly what <span class=SC>scuff-em</span>
does to solve scattering problems, and so *mechanically*
Beyn's algorithm boils down simply to solving integral-equation
scattering problems, at some number of frequencies, with
some number of RHS vectors at each frequency---not different,
in principle, from what we would wind up doing in the brute-force
(frequency-scanning) approach to computing mode frequencies
discussed above. The difference is that Beyn's algorithm
chooses the frequencies at which we calculate, and combines
the results of all those calculations, in clever ways to ensure
that maximal benefit is extracted from the computations,
allowing us to converge quickly to highly accurate values
for the eigenvalues $\omega$ in (1). Of course, in addition
to the eigenvalue $\omega$, Beyn's method also gives us
the eigen*vectors* $\mathbf{v}$, which we can use in
post-processing to do things like visualizing the distribution
of currents and fields in the various eigenmodes we find (see examples
below).

### Beyn's algorithm: Mechanics
Having briefly outlined the idea of Beyn's method---and referring
interested readers to Beyn's paper cited above for more detail on the 
theory---let's now turn to a discussion of the mechanics
of running Beyn-method calculations in [[scuff-spectrum]].

The following figure illustrates a typical elliptical
contour $\mathcal{C}$ in the complex-$\omega$ plane over
which we might wish to execute Beyn's algorithm.
(Contours in [[scuff-spectrum]] are always elliptical.)

![BeynContour](BeynContour.png)

In this figure, blue crosses indicate eigenvalues of (1);
note that there are 3 eigenvalues contained inside the
contour $\mathcal{C}$---labeled $\omega_{1,2,3}$---as
well as several other eigenvalues not enclosed by the 
contour, which we do not label as they will be ignored
by our calculation. The black dots on the contour indicate
quadrature points; the $\omega$ values represented by these points
are the complex frequencies at which [[scuff-scatter]]
will assemble the matrix $\mathbf{M}(\omega)$
and do some linear algebra to compute 
the integrand of the Beyn-method contour integrals.

Items in *red* in the figure above indicate user-specified
inputs to [[scuff-scatter]]: these are

+ the complex-valued frequency $\omega_0$ at which
  the contour $\mathcal{C}$ is centered

+ the horizontal and vertical radii (half-minor axes)
  $R_x$ and $R_y$ of the elliptical contour

+ the number of quadrature points $N$ (12 in this case)

+ an integer $L$ which should be greater than or equal to
  the number of eigenvalues you *expect* to be found
  within the contour; if [[scuff-spectrum]] finds more than
  this number of eigenvalues, Beyn's method breaks down and 
  must be restarted with a larger value of $L$. (You will 
  get a console warning in this case.)

  On the other hand, the algorithm works fine if $L$ winds up
  being *larger* than the number of eigenvalues found inside
  the contour, so you should give yourself some leeway by 
  choosing generous values for $L$ (typical values might 
  be 10 or 20); higher values of $L$ do result in slightly
  greater computation time, but not much. (More specifically, $L$
  is the number of linear-system solves of the form
  $\mathbf{M}(\omega_n)\mathbf{x}=\mathbf{b}$ that must be done
  at each quadrature point $\omega_n$; increasing $L$ from 10 to 20
  or 30 or so does require more solves, but that cost is 
  generally negligible compared to the cost of assembling 
  and factorizing the matrix $\mathbf{M}(\omega)$, so don't
  bother with excessive parsimony in choosing $L$.)

Running [[scuff-spectrum]] with the above inputs as illustrated
in the figure above would yield accurate values for the three 
eigenvalues $\omega_{1,2,3}$ inside the contour---with the
accuracy increasing with the number of quadrature points $N$---as
well as values for the corresponding eigenvectors.

<a name="SphericalDielectricCavity"></a>
# 2. <span class="SC">scuff-spectrum</span> tutorial: Modes of a 
     spherical dielectric cavity

<a name="CommandLineReference"></a>
# 3. <span class="SC">scuff-spectrum</span> command-line reference

[ApplicationSuite]:            ../../index.md/#ApplicationSuite
[SphericalDielectricCavity]:   ../../applications/scuff-spectrum/scuff-spectrum.md/#SphericalDielectricCavity
[BeynMethod]:                  ../../applications/scuff-spectrum/scuff-spectrum.md/#BeynMethod
[PMCHWTSystem]:                ../../forDevelopers/Implementation.md/#PMCHWTSystem
