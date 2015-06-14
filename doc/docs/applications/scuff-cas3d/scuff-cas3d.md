<h1>Computing thermal-equilibrium Casimir energies, forces, 
    and torques with <span class="SC">scuff-cas3d</span>
</h1>

[[scuff-cas3d]] is a command-line application within the [[scuff-em]] suite
for modeling Casimir interactions between compact or extended homogeneous
bodies of arbitary shape and arbitrary (linear, isotropic, piecewise homogeneous)
frequency-dependent permittivity and permeability. [[scuff-cas3d]] implements the 
[*fluctuating-surface current (FSC)* approach][FSCPaper]
to numerical Casimir modeling.]

[[scuff-cas3d]] handles *equilibrium* Casimir interactions, in which
all interacting bodies and the external medium in which they are embedded 
exist at the same temperature (which may be absolute zero). If you
need to model Casimir interactions between bodies at *different* 
temperatures, the tool you want is 
[scuff-neq](../scuff-neq/scuff-neq.md){.SC}.
(However: if your exterior embedding medium is not at zero 
temperature, then the *total* Casimir forces will involve
both non-equilibrium contributions computed with [[scuff-neq]] 
and equilibrium contributions (at the temperature of the 
exterior medium) computed with [[scuff-cas3d]].

The basic flow of a typical [[scuff-cas3d]] run goes something like this:

+ You create a [<span class="SC">scuff-em</span> geometry file][Geometries]
describing the interacting objects or surfaces in your geometry.

+ Optionally, you define a 
[list of geometric transformations][Transformations] 
to be applied to the geometry for Casimir computations.
For example, if your geometry consists of two nanoparticles, you
might ask for the Casimir force between the particles at 10 different
values of the surface--surface separation.

+ You run [[scuff-cas3d]] with various
[command-line options](#CommandLineOptions) specifying the
desired output quantities (energy, $y$-force, etc.), whether
you want frequency-resolved or frequency-summed data, and
other options. This produces various text-based 
[output files](#OutputFiles), which you will typically
plot using [gnuplot](www.gnuplot.info){.SC} or other 
plotting or post-processing tools.

[TOC]  

# 1. What [[scuff-cas3d]] actually computes

## 1a. Compact objects

The Casimir energy $\mathcal{E}$ of a collection of compact bodies,
and the $i$-directed Casimir force $\mathcal{F}_i$ and 
torque $\mathcal{T}_i$ on one of 
those bodies, are computed in the FSC approach as integrals 
over the positive imaginary frequency axis
($\omega=i\xi$) of the form

$$ \begin{array}{ccc}
   \mathcal{E}   &=& \int_0^\infty E(\xi) \, d\xi \\
   \mathcal{F}_i &=& \int_0^\infty F_i(\xi) \, d\xi \\
   \mathcal{T}_i &=& \int_0^\infty T_i(\xi) \, d\xi \\
   \end{array}
$$

(This is for the zero-temperature case; at a finite temperature
$T$ the imaginary-frequency integration is replaced by a 
Matsubara sum according to the prescription
$$ \int_0^\infty F(\xi) \, d\xi 
   \qquad \Longrightarrow \qquad
   (\Delta \xi) \sideset{}{'}{\sum}_{n=0}^\infty F(n\Delta\xi),
   \qquad \Delta\xi =\frac{2\pi k T}{\hbar}
$$
where the primed sum indicates that the $n=0$ term is 
to be weighted by $1/2.$)

The heart of the FSC algorithm implemented by [[scuff-cas3d]]
is an efficient technique for computing the quantities
$\{E,F_i,T_i\}(\xi)$---that is, the contributions of 
individual imaginary frequencies to the total Casimir 
quantities---at arbitrary frequencies $\xi$.
The question of *which* frequencies $\xi$ are sampled
depends on the command-line options you specify:

+ If you use the `--Xi` or `--XiFile` command-line options to 
specify one or more particular values of $\xi$, then [[scuff-cas3d]]
will compute and report just the integrand $F(\xi)$ at those
values. In this case, the code will produce a `.byXi` file,
but no `.out` file.

+ If, instead, you use the `--temperature` command-line option to specify 
a temperature at which to calculate, then [[scuff-cas3d]]
will evaluate the Matsubara sums to compute the full 
Casimir quantities at the given temperature. In this 
case you will get both a `.byXi` and a `.out` file; 
the `.byXi` file will report data on the Casimir
integrands at the Matsubara frequencies
$\xi_n\equiv n\Delta \xi$ for n=0,1,\cdots$

+ If you don''t specify any of the above, then [[scuff-cas3d]
defaults to performing a full numerical frequency integration
to compute zero-temperature Casimir quantities. In 
this case, you will get both a `.byXi` and a `.out`
file; the `.byXi` file reports data on the Casimir
integrands at the frequencies $\{\xi\}$ chosen by 
the built-in integrator.

## 1b. Extended objects

For an extended material configuration described by a
[periodic geometry](../../reference/Geometries.md#Extended)
with Bloch-periodic boundary conditions, the Casimir force 
*density* (that is, the force **per unit length** for a 1D 
extended geometry, or **per unit area** for a 2D extended 
geometry) is computed in the FSC approach as multi-dimensional 
integrals over both imaginary frequencies $\xi$ and Bloch 
vectors $\mathbf{k}$:

$$ \mathcal{F}
   = 
   \int_0^\infty d\xi \,
    \underbrace{\int_{\text{BZ}} f(\xi, \mathbf{k}) d\mathbf{k}}_{\equiv F(\xi)}
$$

(Expressions for the Casimir energy and torque are similar).

<a name="CommandLineOptions"></a>
# 2. <span class="SC">scuff-cas3d</span> command-line options

<a name="OutputFiles"></a>
# 3. <span class="SC">scuff-cas3d</span> output files

<a name="Examples"></a>
# 4. Examples of Casimir calculations using <span class="SC">scuff-cas3d</span>

Hello

[EarlierVersion]: http://homerreid.com/scuff-em/scuff-cas3d
[FSCPaper]: http://dx.doi.org/10.1103/PhysRevA.88.022514
[Geometries]: ../reference/Geometries.md
[Transformations]: ../reference/Transformations.md
