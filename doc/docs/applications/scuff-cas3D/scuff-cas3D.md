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
$\xi_n\equiv n\Delta \xi$ for $n=0,1,\cdots$

+ If you don't specify any of the above, then [[scuff-cas3d]]
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
    \underbrace{\int_{\text{BZ}} 
     f(\xi, \mathbf{k}) 
     \, d\mathbf{k}}_{\equiv F(\xi)}
$$

(Expressions for the Casimir energy and torque are similar).
The Bloch vector $\vb k$ is a one-component vector
for 1D-extended geometries (such as 
[infinite-length cylinders or beams][SiliconBeams],
and a two-component vector
for 2D-extended geometries (such as 
[infinite-area slabs][SiliconSlabs]).
The $\vb k$ integral here ranges over the *Brillouin
zone* (BZ)

[[scuff-cas3d]] uses the FSC algorithm to compute
values of the integrand $f(\xi, \vb k)$ at 
individual $(\xi, \vb k)$ points.

The FSC algorithm implemented by [[scuff-cas3d]]
is an efficient technique for computing the quantity
$ f(\xi, \vb k)$---that is, the contributions of
individual (imaginary frequency, wavevector) pairs
$(\xi,\vb k)$---to the total Casimir quantities.
The question of *which* ($\xi, \vb k)$ points are 
sampled depends on the command-line options you specify:

+ If you use the `--XikBlochFile` command-line option
to specify a list of $(\xi, \vb k)$ points, then 
[[scuff-cas3d]] will compute and report just the 
integrand $f(\xi, \vb k)$ at those
values. 
In this case, the code will produce a `.byXikBloch` file,
but no other output files.

+ If you use the `--Xi` or `--XiFile` command-line options to
specify one or more particular values of $\xi$ (but not
specific values of $\vb k$$), then [[scuff-cas3d]] will
numerically evaluate the Brillouin-zone integral over 
$\vb k$ and will report the resulting value of the 
quantity $F(\xi)$ at each $\xi$ value. 
In this case, the code will produce two output files
**(1)** a `.byXi` file reporting values of $F(\xi)$ at the 
$\xi$ points you specified, and **(2)** a `.byXikBloch`
file reporting values of the $\vb k$ integrand $f(\xi, \vb k)$
at each of the $\vb k$ points sampled by the 
built-in integrator.

+ If you use the `--temperature` command-line option to specify
a temperature at which to calculate, then [[scuff-cas3d]]
will evaluate the Matsubara sums to compute the full 
Casimir quantities at the given temperature. In this 
case you will get three output files: 
**(1)** a `.out` file reporting the full Matsubara-summed
Casimir quantities,
**(2)** a `.byXi` file reporting values of the function
$F(\xi)$ at each Matsubara frequency; and 
**(3)** a `.byXikBloch`
file reporting values of the $\vb k$ integrand $f(\xi, \vb k)$
at each point sampled by the built-in integrator.

+ If you don't specify any of the above, then [[scuff-cas3d]]
defaults to performing a full numerical frequency integration
to compute zero-temperature Casimir quantities. In
this case, you will get the same three output files
as in the case of the previous item (`.out`, `.byXi`, `.byXikBloch`);
the only difference is that the $\xi$ points reported
in the `.byXi` and `.byXikBloch` files are the quadrature 
points chosen by the built-in integrator instead of the 
Matsubara frequencies.

<a name="CommandLineOptions"></a>
# 2. <span class="SC">scuff-cas3d</span> command-line options

## Common options

[[scuff-cas3d]] recognizes the following subset of the 
[list of commonly accepted options to <span class="SC">scuff-em</span> command-line codes][CommonOptions].

+ `--geometry`
+ `--TransFile`
+ ` `
+ `--Xi`
+ `--XiFile`
+ `--XikBlochFile`
+ `--XiQuadrature`
+ `--XiMin`
+ ` ` 
+ `--BZQuadrature`
+ `--BZSymmetry`
+ `--MaxBZSamples`
+ ` `
+ `--AbsTol`
+ `--RelTol`
+ ` ` 
+ `--FileBase`
+ ` `    
+ `--Cache`
+ `--ReadCache`
+ `--WriteCache`
{.toc}

## Options requesting Casimir output quantities

`--Energy`
`--XForce`
`--YForce`
`--ZForce`
`--Torque ax ay az`
{.toc}

Specifies the Casimir quantities in which you are
interested: the energy, the Cartesian components
of the force, or the torque about an axis
passing through the origin and the point with 
Cartesian coordinates `(ax,ay,az).` (Thus, for

You may specify more than one of these options,
but you must specify at least one.

**Note**: [[scuff-cas3d]] always computes the force and torque
on just **one** of the objects or surfaces in your geometry---namely,
the one described by the first `OBJECT` or `SURFACE` 
specification in your `.scuffgeo` file. For geometries
consisting of just two objects or surfaces, the force/torque
on the second object/surface is just the negative of the 
force/torque on the first. 

## Options specifying temperature

`--Temperature 300`
{.toc}

Sets the simulation temperature **in units of Kelvin,**
so that `--temperature 300` requests room-temperature calculations.
This option implies that you are asking [[scuff-cas3d]]
to compute full Matsubara-summed Casimir quantities, so it is
incompatible with options such as `--Xi` or `--XiFile` that
specify particular frequencies at which to compute.

<a name="OutputFiles"></a>
# 3. <span class="SC">scuff-cas3d</span> output files

The base file name of all output files produced by
`scuff-cas3d` may be specified using the `--FileBase`
command-line option; if this option is not specified
then the file base is taken to be the base file
name of the `.scuffgeo` file you specified using
the `--geometry` option.

For all data output files (`.out`, `.byXi`,
`.byXikBloch`), the output file contains a
*header* 

## The `.log` file

Like all command-line codes in the [[scuff-em]] suite,
[[scuff-cas3d]] writes a [`.log` file][LogFiles] that you 
can monitor to keep track of your calculation's progress.

## The `.out` file

If you requested the calculation of
full frequency-integrated or Matsubara-summed
Casimir quantities, these will be written to 
the `.out` file.

## The `.byXi` file

For any problem involving compact geometries,
and for any problem involving extended geometries
in which you requested Brillouin-zone integrations,
the contributions of 
individual imaginary frequencies $\xi$ will be 
written to a file named `.byXi`.

## The `.byXikBloch` file

For any problem involving extended geometries,
the contributions of individual (frequency, Bloch vector)
points $(\xi,\vb k)$ will be written to a file
named `FILEBASE.byXikBloch`.

<a name="Examples"></a>
# 4. Examples of Casimir calculations using <span class="SC">scuff-cas3d</span>

+ Casimir forces in a compact geometry:
  [A cube and a torus immersed in ethanol][CubeTorus]

+ Casimir forces in a 1D extended geometry: 
  [infinite-length silicon beams][SiliconBeams]

+ Casimir forces in a 2D extended geometry: 
  [infinite-area silicon slabs][SiliconSlabs]

[EarlierVersion]: http://homerreid.com/scuff-em/scuff-cas3d
[FSCPaper]: http://dx.doi.org/10.1103/PhysRevA.88.022514
[Geometries]: ../../reference/Geometries.md
[Transformations]: ../../reference/Transformations.md
[CubeTorus]: http://homerreid.dyndns.org/scuff-EM/scuff-cas3D/scuff-cas3D-Tutorial.shtml
[SiliconBeams]: ../../examples/SiliconBeams/SiliconBeams.md
[SiliconSlabs]: ../../examples/SiliconSlabs/SiliconSlabs.md
[LogFiles]:     ../GeneralReference.md#LogFiles
