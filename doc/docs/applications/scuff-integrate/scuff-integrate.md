# Integrating frequency-dependent data with <span class="SC">scuff-integrate</span>

Many application codes in the [[scuff-em]] suite compute
physical quantities defined by definite integrals over real
or imaginary frequencies, with the numerical value of the
integrand at each point obtained by solving individual
[[scuff-em]] scattering problems at that frequency. For example,

+ [ <span class=SC>scuff-cas3d</span> ][scuff-cas3D] and
  [ <span class=SC>scuff-caspol</span> ][scuff-caspol]
  compute zero-temperature Casimir quantities by integrating 
  contributions from imaginary frequencies $\xi$:
<!--%====================================================================%-->
$$ Q = \int_0^\infty I(\xi) \, d\xi \tag{1} $$
<!--%====================================================================%-->
  Here $Q$ is a zero-temperature Casimir energy/force/torque ([[scuff-cas3d]])
  or Casimir-Polder potential ([[scuff-caspol]]) and $I(\xi)$ is the 
  spectral density of contributions to $Q$ from fluctuations at 
  imaginary frequency $\xi$, which may be obtained by solving
  [[scuff-em]] scattering calculations at imaginary frequency $\xi$.

+ [ <span class=SC>scuff-neq</span> ][scuff-neq] computes
    the total rate of energy or momentum transfer
    from a source body $s$ to a destination body $d$ 
    by integrating contributions from real frequencies $\omega$:
<!--%====================================================================%-->
$$ Q_{s\rightarrow d}
   = \int_0^\infty
     \Big[ \Theta(T_s,\omega) - \Theta(T_\text{env},\omega) \Big]
     \Phi_{s\rightarrow d}(\omega) \, d\omega
   \tag{2}
$$
<!--%====================================================================%-->
    Here $Q_{s\rightarrow d}$ is the contribution of body $s$ to
    quantity $Q$ (a heat-transfer rate, force, or torque) for body $d$,
    $\Theta(T,\omega)=\frac{\hbar\omega}{e^{\hbar\omega/kT}-1}$
    is the Bose-Einstein statistical factor, $T_s$ and $T_\text{env}$
    are the temperatures of the source body and the environment,
    and $\Phi(\omega)$ is a "generalized flux" quantity that may be 
    computed by solving [[scuff-em]] scattering calculations
    at frequency $\omega$.
    
The integrals over $\xi$ and $\omega$ are evaluated by
numerical cubature---that is, as weighted sums of integrand
samples. In a perfect world, it would be possible for 
[[scuff-em]] application codes to choose appropriate 
integration strategies automatically, hiding these details
from the user and reporting just the frequency-integrated
quantities $Q$. This is in fact the strategy that was
adopted in early incarnations of the [[scuff-em]] codes.

In the real world, however, the behavior of integrand 
functions like $I(\xi)$ and $\Phi(\omega)$ varies widely
from problem to problem, depending on factors such
as the shapes and materials of bodies in the scattering
geometry and the quantity being computed. For this reason,
it's hard for [[scuff-em]] to make intelligent automatic
choices of integration strategies, and attempts to do so 
without user input may result in misleading or even
flat-out incorrect data.

For this reason, the modern approach to frequency integration
in [[scuff-em]] is to ask users to define a
list of frequencies at which to sample the integrand; this
list is passed to [[scuff-em]] application codes using the 
`--XiFile` or `--OmegaFile` command-line options, and in 
response the code produces an output file reporting values
of the integrand at the specified points. The frequency 
integral may then be calculated as as post-processing step
using the information reported in the frequency-resolved
data files, and this is the task for which [[scuff-integrate]]
exists.

**<span class=SC>scuff-integrate</span> is a general-purpose tool, not limited to <span class=SC>scuff-em</span> applications**

Actually, [[scuff-integrate]] is designed to be a general-purpose
numerical-integration tool, not particularly tied to
[[scuff-em]]; given a data file containing samples of
some function $f(x)$ at sample points $\{x_n\}$
distributed over an interval $[x_a, x_b]$,
[[scuff-integrate]] approximates the integral
$\int_{x_a}^{x_b} f(x) \, dx$. (It does this by
constructing a third-order spline interpolant through
the data, then integrating the interpolant via numerical
quadrature.) The code offers command-line options
to specify how the data file is to be interpreted
(for example, which columns are the $x$ samples and 
which are the $f$ samples). There is no restriction
on the number or spacing of the sample points; the
variable *x* need not correspond to a real or imaginary 
frequency, and the function $f$ need not be a generalized
flux, Casimir integrand, or any other quantity
reported by a [[scuff-em]] calculation.

**But for data files that *do* come from <span class=SC>scuff-em</span> calculations, <span class=SC>scuff-integrate</span> knows what to do automatically with minimal user input**

In addition to the fully general-purpose mode described above,
[[scuff-integrate]] also has built-in automatic support for
certain specialized types of <span class=SC>scuff-em</span>
data files---such as the [`.SIFlux` files produced by <span class=SC>scuff-neq</span>](../scuff-neq/scuff-neq.md#SIFluxFile)---designed
to make it easy to run the tool to get the frequency-integrated data you want with minimal user input
required.

[TOC]

# 1. Using <span class=SC>scuff-integrate</span> as a general-purpose numerical integrator

## Integrating a single function of a single variable

The simplest usage of [[scuff-integrate]] is to integrate
a single function $f$ of a single variable $x$; as noted
above, this is a general-purpose use of the code, independent
of [[scuff-em]], in which $x$ and $f$ need not be a frequency and
a generalized flux, but could have any arbitrary significance.
Suppose we have a file called `fData` in which
are tabulated numerical pairs $(x_n, f_n)$, $n=1,2,\cdots,N$, 
where $f_n=f(x_n)$:
````bash
x1 f1 
x2 f2 
...
xN fN 
````

Then to compute $\int_{x_1}^{x_N} f(x) \, dx$ we can just 
go like this:
````bash 
% scuff-integrate --datafile fData --freqColumn 1 --dataColumn 2
````

This will produce a file named `fData.Integrated` containing
a single line of data: the integrated value of $f$.

If the frequency and/or integrand values are printed in
different columns of the data file, just adjust the 
`--freqColumn` and `--dataColumn` options accordingly. For
example, if the data file looks like this: 

````bash
<stuff> <stuff> x1 <stuff> <stuff> f1 <stuff> ...
<stuff> <stuff> x2 <stuff> <stuff> f2 <stuff> ...
...
<stuff> <stuff> xN <stuff> <stuff> fN <stuff> ...
````

you would use `--freqcolumn 3 --datacolumn 6.` 
In this case, the content of the other columns (the `<stuff>`
in the above snippet) is ignored.

### Integrating multiple functions of frequency

More generally, you may have multiple integrand functions 
$f_1, \cdots, f_N$, each sampled at the same set of $x$ values,
You can integrate all of these
at once simply by specifying multiple `--dataColumn`
options. For example, if you have functions $f$ and $g$
and you have a file named `fgData` with the format

````bash
x1 f1 g1 
x2 f2 g2 
...
xN fN gN 
````

then you can say
````bash 
% scuff-integrate --datafile fgData --freqColumn 1 --dataColumn 2 --datacolumn 3
````

and the resulting output file `fgData.Integrated`
will report values for both $\int f(x)\,dx$ and $\int g(x)\,dx$.

### Giving names to data columns

In the legend at the top of the `.Integrated` output file,
the values of the various integrated functions will by default
be labeled `data 0`, `data 1`, etc. If you want to give more
descriptive names, just follow each `--dataColumn` option
with a `--dataName` option. 

For example, if 
columns 8 and 11 of your data file respectively
report values of force and torque integrands, you might say
, say 
`--dataColumn 8 --dataName Force --dataColumn 11 --dataName Torque`.

### Integrating functions of frequency and other parameters

In many cases we will have functions that depend on various
parameters beside the integration variable, i.e. functions of the form
$f(x;p_1, \cdots, p_N)$ where the data file includes integrand
samples for multiple sets of values of the $(p_1,\cdots,p_n)$ 
parameters. (In [[scuff-cas3d]], for example,
we might compute Casimir forces between particles separated
by various distances $d$, so the integrand function may 
be thought of as a function $f(x,d)$ of both distance and frequency.)
In these cases you will generally want to evaluate *separate* integrals
over $x$ for each set of parameter values represented in your
data file.

For example, suppose your data file is called `pxfgData` and 
looks something like
````
p1 x1 f11 g11
p1 x2 f12 g12
....
p1 xN f1N g1N
p2 x1 f21 g21
p2 x2 f22 g22
....
p2 xN f2N g2N
....
pM xN fMN gMN
````

where `p1,` `p2,` ..., `pM` denote $M$ distinct values
of some parameter $p$ and `fmn,gmn` are the numerical
integrand values $f(p_m, x_n), g(p_m,x_n)$. 
In this case you can't simply say
`--freqColumn2 --dataColumn 3 --dataColumn 4,` because then
data for all parameter values will be mashed all together and 
integrated as a single function of frequency, yielding nonsense.

Instead, you handle this situation by specifying the additional
command-line parameter `--tagColumn 1` to tell `scuff-integrate`
to interpret data lines with different values in column 1 
as samples of different functions:
````bash 
% scuff-integrate --datafile pxfgData --tagcolumn 1 --freqColumn 2 --dataColumn 3 --dataColumn 4
````

In this case, the output file `pxfgData.Integrated` will report
$x$-integrated values of $f$ and $g$ separately for each value of
$p$.

If your integrands depend on multiple parameters $(p,q,\cdots)$,
you may specify multiple `--tagColumn` options to specify
the columns in which values of the various parameters live.
Then each line of the `.Integrated` output file will report
$x$-integrated values of all functions for a single tuple of 
parameter values $(p,q,\cdots).$

# 2. Using <span class=SC>scuff-integrate</span> as a specialized integrator for <span class=SC>scuff-em</span> data

As noted above, [[scuff-integrate]] has built-in knowledge
of the structure of various data files produced by
<span class=SC>scuff-em</span> calculations. This streamlines
calculations by allowing you
you to feed those data files directly into [[scuff-integrate]]
without needing to tell the code which column of the data file
is which.
 
## 2A. Integrating [[scuff-neq]] data

As discussed at the top of the page, the frequency-resolved
data reported by [<span class=SC>scuff-neq</span>][scuff-neq]
are generalized fluxes describing temperature-independent
rates of energy and momentum transfer from each body
to each other body in your geometry. To turn these
into total heat-transfer rates, forces, and torques
for a given set of object and environment temperatures,
we must **(a)** integrate over frequency with
Bose-Einstein thermal weight factors appropriate for
the given temperatures, and **(b)** sum the
contributions of multiple source bodies to yield
the total rates of power and momentum transfer to
each destination body. 
<span class=SC>scuff-integrate</span> already knows the
file format of the `.SIFlux` and `.SRFlux` frequency-resolved
data files produced by [<span class=SC>scuff-neq</span>][scuff-neq],
so there is no need to specify `--FreqColumn` and `--DataColumn`
options. 

### Specifiying temperatures

Instead, the only input you need to supply (besides
the `.SIFlux` or `.SRFlux` data file) is a specification
of the temperatures of the bodies in your geometry,
and of the surrounding environment. By default, all bodies
and the environment are at absolute zero, $T=0$ K; if
you do not modify this situation by setting a nonzero temperature
for at least one body (or the environment), all heat-transfer
rates and forces/torques will be zero. The command-line option
for setting 

+ `--Temperature N T`

    If $N\ge 1$, this sets the temperature of the $N$th object/surface in the geometry to `T`
    Kelvin.
    (Here objects/surfaces are indexed using a 1-based convention; to set the temperature
    of the first object/surface specified in the `.scuffgeo` file to room temperature
    you would say `--Temperature 1 300.`)

    If $N=0$, this instead sets the temperature of the environment to `T`.

### Specifiying multiple temperature sets

A bonus feature
of the separation between [[scuff-neq]] and [[scuff-integrate]]
is that a single set of frequency-resolved flux data (in
an `.SIFlux` or `.SRFlux` file) can be used to obtain data on
heat-transfer rates and forces and torques at multiple sets
of temperatures for the bodies and the environment, simply by 
evaluating integral (2) multiple times with different temperatures
but the same flux data.


### Spatially-integrated flux data (`.SIFlux` files)

knows how to do this automatically 

In this case, for a geometry containing $N$ bodies,
each line of the `.SIFlux` output file is tagged with
a data field of the form $sd$ (where $s$ and $d$ are
integers between 1 and $N$) to label the contributions
of sources in body $s$ to the power, force, and/or torque
(PFT) on body $d$. (For example, lines for which this field reads
`13` give contributions of body 1 to the PFT for body 3).
The actual data quantities reported in the
`.SIFlux` file are the generalized fluxes $\Phi_{s\rightarrow d}$
in equation (2) above, and to evaluate the $\omega$ integral
here we need to know the temperatures of the environment
and of all bodies in the geometry, which enter through 
the Bose-Einstein factors in (2).

To handle these complications, [[scuff-integrate]]
supports the following additional command-line options:

+ `--sdColumn xx` 

    Specifies that the $sd$ indicator field appears on
    column `xx` of the data file. (The default is `--sdColumn 3`,
    matching the default file format of the `.SIFlux` files
    produced by [[scuff-neq]], so for those files this option
    may be omitted.)

&nbsp;

&nbsp;

+ `--TemperatureFile TFile`

    Specifies a file containing multiple temperature configurations
    at which to compute total PFTs. For an $N$-body geometry,
    each line of `TFile` should contain $N+1$ space-separated
    numbers in the same format as the arguments to the `--Temperature`
    option, i.e. `TEnv T1 ... TN.`
    
    For example, to compute PFTs in a two-body geometry with the
    temperature of body 1 scanned from 10 to 300 Kelvin, the
    temperature of body 2 held fixed at room temperature, and the
    environment temperature fixed at 0, `TFile` would look like

````bash
0 10  300
0 20  300
...
0 300 300
````

# 3. <span class=SC>scuff-integrate</span> Command-Line Reference

# 4. Miscellaneous notes

## Only numerical data columns are counted as columns

There is one potentially confusing aspect of the way
[[scuff-integrate]] interprets column indices as specified
by command-line arguments such as `--FreqColumn` or `DataColumn`.
This is that *[[scuff-integrate]] treats non-numerical data columns
as white space*, and in particular does not include data columns
containing text strings when counting column indices.

Thus, for example, if your data file contains frequency and
integrand data in the second and third columns, with the first
column containing a character string, like this:

````bash
DEFAULT 0.1 3.45e-5
DEFAULT 0.2 7.82e-5
DEFAULT 0.3 1.10e-4
...
````

then [[scuff-integrate]] ignores the `DEFAULT` column and considers
the first column with numerical data to be column 1, so here
you would say `--freqColumn 1 --dataColumn 2.`

In contrast, if your data file looks instead like this: 

````bash
4.00000 0.1 3.45e-5
4.00000 0.2 7.82e-5
4.00000 0.3 1.10e-4
...
````

you would want to say `--freqColumn 2 --dataColumn 3.`

{!Links.md!}
