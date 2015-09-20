# General reference for [[scuff-em]] command-line applications

This page collects some general information that applies
to many or all of the standalone command-line applications
in the [[scuff-em]] suite.

[TOC]

<a name="CommonOptions"></a>
# 1. Common command-line arguments

The various standalone applications in the [[scuff-em]] suite
share a number of command-line arguments in common,
as described below.
Not all codes accept all arguments (for example, 
[[scuff-transmission]] does not accept ``--TransFile``), but the 
format of each arguments is standardized among all codes that *do* 
accept that argument.

### *Options specifying geometry inputs*

`--geometry MyGeometry.scuffgeo`
{.toc}

> Specifies the 
> [<span class="SC">scuff-em</span> geometry file](../reference/Geometries.md)
> describing your geometry. This option is always mandatory.

`--TransFile MyTransformations.trans`
{.toc}

> Specifies an optional file describing one or more
> [geometrical transformations](../reference/Transformations.md)
> to be applied to your geometry. Omitting the `--TransFile` option
> (when running a code that accepts it) is equivalent to specifying
> an empty transformation named `DEFAULT` that leaves the geometry
> unchanged from the configuration described by the `.scuffgeo` file.

### *Options specifying individual frequencies and Bloch vectors at which to calculate*

  ````
--Omega 3.34
--Omega 4.25+0.9i
--Lambda 2.3
--OmegaFile MyOmegaFile
--LambdaFile MyLambdaFile
  ````
{.toc} 

> Specifies one or more frequencies at which to perform calculations.
>
> `--Omega` specifies the angular frequency in units of 
> $\omega_0=3\cdot 10^{14}$ rad/sec. The argument of `--Omega` may be a
> [complex number](../applications/GeneralReference.md#Complex).
>
> The alternative option `--Lambda` instead specifies the frequency
> in terms of the corresponding free-space wavelength $\lambda=c/\omega$ 
> (where $c$ is the *vacuum* speed of light, irrespective of the material
> properties of your geometry.)
> `--Lambda` values are interpreted in units of microns ($\mu$m).
>
> The options `--OmegaFile` or `--LambdaFile` may be used to specify a 
> file containing one or more `--Omega` or `--Lambda` values, one 
> per line (blank lines and comments are ignored.)


  ````
--OmegakBlochFile MyOkBFile
  ````
{.toc}

> Similar to `--OmegaFile`, but specifies a list of
> (frequency, Bloch vector) points at which to perform calculations.
> This option only makes sense when used with 
> [extended geometries][ExtendedGeometries]
>
> The argument specified for `--OmegakBlochFile` should
> be a file containing two numbers on each line (for 1D extended 
> geometries) or three numbers on each line (for 2D extended 
> geometries). (Blank lines and comment lines beginning with `#`
> are ignored.) The first number on each line is the `--Omega`
> value; the next one or two numbers are the components of the 
> 1D or 2D Bloch wavevector, measured in units of ($\mu$m)$^{-1}$.


  ````
--Xi 0.39 
--XiFile MyXiFile
--XikBlochFile MyXkBFile
  ````
{.toc} 

> Similar to `--Omega`, `--OmegaFile`, and `--OmegakBlochFile,` 
> but used for codes ([scuff-cas3d][scuff-cas3D]{.SC} and
> [scuff-caspol][scuff-caspol]{.SC}) that perform calculations
> at pure imaginary frequencies, $\omega=i\xi$.
> Values specified for `--Xi` should be positive real numbers.

<a name="EvaluationPoints"></a>
### *Options specifying evaluation points*

`--EPFile MyEPFile`
{.toc}

> For codes that compute spatially-resolved output quantities,
> this option specifies a file describing a list of spatial 
> evaluation points. (*Which* output quantity depends
> on the code you are running; for example,
> [scuff-scatter](../applications/scuff-scatter/scuff-scatter.md){.SC}
> will report components of the scattered and total fields at the 
> evaluation points, while
> [scuff-neq](../applications/scuff-neq/scuff-neq.md){.SC}
> will report values of the thermally-averaged fluxes of
> energy and momentum at the evaluation points.)
>
> The argument to `--EPFile` should be a file containing one or more
> lines, each of which contains three space-separated numbers
> (the Cartesian coordinates of the evaluation point). 
> Blank lines and comments (lines beginning with `#`) are ignored.
> For example,
> to ask [[scuff-scatter]] for the scattered field components
> at points on the *z* axis, the file might look like this:

> ````
>     # evaluation points
> 0.0 0.0 -2.0
> 0.0 0.0 -1.9
> ...
> 0.0 0.0  1.9
> 0.0 0.0  2.0
> ````
{.hljs}

### *Options controlling frequency integrations*

  ````
--OmegaQuadrature [adaptive | cliff]
--OmegaMin 0.01
--XiQuadrature [adaptive | cliff]
--XiMin 0.01
--AbsTol
--RelTol
  ````
{.toc}


> These arguments affect the behavior of application codes
> that compute output quantities by performing numerical
> integrations over angular frequencies---either real
> angular frequencies $\omega$ or imaginary angular
> frequencies $\xi$. (More specifically, [[scuff-neq]]
> performs $\omega$ integrations, while [[scuff-cas3d]] and
> [[scuff-caspol]] perform $\xi$ integrations.)
>
> `--OmegaQuadrature` or `--XiQuadrature` specify the
> numerical quadrature algorithm. If these
> options are left unspecified, an appropriate algorithm 
> is chosen automatically.
>
> `--OmegaMin` or `--XiMin` specify the minimum angular
> frequency at which numerical [[scuff-em]] calculations 
> are performed. For cases in which the lower limit of the 
> $\omega$ or $\xi$ integration is 0, the integrand is 
> assumed to be constant between 0 and the value specified here.
> These options are interpreted in the usual [[scuff-em]]
> frequency units of $3\times 10^{14}$ rad/sec, so typical
> values will be something like `0.001.`
>
> `--AbsTol` and `--RelTol` may be used to specify
> absolute and relative error tolerances for adaptive
> quadrature algorithms. If an adaptive cubature method
> seems be spending too much time attempting to achieve
> high accuracy in a frequency quadrature, try increasing
> `--RelTol` to something like `0.1` or even `0.5.`

### *Options controlling Brillouin-zone integrations*

 ````
--BZQuadrature [adaptive | cliff]
--BZSymmetry [adaptive | cliff]
--MaxBZSamples 1000;
 ````
{.toc}

> These arguments affect the behavior of application codes
> that compute output quantities for periodic geometries
> by performing numerical integrations over the Brillouin
> zone. Such codes include [[scuff-ldos]], [[scuff-neq]],
> [[scuff-cas3d]], and [[scuff-caspol]].
>
> `--BZSymmetry` may be used for 2D periodic geometries to
> declare that the Brillouin-zone integrand $f(k_x, k_y)$
> is symmetric under the interchange $k_x \leftrightarrow k_y$.
>
> `--MaxBZSamples 1000` may be used to restrict adaptive
> integration algorithms to a maximum of 1,000 evaluations
> of the Brillouin-zone integrand.
 
### *Miscellaneous options*

`--FileBase MyFileBase`
{.toc}

> Specifies the base file name for output files (so that, for example,
> the frequency-resolved output file written by [[scuff-cas3d]]
> will be `MyFileBase.byXi`, while the frequency-integrated 
> output file will be `MyFileBase.out`). If this option is not 
> specified, the file base is taken to be the base filename of the 
> `.scuffgeo` file.

# 2. Passing command-line options via text file

All of the standalone applications in the [[scuff-em]] suite allow 
their command-line options to be passed via a text file fed into 
standard input.

Each line of this text file should consist of a single 
command-line option (minus the -- at the beginning) followed by any 
arguments the option might take.

For example, running [[scuff-scatter]] with the command-line options

````bash
% scuff-scatter --geometry Spheres.scuffgeo --omega 1.0 --pwPolarization 1 0 0 --pwDirection 0 0 1 --EPFile MyEPFile 
````
    
is equivalent to running

````bash
% scuff-scatter < MyOptionsFile
````
    
where the file ``MyOptionsFile`` looks like this:

````
# options for scuff-scatter 
geometry Spheres.scuffgeo
omega 1.0

pwPolarization 1 0 0 
pwDirection 0 0 1

EPFile MyEPFile
````
    
Note that blank lines and comments (lines starting with #) are ignored.

You may also combine the two methods of specifying options by passing 
some options via text file and others on the command line. If there are 
any conflicts, the values specified on the command line take precedence. 
For instance, to re-run the example above at a new frequency with 
everything else unchanged, you could say

````
 % scuff-scatter --Omega 2.0 < MyOptionsFile
````

<a name="Complex"></a>
# 3. Complex numbers

Many of the standalone programs in the [[scuff-em]] suite have 
options for which you may specify complex numbers. (An example 
is the --omega option accepted by [[scuff-scatter]] and other 
codes, for which you may specify complex or even pure imaginary 
numbers to do calculations at complex frequencies.)

To specify a complex number as a parameter value, write both the 
real and imaginary parts together as a single string (no spaces), 
separated by ``+`` or ``-``, with the imaginary part terminated by 
``i`` or ``I`` (you may also use ``j`` or ``J``). For example, all 
of the following are valid frequency specifications:

````
 --omega 2.3+4.5i
 --omega 2.3
 --omega 4.5j
 --omega 12.3e2+45.4e2I
````

<a name="LogFiles"></a> 
# 4. Log files 

All command-line codes in the [[scuff-em]] suite
write logging information to text-based logfiles
with extension `.log.` You can monitor these
files to follow the progress of your calculations.

For example, after launching 
[this sample <span class="SC">scuff-cas3d</span> run][SiliconSlabs],
type the following at a terminal window:

````bash
% tail -f SiliconSlabs_L2_40.log
````
This produces a running list of log messages, something like this:
````
06/15/15::10:13:40: scuff-cas3D running on superhr1
06/15/15::10:13:40: Added /home/homer/work/scuff-em-sandbox/mshFiles to mesh search path.
06/15/15::10:13:40: Adding lattice basis vector (2,0).
06/15/15::10:13:40: Adding lattice basis vector (0,2).
06/15/15::10:13:40: Computing Casimir integrand at (Xi,kx,ky)=(0.5,0.19635,0.19635)
06/15/15::10:13:49: Computing Casimir integrand at (Xi,kx,ky)=(0.5,0.217586,0.172518)
06/15/15::10:13:51: Computing Casimir integrand at (Xi,kx,ky)=(0.5,0.256543,0.106264)
06/15/15::10:13:52: Computing Casimir integrand at (Xi,kx,ky)=(0.5,0.275845,0.0318681) 
````

<a name="OutputFiles"></a>
# 5. Output Files

The [[scuff-em]] application codes generally produce
output in the form of human-readable text-based
data files. The typical naming convention for 
these files is `Geometry.Extension`, where `Geometry.scuffgeo`
is the name of the [[scuff-em]] geometry file 
on which the calculation was run, and where `Extension`
is an application-specific file extension that attempts 
to describe the content of the file; for example, 
[[scuff-scatter]] produces files named `Geometry.PFT`
to report the power, force, and torque on bodies 
irradiated by external fields. You can use the 
`--filebase` command-line option to [[scuff-em]]
codes to select a base filename other than 
`Geometry.`

All text-based data files produced by [[scuff-em]]
codes should contain a human-readable *header* at the 
top of the file explaining how to interpret its content.
For example, the first few lines of the 
`.PFT` output file produced by [[scuff-scatter]]
look like this:

````
# scuff-scatter running on hikari 08/03/15::23:15:45
# data file columns: 
# 1 omega 
# 2 surface label 
# 3 absorbed power (watts)
# 4 scattered power (watts)
# 5 x-force (nanoNewtons)
# 6 y-force (nanoNewtons)
# 7 z-force (nanoNewtons)
# 8 x-torque (nanoNewtons*microns)
# 9 y-torque (nanoNewtons*microns)
# 10 z-torque (nanoNewtons*microns)
0.05 Particle 6.994160e-04 1.562608e-03 -1.120907e-04 7.169210e-03 7.128598e-03 -3.622586e-03 1.005173e-03 9.428114e-02 
````

If you encounter a situation in which a [[scuff-em]]
application code fails to write an appropriate file
header to an output file, please file an issue
on the 
[<span class="CodeName">scuff-em</span> GitHub page][GitHub].


Other types of output files produced by [[scuff-em]]
application codes include `.hdf5` binary data
files and `.pp` files containing visualization data
that may be viewed in [[gmsh]].

# 6. Environment variables

Here are some environment-variable settings that 
affect the behavior of [[scuff-em]].

````bash
% export SCUFF_MESH_PATH=/path/to/msh/files
````

> Specifies a directory in which to look for mesh files
> (such as `.msh` files produced by [[gmsh]]) 
> referred to by `.scuffgeo` files.

````bash
% export SCUFF_LOGLEVEL="NONE"
% export SCUFF_LOGLEVEL="TERSE"
% export SCUFF_LOGLEVEL="VERBOSE"
% export SCUFF_LOGLEVEL="VERBOSE2"
````

> Sets the verbosity of messages written to the `.log` file.

````bash
% export SCUFF_INTERPOLATION_TOLERANCE=1.0e-3
````

> This option, which is only relevant for
> [extended geometries][ExtendedGeometries],
> set an internal tolerance parameter that controls
> the accuracy with which BEM matrix elements and
> scattered fields are computed. Its default value
> is `1.0e-6`, but this is probably overly stringent,
> and may generally be relaxed to `1e-3` or so
> to reduce memory usage and CPU time without 
> significant accuracy penalties. If your calculation
> is running out of memory or taking too long to run,
> try setting it to `1.0e-4` or `1.0e-3.`
> Please tell us about your experiences 
> with this parameter!
>
> (More specifically: [[scuff-em]] uses Ewald summation
> to accelerate the calculation of the periodic Green's
> function, but even this accelerated calculation is 
> not fast enough to handle the many millions of 
> evaluations needed to assemble the full BEM matrix.
> For this reason, when assembling the BEM matrix at
> a given frequency and Bloch vector, [[scuff-em]] first
> precomputes Ewald-summed values of the periodic DGF 
> at grid points of an interpolation grid, after which
> values are obtained by interpolation (bypassing 
> Ewald summation). The spacing of the grid points is 
> chosen automatically to ensure that the maximum relative 
> error between the interpolated and exact values at any 
> point within the grid boundaries is less than 
> `SCUFF_INTERPOLATION_TOLERANCE.`


````bash
% export OMP_NUM_THREADS="8"
% export GOMP_CPU_AFFINITY="0-7"
````

> These options **should not** be necessary, but **may** be
> needed to ensure that [[scuff-em]] takes advantage
> of all available CPU cores on your system. The former 
> option says that you want to use 8 cores, and the latter 
> option says that you want these 8 cores to be the first
> 8 available (as opposed to, say, the second set of 8
> available cores on a 16-core machine).

[scuff-cas3D]:  scuff-cas3D/scuff-cas3D.md
[scuff-caspol]: scuff-caspol/scuff-caspol.md
[SiliconBeams]: ../../examples/SiliconBeams/SiliconBeams.md
[SiliconSlabs]: ../../examples/SiliconSlabs/SiliconSlabs.md
[ExtendedGeometries]: ../reference/Geometries.md#Extended
[GitHub]:       https://github.com/HomerReid/scuff-em.git
