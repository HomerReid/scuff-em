# General reference for [[scuff-em]] command-line applications

This page collects some general information that applies
to many or all of the standalone command-line applications
in the [[scuff-em]] suite.

<a name="CommonOptions"></a>
# Common command-line arguments

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
> [extended geometries](../reference/Geometries.md#Extended).
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
--OmegaQuadrature [adaptive | cliff ]
--OmegaMin 0.01
--XiQuadrature [adaptive | cliff ]
--XiMin 0.01
--AbsTol
--RelTol
  ````
{.toc}

### *Options controlling Brillouin-zone integrations*

 ````
--BZQuadrature [adaptive | cliff ]
--BZSymmetry [adaptive | cliff ]
--MaxBZSamples 1000;
 ````
{.toc}
 
### *Miscellaneous options*

`--FileBase MyFileBase`
{.toc}

> Specifies the base file name for output files (so that, for example,
> the frequency-resolved output file written by [[scuff-cas3d]]
> will be `MyFileBase.byXi`, while the frequency-integrated 
> output file will be `MyFileBase.out`). If this option is not 
> specified, the file base is taken to be the base filename of the 
> `.scuffgeo` file.

# Passing command-line options via text file

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
# Complex numbers

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
 
# Log files 

All command-line codes in the [[scuff-em]] suite
write logging information to text-based logfiles
with extension `.log.` You can monitor these
files to follow the progress of your calculations.

For example, after launching 
[ this sample <span class="SC">scuff-cas3d</span> run ][SiliconSlabs],

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

# Environment variables

Here are some environment-variable settings that 
affect the behavior of [[scuff-em]].

````bash
% export SCUFF_MESH_PATH=/path/to/msh/files`
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
