# Computing T-Matrices of arbitrary objects with scuff-tmatrix

The well-known [*T-matrix method*](http://en.wikipedia.org/wiki/T-matrix_method) is a widely used
technique for solving problems involving electromagnetic scattering from
compact objects. In this method, the scattering properties of a compact
body are encapsulated in its *T-matrix*, whose entries give the
amplitudes of the *outgoing* spherical waves that arise from irradiating
the object with a single *regular* spherical wave.

<span class="SC">scuff-tmatrix</span> is a specialized tool within
the <span class="SC">scuff-em</span> code suite for computing the
*T*-matrix of arbitrarily-shaped objects with arbitrary
frequency-dependent material properties.

To compute the *T-*matrix of an object (or a collection of objects)
using <span class="SC">scuff-tmatrix,</span> you will

1.  create a [geometry file](../../reference/scuffEMGeometries.shtml)
    describing the shapes and material properties of the scattering
    objects in your geometry;

2.  run <span class="SC">scuff-tmatrix</span> with command-line
    options specifying the geometry, the maximum $\ell$ value of the
    spherical wave to consider (which determines the dimension of the 
    computed *T-*matrix), and the frequency range over which to run computations.

You will get back

1.  a text file listing all *T-*matrix elements (with $\ell$-values up to
    the maximum $\ell-$value you specified) at all frequencies you
    requested, and

2.  optionally, a binary `.HDF5` file containing the *T-*matrix data
    together with simple scripts for importing the data into 
    <span class="SC">julia.</span>

+--------------------------------------------------------------------------+
| Table Of Contents                                                        |
+==========================================================================+
| [1. <span class="SC">scuff-tmatrix</span> Command-Line             |
| Options](scuff-tmatrix.md#Options)                                       |
+--------------------------------------------------------------------------+
| [2. <span class="SC">scuff-tmatrix</span> Output                   |
| Files](scuff-EM/scuff-tmatrix/index.shtml#OutputFiles)                   |
+--------------------------------------------------------------------------+
| [3. <span class="SC">scuff-tmatrix</span>                          |
| Examples](scuff-EM/scuff-tmatrix/index.shtml#Examples)                   |
| >   -------------------------------------------------------------------- |
| >   [4a. Dielectric Sphere](scuff-EM/scuff-tmatrix/index.shtml#Sphere)   |
| >   [4b. Dielectic Cube](scuff-EM/scuff-tmatrix/index.shtml#Cube)        |
| >   -------------------------------------------------------------------- |
+--------------------------------------------------------------------------+
| [4. Reading T-Matrix data into <span                                     |
| class="SC">julia</span>](scuff-EM/scuff-tmatrix/index.shtml#Julia) |
+--------------------------------------------------------------------------+

# 1. <span class="SC">scuff-tmatrix</span> Command-Line Options
------------------------------------------------------------------------

------------------------------------------------------------------------

The following table summarizes all command-line options currently
available in <span class="SC">scuff-tmatrix.</span>

As is true for all programs in the <span
class="SC">scuff-em</span> suite, command-line options may be
specified in a text file catted to standard input; see
[here](scuff-EM/reference/scuffEMMisc.shtml#Options) for an example of
how this works.

\

Option

Description

*Options controlling the scattering geometry*

    --geometry MyGeometry.scuffgeo

Specifies the [geometry
file](scuff-EM/reference/scuffEMGeometries.shtml) describing the
scattering geometry. This option is always mandatory.

*Options specifying the frequencies considered*

     --Omega xx 

Specifies an angular frequency at which to run computations, in units of
3•10^14^ rad/s (=c / 1 μm).

The value specified for `--Omega` may be a [complex
number.](scuff-em/reference/scuffEMMisc.shtml#Complex)

You may request computations at more than one frequency by using the
`--Omega` option more than once, i.e. you may say

``` {.listing}
            --Omega 1e-2 --Omega 1e-1 --Omega 1 --Omega 10
           
```

(However, for more than a few frequencies it is more convenient to use
the `--OmegaFile` option discussed below.)

     --OmegaFile xx 

Specifies the name of a file containing a list of frequencies at which
to do computations.

The file should contain one frequency per line; blank lines and comments
(lines beginning with `#`) are ignored.

*Options describing the range of *l*-values considered*

+--------------------------------------------------------------------------+
|      --lMax xx                                                           |
+--------------------------------------------------------------------------+

Request calculation of *T-*matrix entries with spherical wave indices up
to *l=*`xx`.

The dimension of the *T-*matrix computed at each frequency is `DxD`,
where `D=(lMax+1)(lMax+1)-1`.

For example, if you specify `--lMax 3`, then the *T-*matrices computed
will be 15x15 matrices, with entries indexed by spherical-wave indices
*(l,m)*=(1,-1), (1,0), (1,0), (2,-2), ..., (3,3).

If you do not specify this option, the default value is `lMax=5`.

*Options controlling output files*

+--------------------------------------------------------------------------+
|      --OutputFile MyOutputFile                                           |
+--------------------------------------------------------------------------+

Requests that the text-based output data be written to the file
`MyPowerFile.dat.`

If this file is not specified, the text-based output data are written by
default to a file named `Geometry.TMatrix`, where `Geometry.scuffgeo`
was the name of the <span class="SC">scuff-em</span> geometry file
you specified with the `--geometry` command-line option.

+--------------------------------------------------------------------------+
|      --HDF5File MyFile.hdf5                                              |
+--------------------------------------------------------------------------+

Requests that binary *T*-matrix data be written to the file
`MyFile.hdf5.` (Specifying this option will automatically enable the
generation of a <span class="codename">matlab</span> import script named
`MyFile.m` which you can execute at the <span
class="codename">matlab</span> command line to import *T-*matrix data;
see the examples below.)

Note that, in contrast to text-based output, binary data output is
disabled by default; you must specify this option to enable binary data
output.

*Cache options*

+--------------------------------------------------------------------------+
|      --Cache MyCache.scuffcache                                          |
+--------------------------------------------------------------------------+
|                                                                          |
+--------------------------------------------------------------------------+
|      --ReadCache InCache1.scuffcache                                     |
+--------------------------------------------------------------------------+
|      --ReadCache InCache2.scuffcache                                     |
+--------------------------------------------------------------------------+
|      --WriteCache OutCache.scuffcache                                    |
+--------------------------------------------------------------------------+

Specifies the names of [cache files for geometric
data.](scuff-em/reference/scuffEMMisc.shtml#Caching)

The `--ReadCache` option allows you to specify a file from which
geometric data will be preloaded before the calculation begins. This
option may be specified any number of times. If the specified file does
not exist, the option is silently ignored.

The `--WriteCache ` option allows you to specify a file to which
geometric data will be written after the calculation has completed. The
resulting cache file will contain any data that were preloaded from
`--ReadCache` files, plus any data that were newly generated during the
course of the <span class="SC">scuff-tmatrix</span> run.

The `--Cache XX` option is equivalent to saying
`--ReadCache XX --WriteCache XX.`

For more information on geometric data caching in <span
class="SC">scuff-em</span>, see
[here.](scuff-em/reference/scuffEMMisc.shtml#Caching)

\

[]() 2. <span class="SC">scuff-tmatrix</span> Output Files
----------------------------------------------------------------

------------------------------------------------------------------------

### 1. *T-*matrix data in text form

<span class="SC">scuff-tmatrix</span> always writes *T-*matrix
data to a text-based output file. (By default this file is named
`MyGeometry.TMatrix`, where `MyGeometry.scuffgeo` was the name of the
geometry file you specified with the `--geometry` option; alternatively,
you can specify your own output file name using the `--OutputFile`
option.)

Each line of the text-based output file contains a single *T-*matrix
element at a single frequency. The format of each line is this:

  -----------------------------------------------------
  *Omega P L M PPrime LPrime MPrime real(T) imag (T)*
  -----------------------------------------------------

The triple *(P,L,M)* labels the spherical wave that constitutes the
**row index** for the *T-*matrix entry. (Here *P* identifies the
polarization and is either `M` or `E` for *M-* or *E-*type waves,
respectively; *L* is an integer between 1 and `lMax` inclusive; and *M*
is an integer between *-L* and *L* inclusive.

Similarly, the triple *(PPrime,LPrime,MPrime)* labels the spherical wave
that constitutes the **column index** into the *T-*matrix.

The final two entries on the line are the real and imaginary parts of
the *T-*matrix entry for the given pair of spherical waves at the given
frequency.

Physically, the *T-*matrix element with row index *(P,L,M)* and column
index *(PPrime, LPrime, MPrime)* is the amplitude of the outgoing
spherical wave with indices *(P,L,M)* that arises from illuminating your
object with a regular spherical wave with indices
*(PPrime,LPrime,MPrime.)*

See below for a simple shell script you can use to extract a particular
*T-*matrix entry vs. frequency from the `.TMatrix` file.

Also see below for a [julia](http://julialang.org){.SC} code that
you can use to import T-matrix data into a
[julia](http://julialang.org){.SC} session.

### 2. *T-*matrix data in binary form

If you specify the `--HDF5File` command-line argument, then *T-*matrix
data will be written in binary `HDF5` format to the specified file. The
*T-*matrix for each frequency will be stored as a separate dataset with
a title like `T_1.23456` where the number indicates the frequency.

Also, if you say `--HDFFile MyFile.hdf5`, then, in addition to the
`MyFile.hdf5` file, you will get a file `MyFile.m`, which is a simple
script for importing the *T-*matrix data in `MyFile.hdf5` into a <span
class="SC">matlab</span> session. See below for examples of how
this works.

\

[]() 3. <span class="SC">scuff-tmatrix</span> Examples
------------------------------------------------------------

------------------------------------------------------------------------

Here are some examples of calculations you can do with <span
class="SC">scuff-tmatrix</span>. Input files and command-line
runscripts for all these examples are included in the
`share/scuff-em/examples` subdirectory of the <span
class="SC">scuff-em</span> installation.

------------------------------------------------------------------------

### []() 4a. Dielectric sphere

We start with the canonical textbook stalwart of scattering from a
dielectric sphere. This is an example (indeed, the only example) of a
dielectric object whose *T-*matrix may be computed analytically, making
it a useful benchmark for our numerical computation.

The first step is to create a geometry mesh and a `.scuffgeo` file for
the sphere. This process is described in detail
[here.](scuff-EM/scuff-scatter/index.shtml#Mie) For simplicity in this
case we will take the sphere to have a constant relative dielectric
permittivity of ε~r~=10 at all frequencies. We will start with the
sphere mesh containing 327 internal edges; the file
`E10Sphere_327.scuffgeo` looks like this:

``` {.listing}
OBJECT E10Sphere
        MESHFILE Sphere.msh
        MATERIAL CONST_EPS_10
ENDOBJECT
 
```

We create a simple file called
[`OmegaValues.dat`](scuff-em/scuff-scatter/OmegaValues.dat) containing a
list of angular frequencies (in units of 3•10^14^ rad/s) at which to
compute *T-*matrices:

``` {.listing}
    0.010
    0.013
    ...
    5.0
    
```

And now we launch <span class="codename">scuff-tmatrix: </span>

``` {.listing}
    % scuff-tmatrix --geometry E10Sphere_654.scuffgeo --omegafile OmegaValues.dat --cache Sphere.cache 
    
```

This produces the file `E10Sphere_654.TMatrix`, which contains one
*T-*matrix entry per line:

``` {.listing}
    
```

To assess the impact of meshing fineness, let's re-run the example with
a finer mesh. We will use the sphere mesh with 1362 interior edges:

``` {.listing}
   % scuff-tmatrix --geometry E10Sphere_1362.scuffgeo --omegafile OmegaValues.dat --cache Sphere.cache 
    
```

This produces the file `E10Sphere_1362.TMatrix`, with file format
similar to the above.

To look at individual *T-*matrix entries as a function of frequency, it
is convenient to extract those entries from the `.TMatrix` file into
their own little file. [Here](scuff-em/scuff-tmatrix/Filter) is a simple
shell script, named simply `Filter`, which takes either 4 or 7
command-line arguments, as follows:

``` {.listing}
 
     Filter MyFile.TMatrix P L M 

     Filter MyFile.TMatrix P L M PPrime LPrime MPrime
    
```

The first invocation will filter out only the diagonal *T-*matrix
entries with row and column index equal to *(P,L,M).* The second will
filter out off-diagonal *T-*matrix entries with given row and column
indices. In both cases the filtered data are written to a file named
something like `MyFile.M11` or `MyFile.M11.M22` with three numbers per
line:

  --------------------------
  *Omega real(T) imag (T)*
  --------------------------

Let's extract a few diagonal matrix entries:

``` {.listing}
 
     Filter Sphere_327.TMatrix E 1 1 # get the (M11, M11) diagonal entry; produces file Sphere327.M11
     Filter Sphere_327.TMatrix M 2 0 
    
```

and a couple of off-diagonal matrix entries:

``` {.listing}
 
     Filter Sphere_327.TMatrix M 1 1 E 1 1 # produces file Sphere327.M11.E11
     Filter Sphere_327.TMatrix M 2 1 M 1 1
    
```

![](scuff-em/scuff-tmatrix/TMatrix.png)

\

[]() 4. Reading T-Matrix data into <span class="SC">julia</span>
----------------------------------------------------------------------

------------------------------------------------------------------------

Here's a [julia](http://julialang.org){.SC} function that you can
use to read T-matrix data from the `.TMatrix` output file produced by
<span class="SC">scuff-tmatrix</span>:
[`ReadTMatrix.jl`](scuff-em/scuff-tmatrix/ReadTMatrix.jl).

For example, after following the example outlined above to produce a
data file named `E10Sphere_654.TMatrix`, you could go like this to
import the T-matrix data into <span class="SC">julia</span>:


``` {.listing}
 
     % julia
     julia> reload("ReadTMatrix.jl");
     julia> (OmegaList, TList) = ReadTMatrix("E10Sphere_654.TMatrix");
      11 frequencies, LMax=2, TMatrix dimension=256
    
```

Now the vector `OmegaList` contains the unique values of the angular
frequency at which you requested T-matrix data, and `TList[n,:,:]` is
the T-matrix corresponding to angular frequency `OmegaList[n].`
