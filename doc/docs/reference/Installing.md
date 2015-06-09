## Installing [[scuff-em]]

## 1. External packages

[[scuff-em]] relies on a small number of well-established free 
software packages to implement certain non-essential functionality. 
[[scuff-em]] can be compiled and installed without any of these packages, 
but in this case the code will be somewhat crippled.

+ If you actually want to solve scattering problems (instead of 
  just setting them up), you will need
  [LAPACK/BLAS](http://www.netlib.org/lapack).

+ If you want the capacity to write output files in the 
  standard HDF5 binary format, you will need
  [HDF5](http://www.hdfgroup.org/HDF5).

+ If you want to compile the [[python]] interface, you will need the
  [[python]] development files.

+ Although not required to install, compile, or use
  [[scuff-em]],
  [<span class="SC">gmsh</span>](http://geuz.org/gmsh)
  is an extremely valuable open-source meshing and visualization
  tool that is used throughout the
  [[scuff-em]] documentation.

On Debian/Ubuntu Linux systems, you can fetch all of these packages by doing a 

````bash
% sudo apt-get install liblapack-dev libblas-dev libhdf5-serial-dev python-dev gmsh
````

> Note: In some cases it seems the ``gmsh`` package conflicts with 
> the ``libhdf5-serial-dev`` package. In this case, just 
> remove ``gmsh`` from the above ``apt-get`` statement; you can 
> install it by hand following the instructions
> on the [GMSH website](http://geuz.org/gmsh).
> (Note that [[gmsh]], though very useful, 
> is not necessary to compile or run [[scuff-em]].)

## 2. Cloning the GitHub repository and building the code

[[scuff-em]] is hosted on [GitHub][GitHub].
To fetch and install the latest version of the 
code, execute the following steps. (Replace the string
``/path/to/scuff-em-installation-directory``
with your desired installation directory.)

````bash
% git clone https://homerreid@github.com/HomerReid/scuff-em.git
% cd scuff-em
% sh autogen.sh --prefix=/path/to/scuff-em-installation-directory
% make install
````

If this succeeds, the executable versions of the application
programs (such as ``scuff-scatter``, ``scuff-rf``, etc.) will be 
installed in the directory ``PREFIX/bin/`` 
and the demonstration examples for the various application programs 
will be in ``PREFIX/share/scuff-em/examples``
(where ``PREFIX`` is the directory you specified using the 
``--prefix`` option above).

If you have trouble installing [[scuff-em]],
please file an issue on the 
[<span class="SC">scuff-em</span> GitHub page](GitHub).

### Build options

You may specify options to the ``autogen.sh``
(or ``configure``) command to guide the compilation process. 
For a full list of available options,
type ``configure --help.`` Here we summarize some of the
more salient possibilities.

<a name="Threading"></a>
#### Multithreading: [[openmp]] vs. [[pthreads]] 

The [[scuff-em]] core library 
uses multithreading for all steps in the
<a href="scuff-em/libscuff/MainFlow.shtml">main flow</a> of
the BEM scattering procedure. You can use ``configure``
options to select whether this multithreading is implemented
using [[openmp]] or [[pthreads]]. The former is the default,
while the latter may be enabled like this:

````bash
% ./configure --without-openmp --with-pthreads
````

The default [[openmp]]
multithreading tends to play better with other multithreaded 
software packages and is definitely the right choice if you 
will be operating on a shared machine with some CPU cores 
occupied by other users. 
(You can monitor performance by inspecting
<a href="scuff-em/reference/scuffEMMisc.shtml#LogFiles"``.log`` files.</a>)

Support for [[pthreads]]
is a legacy feature that will be discontinued in future versions of 
[[scuff-em]].

Note: in some cases you may need to tweak certain environment 
variables to achieve maximal 
[[openmp]] performance.
For example, on my workstation (which has 8 CPU cores),
in order to get [[openmp]] codes
to use all 8 cores I need to set the following environment
variable:

````bash
% export GOMP_CPU_AFFINITY=0-7
````
