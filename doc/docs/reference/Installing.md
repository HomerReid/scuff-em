## Installing [[scuff-em]]

## 0. Try <span class="SC">scuff-em</span> online before installing anything!

[Johannes Feist](http://www.johannesfeist.eu/) has set up
an interactive online notebook using [MyBinder](http://mybinder.org/)
that you can use to run some simple [[scuff-em]] calculations on
the cloud without having to install anything. This free service
offers limited computing power, so you won't want to attempt any
heavy calculations here, but you can follow Johannes' tutorial
to observe the flow of a typical [[scuff-em]] calculation,
play around with your own calculations, and get a feel for what 
using [[scuff-em]] will like once you have installed it on your 
machine.

Here's the link:
[http://mybinder.org/repo/jfeist/scuff-em-mybinder](http://mybinder.org/repo/jfeist/scuff-em-mybinder).

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
The current build status of the [[scuff-em]] master branch is:

![Build Status](https://travis-ci.org/HomerReid/scuff-em.svg?branch=master)

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
[<span class="SC">scuff-em</span> GitHub page][GitHub].

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
[`.log` files][LogFiles].

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

<a name="Disabling Python"></a>
#### Disabling the python interface to speed the build process

Compiling the python interface is slow---it accounts for
more than half of the build time on some systems.
If you don't need the python interface to [[scuff-em]],
use the option `--without-python` when running `configure`
to accelerate the build process.

<a name="Debugging"></a>
#### Building for debugging

If you would like to run [[scuff-em]] API codes in a debugger
like [<span class="SC">gdb</sc>](https://www.gnu.org/software/gdb),
you will want to modify the build options to **(a)** include
debugging symbols, **(b)** turn off optimization, **(c)**
disable [[openmp]] multithreading, which does not play well
with [[GDB]]. 

Here is the script that works for me to achieve these goals:

````bash
#!/bin/bash

CC="gcc -ggdb -O0"
CXX="g++ -ggdb -O0"
export CFLAGS="-O0"
export CXXFLAGS="-O0"
sh autogen.sh --enable-debug --without-openmp --disable-shared
````

(It shouldn't be necessary to have to add `-O0` to both the
environment variables and the compiler command lines, but
this seems to be the only way things work for me.)

After running this script to reconfigure for building with
compiling support, you will want to `make clean; make`
to rebuild everything with debugging support. Of course,
after you are finished debugging you will need to reconfigure
with debugging support removed and then re-do the 
clean build. Because this is time-consuming, I typically
maintain two separate copies of the code base, one for
debugging and one for running at full speed.

Once debugging support has been added, you can run 
the code in [[gdb]]. Here's a sample session:

````bash
 % gdb /path/to/debug/repository/bin/scuff-ldos

 (gdb) set args < scuff-ldos.args

 (gdb) break GetHalfSpaceDGFs
Breakpoint 1 at 0x409626: file AnalyticalDGFs.cc, line 217.

 (gdb) run

Breakpoint 1, GetHalfSpaceDGFs (Omega=..., kBloch=0x7fffffffd5c0, 
    zp=0.10000000000000001, LBasis=0x7fffffffd380, MP=0x130a7c0, RelTol=0.01, 
    AbsTol=1e-10, MaxCells=1000, GE=0x7fffffffd3a0, GM=0x7fffffffd430)
    at AnalyticalDGFs.cc:217

217	  double BZVolume=GetRLBasis(LBasis, RLBasis);

(gdb) print Omega
$1 = {_M_value = 0.10000000000000001 + 0 * I}

(gdb) u 230

GetHalfSpaceDGFs (Omega=..., kBloch=0x7fffffffd5c0, zp=0.10000000000000001, 
    LBasis=0x7fffffffd380, MP=0x130a7c0, RelTol=0.01, AbsTol=1e-10, 
    MaxCells=1000, GE=0x7fffffffd3a0, GM=0x7fffffffd430)
    at AnalyticalDGFs.cc:230

230	  cdouble Sum[18];

(gdb) print Data->Epsilon
$2 = {_M_value = -45765.335680436379 + 115960.58424811834 * I}
````

**Note**: If you have a better debugging solution that 
does not require steps **(b)** and/or **(c)** above , please
let me know about it. It stinks to have to run the codes
at greatly reduced speed when debugging, because often the
problem spots lie after expensive sections like BEM matrix
assembly, and then it takes forever for the code to run
in the debugger to get there.

[GitHub]:                      https://github.com/HomerReid/scuff-em/
[LogFiles]:                    ../applications/GeneralReference.md#LogFiles
