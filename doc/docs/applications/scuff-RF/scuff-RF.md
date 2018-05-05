# <span class=SC>scuff-rf</span> The <span class=SC>scuff-em</span> RF module

The <span class=SC>scuff-em</span> RF module is an extension of the
[<span class=SC>scuff-em</span> core library][libscuff] to enable
RF device modeling within the framework of the surface-integral-equation (SIE)
approach to electromagnetic scattering implemented by <span class=SC>scuff-em</span>.
More specifically, the RF module extends <span class=SC>scuff-em</span>'s
core SIE solver in two key ways:

+ As a new type of incident-field source, <span class=SC>scuff-rf</span>
introduces the notion of *RF port.* An RF port is simply
a localized region of a material body into which an RF current
may be injected. The fields radiated by currents forced into RF ports
are taken to be the incident fields in an SIE scattering problem, which 
<span class=SC>scuff-em</span> solves like any other scattering problem
 to determine the induced currents everywhere *else* (outside the ports) 
on the surfaces in a geometry. 

Having computed the induced currents in response to an RF port excitation,
we can now make use of the full panoply of existing post-processing options
(calculation of scattered and total field components, induced moments, 
absorbed and scattered power, visualization of fields and currents, etc)
provided by the  [<span class=SC>scuff-em</span> core library](libscuff).
However, 

+ As a specific new type of post-processing calculation, <span class=SC>scuff-rf</span>
implements an algorithm for computing the matrix of *impedance parameters*
for a multiport RF structure.

Technical details regarding the implementation of both of these features
are discussed in the [<span class=SC>scuff-rf</span> technical memo][scuffRFMemo].

## Interfaces to the <span class=SC>scuff-em</span> RF module: C++, Python, command-line

The features described above are implemented by a library named `libRFSolver`
packaged within the <span class=SC>scuff-em</span> source distribution. There
are multiple ways to access the functionality of `libRFSolver`

+ You can write C++ or python codes that make calls to 'libRFSolver' API routines.

+ You can use the `scuff-rf` command-line module (built and installed by the standard
<span class=SC>scuff-em</span> installation process).

## Interfaces to the <span class=SC>scuff-em</span> RF module: C++, Python, command-line

The documentation for [[scuff-RF]] has not yet been
ported from its earlier version. For the time being, please
[access the earlier version of the documentation.][EarlierVersion]

[EarlierVersion]: http://homerreid.com/scuff-em/scuff-RF

[GMSH]:                 http://www.geuz.org/gmsh
[scuffGeometries]:      ../../reference/Geometries
[scuff-rf]              ../../applications/scuff-RF/scuff-RF.md
[libscuff]              ../../API/libscuff.md
[libSubstrate]          ../../tex/FullWaveSubstrate.pdf
[PCBTraceMemo]:		http://www.analog.com/media/en/training-seminars/tutorials/MT-094.pdf
[scuffRFMemo]:		../../tex/scuff-rf.pdf
