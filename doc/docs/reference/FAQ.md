# Frequently asked questions about [[scuff-em]]

<a name="Units">
## What units does [[scuff-em]] use for physical quantities like length, frequency, field strengths, power, force, torque, etc?

### Length and frequency 

Short answer: The default units are $L_0=1\, \mu\text{m}$ for length
and $\omega_0=3\cdot 10^{14}$ rad/sec (=$c/L_0)$ for angular frequency.

Thus, if you write a 
[<span class="SC">gmsh</sc>][gmsh]
geometry (`.geo`) file describing a sphere
of radius `1.3` and use the resulting surface mesh in a
[<span class="SC">scuff-scatter</sc>][scuff-scatter] calculation
with an angular-frequency specification of `--omega 2.0`,
then you will be studying a sphere of radius 1.3 microns
at an angular frequency of $6\cdot 10^{14}$ rad/sec 
(corresponding to a free-space wavelength of 
$\lambda=\frac{2\pi c}{\omega} = 3.1415\,\mu$m,
which could alternatively be specified by saying
`--lambda 3.1415` instead of `--omega 2.0`).
[The dimensionless quantity $a=\frac{\omega R}{c}$
 (the "size parameter" in Mie theory) 
  is just the product of the numerical values specified
  for the radius and `--omega`, i.e. $a=1.3 \cdot 2.0=2.6$
  in this case.]

Longer answer: 
For a problem involving only bodies with
[frequency-independent material properties][Materials]
(permeability $\epsilon$ and permittivity $\mu$),
including perfectly-conducting (`PEC`) bodies and
frequency-independent dielectrics
(such as `CONST_EPS_10+1i`), the scale invariance
of Maxwell's equations means that the same
computational results can be interpreted on different
length scales. For example, if your mesh file
describes a sphere of radius `0.9' and you run
a [[scuff-em]] calculation at angular frequency `1.2`, 
then the results 

+ 

However, this scale invariance is broken by
frequency-dependent dielectric functions, such
as the `SILICON` material definition in 
[this example][SiliconSlabs].

### Field strengths

The only code in the [[scuff-em]] suite 

### Power, force, torque


## I'm working on a big project involving multiple calculations on a geometry with various different values of geometric parameters, material properties, and meshing fineness. It's getting unwieldy to have all of these `.geo` and `.msh` and `.scuffgeo` files cluttering up my project directory. How would you suggest organizing things?

Here is what I typically do:

+ Within your top-level project directory (I'll call it
      `~/myProject/`), create the following subdirectories:
    + `~/myProject/geoFiles` (for [<span class="SC">gmsh</span>][gmsh] geometry files)
    + `~/myProject/mshFiles` (for surface mesh files)
    + `~/myProject/scuffgeoFiles` (for [<span class="SC">scuff-em</span> geometry files][Geometries]).

+ Before running any [[scuff-em]] calculations, set the
    environment variable `SCUFF_MESH_PATH` to the directory
    you created for mesh files:

````bash
% export SCUFF_MESH_PATH=~/myProject/mshFiles
````

(Or just include this line in any scripts you write to launch jobs; see below).

+ From the top-level project directory, create subdirectories for
    each separate run you plan to do. For example, to do separate
    runs for PEC and real gold, with coarse and fine resolutions 
    in both cases, I might create four directories called 
    `PEC_Coarse`, `PEC_Fine`, `Gold_Coarse`, and `Gold_Fine.`

+ Assuming you want to look at the same frequencies, evaluation points,
    incident fields, geometrical transformations, etc. in each case, 
    put files like `OmegaFile`, `EPFile`, `IFFile`, and `TransFile` 
    in the top-level directory (`~/myProject/`).

+ Within e.g. the `PEC_Fine` subdirectory, create a run script
    that explicitly specifies the locations in which it expects
    to find files, something like this:

````bash
#!/bin/bash

export FILEBASE=${HOME}/myProject
export SCUFF_MESH_PATH=${FILEBASE}/mshFiles
export SCUFF_GEO_PATH=${FILEBASE}/scuffgeoFiles

ARGS=""
ARGS="${ARGS} --geometry  ${SCUFF_GEO_PATH}/PEC_Fine.scuffgeo"
ARGS="${ARGS} --OmegaFile ${FILEBASE}/OmegaFile"
ARGS="${ARGS} --EPFile    ${FILEBASE}/EPFile"
ARGS="${ARGS} --IFFile    ${FILEBASE}/IFFile"  

scuff-scatter ${ARGS}

````

+ Now you can copy this run script to each new run directory 
    and make only minor changes (i.e. specify different `.scuffgeo`)
    files to launch the new job.

+ When running multiple jobs simultaneously on a single
    multi-core workstation, I usually use the environment variables
    `OMP_NUM_THREADS` and `GOMP_CPU_AFFINITY` to specify an 
    explicit divvying up of the available CPU cores so that 
    the various jobs don't step on each others' toes. (The
    OS scheduler should be able to do this automatically, but 
    I haven't had good luck with that.)
    
    For example, suppose I'm on a workstation that has 24 CPU cores,
    and I want to run 3 simultaneous [[scuff-em]] jobs.
    Then in the run scripts for the three jobs I will include
    the following lines:

````bash
...
export OMP_NUM_THREADS=8
export GOMP_CPU_AFFINITY="0-7"
...
````
for the first run script,

````bash
...
export OMP_NUM_THREADS=8
export GOMP_CPU_AFFINITY="8-15"
...
````

for the second run script, and 

````bash
...
export OMP_NUM_THREADS=8
export GOMP_CPU_AFFINITY="16-23"
...
````
for the third run script. Once all three jobs
are running, you can use 
[<span class="SC">htop</sc>](http://hisham.hm/htop/)
or similar utilities to double check that each 
job is running on its own set of 8 cores and 
not interfering with the other jobs.

[GMSH]:                        http://www.geuz.org/gmsh
[scuff-scatter]:               ../applications/scuff-scatter/scuff-scatter.md
[Geometries]:                  ../reference/Geometries.md
[Materials]:                   ../reference/Materials.md
[SiliconSlabs]:                ../examples/SiliconSlabs/SiliconSlabs.md
