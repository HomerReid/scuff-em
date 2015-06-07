# General reference material for [[scuff-em]] command-line applications

This page collects some general information that applies to
most or all of the standalone command-line 

## Commonly accepted command-line arguments

Note that not all codes accept all arguments 
(for example, [[scuff-transmission]] does not accept ``--TransFile``),
but the format of each option is standardized among all codes that 
*do* accept that option.

** Arguments specifying inputs to a calculation **

+ ``--geometry  MyGeometry.scuffgeo``

> Specifies the [[scuff-em]] geometry file. 

+ ``--TransFile MyTransformations.trans`` 

> Specifies a list of geometrical transformations.

** Arguments specifying outputs from a calculation**

+ ``--FileBase MyFileBase`` 

> Specifies the base file name for output files (so that, for example,
> the frequency-resolved output file written by [[scuff-cas3d]]
> will be ``MyFileBase.byXi``, while the frequency-integrated 
> output file will be ``MyFileBase.out.`` If this option is not 
> specified, the file base is taken to be the base filename of the 
> ``.scuffgeo`` file.

## Passing command-line options via text file

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
## Complex numbers

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
