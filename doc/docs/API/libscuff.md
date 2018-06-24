<h1> <span class=SC>libscuff</span>: The core library
   underlying the <span class=SC>scuff-em</span> suite
</h1>

The back engine implementing the computational algorithms that underlie
the <span class=SC>scuff-em</span> suite is a C++ library called
<span class=SC>libscuff.</span> If you find that the 
[standalone application programs][applications]
distributed with <span class=SC>scuff-em</span> are not sufficiently
flexible for your computational needs, you may wish to write your own
codes that access the core library directly.

The most direct way to access the core library is of course to write 
your own C/C++ codes that link against the 
<span class=SC>libscuff</span> library binary
(`libscuff.a` or `libscuff.so`), 
but the <span class=SC>scuff-em</span> distribution also comes
with wrappers that allow you to access the core library 
from <span class="SC">python</span>.

The primary entity in <span class=SC>libscuff</span> 
is a C++ class named `RWGGeometry.` (The name 
derives from the fact that <span class=SC>libscuff</span>
uses a set of basis functions known as 
["RWG functions"](http://dx.doi.org/10.1109/TAP.1982.1142818)
for representing surface currents.) An `RWGGeometry`
stores all necessary data on a scattering geometry described by a
[`.scuffgeo` file](scuff-EM/reference/scuffEMGeometries.shtml),
and provides class methods for setting up, solving, and 
procssing the solutions of scattering problems in that
geometry.

The functions provided by <span class=SC>libscuff</span> 
may be broadly classified into three categories:
*main flow* routines, which handle each of the major steps 
in the process of setting up and solving electromagnetic scattering 
problems; *ancillary* routines on `RWGGeometry` 
structures, which offer various additional functionality such as 
visualization; and lower-level *support* routines, which 
allow you to manipulate some of the constituent entities 
(such as matrices and vectors) involved in scattering problems.

The <span class=SC>libscuff</span> API may be accessed from C++ or 
<span class=SC>python</span> 
programs. In the former case, source files should 
include the lines 

```C++
#include <libscuff.h>
using namespace scuff;
```

and should be linked with `-lscuff.`

Python code should include the line

 
```python
  import scuff;
```

+ [1. Core Library Reference: Main flow routines][MainFlowAPI]
+ [2. Core Library Reference: Describing incident fields][IncFields]
+ [3. Core Library Reference: Ancillary routines][Ancillary]
+ [4. Core Library Reference: Utility sublibraries][SublibraryAPI]
+ [5. Core Library Reference: C++ examples][CXXExamples]
+ [6. Core Library Reference: Python examples][PythonExamples]

{!Links.md!}
