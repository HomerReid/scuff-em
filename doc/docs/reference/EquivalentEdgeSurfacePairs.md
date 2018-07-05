# Automatic detection of equivalent edges and surfaces in <span class="SC">scuff-em</span>

<span class=SC>scuff-em</span> versions 0.96 (June 2018) and later
include a new optimization that automatically detects and exploits
two common symmetries in
[<span class=SC>scuff-em</span> surface-mesh geometries][scuffEMGeometries],
namely

+ For surface meshes that are *highly structured* (contain many pairs of
    identical triangles), <span class=SC>scuff-em</span>
    implements *automatic detection of equivalent edge pairs.*
    This can yield an *enormous* reduction (as much as 90% or more)
    in the cost of assembling the system matrix, as well as
    computing power, force, and torque via the EMTPFT algorithm.

+ For geometries containing pairs of identical surface meshes,
    <span class=SC>scuff-em</span> implements
    implements *automatic detection of equivalent surface pairs.*
    This allows the assembly of entire blocks of the system matrix
    to be bypassed entirely.

[TOC]

<a name="EquivalentEdgePairs"></a>
## Equivalent edge-pair detection

Insert documentation here!

Note: In <span class=SC>scuff-em</span> version 0.96, 
equivalent edge-pair detection is *disabled* by default,
but may be enabled by setting the environment variable
`SCUFF_MATRIX_2018=1`.


<a name="EquivalentSurfacePairs"></a>
## Equivalent surface-pair detection

{!Links.md!}
