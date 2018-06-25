<h1> <span class=SC>scuff-solver</span>: A high-level interface to <span class=SC>scuff-em</span>
</h1>

Since its inception, <span class=SC>scuff-em</span> has offered two distinct
portals for accessing the underlying computational functionality of the core library:

+  The [suite of command-line tools][CommandLineSuite] such as
   [<span class=SC>scuff-ldos</span>][scuff-ldos] and
   [<span class=SC>scuff-neq</span>][scuff-neq] provide pre-implemented
   modules targeting specific physical problems such as scattering
   or heat transfer. 
   These tools are designed to be used off-the-shelf by users interested in 
   (say) Casimir forces but no urge to delve into the
   [integral-equation formalism][Implementation] underlying the computation, 
   and they have the advantage of offering access to
   <span class=SC>scuff-em</span> algorithms with no need to understand
   anything about what's going on under the hood.
   The difficulty is that, although tools such as 
   [<span class=SC>scuff-scatter</span>][scuff-scatter] aspire to offer a
   [rich smorgasbord of command-line options](/applications/scuff-scatter/scuff-scatter.md#Options)
   to encompass as wide range of possible use cases, the command-line interface 
   unavoidably suffers somewhat in the departments of customizability and interactivity.

+  Meanwhile, the [low-level API][libscuff] offers ultimate customizability
   to users willing to write C++ or python codes, but this requires some
   knowledge of integral-equation methods in general and
   <span class=SC>scuff-em</span> in particular; moreover, 

As a third portal aspiring to capture the best of both worlds,
<span class=SC>scuff-em</span> version 0.99
introduces a preliminary version of a *high-level interface.*
This is a glue layer, designed to sit between the low-level C++ solver
and a high-level scripting language such as python, that facilitates
interactive sessions for exploring <span class=SC>scuff-em</span>
functionality, as well as complex sequences of calculations driven by
scripts.

This page offers an initial peek at the
<span class=SC>scuff-em</span> high-level interface
and a flavor for the computational methodology it seeks to enable.
Readers should be aware that 
**the <span class=SC>scuff-em</span> high-level interface is a
work in progress, and we welcome suggestions for how it could be made
more convenient,** which may be contributed via 
the [GitHub issues page](https://github.com/HomerReid/scuff-em/issues).

Note: Fans of the finite-difference solver [<span class=SC>meep</span>](https://meep.readthedocs.io)
may note some parallels to the 
[high-level python interface](http://meep.readthedocs.io/en/latest/Python_Tutorials/Basics/)
recently added to that tool, which was one of the inspirations for the
<span class=SC>scuff-em</span> high-level interface.

[TOC]

## Objectives of the high-level interface

### Streamline the process of solving scattering problems

### Exploit symmetries

### Access application-specific modules

To avoid blowing the core library up to an unwieldy size,
certain application-specific algorithms---in particular,
the [electrostatics][scuff-static] and [RF][scuff-rf] modules---
are implemented as separate modules lying above <span class=SC>libscuff</span>.
The high-level interface makes all of this functionality available
on the same footing as core-library functionality.

## Example

Here's an example of a python script that uses the high-level interface
to drive a calculation involving microstrip devices.

{!Links.md!}
