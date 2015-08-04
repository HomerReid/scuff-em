Computing the photonic local density of states with [[scuff-ldos]]
===========================================================

[[scuff-ldos]]
is a tool for computing the electromagnetic local density 
of states (LDOS) at points inside or outside compact or 
extended material bodies.

[TOC]

--------------------------------------------------
# 1. What <span class="SC">scuff-ldos</span> actually computes

--------------------------------------------------
<a name="CommandLineOptions"></a>
# 2. <span class="SC">scuff-ldos</span> command-line options

### Common options

[[scuff-ldos]] recognizes the following subset of the 
[list of commonly accepted options to <span class="SC">scuff-em</span> command-line codes][CommonOptions].

 
  ````
--geometry
--EPFile
--Omega
--OmegaFile
--OmegakBlochFile
--BZSymmetry
--AbsTol
--RelTol
--FileBase
--Cache
--ReadCache
--WriteCache
  ````
{.toc}

### Options requesting special calculations

  ````
--HalfSpace Aluminum
  ````
{.toc}

--------------------------------------------------
# 3. <span class="SC">scuff-ldos</span> output files

--------------------------------------------------
<a name="Examples"></a>
# 4. Examples of calculations using <span class="SC">scuff-ldos</span>

+ [LDOS above a dielectric half-space][HalfSpaceLDOS]

[HalfSpaceLDOS]:    ../../examples/HalfSpaceLDOS/HalfSpaceLDOS.md
