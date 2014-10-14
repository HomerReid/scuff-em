scuff-em
========

A comprehensive and full-featured computational physics suite for boundary-element 
analysis of electromagnetic scattering, fluctuation-induced phenomena (Casimir forces 
and radiative heat transfer), nanophotonics, RF device engineering, electrostatics, 
and more. Includes a core library with C++ and python APIs as well as many command-line 
applications.

SCUFF-EM stands for

 Surface-CUrrent-Field Formulation of ElectroMagnetism

For documentation and further information on SCUFF-EM, visit the SCUFF-EM homepage:

http://www.homerreid.com/scuff-em

The SCUFF-EM distribution includes the following components:

 -- libscuff, a C++ library that implements the basic boundary-element
    solver functionality

 -- a number of ancillary C++ libraries implementing utility functions

 -- a collection of standalone application programs, built atop
    libscuff and its associated libraries, for solving specialized
    problems in computational electromagnetism and related fields.

 -- a collection of simple C++ code examples demonstrating how to use
    the libscuff API

The standalone application programs distributed with SCUFF-EM
include the following:

 -- scuff-scatter: A comprehensive tool for computing the fields
                   scattered from arbitrary geometries irradiated
                   by arbitrary user-specified incident fields


 -- scuff-RF:      A tool for analyzing passive RF and microwave 
                   devices: RF antennas, coaxial cables, MRI coils,
                   lumped-element passives, etc. 
                   scuff-RF can compute 

 -- scuff-cas2D
 -- scuff-cas3D
 -- scuff-caspol:  Tools for computing Casimir forces between 2D and 
                   3D objects of arbitrary geometries and material 
                   properties, as well as Casimir-Polder potentials for
                   polarizable molecules in the vicinity of complex
                   material surfaces.

 -- scuff-neq      A tool for modeling non-equilibrium fluctuation-induced
                   electromagnetic phelomena such as radiative heat transfer,
                   thermal self-propulsion and self-rotation, and 
                   non-equilibrium Casimir forces. 

 -- scuff-static   A tool for computational electrostatics: Capacitance
                   matrices of conductor arrays, solutions of Laplace's
                   equation in the presence of dielectrics and conductors 
                   held at fixed potentials, etc.

