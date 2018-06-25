# SCUFF-EM Release Notes

## SCUFF-EM 0.96 6/25/2018


  * Pre-release in anticipation of first official version 1.0 release

  * Added unit test suite, installed in $(PREFIX)/share/scuff-em/tests.
      --Tests currently implemented, and what they test:
          --Mie
              -- Core solver functionality (BEM matrix and RHS vector) for scattering from compact dielectric body
              -- Ccattered field components
              -- Induced dipole moments
              -- PFT algorithms: DSIPFT and EMTPFT
          --Fresnel
              -- Periodic boundary conditions, Ewald summation, acceleration via interpolation
              -- Direct computation of transmission/reflection coefficients from surface currents

  * Added high-level interface (libSCUFFSolver)
          -- revamped and modernized support for RF calculations
          -- More efficient handling of ports
          -- New algorithm for computing impedance (Z-) parameters
          -- Support for microstrip geometries via `libSubstrate` integration
          -- New tutorial examples of microstrip geometries, including python examples

  * Automatic detection of equivalent edge pairs, offering enormous speedup for highly structured meshes

  * Initial support for implicit dielectric substrates in full-wave calculations.
          -- Implemented by new `lib/libSubstrate` library

  * Support for geometries specified by GDSII files

  * New `scuff-spectrum` command-line application for computing resonant modes
          -- Uses Beyn algorithm for nonlinear eigenproblems, implemented by `lib/libBeyn` library
