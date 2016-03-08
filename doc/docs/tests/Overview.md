# Overview of the <span class="SC">scuff-em</span> Validation Test Suite

The [[scuff-em]] distribution includes a number of validation tests
that use the various dedicated
[application modules][Applications]
to solve physics problems with known analytical solutions.

In addition to demonstrating the core functionality of
the application modules in simple cases,
these to 
these are used as unit tests to catch bugs and regressions
in the [[scuff-em]] development process; you can
also use them to verify the correct performance of your 
[[scuff-em]] installation.

## Descriptions of the individual tests in the validation suite

The various tests in the test suite are described in detail
on the following pages, which also present comparisons
of [[scuff-em]] results to known analytical solutions.

* [Mie scattering](MieScattering/MieScattering.md)
* [Fresnel scattering](FresnelScattering/FresnelScattering.md)
* [Equilibrium Casimir forces between spheres](CasimirSpheres/CasimirSpheres.md)
* [Equilibrium Casimir forces between plates](CasimirPlates/CasimirPlates.md)
* [Equilibrium Casimir-Polder potential near a sphere](CPSphere/CPSphere.md)
* [Equilibrium Casimir-Polder potential near a plate](CPPlate/CPPlate.md)
* [Heat transfer and non-equilibrium Casimir forces between spheres](NEQSpheres/NEQSpheres.md)
* [Low-level tests of the <span class="CodeName">scuff-em</span> core library](libscuff/libscuff.md)

## Running the [[scuff-em]] tests

The [[scuff-em]] distribution includes a command-line test harness
application that you may use to run individual tests, or the entire
test suite, to verify correct performance of your [[scuff-em]]
installation.

The test harness is named ``scuff-test-harness,`` and it lives
in the `tests` subdirectory of the [[scuff-em]] repository.
The `tests` folder also contains

* a collection of text files with file extension `.scuffTest`,
  each of which describes a single validation test in a format
  understood by `scuff-test-harness`

* various input files (geometry files, surface meshes, etc.) 
  needed to run the tests

* a collection of data files with file extension `.ref`,
  which contain the correct results of the validation tests
  and are used by `scuff-test-harness` to determine whether 
  or not the test succeeded.

### Running individual tests

To run an individual validation test, use the `--test` option
to `scuff-test-harness:`

```
 % scuff-test-harness --test MieScattering
```

The argument passed to `--test` should be the base
file name of one of the `.scuffTest` files.

### Running the full test suite or a subset

To run the entire suite of validation tests, use the `--allTests`
option:

```
 % scuff-test-harness --allTests
```

## Checking the results of validation tests
 
[Applications]:        ../reference/TopLevel.md#AvailableApplications
