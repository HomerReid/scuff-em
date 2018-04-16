# Overview of the <span class="SC">scuff-em</span> Validation Test Suite

The [[scuff-em]] distribution includes a number of validation tests
that use both the
[core library][CoreLibrary] and the
various dedicated
[application modules][Applications]
to solve physics problems with known (usuallly analytical) solutions.

One purpose of the test suite is to serve as battery of unit tests
to catch bugs and regressions in the [[scuff-em]] development process.
You can run the full unit-test suite in your local copy of the repository
by saying `% make check` from the top of the source tree.
The unit-tests are also incorporated into
[Travis continuous-integration framework](https://travis-ci.org/)
to run the test suite automatically whenever any new code is
committed to the 
[<span class=SC>scuff-em</span> Github repository](http://github.com/HomerReid/scuff-em);
you can check the status of the test suite in the latest repository versions
[here](https://travis-ci.org/HomerReid/scuff-em) or from the
icon in the README section of the 
[GitHub page](http://github.com/HomerReid/scuff-em).

Beyond unit testing---which is really intended to catch minor things like
inadvertent programming errors and missing files---the test suite is
also useful for benchmarking and validating the performance of
<span class=SC>scuff-em</span>
on nontrivial physics and engineering problems.

## Structure of the test suite

The various individual tests live in subfolders of the `tests` subdirectory
of the <span class=SC>scuff-em</span> source repository; after
`make install` they are copied to the directory `$(prefix)/share/scuff-em/tests.`
The tests come in various formats:
    --C++ codes to be compiled and linked against `libscuff`
    --python scripts to be run with the <span class=SC>scuff-em</span> python module, and
    --shell scripts that invoke [`scuff-scatter`](../applications/scuff-scatter/scuff-scatter.md) or other
     [<span class=SC>scuff-em</span> command-line applications](../applications)
     and compare the resulting output files to reference files.

### `CheckSCUFFData`

Upon a successful <span class=SC>scuff-em</span> build and install,
a binary utility code called `CheckSCUFFData` will be placed
in `$(prefix)/share/scuff-em/tests`. This is a simple utility
used by some tests to compare data files produced by 
[<span class=SC>scuff-em</span> command-line applications](../applications)
to reference data files. See the examples below for
information on how this tool is used.

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

## Checking the results of validation tests
 
[CoreLibrary]:         ../API/libscuff.md
[Applications]:        ../reference/TopLevel.md#AvailableApplications../API
