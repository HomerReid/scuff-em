# Power, force, and torque in Mie scattering

This test validates the algorithms implemented in
[[scuff-em]] for computing the power, force, and torque (PFT)
on bodies irradiated by external fields. The test 
uses the 
[<span class="SC">scuff-scatter</span>][scuff-scatter] application module
to compute the PFT for a lossy dielectric sphere
irradiated by a circularly-polarized plane wave
and compares the results to the predictions of 
Mie scattering theory.

[scuffScatter]:               ../../applications/scuff-scatter/index.md
[scuffEMGeometries]:          ../../reference/Geometries.md
[scuffEMTransformations]:     ../../reference/Transformations.md
[scuffEMMaterials]:           ../../reference/Materials.md
[scuffEMInstallation]:        ../../reference/Installation.md
[EmigPaper]:                  http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.170403
[gnuplot]:                    https://www.gnuplot.info
[BohrenHuffman]:              http://onlinelibrary.wiley.com/book/10.1002/9783527618156
[MarstonCrichton]:            http://journals.aps.org/pra/abstract/10.1103/PhysRevA.30.2508
