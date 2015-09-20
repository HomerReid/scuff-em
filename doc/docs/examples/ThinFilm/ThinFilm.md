# 4c. Spatially-resolved study of plane-wave transmission through an infinite-area thin dielectric film

The previous examples dealt with *compact* scatterers. We'll next consider an [*extended*](scuff-em/reference/scuffEMGeometries.shtml#ExtendedGeometries) geometry -- namely, a thin dielectric film of finite thickness in the *z* direction but infinitely extended in the *x* and *y* directions. This is the same geometry for which we used scuff-transmission to look at the plane-wave transmission and reflection coefficients as a function of frequency in [this example,](scuff-em/scuff-transmission/index.shtml#ThinFilm) but here we'll do a different calculation -- namely, we'll pick a single frequency and look at how the electric and magnetic fields vary in space, both inside and outside the thin film. (The files for this example may be found in the `share/scuff-em/examples/ThinFilm` directory of the scuff-em installation.)

The mesh file and `.scuffgeo` file for this geometry are discussed in the [documentation for scuff-transmission.](scuff-em/scuff-transmission/index.shtml#ThinFilm) The geometry consists of a film of thickness *T*=1μm, with relative dielectric constant ε<sup>*r*</sup>=100, illuminated from below by a plane wave at normal incidence. (We'll take the incident field to be linearly polarized with **E** field pointing in the *x* direction.) The lower and upper surfaces of the film are at *z=0* and *z=T.* For this geometry it is easy to solve Maxwell's equation directly to obtain the **E** and **H** fields directly at points below, within, and above the film:

![](scuff-em/scuff-scatter/EHvsz.png)

We will try to reproduce this behavior using scuff-scatter. First create a little text file ([ThinFilm.EvalPoints](scuff-em/scuff-scatter/ThinFilm.EvalPoints)) containing the coordinates of a bunch of points on a straight line passing from below the film to above the film. Then put the following command-line arguments into a file called `ThinFilm_58.args:`

~~~~ {.Listing}
 geometry ThinFilm_58.scuffgeo
 cache ThinFilm_58.scuffcache
 omega 1.0
 EPFile ThinFilm.EvalPoints
 pwDirection 0 0 1
 pwPolarization 1 0 0
    
~~~~

and pipe it into scuff-scatter:

~~~~ {.Listing}
 scuff-scatter < scuff-scatter.args
    
~~~~

This produces files named `ThinFilm.scattered` and `ThinFilm.total`. Plotting the 4th vs the 3rd column [(Re *E<sub>x</sub>* vs. *z*)](scuff-em/scuff-scatter/scuffScatterFiles.shtml#FieldOutput) as well as the the 12th vs the 3rd column [(Re *H<sub>y</sub>* vs. *z*)](scuff-em/scuff-scatter/scuffScatterFiles.shtml#FieldOutput) of the `.total` file yields excellent agreement with theory:

![](scuff-em/scuff-scatter/ThinFilmEField.png)

![](scuff-em/scuff-scatter/ThinFilmHField.png)
