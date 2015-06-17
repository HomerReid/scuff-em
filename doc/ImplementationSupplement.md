## The standalone scattering code

The standalone scattering code distributed with the [[scuff-em]] suite 
s 
[<span class="SmallCaps">scuff-scatter</span>][scuff-scatter].

This is essentially just a command-line interface to the core library for
solving scattering problems involving compact or extended 3D 
geometries using the algorithms outlined above.

(For extended geometries, the method described above is modified 
somwhat. In this case, the usual 3D Helmholtz kernels 
that enter into the computation of the BEM matrix elements
must be replaced by their periodic versions, which include contributions 
from all cells in the infinite lattice weighted by Bloch phase factors. 
To accelerate the computation of these lattice sums,
[[scuff-em]] implements the Ewald-summation method,
discussed (among other places) in  this paper: 

> "An Efficient Numerical Evaluation of the Green's Function for the
> Helmholtz Operator on Periodic Structures," by Kirk E. Jordan, 
> Gerard R. Richter, and Ping Sheng, 
> *Journal of Computational Physics* **63** 222 (1986)
> [http://dx.doi.org/10.1016/0021-9991(86)90093-8](http://dx.doi.org/10.1016/0021-9991(86)90093-8)

    <!---------------------------------------------------------------------->
    <!---------------------------------------------------------------------->
    <!---------------------------------------------------------------------->
    <p>
    <h2> The Casimir codes </h2>
    <hr> 

    <p>  
    The [[scuff-em]] suite includes three 
    standalone codes for Casimir computations.

    <p>
    <a href="scuff-EM/scuff-cas3D" class="CodeName">scuff-cas3d</a>
    is an implementation of the "fluctuating-surface-current" (FSC) 
    approach to Casimir computations for compact 3D objects. In this 
    case, the Casimir energy of a configuration of objects, and the 
    force on one of the objects, are given by the expressions
 
    <p align="center">
    <img align="center" src="scuff-em/reference/Equations/FSCEquations.png">

    <p>
    where the matrix **M** is simply the BEM matrix for the given 
    configuration of objects, evaluated at imaginary frequency 
    *&omega;=i&xi;*. (The matrix in the denominator of the first
    equation is just the BEM matrix for the configuration in which all
    objects are separated from each other by infinite distances, which
    amounts to just zeroing out the off-diagonal blocks of the original
    matrix. The derivative of the matrix in the second equation is 
    taken with respect to a rigid displacement of the object on which
    we are computing the force.)

    <p>
    In its default invocation, what [[scuff-cas3d]] 
    does is to evaluate the imaginary-frequency integrals in the above
    equations (or the Matsubara sums that replace them for finite-temperature
    calculations), with the BEM matrix and its derivatives computed at 
    each frequency by the standard [[libscuff]]
    routine for assembling the BEM matrix.
 
    <p>
    The FSC approach to Casimir computations, and in particular the fact 
    that the Casimir energy and force can be computed directly from the 
    BEM matrix, was first noted in my 
    <a href="research/memos.shtml#thesis">PhD thesis</a>
    (see also my <a href="research/talks.shtml#defense">thesis defense presentation</a>)
    and in these references:

    <ul>
     <p><li> "Fluctuating Surface Currents: A New Algorithm for 
              Efficient Prediction of Casimir Interactions among 
              Arbitrary Materials in Arbitrary Geometries. I. Theory,"
              by M. T. Homer Reid, Jacob White, and Steven G. Johnson
              (<a href="http://arxiv.org/abs/1203.0075v1">http://arxiv.org/abs/1203.0075v1</a>)
     <p><li> "Computation of Casimir Interactions between Arbitrary 
              3D Objects with Arbitrary Material Properties," by
              M. T. Homer Reid, Jacob White, and Steven G. Johnson, 
              *Physical Review A* **84** 010503(R) (2011)
              (<a href="http://dx.doi.org/10.1103/PhysRevA.84.010503">http://dx.doi.org/10.1103/PhysRevA.84.010503"</a>)
     <p><li> "Efficient Computation of Casimir Interactions between 
              Arbitrary 3D Objects," by M. T. Homer Reid, 
              Alejandro W. Rodriguez, Jacob White, Steven G. Johnson,
              *Physical Review Letters* **103** 040401 (2009)
              (<a href="http://dx.doi.org/10.1103/PhysRevLett.103.040401">http://dx.doi.org/10.1103/PhysRevLett.103.040401</a>)
    </ul>

    <p>  
    <a href="scuff-EM/scuff-cas2D" class="CodeName">scuff-cas2D</a>
    is an implementation of the FSC approach approach to Casimir 
    computations for *quasi-2D* objects---that is, 3D objects 
    of infinite spatial extent in one direction and constant 
    cross-sectional shape in the transverse directions. 

    <p>  
    <a href="scuff-EM/scuff-caspol" class="CodeName">scuff-caspol</a>
    computes the Casimir-Polder potential for a polarizable molecule 
    near a surface or a collection of surfaces. The method used 
    is essentially that discussed in Section 5.2 (equation 401) 
    of this paper:

    <p>
    <ul>
     <p><li> "Macroscopic QED - Concepts and Applications,"
              by Stefan Scheel and Stefan Yoshi Buhmann,
              *Acta Physica Slovaca* 58 (5), 675 (2008)
              (<a href="http://www.physics.sk/aps/pub.php?y=2008&pub=aps-08-05">http://www.physics.sk/aps/pub.php?y=2008&pub=aps-08-05</a>)
    </ul>

# The RF / microwave code </h2>
    <hr> 

    <p>  
    The microwave engineering code distributed with the 
    [[scuff-em]] suite is
    <a href="scuff-EM/scuff-RF" class="CodeName">scuff-rf.</a>

    <p> 
    The functionality implemented by [[scuff-rf]]
    is akin to that of [[scuff-scatter]], but with
    a couple of important extensions that are not present in the 
    [[scuff-em]] core library. These are:

    <p> 
    <ul>
      <p>
      <li> Support for RF/microwave ports in meshed geometric 
           structures. A *port* is a region of the structure
           through which the simulator forces a fixed current;
           the fields radiated by the currents forced through the
           port regions of the geometry are used as the incident 
           fields in a scattering calculation.
      <p>
      <li> Functionality to compute port "voltages" as the sum of 
           in the "scattered" scalar potential difference between 
           two points in a structure plus the line integral of the 
           vector potential between those points.
      <p>
      <li> Functionality to compute *z-parameters*  
           (*impedance parameters*) for
           a multiport network by driving a single port with a 
           unit-strength current and computing all the 
           resulting port voltages. (After we have computed the 
           matrix of *z-*parameters it is a simple transformation
           to convert them into the more commonly encountered
           *s-*parameters.)
    </ul>


<!----------------------------------------------------------------------->
<!-- end main page table                                               -->
<!----------------------------------------------------------------------->
</tr></table>
   
<!----------------------------------------------------------------------->
<!-- end main page body ------------------------------------------------->
<!----------------------------------------------------------------------->
