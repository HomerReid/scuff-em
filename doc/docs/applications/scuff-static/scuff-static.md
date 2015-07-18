<h1> Solving electrostatics problems with
     <span class="SC">scuff-static</span>
</h1>

[[scuff-static]] is a tool within the 
[[scuff-em]] code suite for solving 
a broad class of electrostatics problems.


The calculations that [[scuff-static]] can 
perform include the following:

+Compute the capacitance matrix (i.e. the self- and mutual-
capacitances) for a collection of conductors.

+Compute the DC polarizability of a conducting or 
dielectric body.

+Compute the electrostatic potential and field
at arbitrary user-specified points in the vicinity 
of conducting or dielectric bodies, with the  
conductors maintained at arbitrary user-specified 
potentials and (optionally) an arbitrary user-specified
external forcing field.

+Compute the *C-matrix*, a sort of electrostatic
version of the 
["T-matrix"](../scuff-tmatrix/scuff-tmatrix.md)
used to characterize the scattering properties
of bodies at nonzero frequencies.

(As a technical detail, we note that the implementation of 
[[scuff-static]] actually differs in some 
significant ways from the other codes in the 
[[scuff-em]] suite; in particular,
as compared to the 
[[scuff-em]] core library,
[[scuff-static]] uses different basis 
functions and a fundamentally different formulation of the 
boundary-element method, as appropriate for zero-frequency 
problems. However, it turns out that the calculations 
needed to implement the electrostatics calculations in 
[[scuff-static]]
are, for the most part, a *subset* of the calculations already 
implemented in [[scuff-em]], which
is why it makes sense to package these codes together.)


Here is a brief 
[technical memo](scuff-static.pdf)
discussing the implementation of [[scuff-static]],
including both the underlying BEM electrostatics formulation
and the execution of the various types of calculation
(capacitance, polarizability, etc.) that the code can do.

<!---------------------------------------------------->
<!---------------------------------------------------->
<!---------------------------------------------------->
<p width="75">
<table class="TOC" cellpadding="5" cellspacing="5">

<tr> <th> Table Of Contents </th></tr>

<tr> <td>
<a href="scuff-EM/scuff-static/scuffStaticOptions.shtml">
1. [[scuff-static]] Command-Line Options 
</a>
</td>
</tr>

<tr> <td>
<a href="scuff-EM/scuff-static/scuffStaticFiles.shtml">
2. [[scuff-static]] Output Files
</a>
</td>

<tr> <td>
<a href="scuff-EM/scuff-static/scuffStaticExamples.shtml">
3. [[scuff-static]] Examples 
</a>
<blockquote>
<table cellpadding="5" cellspacing="5">

<tr><td>
<a href="scuff-EM/scuff-static/scuffStaticExamples.shtml#Capacitance">
3a. Self- and mutual-capacitance of irregularly shaped 
conductors
</a>
</td></tr>

<tr><td>
<a href="scuff-EM/scuff-static/scuffStaticExamples.shtml#Platonic">
3b.  Polarizability of the platonic solids 
</a>
</td></tr>

<tr><td>
<a href="scuff-EM/scuff-static/scuffStaticExamples.shtml#Gates">
3c. Electrostatic fields in the vicinity of a complicated
gate array
</a>
</td></tr>
</table>
</blockquote>
</td>
</tr>
</table>
