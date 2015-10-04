<!--#set 
    var="title" 
    value="Solving electrostatics problems with scuff-static"
  -->
<!--#include virtual="/pagetop.shtml"-->

<!-- begin main page body -->

<!----------------------------------------------------------------------->
<!-- main page table, with one row and three columns:                  -->
<!--  navbar, trough, body.                                            -->
<!----------------------------------------------------------------------->
<table cellspacing="0" cellpadding="0" width="100%"><tr>

<!----------------------------------------------------------------------->
<!-- left column of main page table: nav bar. --------------------------->
<!----------------------------------------------------------------------->
<td valign="top" width="180"> 
  <!--#include virtual="/research/navbar.shtml"-->
  </td> 

<!----------------------------------------------------------------------->
<!-- central column of main page table: separation trough --------------->
<!----------------------------------------------------------------------->
<td width="5%"></td>
 
<!----------------------------------------------------------------------->
<!-- right column of main page table: content of page. ------------------>
<!----------------------------------------------------------------------->
<td valign="top">
    
<p>
<br>
<p align="center">
   <table align="center"> 
     <tr> 
     <td>
<img border="0" width="300" height="150" 
              src="scuff-em/scuff-static/scuff-static-thumb.png">
</td>
     <td width="5%"> </td>
     <td> <h1>
          
Solving electrostatics problems with <span
class="CodeName">scuff-static</span>
</td>
     </tr>
   </table> 
    
    <!----------------------------------------------------------------------->
    <!-- big table encompassing the entire page ----------------------------->
    <!----------------------------------------------------------------------->
    <table align="center" width="90%"><tr><td>

    <p>
    <br>
    <p>
    <span class=CodeName>scuff-static</span> is a tool within the 
    <span class=CodeName>scuff-em</span> code suite for solving 
    a broad class of electrostatics problems.

    <p>
    The calculations that <span class=CodeName>scuff-static</span> can 
    perform include the following:
    <ul>

      <p>
      <li> Compute the capacitance matrix (i.e. the self- and mutual-
           capacitances) for a collection of conductors.
      </li>

      <p>
      <li> Compute the DC polarizability of a conducting or 
           dielectric body.
      </li>

      <p>
      <li> Compute the electrostatic potential and field
           at arbitrary user-specified points in the vicinity 
           of conducting or dielectric bodies, with the  
           conductors maintained at arbitrary user-specified 
           potentials and (optionally) an arbitrary user-specified
           external forcing field.
      </li>

      <p>
      <li> Compute the <i>C-matrix</i>, a sort of electrostatic
           version of the 
           <A href="scuff-em/scuff-tmatrix">``T-matrix''</a> 
           used to characterize the scattering properties
           of bodies at nonzero frequencies.
      </li>
    </ul>

    <p>
    (As a technical detail, we note that the implementation of 
     <span class=CodeName>scuff-static</span> actually differs in some 
     significant ways from the other codes in the 
     <span class=CodeName>scuff-em</span> suite; in particular,
     as compared to the 
     <span class=CodeName>scuff-em</span> core library,
     <span class=CodeName>scuff-static</span> uses different basis 
     functions and a fundamentally different formulation of the 
     boundary-element method, as appropriate for zero-frequency 
     problems. However, it turns out that the calculations 
     needed to implement the electrostatics calculations in 
     <span class=CodeName>scuff-static</span>
     are, for the most part, a <i>subset</i> of the calculations already 
     implemented in <span class=CodeName>scuff-em</span>, which
     is why it makes sense to package these codes together.)
     
    <p>
    Here is a brief 
    <a href="scuff-EM/scuff-static/scuff-static.pdf">technical memo</a>
    discussing the implementation of <span class=CodeName>scuff-static</span>,
    including both the underlying BEM electrostatics formulation
    and the execution of the various types of calculation
    (capacitance, polarizability, etc.) that the code can do.

    <!---------------------------------------------------->
    <!---------------------------------------------------->
    <!---------------------------------------------------->
    <p width="75%">
    <table class="TOC" cellpadding="5" cellspacing="5">

      <tr> <th> Table Of Contents </th></tr>

      <tr> <td>
           <a href="scuff-EM/scuff-static/scuffStaticOptions.shtml">
            1. <span class=CodeName>scuff-static</span> Command-Line Options 
           </a>
           </td>
      </tr>

      <tr> <td>
           <a href="scuff-EM/scuff-static/scuffStaticFiles.shtml">
            2. <span class=CodeName>scuff-static</span> Output Files
           </a>
           </td>

      <tr> <td>
           <a href="scuff-EM/scuff-static/scuffStaticExamples.shtml">
            3. <span class=CodeName>scuff-static</span> Examples 
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
    <!---------------------------------------------------->
    <!---------------------------------------------------->
    <!---------------------------------------------------->

    <!----------------------------------------------------------------------->
    <!----------------------------------------------------------------------->
    <!----------------------------------------------------------------------->
    <!--#include virtual="/scuff-EM/scuffEMFooter.shtml"-->

    <!----------------------------------------------------------------------->
    <!-- end big table encompassing the entire page ------------------------->
    <!----------------------------------------------------------------------->
    </td></tr></table>

<!----------------------------------------------------------------------->
<!-- end main page table                                               -->
<!----------------------------------------------------------------------->
</tr></table>
   
<!----------------------------------------------------------------------->
<!-- end main page body ------------------------------------------------->
<!----------------------------------------------------------------------->
   
<!--#include virtual="/pageend.shtml"-->



