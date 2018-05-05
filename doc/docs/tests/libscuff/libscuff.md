<h1> Tests of the <span class=SC>scuff-em</span> core library </h1>

The `libscuff` unit of the [<span class=SC>scuff-em</span> test suite][TestSuite]
tests several elemental calculations implemented by `libscuff` routines
for assembling and post-processing the linear  SIE system.

[TOC]

{!LaTeXCommands.md!}

# SIE matrix elements: Singular and non-singular 4-dimensional integrals

This subunit tests `libscuff` routines for computing element of the SIE
system matrix, which involve singular and nonsingular 4-dimensional
integrals over pairs of triangular panels.

## SIE matrix elements in <span class=SC>scuff-em</span>

First we review the computation of SIE matrix elements in <span class=SC>scuff-em</span>.

### Matrix elements

In the most general case, each pair of RWG basis functions $\{\vb b_{\alpha}, \vb b_{\beta}\}$ 
contributes a $2\times 2$ block of entries to the system matrix:
<!----------------------------------------------------->
$$ \Delta \vb M(\vb b_{\alpha},  \vb b_{\beta})
   =\sum_r
    \left(\begin{array}{cc}
              i\omega \mu_r      \Vmv{\vb b_{\alpha}}{\mb G(k_r)}{\vb b_{\beta}} &
                      -ik_r      \Vmv{\vb b_{\alpha}}{\mb C(k_r)}{\vb b_{\beta}} \\
                      -ik_r      \Vmv{\vb b_{\alpha}}{\mb C(k_r)}{\vb b_{\beta}} &
             -i\omega \epsilon_r \Vmv{\vb b_{\alpha}}{\mb G(k_r)}{\vb b_{\beta}} \\
    \end{array}\right)
$$
<!----------------------------------------------------->
where $r$ runs over the indices of the (zero, one, or two) regions through which
the basis functions interact, $\{\epsilon_r, \mu_r, k_r\}$ are the material properties
and wavevector magnitude in region $r$, and $\mb G, \mb C\propto \nabla \times \mb G$
are the homogeneous dyadic Green's functions of Maxwell's equations.

If one or both basis functions lie on PEC surfaces then
$\Delta\vb M$ reduces to a $2\times 1$, $1\times 2$, or $1\times 1$ matrix.


### Edge-edge interactions

Matrix elements of the $\mb G$ and $\mb C$ operators are 4-dimensional integrals
over the supports of the basis functions that may be written in terms of the
scalar Greens' function:
<!----------------------------------------------------->
\begin{align}
 \mb G_{\alpha\beta} \equiv
  \Vmv{\vb b_{\alpha}}{\mb G}{\vb b_{\beta}}
&= \mc I_{\bullet \alpha\beta} -\frac{4}{k^2} \mc I_{\nabla \alpha\beta},
\qquad 
 \mb C_{\alpha\beta} \equiv
  \Vmv{\vb b_{\alpha}}{\mb C}{\vb b_{\beta}}
   \mc I_{\times \alpha\beta}
\\
\mc I_{\bullet\alpha\beta}
&\equiv \iint \Big( \vb b_{\alpha} \cdot \vb  b_{\beta}\Big) \Phi(r_{\alpha\beta}) \, d\vb x_{\alpha} d\vb x_{\beta}
\\
\mc I_{\nabla\alpha\beta}
&\equiv \iint \Phi(r_{\alpha\beta}) \, d\vb x_\alpha d\vb x_\beta
\\
\mc I_{\times\alpha\beta}
&\equiv \iint \Big[ \big(\vb b_{\alpha} \cdot \vb b_{\beta}\big)\cdot \vb r_{\alpha\beta}\Big]
 \Psi(r_{\alpha\beta}) \, d\vb x_{\alpha} d\vb x_{\beta}
\end{align}
<!----------------------------------------------------->
with
<!----------------------------------------------------->
$$ \vb r_{\alpha\beta} \equiv \vb x_{\alpha} - \vb x_{\beta},
   \qquad
   r_{\alpha\beta} \equiv |\vb r_{\alpha\beta}|
$$
$$ \Phi(r) \equiv \frac{e^{ik r}}{4\pi r},
   \Psi(r) \equiv (ikr-1)\frac{e^{ikr}}{4\pi r^3}.
$$
<!----------------------------------------------------->
<!----------------------------------------------------->

### Panel-panel interactions

Each edge-edge integral is a sum of 4 panel-panel integrals, e.g.

<!----------------------------------------------------->
$$ \mc I_{\bullet \alpha\beta}
=  \mc I^{++}_{\bullet\alpha\beta}
  -\mc I^{+-}_{\bullet\alpha\beta}
  -\mc I^{-+}_{\bullet\alpha\beta}
  +\mc I^{--}_{\bullet\alpha\beta}
$$
<!----------------------------------------------------->

(This is for the general case in which both basis functions are *full*
RWG basis functions; if one or both are *half*-RWG functions then the
number of contributing panel pairs reduces to 2 or 1.) 

The panel-panel integrals are
<!----------------------------------------------------->
\begin{align*}
 I^{\sigma\tau}_{\bullet\alpha\beta}
  &=\frac{\ell_{\alpha} \ell_{\beta}}{4A_{\alpha}^{\sigma} A_{\beta}^{\tau}}
     \int_{{\mc P}_{\alpha}^{\beta}} d \vb x_{\alpha}
     \int_{{\mc P}_{\beta}^{\tau}}   d \vb x_{\beta}
     \big(\vb x_{\alpha} - \vb Q_{\alpha}^{\sigma}\big) \cdot
     \big(\vb x_{\beta} - \vb Q_{\beta}^{\tau}\big) \Phi(r_{\alpha\beta})
\\
 I^{\sigma\tau}_{\nabla\alpha\beta}
  &=\frac{\ell_{\alpha} \ell_{\beta}}{4A_{\alpha}^{\sigma} A_{\beta}^{\tau}}
     \int_{{\mc P}_{\alpha}^{\beta}} d \vb x_{\alpha}
     \int_{{\mc P}_{\beta}^{\tau}} d \vb x_{\beta}
     \Phi(r_{\alpha\beta})
\\
 I^{\sigma\tau}_{\times\alpha\beta}
  &=\frac{\ell_{\alpha} \ell_{\beta}}{4A_{\alpha}^{\sigma} A_{\beta}^{\tau}}
     \int_{{\mc P}_{\alpha}^{\beta}} d \vb x_{\alpha}
     \int_{{\mc P}_{\beta}^{\tau}}   d \vb x_{\beta}
     \bigg[ \big(\vb x_{\alpha} - \vb Q_{\alpha}^{\sigma}\big)
             \times
             \big(\vb x_{\beta} - \vb Q_{\beta}^{\tau}\big)
             \cdot \vb r_{\alpha\beta}
     \bigg]
     \Psi(r_{\alpha\beta})
\end{align*}
<!----------------------------------------------------->

### Desingularization

For pairs of basis functions with 1 or more common vertices I
desingularize by subtracting off the first 3 nonvanishing terms
in the small-$r$ series for $\Phi(r)$ and $\Psi(r)$ and integrating
these terms separately.
Thus I put
<!----------------------------------------------------->
\begin{align*}
 \Phi(r) \equiv C_{-1} r^{-1} + C_0 + C_1 r + C_2 r^2  + \Phi^{\text{DS}{(r) \\
 \Psi(r) \equiv D_{-3} r^{-3} + D_{-1} r^{-1} + D_0 + D_1 r  + \Psi^{\text{DS}}(r)
\end{align*}
<!----------------------------------------------------->
and e.g.
<!----------------------------------------------------->
$$  \mc{I}_{\bullet\alpha\beta}
  = \sum_{n=-1}^2 C_n \mc{I}^{(n)}_{\bullet \alpha \beta}
                     +\mc{I}^{\text{DS}}_{\bullet \alpha \beta}
$$
<!----------------------------------------------------->

## Reference values

As reference values for testing `libscuff` routines for evaluating RWG integrals, I consider matrix elements between
RWG functions on the following surface mesh, in which the triangles are sufficiently regular that it's
easy to write down and evaluate the required 4-dimensional integrals .

For example, the positive-positive-panel contribution to $\mc{I}_{\nabla 0 0}$ reads
<!----------------------------------------------------->
$$ \mc{I}_{\nabla 0 0}^{++} = \frac{\ell^2}{4A_0^2}
   \iint_{A_0\times A_0} \Phi(r_{\alpha\beta})d\vb x_{\alpha} d\vb x_{\beta}
$$
<!----------------------------------------------------->
Parameterize the integral according to 
<!----------------------------------------------------->
$$ \vb x_{\alpha} = u_1 L \hat{\vb x} + u_2 L \hat{\vb y},
   \quad
   \vb x_{\beta } = u_3 L \hat{\vb x} + u_4 L \hat{\vb y}, 
   \quad
   r_{\alpha\beta} = L\sqrt{ (u_1-u_3)^3 + (u_2-u_4)^2}
$$
<!----------------------------------------------------->
to yield an integral over the unit 4-dimensional hypercube:
<!----------------------------------------------------->
$$
\mc I_{\nabla 0 0}^{++} = \ell^2 \int_{[0,1]^4} u_1 u_3 \Phi(r_{\alpha\beta}(\vb u)  )d\vb u
$$
<!----------------------------------------------------->

{!Links.md!}
