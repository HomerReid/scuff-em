/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * PanelCubature.h -- header file for PanelCubature and PanelPanelCubature
 *                 -- routines for integrating over panels in libscuff
 *                 -- geometries.
 *
 * homer reid      -- 5/2013
 */
#ifndef PANELCUBATURE_H
#define PANELCUBATURE_H

#include <libscuff.h>

namespace scuff{

/***************************************************************/
/* PCData is the data structure passed to the user's PCFunction*/
/* routine; it contains data on the panel and (optionally) the */
/* surface-current distribution on the panel.                  */
/***************************************************************/
typedef struct PCData
 { 
   double *nHat;       // panel normal
   cdouble Omega;      // angular frequency
   cdouble *K, *N;     // surface currents at eval point
   cdouble DivK, DivN; // divergences of surface currents

 } PCData;

/***************************************************************/
/* This is the prototype for the integrand routine to be       */
/* written by the user and passed as a parameter to            */
/* GetPanelCubature.                                           */
/***************************************************************/
typedef void (*PCFunction)(double *x, PCData *PCD,
                           void *UserData, double *Result);

/***************************************************************/
/* Routine to evaluate a two-dimensional numerical integral    */
/* over a single panel.                                        */
/***************************************************************/
void GetPanelCubature(RWGGeometry *G, int ns, int np,
                      PCFunction Integrand, void *UserData, int IDim,
                      int MaxEvals, double RelTol, double AbsTol,
                      cdouble Omega, HVector *KN, int iQ, double BFSign,
                      double *Result);

/***************************************************************/
/* Similar to GetPanelCubature, except the integral is taken   */
/* over the support of an RWG basis function, which amounts to */
/* computing two panel cubatures.                              */
/***************************************************************/
void GetBFCubature(RWGGeometry *G, int ns, int ne,
                   PCFunction Integrand, void *UserData, int IDim,
                   int MaxEvals, double RelTol, double AbsTol,
                   cdouble Omega, HVector *KN,
                   double *Result);

/***************************************************************/
/* PPCData is the structure passed to the user's PPCFunction   */
/* routine; it contains data on the two panels and (optionally)*/
/* the surface-current distributions on the two panels.        */
/***************************************************************/
typedef struct PPCData
 { 
   double *nHat1, *nHat2;       // panel normals
   cdouble Omega;               // angular frequency
   cdouble *K1, *N1, *K2, *N2;  // surface currents at eval points
   cdouble DivK1, DivN1, DivK2, DivN2; // divergences of surface currents

 } PPCData;

/***************************************************************/
/* This is the prototype for the integrand routine to be       */
/* written by the user and passed as a parameter to            */
/* GetPanelPanelCubature.                                      */
/***************************************************************/
typedef void (*PPCFunction)(double *x, double *xp, PPCData *PPCD,
                            void *UserData, double *Result);

/***************************************************************/
/* Routine to evaluate a four-dimensional numerical cubature   */
/* over a pair of panels.                                      */
/***************************************************************/
void GetPanelPanelCubature(RWGGeometry *G, int ns1, int np1, int ns2, int np2,
                           double *Displacement,
                           PPCFunction Integrand, void *UserData, int IDim,
                           int MaxEvals, double RelTol, double AbsTol,
                           cdouble Omega, HVector *KN, 
                           int iQ1, double BFSign1, int iQ2, double BFSign2,
                           double *Result);

/***************************************************************/
/* Similar to GetPanelPanelCubature, except the integral is    */
/* over the supports of two RWG basis functions, which amounts */
/* to computing four panel-panel cubatures.                    */
/***************************************************************/
void GetBFBFCubature(RWGGeometry *G, int ns1, int ne1, int ns2, int ne2,
                     double *Displacement,
                     PPCFunction Integrand, void *UserData, int IDim,
                     int MaxEvals, double RelTol, double AbsTol,
                     cdouble Omega, HVector *KN, 
                     double *Result);


/***************************************************************/
/* 20151118 streamlined implementation of GetBFCubature and    */
/* GetBFBFCubature that are 2x and 4x faster respectively.     */
/*                                                             */
/* In these versions, the user's integrand function must       */
/* *accumulate* contributions to Integral with weight Weight.  */
/***************************************************************/
typedef void PCFunction2(double x[3], double b[3], double Divb,
                         void *UserData, double Weight, double *Integral);


typedef void PPCFunction2(double x[3], double b[3], double Divb,
                          double xp[3], double bp[3], double Divbp,
                          void *UserData,
                          double Weight, double *Integral);

void GetBFCubature2(RWGGeometry *G, int ns, int ne,
                    PCFunction2 Integrand, void *UserData, int IDim,
                    int Order, double *Integral, int PanelOnly=-1);

void GetBFBFCubature2(RWGSurface *S, int ne, 
                      RWGSurface *SP, int neP,
                      PPCFunction2 Integrand, void *UserData, int IDim,
                      int Order, double *Integral, 
                      int PanelOnlyA=-1, int PanelOnlyB=-1);

void GetBFBFCubature2(RWGGeometry *G,
                      int ns, int ne, int nsP, int neP,
                      PPCFunction2 Integrand, void *UserData, int IDim,
                      int Order, double *Integral, 
                      int PanelOnlyA=-1, int PanelOnlyB=-1);

void GetBFBFCubatureTD(RWGSurface *SA, int neA,
                       RWGSurface *SB, int neB,
                       PPCFunction2 Integrand, void *UserData,
                       int IDim, double *Integral, double *Error,
                       int MaxEvals=0, double AbsTol=0.0, double RelTol=1.0e-4);

void GetBFBFCubatureTD(RWGGeometry *G,
                       int nsA, int neA, int nsB, int neB,
                       PPCFunction2 Integrand, void *UserData,
                       int IDim, double *Integral, double *Error,
                       int MaxEvals=0, double AbsTol=0.0, double RelTol=1.0e-4);

} // namespace scuff{


#endif // #ifndef PANELCUBATURE_H
