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
 * TaylorDuffy.h -- header file for implementation of the taylor-duffy
 *               -- scheme for computing singular panel-panel integrals
 *
 * homer reid    -- 5/2009 -- 2/2012
 */
#ifndef TAYLORDUFFY_H
#define TAYLORDUFFY_H

#include <libhrutil.h>

namespace scuff {

#define TM_INTERVALS 10000

/***************************************************************/
/* constants ***************************************************/
/***************************************************************/
//values for the WhichG parameter to the Taylor_xx routines 
#define TM_RP                   0
#define TM_HELMHOLTZ            1
#define TM_GRADHELMHOLTZ        2
#define TM_HIGHK_HELMHOLTZ      3
#define TM_HIGHK_GRADHELMHOLTZ  4

//values for the WhichH parameter to the Taylor_xx routines
#define TM_ONE                 0
#define TM_DOT                 1
#define TM_DOTPLUS             2
#define TM_CROSS               3

// values for the WhichCase parameter. note the values correspond to
// the numbers of common vertices in each case.
#define TM_COMMONVERTEX        1
#define TM_COMMONEDGE          2
#define TM_COMMONTRIANGLE      3

/***************************************************************/
/* Data structure containing various data passed back and      */
/* forth among taylor-method routines.                         */
/***************************************************************/
typedef struct TMWorkspace
 {
   /* geometric data on triangles */
   double A2, B2, AP2, BP2, L2;
   double AdB, AdAP, AdBP, AdL, AdD, AdDP;
   double BdAP, BdBP, BdDP;
   double APdBP, APdD;
   double BPdD, BPdL;
   double DdDP;

   double AdQxQP, APdQxQP, BdQxQP, BPdQxQP, LdQxQP; 
   double V1xAdQmQP, V1xBdQmQP, V1xAPdQmQP, V1xBPdQmQP; 
   double AxBdQmQP, AxAPdQmQP, AxBPdQmQP, BxAPdQmQP, BxBPdQmQP;

   /* */
   int WhichCase;
   int NumPKs;
   int *PIndex;
   int *KIndex;
   cdouble *KParam;
   cdouble CurrentKParam;
   int TwiceIntegrable;

   int nCalls;

 } TMWorkspace;


/*--------------------------------------------------------------*/
/*- TaylorDuffy() ---------------------------------------------*/
/*--------------------------------------------------------------*/
typedef struct TaylorDuffyArgStruct
 { 
    int WhichCase;

    int NumPKs;
    int *PIndex;
    int *KIndex;
    cdouble *KParam;

    double *V1, *V2, *V3;
    double *V2P, *V3P;
    double *Q, *QP;

    double AbsTol, RelTol;
    int MaxEval;
    int ForceOnceIntegrable;

    int nCalls;

    cdouble *Result, *Error;

 } TaylorDuffyArgStruct;

void TaylorDuffy(TaylorDuffyArgStruct *Args);
void InitTaylorDuffyArgs(TaylorDuffyArgStruct *Args);

} // namespace scuff

#endif
