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
 * libscuffInternals.h -- some definitions and prototypes that are used 
 *                     -- internally within libscuff code files, but which
 *                     -- probably don't need to be exported as part of the API
 *
 * homer reid          -- 10/2005 -- 10/2011
 */

#ifndef LIBSCUFFINTERNALS_H 
#define LIBSCUFFINTERNALS_H

#include "libscuff.h"
#include "rwlock.h"
#include "GBarAccelerator.h"

namespace scuff {

// values for the WhichAlgorithm field in the PanelPanelInteractions structure
// (low-order cubature, high-order cubature, taylor-duffy, high-k taylor-duffy, 
//  or desingularization)
#define PPIALG_LOCUBATURE    0
#define PPIALG_HOCUBATURE    1
#define PPIALG_TD            2
#define PPIALG_HKTD          3
#define PPIALG_DESING        4
#define NUMPPIALGORITHMS     5

/***************************************************************/ 
/* 1. argument structures for routines whose input/output      */
/*    interface is so complicated that an ordinary C++         */
/*    function prototype would be too unwieldy.                */
/*                                                             */
/*    note: in each case, we provide an Init() function that   */
/*    fills in default values for the lesser-used fields       */
/*    in the argument structure; however, this leaves many     */
/*    fields uninitialized, and it is up to the caller to      */
/*    fill those in.                                           */
/***************************************************************/

/*--------------------------------------------------------------*/
/*- TaylorDuffy() ----------------------------------------------*/
/*- 20121101 this content has been moved to TaylorDuffy.h       */
/*--------------------------------------------------------------*/

/*--------------------------------------------------------------*/
/*- GetPanelPanelInteractions() --------------------------------*/
/*--------------------------------------------------------------*/
typedef struct GetPPIArgStruct
 { 
   // input fields to be filled in by caller
   RWGSurface *Sa, *Sb;
   int npa, npb;
   int iQa, iQb;
   cdouble k;

   int NumGradientComponents;
   int NumTorqueAxes; 
   int ForceTaylorDuffy;
   double *GammaMatrix;
   void *opFC; // 'opaque pointer to FIPPI cache'

   // this is an optional 3-vector displacement applied to object b
   double *Displacement;

   // this field is filled in by the PanelPanelInteractions() routine
   // to indicate which of the various computational algorithms was   
   // used to compute the panel-panel integrals
   int WhichAlgorithm;

   // if this object is nonzero, it is used to compute the 
   // periodic kernel; otherwise we use the usual Helmholtz kernel.
   GBarAccelerator *GBA;
   bool ForceFullEwald;

   // output fields filled in by routine
   // note: H[0] = HPlus ( = HDot + (1/(ik)^2) * HNabla )
   // note: H[1] = HTimes
   // note: GradH[3*Mu + 0 ] = dHPlus/dR_\Mu
   // note: GradH[3*Mu + 1 ] = dHTimes/dR_\Mu
   // note: dHdT[3*Mu + 0 ] = dHPlus/dTheta_\Mu
   // note: dHdT[3*Mu + 1 ] = dHTimes/dTheta_\Mu
   cdouble H[2];
   cdouble GradH[6];
   cdouble dHdT[6];

 } GetPPIArgStruct;

void InitGetPPIArgs(GetPPIArgStruct *Args);
void GetPanelPanelInteractions(GetPPIArgStruct *Args);
void GetPanelPanelInteractions(GetPPIArgStruct *Args,
                               cdouble *H,
                               cdouble *GradH, 
                               cdouble *dHdT);

/*--------------------------------------------------------------*/
/*- GetEdgeEdgeInteractions() ----------------------------------*/
/*--------------------------------------------------------------*/
#define EEI_NOFORCE  0       // values for the 'Force' field
#define EEI_FORCE_PP 12345   //  of the GetEEIArgStruct
#define EEI_FORCE_SM 23456
typedef struct GetEEIArgStruct
 { 
   // input fields to be filled in by caller
   RWGSurface *Sa, *Sb;
   int nea, neb;
   cdouble k;

   int NumGradientComponents;
   int NumTorqueAxes; 
   double *GammaMatrix;

   // if this object is nonzero, it is used to compute the 
   // periodic kernel; otherwise we use the usual Helmholtz kernel.
   GBarAccelerator *GBA;
   bool ForceFullEwald;
   
   // this is an optional 3-vector displacement applied to object b
   double *Displacement;

   void *opFC; // 'opaque pointer to FIPPI cache'

   // this is used to force the code to use a specific
   // panel-integration algorithm; for diagnostic purposes only
   int Force;
   
   // this field provides diagnostic information on 
   // how many times the various panel-panel integral 
   // algorithms were invoked 
   unsigned PPIAlgorithmCount[NUMPPIALGORITHMS];

   // output fields filled in by routine
   // note: GC[0] = <f_a|G|f_b>
   // note: GC[1] = <f_a|C|f_b>
   // note: GradGC[3*Mu + 0 ] d/dR_\Mu (<f_a|G|f_b>)
   // note: GradGC[3*Mu + 1 ] d/dR_\Mu (<f_a|C|f_b>)
   // note: dGCdT[3*Mu + 0 ] d/dTheta_\Mu (<f_a|G|f_b>)
   // note: dGCdT[3*Mu + 1 ] d/dTheta_\Mu (<f_a|C|f_b>)
   cdouble GC[2];
   cdouble GradGC[6];
   cdouble dGCdT[6];

 } GetEEIArgStruct;

void InitGetEEIArgs(GetEEIArgStruct *Args);
void GetEdgeEdgeInteractions(GetEEIArgStruct *Args);

/*--------------------------------------------------------------*/
/*- GetSurfaceSurfaceInteractions() ----------------------------*/
/*--------------------------------------------------------------*/
typedef struct GetSSIArgStruct
 {
   // input fields to be filled in by caller
   RWGGeometry *G;
   RWGSurface *Sa, *Sb;
   cdouble Omega;    

   int NumTorqueAxes;
   double *GammaMatrix;

   int RowOffset, ColOffset;

   // if this flag is true, then the routine only
   // computes the upper triangle of the matrix, then
   // fills in the lower triangle assuming the matrix
   // block is symmetric
   bool Symmetric;

   // if these are nonzero, then the usual Helmholtz  
   // Green's function is replaced with its periodic
   // equivalent as computed using GBA1/2
   GBarAccelerator *GBA1, *GBA2;

   // this is an optional 3-vector displacement applied to object b
   double *Displacement;

   // these flags allow the caller to request the omission of 
   // contributions from one of the regions through which
   // the surfaces interact
   bool OmitRegion1, OmitRegion2;

   // if this flag is true, the call to GetSurfaceSurfaceInteraction
   // augments (does not overwrite) the matrix entries
   bool Accumulate;

   // output fields filled in by routine
   HMatrix *B;
   HMatrix **GradB;
   HMatrix **dBdTheta;

   // additional fields used internally that may be ignored by 
   // the caller both before and after the call
   double SignA, SignB;
   cdouble EpsA, EpsB;
   cdouble MuA, MuB;
   bool SaIsPEC, SbIsPEC;

 } GetSSIArgStruct;

void InitGetSSIArgs(GetSSIArgStruct *Args);
void GetSurfaceSurfaceInteractions(GetSSIArgStruct *Args);
void AddSurfaceSigmaContributionToBEMMatrix(GetSSIArgStruct *Args);

/***************************************************************/
/* 2. definition of data structures and methods for working    */
/*    with frequency-independent panel-panel integrals (FIPPIs)*/
/*                                                             */
/* note:                                                       */
/*  'FIPPI'    = 'frequency-independent panel-panel integral'  */
/*  'QIFIPPID' = 'Q-independent FIPPI data'                    */
/*  'QDFIPPID' = 'Q-dependent FIPPI data'                      */
/***************************************************************/

// a 'QIFIPPIData' is the minimal chunk of data that needs
// to be stored for each panel-panel pair.
typedef struct QIFIPPIData
 { 
   double xMxpRM3[3], xXxpRM3[3];
   double uvupvpRM1[9];
   double uvupvpR1[9];
   double uvupvpR2[9];
#define FIPPIFIX_20150306
#ifdef FIPPIFIX_20150306
   double x0[3];
#endif 
 } QIFIPPIData;

// a 'QDFIPPIData' is ...
//
typedef struct QDFIPPIData
 { 
   double hTimesRM3;
   double hDotRM1, hNablaRM1, hTimesRM1;
   double hDotR0,  hNablaR0,  hTimesR0;
   double hDotR1,  hNablaR1,  hTimesR1;
   double hDotR2,  hNablaR2;
 } QDFIPPIData;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void ComputeQIFIPPIData(double **Va, double **Vb, int ncv, QIFIPPIData *QIFD);

/*--------------------------------------------------------------*/
/*- GetQDFIPPIData is the basic routine that is exported to the */
/*- outside world.                                              */
/*--------------------------------------------------------------*/
void GetQDFIPPIData(double **Va, double *Qa, double **Vb, double *Qb, 
                    int ncv, void *opFC, QDFIPPIData *QDFD);

/*--------------------------------------------------------------*/
/* 'FIPPICache' is a class that implements efficient storage    */
/* and retrieval of QIFIPPIData structures for many panel pairs.*/
/* i am encapsulating this as its own separate class to allow   */
/* easy experimentation with various implementations.           */
/*--------------------------------------------------------------*/
class FIPPICache
 { 
  public:

    // constructor, destructor 
    FIPPICache();
    ~FIPPICache();

    // store/retrieve cache to/from binary file
    void Store(const char *FileName);
    void PreLoad(const char *FileName);
    
    // look up an entry 
    QIFIPPIData *GetQIFIPPIData(double **OVa, double **OVb, int ncv);

    int Hits, Misses;

  private:

    // any implementation of this class will have some kind of 
    // storage table, but to allow maximal flexibility in implementation
    // i am just going to store an opaque pointer to this table
    // in the class body, with all the details left up to the 
    // implementation 
    void *opTable;

    rwlock FCLock;

    char *PreloadFileName;
    unsigned int RecordsPreloaded;

 };

/***************************************************************/   
/* single global instance of FIPPICache that is used whenever  */   
/* no alternative is provided                                  */   
/***************************************************************/   
extern FIPPICache GlobalFIPPICache;

/****************************************************************/   
/*- AssessPanelPair counts common vertices in a pair of panels, */
/*- and puts arrays of panel vertices into certain orders that  */
/*- are expected by subsequent algorithms that work on the      */
/*- panel pairs.                                                */
/*- CanonicallyOrderVertices is an optional follow-up routine   */
/*- to AssessPanelPair that further orders the vertices in a    */
/*- canonical way for use in FIPPI cache lookups.               */
/****************************************************************/   
int AssessPanelPair(double **Va, double **Vb, double rMax);
int AssessPanelPair(double **Va, double **Vb);

int AssessPanelPair(RWGSurface *Sa, int npa, 
                    RWGSurface *Sb, int npb,
                    double *rRel, 
                    double **Va, double **Vb);

int AssessPanelPair(RWGSurface *Sa, int npa, 
                    RWGSurface *Sb, int npb,
                    double *rRel);

int NumCommonVertices(RWGSurface *Sa, int npa, RWGSurface *Sb, int npb);

int NumCommonBFVertices(RWGSurface *Sa, int nea, RWGSurface *Sb, int neb);

int CanonicallyOrderVertices(double **Va, double **Vb, int ncv,
                             double **OVa, double **OVb);

} // namespace scuff

#endif //LIBSCUFFINTERNALS_H
