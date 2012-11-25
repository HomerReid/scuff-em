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

namespace scuff {

// values for the WhichAlgorithm field in the PanelPanelInteractions structure
// (straight-up cubature, taylor-duffy, high-k taylor-duffy, or desingularization)
#define PPIALG_CUBATURE    0
#define PPIALG_TD          1
#define PPIALG_HKTD        2
#define PPIALG_DESING      3
#define NUMPPIALGORITHMS   4

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

   // if this field is nonzero, it points to an Interp3D object
   // for the kernel function; otherwise the kernel function 
   // is the usual Helmholtz kernel, possibly desingularized 
   Interp3D *GInterp;

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

   // if this object is nonzero, it is used as an interpolation 
   // table to compute values of the kernel (otherwise, the 
   // usual Helmholtz kernel is used) 
   Interp3D *GInterp;
   
   // this is an optional 3-vector displacement applied to object b
   double *Displacement;

   void *opFC; // 'opaque pointer to FIPPI cache'

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

   int Symmetric;

   // if this flag is true, it means the caller wants the interaction 
   // of the two surfaces as mediated by the periodic Green's function 
   // with the innermost 9 cell contributions omitted (the 'all-but-9' 
   // kernel). otherwise, we use the usual (direct) Helmholtz kernel.
   bool UseAB9Kernel;

   // this is an optional 3-vector displacement applied to object b
   double *Displacement;

   // these flags allow the caller to request the omission of 
   // contributions from one or both of the regions through which
   // the surfaces interact
   bool OmitRegion1, OmitRegion2;

   // if this field is true, the call to GetSurfaceSurfaceInteraction
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
   int SaIsPEC, SbIsPEC;
   Interp3D *GInterpA, *GInterpB;

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

int NumCommonVertices(RWGSurface *Sa, int npa, 
                      RWGSurface *Sb, int npb);

int CanonicallyOrderVertices(double **Va, double **Vb, int ncv,
                             double **OVa, double **OVb);

/***************************************************************/
/* routine for computing the periodic green's function via     */
/* ewald summation                                             */
/***************************************************************/
void GBarVDEwald(double *R, cdouble k, double *kBloch, double **LBV,
                 double E, int ExcludeFirst9, cdouble *GBarVD);

/***************************************************************/
/* this is an alternative interface to GBarVDEwald that has the*/
/* proper prototype for passage to my Interp3D class routines  */
/***************************************************************/
typedef struct GBarData 
 { 
   cdouble k;           // wavenumber 
   double *kBloch;      // bloch vector 
   double *LBV[2];      // lattice basis vectors 
   double E;            // ewald separation parameter
   bool ExcludeInner9;  
 
 } GBarData;

void GBarVDPhi3D(double X1, double X2, double X3, 
                 void *UserData, double *PhiVD);

} // namespace scuff

#endif //LIBSCUFFINTERNALS_H
