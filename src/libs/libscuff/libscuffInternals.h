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
/*- GetPanelPanelInteractions() --------------------------------*/
/*--------------------------------------------------------------*/
typedef struct GetPPIArgStruct
 { 
   // input fields to be filled in by caller
   RWGObject *Oa, *Ob;
   int npa, npb;
   int iQa, iQb;
   cdouble k;

   int NumGradientComponents;
   int NumTorqueAxes; 
   int ForceTaylorDuffy;
   double *GammaMatrix;
   void *opFT;

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
   RWGObject *Oa, *Ob;
   int nea, neb;
   cdouble k;

   int NumGradientComponents;
   int NumTorqueAxes; 
   double *GammaMatrix;
   void *opFT;

   int Force;

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
/*- AssembleBEMMatrixBlock() -----------------------------------*/
/*--------------------------------------------------------------*/
typedef struct ABMBArgStruct
 {
   // input fields to be filled in by caller
   RWGGeometry *G;
   RWGObject *Oa, *Ob;
   cdouble Frequency;
   int nThread;

   int NumTorqueAxes;
   double *GammaMatrix;
  
   int RowOffset, ColOffset;

   int Symmetric;

   // output fields filled in by routine
   HMatrix *B;
   HMatrix **GradB;
   HMatrix **dBdTheta;

   // additional fields used internally that may be ignored by 
   // the called both before and after the call
   double Sign;
   cdouble EpsA, EpsB; 
   double MuA, MuB;
   int OaIsPEC, ObIsPEC;

 } ABMBArgStruct;

void InitABMBArgs(ABMBArgStruct *Args);
void AssembleBEMMatrixBlock(ABMBArgStruct *Args);

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
                    int ncv, void *opFT, QDFIPPIData *QDFD);

/*--------------------------------------------------------------*/
/* 'FIPPITable' is a class that implements efficient storage    */
/* and retrieval of QIFIPPIData structures for many panel pairs.*/
/* i am encapsulating this as its own separate class to allow   */
/* me easily to experiment with various implementations.        */
/*--------------------------------------------------------------*/
class FIPPITable
 { 
  public:
    FIPPITable();
    ~FIPPITable();
    QIFIPPIData *GetQIFIPPIData(double **OVa, double **OVb, int ncv);

    int DoNotCompute;
 
 };

/***************************************************************/   
/* 3. some additional non-class-method routines used internally*/
/***************************************************************/
cdouble TaylorMaster(int WhichCase, int WhichG, int WhichH, cdouble GParam,
                     double *V1, double *V2, double *V3,
                     double *V2P, double *V3P, double *Q, double *QP);

/*--------------------------------------------------------------*/
/*- AssessPanelPair counts common vertices in a pair of panels, */
/*- and puts arrays of panel vertices into certain orders that  */
/*- are expected by subsequent algorithms that work on the      */
/*- panel pairs.                                                */
/*--------------------------------------------------------------*/
int AssessPanelPair(double **Va, double **Vb, double rMax);
int AssessPanelPair(double **Va, double **Vb);

int AssessPanelPair(RWGObject *Oa, int npa, 
                    RWGObject *Ob, int npb,
                    double *rRel, 
                    double **Va, double **Vb);

int AssessPanelPair(RWGObject *Oa, int npa, 
                    RWGObject *Ob, int npb,
                    double *rRel);

int NumCommonVertices(RWGObject *Oa, int npa, 
                      RWGObject *Ob, int npb);

int CanonicallyOrderVertices(double **Va, double **Vb, int ncv,
                             double **OVa, double **OVb);

#endif //LIBSCUFFINTERNALS_H
