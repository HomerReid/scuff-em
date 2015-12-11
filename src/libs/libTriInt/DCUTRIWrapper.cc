/*****************************************************/
/* libDCUTRI.cc -- wrapper around netlib DCUTRI      */
/*              -- function for adaptive cubature    */
/*              -- over triangles                    */
/* homer reid   -- 7/2007                            */
/*****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libTriInt.h>

#ifndef HAVE_DCUTRI

void *CreateDCUTRIWorkspace(int numfun, int maxpts)
{ 
  (void) numfun;
  (void) maxpts;

  ErrExit("SCUFF-EM was not compiled with DCUTRI; configure --with-dcutri");
  return 0;
}

int DCUTRI(void *opW, double **Vertices, TriIntFun func,
           void *parms, double epsabs, double epsrel,
           double *I, double *abserr)
{ 
  (void)opW; 
  (void)Vertices;
  (void)func;
  (void)parms;
  (void)epsabs;
  (void)epsrel;
  (void)I;
  (void)abserr;

  CreateDCUTRIWorkspace(0,0);
  return 0;
}
#else

/***************************************************************/
/* data structure used to store data for calls to DCUTRI       */
/***************************************************************/
typedef struct DCUTRIWorkspace
 { double V0[3], V1mV0[3], V2mV0[3]; /* vertices  */
   TriIntFun func;                   /* user's function */
   void *parms;                      /* user's data */
   int numfun;                       /* number of components of integrand */
   int maxpts;                       /* max allowed number of func evals */
   int neval;                        /* actual number of func evals */
   double *ver, *work;               /* dcutri internal workspace */
   int nw, lenver, *iwork;
 } DCUTRIWorkspace;

/****************************************************************/
/* prototype for fortran dcutri routine *************************/
/****************************************************************/
typedef void (*funsub_t)(double *, int *,  double *, double *);
extern "C" {
void dcutri_(funsub_t funsub, double *parms, int *numfun, double *ver,
int *numtri, int *minpts, int *maxpts, double *epsabs,
double *epsrel, int *lenver, int *nw, int *restar,
double *result, double *abserr, int *neval, int *ifail,
double *work, int *iwork);
}

/***************************************************************/
/* funsub routine passed to dcutri_  ***************************/
/***************************************************************/
void funsub(double *x, int *nfun, double *parms, double *f)
{ 
  double X[3];
  int i;
  DCUTRIWorkspace *W=(DCUTRIWorkspace *)parms;

  for(i=0; i<3; i++)
   X[i]=W->V0[i] + x[0]*W->V1mV0[i] + x[1]*W->V2mV0[i];

  W->func(X, W->parms, f);
   
} 

/**************************************************************/
/* allocate and return an opaque pointer to a new             */
/* DCUTRIWorkspace structure for passage to DCUTRI            */
/**************************************************************/
void *CreateDCUTRIWorkspace(int numfun, int maxpts)
{
  DCUTRIWorkspace *W;
  int numtri=1, lenver;

  W=(DCUTRIWorkspace *)malloc(sizeof(DCUTRIWorkspace));
  W->numfun=numfun;
  W->maxpts=maxpts;
  W->lenver=lenver=3*(maxpts-37)/(4*37) + numtri;
  W->ver=(double *)malloc(6*lenver*sizeof(double));
  W->nw=lenver*(2*numfun+1) + 32*numfun + 1;
  W->work=(double *)malloc( (W->nw)*sizeof(double));
  W->iwork=(int *)malloc( (lenver+1)*sizeof(int));

  return (void *)W;
  
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
int DCUTRI(void *opDCW, double **Vertices, TriIntFun func, 
           void *parms, double epsabs, double epsrel, 
           double *I, double *abserr)
{ 
  DCUTRIWorkspace *W=(DCUTRIWorkspace *)opDCW;
  int i, restar, numtri=1, minpts=1, ifail;
  double VC[3], Jacobian;

  /* fill in some fields in the DCW structure */
  memcpy(W->V0,Vertices[0],3*sizeof(double));
  for(i=0; i<3; i++)
   { W->V1mV0[i]=Vertices[1][i]-Vertices[0][i];
     W->V2mV0[i]=Vertices[2][i]-Vertices[0][i];
   };
  W->parms=parms;
  W->func=func;
  
  /* call dcutri_ to evaluate the integral over the triangle */
  restar=0;
  numtri=1;
  W->ver[0]=0.0; W->ver[1]=0.0;
  W->ver[2]=1.0; W->ver[3]=0.0;
  W->ver[4]=0.0; W->ver[5]=1.0;
  W->neval=0;

  dcutri_(funsub, (double *)W, &(W->numfun), W->ver, &numtri, &minpts, 
          &(W->maxpts), &epsabs, &epsrel, &(W->lenver), &(W->nw), 
          &restar, I, abserr, &(W->neval), &ifail, W->work, W->iwork);


  /* put in the jacobian factor. without this, the answer returned
     is the value of the integral/(2*Area). in other words, you need
     to multiply the return value by 2*Area to get the full integral.
  
     20091121 i think i am going to disable this by default.
     20100109 no, i am going to restore it.                  
   */
//#undef PUTINJACOBIAN
#define PUTINJACOBIAN
#ifdef PUTINJACOBIAN
  VC[0]=W->V1mV0[1]*W->V2mV0[2] - W->V1mV0[2]*W->V2mV0[1];
  VC[1]=W->V1mV0[2]*W->V2mV0[0] - W->V1mV0[0]*W->V2mV0[2];
  VC[2]=W->V1mV0[0]*W->V2mV0[1] - W->V1mV0[1]*W->V2mV0[0];
  Jacobian=sqrt(VC[0]*VC[0] + VC[1]*VC[1] + VC[2]*VC[2]);
  for(i=0; i<W->numfun; i++)
   { I[i]*=Jacobian;
     abserr[i]*=Jacobian;
   };
#endif
  
  return W->neval;
  
}  

#endif // #HAVE_DCUTRI
