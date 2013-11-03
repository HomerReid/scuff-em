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
 * scuff-static.cc   -- a standalone code within the scuff-EM suite for 
 *                   -- solving electrostatics problems 
 *
 * homer reid        -- 10/2006 -- 5/2013
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <libscuff.h>
#include <SSSolver.h>
#include <cmatheval.h>

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXEPF   10    // max number of evaluation-point files
#define MAXCACHE 10    // max number of cache files for preload

#define MAXSTR   1000

/***************************************************************/
/* in this case *UserData is a integer: 0, 1, 2 for x,y,z ******/
/***************************************************************/
void ConstantField(double *x, void *UserData, double PhiE[4])
{ 
  int WhichDirection = *(int *)UserData;
  memset(PhiE,0,4*sizeof(double));
  PhiE[0] = -1.0*x[WhichDirection];
  PhiE[1+WhichDirection] = 1.0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct USFData 
 {
   void *PhiEvaluator;
   void *EEvaluator[3];
 } USFData;

// 2.0e-20
//#define DELTA 9.5367431640625e-7
#define DELTA 1.0e-6
void UserStaticField(double *x, void *UserData, double PhiE[4])
{
  USFData *SFD = (USFData *)UserData;

  static const char *VariableNames[3] = { "x", "y", "z" };
  cdouble VariableValues[3];
  VariableValues[0] = cdouble(x[0], 0.0 );
  VariableValues[1] = cdouble(x[1], 0.0 );
  VariableValues[2] = cdouble(x[2], 0.0 );

  memset(PhiE, 0, 4*sizeof(double));
  if ( SFD==0 || SFD->PhiEvaluator==0 ) return;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  PhiE[0] = real( cevaluator_evaluate(SFD->PhiEvaluator, 3, 
                                      const_cast<char **>(VariableNames), 
                                      VariableValues) );

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int Mu=0; Mu<3; Mu++)
   { if (SFD->EEvaluator[Mu])
      { PhiE[1+Mu] = real( cevaluator_evaluate(SFD->EEvaluator[Mu], 3, 
                                               const_cast<char **>(VariableNames), 
                                               VariableValues) 
                         );
      }
     else
      { VariableValues[Mu] *= (1.0 + DELTA);
        double Temp = real( cevaluator_evaluate(SFD->PhiEvaluator, 3,
                                                const_cast<char **>(VariableNames), 
                                                VariableValues) 
                          );
        VariableValues[Mu] /= (1.0 + DELTA);
        PhiE[1+Mu] = -(Temp - PhiE[0]) / DELTA; 
      };
   };
   
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PhiEConstant(double *x, void *UserData, double PhiE[4])
{ 
  memset(PhiE, 0, 4.0*sizeof(double));
  int Direction = *((int *)UserData);
  PhiE[0]             = -x[Direction];
  PhiE[1 + Direction] = 1.0;
 
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPolarizabilities(SSSolver *SSS, HMatrix *M, HVector *Sigma, char *FileName)
{
  RWGGeometry *G = SSS->G;
  int NS = G->NumSurfaces;

  HMatrix *PolMatrix = new HMatrix(NS, 3);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *QP        = new HMatrix(NS, 4);
  for(int Mu=0; Mu<3; Mu++)
   { 
     SSS->AssembleRHSVector(0, PhiEConstant, (void *)(&Mu), Sigma);
     M->LUSolve(Sigma);
     SSS->GetCartesianMoments(Sigma, QP);
     for(int ns=0; ns<NS; ns++)
      PolMatrix->SetEntry(ns, Mu, QP->GetEntryD(ns,Mu+1));
   };
  delete QP;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=fopen(FileName,"w");
  for(int ns=0; ns<NS; ns++)
   fprintf(f,"%s %e %e %e \n",G->Surfaces[ns]->Label,
              PolMatrix->GetEntryD(ns,0),
              PolMatrix->GetEntryD(ns,1),
              PolMatrix->GetEntryD(ns,2));

  delete PolMatrix;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetCapacitanceMatrix(SSSolver *SSS, HMatrix *M, HVector *Sigma, HMatrix *C)
{

  RWGGeometry *G = SSS->G;
  int NS = G->NumSurfaces;
  double *Potentials = new double[NS];

  /*--------------------------------------------------------------*/
  /*- (re)allocate matrix as necessary ---------------------------*/
  /*--------------------------------------------------------------*/
  if (C)
   { if ( C->NR!=NS || C->NC!=NS )
      Warn("%s:%i: incorrect matrix passed to GetCapacitanceMatrix (reallocating)...",__FILE__,__LINE__);
      delete C;
      C=0;
   };
  if (!C)
   C=new HMatrix(NS, NS);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *QP = new HMatrix(NS, 4);
  for(int ns=0; ns<NS; ns++)
   { 
     memset(Potentials, 0, NS*sizeof(double));
     Potentials[ns]=1.0;
     SSS->AssembleRHSVector(Potentials, 0, 0, Sigma);
     M->LUSolve(Sigma);
     SSS->GetCartesianMoments(Sigma, QP);
     for(int nsp=0; nsp<NS; nsp++)
      C->SetEntry(nsp, ns, QP->GetEntry(nsp,0));
   };
  delete QP;

  return C;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetCapacitanceMatrices(SSSolver *SSS, HMatrix *M, HVector *Sigma, char *FileName)
{ 
  HMatrix *QP = GetCapacitanceMatrix(SSS, M, Sigma, 0);
  QP->ExportToText(FileName);
  delete QP;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/  
  /* process options *********************************************/
  /***************************************************************/
  InstallHRSignalHandler();
  char *GeoFile = 0;
  char *PolFile = 0;
  char *CapFile = 0;
  char *PotFile = 0;
  char *PhiExt  = 0;
  char *EPFiles[MAXEPF];             int nEPFiles;
  char *PlotFile = 0;
  char *Cache=0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache=0;
  char *ConstField=0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
/**/
     {"polfile",        PA_STRING,  1, 1,       (void *)&PolFile,    0,             "polarizability output file"},
/**/
     {"capfile",        PA_STRING,  1, 1,       (void *)&CapFile,    0,             "capacitance output file"},
/**/
     {"PotentialFile",  PA_STRING,  1, 1,       (void *)&PotFile,    0,             "list of conductor potentials"},
/**/
     {"PhiExternal",    PA_STRING,  1, 1,       (void *)&PhiExt,     0,             "functional form of external potential"},
     {"ConstField",     PA_STRING,  1, 1,       (void *)&ConstField, 0,             "direction of constant unit-strength E field (x,y,z) "},
/**/
     {"EPFile",         PA_STRING,  1, MAXEPF,  (void *)EPFiles,     &nEPFiles,     "list of evaluation points"},
     {"PlotFile",       PA_STRING,  1, 1,       (void *)&PlotFile,   0,             "surface charge visualization output file"},
/**/
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,      0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,   &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache, 0,             "write cache"},
/**/
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");

  /*******************************************************************/
  /* sanity check on input arguments *********************************/
  /*******************************************************************/
  if ( (PolFile || CapFile) && (nEPFiles>0 || PotFile!=0 || PhiExt!=0 ) )
   ErrExit("(--EPFile,--PotFile,--PhiExternal) may not be used with (--polfile, --capfile)");
  if ( nEPFiles==0 && !PolFile && !CapFile && !PlotFile )
   OSUsage(argv[0], OSArray, "you have not selected any type of calculation");

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  int ConstFieldDirection=-1;
  if (ConstField)
   { switch(tolower(ConstField[0]))
      { case 'x': ConstFieldDirection=0; break;
        case 'y': ConstFieldDirection=1; break;
        case 'z': ConstFieldDirection=2; break;
        default:  ErrExit("invalid --ConstField specification");
      };
   };

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  SetLogFileName("scuff-static.log");
  Log("scuff-static running on %s",GetHostName());

  /*******************************************************************/
  /* create the ScuffStaticGeometry **********************************/
  /*******************************************************************/
  SSSolver *SSS   = new SSSolver(GeoFile);
  RWGGeometry *G  = SSS->G;
  HMatrix *M      = SSS->AllocateBEMMatrix();
  HVector *Sigma  = SSS->AllocateRHSVector();

  /*******************************************************************/
  /* preload the scuff cache with any cache preload files the user   */
  /* may have specified                                              */
  /*******************************************************************/
  if ( Cache!=0 && WriteCache!=0 )
   ErrExit("--cache and --writecache options are mutually exclusive");
  if (Cache) 
   WriteCache=Cache;
  for (int nrc=0; nrc<nReadCache; nrc++)
   PreloadCache( ReadCache[nrc] );
  if (Cache)
   PreloadCache( Cache );

  /*******************************************************************/
  /* assemble and factorize the BEM matrix                           */
  /*******************************************************************/
  SSS->AssembleBEMMatrix(M);
  M->LUFactorize();

  /*******************************************************************/
  /* now switch off depending on what the user requested             */
  /*******************************************************************/
  if (PolFile)
   GetPolarizabilities(SSS, M, Sigma, PolFile);
  if (CapFile)
   GetCapacitanceMatrices(SSS, M, Sigma, CapFile);
  if ( nEPFiles>0 || PlotFile )
   { 
     /***************************************************************/
     /* process user's conductor potential file if present          */
     /***************************************************************/
     double *Potentials = new double[G->NumSurfaces];
     if (PotFile)
      { 
        ErrExit("--potfile option not yet supported");
      }
     else
      memset( Potentials, 0.0, G->NumSurfaces * sizeof(double) );

     /***************************************************************/
     /* process user's external-field specification if present      */
     /***************************************************************/
     if (PhiExt)
      { 
        ErrExit("--phiexternal option not yet supported");
        /*
        USFData MyData, *D = &MyData;
        D->PhiEvaluator = cevaluator_create(PhiExt);
        D->EEvaluator[0] = D->EEvaluator[1] = D->EEvaluator[2] = 0;
        SSS->AssembleRHSVector(Potentials, UserStaticField, (void *)D, Sigma);
        */
      }
     else if ( ConstFieldDirection!=1 )
      SSS->AssembleRHSVector(Potentials, ConstantField, &ConstFieldDirection, Sigma);
     else
      SSS->AssembleRHSVector(Potentials, 0, 0, Sigma);

     SSS->PlotChargeDensity(Sigma, "%s.RHS", PlotFile);

     /***************************************************************/
     /* solve the problem *******************************************/
     /***************************************************************/
     M->LUSolve(Sigma);

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/ 
     if (PlotFile)
      SSS->PlotChargeDensity(Sigma, PlotFile, 0);

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     for(int nepf=0; nepf<nEPFiles; nepf++)
      {
         HMatrix *X = new HMatrix(EPFiles[nepf]);
         if (X->ErrMsg)
          ErrExit(X->ErrMsg);

         HMatrix *PhiE;
         if (PhiExt)
          ErrExit("%s:%i: internal error ",__FILE__,__LINE__);
         if ( ConstFieldDirection!=1 )
          PhiE = SSS->GetFields(ConstantField, &ConstFieldDirection, Sigma, X, 0);
         else
          PhiE = SSS->GetFields(0, 0, Sigma, X, 0);

         FILE *f=vfopen("%s.out","w",GetFileBase(EPFiles[nepf]));
         for(int nr=0; nr<X->NR; nr++)
          fprintf(f,"%e %e %e %e %e %e %e\n",
                     X->GetEntryD(nr,0), X->GetEntryD(nr,1), X->GetEntryD(nr,2),
                     PhiE->GetEntryD(nr,0), PhiE->GetEntryD(nr,1), 
                     PhiE->GetEntryD(nr,2), PhiE->GetEntryD(nr,3));
         fclose(f);
      };

   };

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  printf("Thank you for your support.\n");

}
