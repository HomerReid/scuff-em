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
 * scuff-integrate -- simple utility for numerically integrating functions
 *                 -- from samples tabulated in data files, possibly with
 *                 -- additional integrand factors inserted
 *
 * homer reid  -- 5/2012
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include <string.h> 
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fenv.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <libhrutil.h>
#include <libSGJC.h>
#include <libMDInterp.h>
#include <libscuff.h>

#define II cdouble(0.0,1.0)

using namespace scuff;

#define MAXDATA     21
#define MAXPARMS    20
#define MAXOBJ      9

typedef std::vector<int>     ivec;
typedef std::vector<char>    cvec;
typedef std::vector<double>  dvec;
typedef std::vector<cdouble> zvec;

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define HBAROMEGA02 9.491145534e-06
#define BOLTZMANNK 4.36763e-4
double GetThetaFactor(double Omega, double T)
{ 
  if (T==0.0)
   return 0.0;

  double ExpArg = Omega/(BOLTZMANNK*T);
  if (ExpArg>100.0)
   return 0.0;

  return Omega / ( exp(ExpArg) - 1.0 );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct IntegrandData
{
  Interp1D *I1D;
  double Temperature;
  double TEnvironment;
  double PreFactor;
  FILE *LogFile;

} IntegrandData;

int Integrand(unsigned ndim, const double *x, void *UserData,
              unsigned fdim, double *fval)
{
  (void )ndim; // unused

  IntegrandData *Data = (IntegrandData *)UserData;
  Interp1D *I1D       = Data->I1D;
  double Temperature  = Data->Temperature;
  double TEnvironment = Data->TEnvironment;
  double PreFactor    = Data->PreFactor;
  FILE *LogFile       = Data->LogFile;

  double Omega=x[0];
  I1D->Evaluate(Omega, fval);

  for(unsigned nf=0; nf<fdim; nf++)
   fval[nf]*=PreFactor;

  if (Temperature>=0.0)
   { 
     double DeltaTheta =
      GetThetaFactor(Omega, Temperature)
       -GetThetaFactor(Omega, TEnvironment);

     for(unsigned nf=0; nf<fdim; nf++)
      fval[nf]*=DeltaTheta;
   };

  if (LogFile)
   { fprintf(LogFile,"%e ",Omega);
     for(unsigned nf=0; nf<fdim; nf++)
      fprintf(LogFile,"%e ",fval[nf]);
     fprintf(LogFile,"\n");
   };

  return 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool Equal(double X, double Y, double RelTol=1.0e-7, double AbsTol=0.0)
{
  if ( fabs(X)<=AbsTol && fabs(Y)<=AbsTol )
   return true;
  double Delta = fabs(X-Y);
  if ( Delta<=AbsTol || Delta<=RelTol*fmax(fabs(X),fabs(Y)) )
   return true;
  return false;
}

/***************************************************************/
/* return true if the values of all parameter columns in       */
/* row nrA equal their counterparts in nrB                     */
/***************************************************************/
bool ParmsMatch(HMatrix *M, int nrA, int nrB,
                int *ParmColumns,int NumParms,
                double RelTol=1.0e-4, double AbsTol=1.0e-6)
{
  for(int np=0; np<NumParms; np++)
   { int nc = ParmColumns[np]-1;
     cdouble zA = M->GetEntry(nrA, nc), zB = M->GetEntry(nrB, nc);
     double dz = abs(zA-zB), Scale=fmin(abs(zA),abs(zB));
     if ( (dz>AbsTol) && (dz>Scale*RelTol) )
      return false;
   };
  return true;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
#ifndef __APPLE__
  feenableexcept(FE_INVALID | FE_OVERFLOW);
#endif

  /***************************************************************/
  /* process command-line arguments.                             */
  /* Note: All column-specification arguments are interpreted as */
  /* ONES-BASED, so leftmost column in a data file is column #1  */
  /* and there is no column 0.                                   */
  /***************************************************************/
  char *DataFileName=0;
//
  int FreqColumn=1;                 int nFreqColumn=0;
//
  int DataColumns[MAXDATA];         int NumData=0;
  char *DataNames[MAXDATA];         int nDataNames=0;
//
  int ParmColumns[MAXPARMS];        int NumParms=0;
  char *ParmNames[MAXPARMS];        int nParmNames=0;
//
  int SDColumn=0;
//
  char *TemperatureFile=0;
  char *TempStrings[2*MAXOBJ];      int nTempStrings;
//
  double PreFactor=0.0;
//
  double AbsTol=0.0;
  double RelTol=1.0e-4;
  double DupAbsTol=0.0;
  double DupRelTol=1.0e-6;
  double FreqRelTol=1.0e-6;
//
  char *IntegrandFileName=0;
  char *OutFileName=0;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"DataFile",        PA_STRING,  1, 1,       (void *)&DataFileName,       0,   "data file"},
//
     {"FreqColumn",      PA_INT,     1, 1,       (void *)&FreqColumn,         &nFreqColumn,   "frequency column index in data file"},
     {"xColumn",         PA_INT,     1, 1,       (void *)&FreqColumn,         0,   "(synonym for --FreqColumn)"},
//
     {"DataColumn",      PA_INT,     1, MAXDATA, (void *)DataColumns,   &NumData,        "data column index in data file"},
     {"DataName",        PA_STRING,  1, MAXDATA, (void *)DataNames,     &nDataNames,     "name for data column"},
//
     {"ParmColumn",      PA_INT,     1, MAXPARMS,(void *)ParmColumns,   &NumParms,       "parameter column index in data file"},
     {"ParmName"  ,      PA_STRING,  1, MAXPARMS,(void *)ParmNames,    &nParmNames,     "name of parameter"},
//
     {"SDColumn",        PA_INT,     1, 1,       (void *)&SDColumn,           0,   "(source object, dest object) column index in data file"},
//
     {"TemperatureFile", PA_STRING,  1, 1,       (void *)&TemperatureFile,    0,   "list of object temperatures"},
     {"Temperature",     PA_STRING,  2, MAXOBJ,  (void *)TempStrings,         &nTempStrings,   "set object temperature"},
//
     {"PreFactor",       PA_DOUBLE,  1, 1,       (void *)&PreFactor,          0,   "overall multiplicative prefactor"},
     {"IntegrandFile",   PA_STRING,  1, 1,       (void *)&IntegrandFileName,  0,   "name of file for integrand data"},
     {"OutFile",         PA_STRING,  1, 1,       (void *)&OutFileName,        0,   "name of output file"},
     {"DupRelTol",       PA_DOUBLE,  1, 1,       (void *)&DupRelTol,          0,   "relative tolerance for identifying equivalent data points"},
     {"DupAbsTol",       PA_DOUBLE,  1, 1,       (void *)&DupAbsTol,          0,   "absolute tolerance for identifying equivalent data points"},
     {"FreqRelTol",      PA_DOUBLE,  1, 1,       (void *)&FreqRelTol,         0,   "relative tolerance for identifying equivalent frequency points"},
//
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (DataFileName==0)
   OSUsage(argv[0], VERSION,OSArray,"--DataFile option is mandatory");

  /***************************************************************/
  /* auto-detect special known file types and autoset values of  */
  /* input parameters                                            */
  /***************************************************************/
  if( strcasestr(DataFileName,".SIFlux") )
   { 
     Log("Autodetecting .SIFlux data file...");

     ParmColumns[0] = 1;    ParmNames[0]=strdup("Geometrical transform");
     FreqColumn     = 2;
     // FUTURE: handle kBloch 
     int nc = 3;
     ParmColumns[1] = nc++; ParmNames[1]=strdup("(source, dest)");
     SDColumn       = ParmColumns[1];
     NumParms=nParmNames=2;

     DataColumns[0] = nc++; DataNames[0]=strdup("PAbs");
     DataColumns[1] = nc++; DataNames[1]=strdup("PRad");
     DataColumns[2] = nc++; DataNames[2]=strdup("XForce");
     DataColumns[3] = nc++; DataNames[3]=strdup("YForce");
     DataColumns[4] = nc++; DataNames[4]=strdup("ZForce");
     DataColumns[5] = nc++; DataNames[5]=strdup("XTorque");
     DataColumns[6] = nc++; DataNames[6]=strdup("YTorque");
     DataColumns[7] = nc++; DataNames[7]=strdup("ZTorque");
     NumData=nDataNames=8;
   }
  else if( strcasestr(DataFileName,".SRFlux") )
   { 
     Log("Autodetecting .SRFlux data file");
     ParmColumns[0] = 1;    ParmNames[0]=strdup("Geometrical transform");
     FreqColumn     = 2;
     // FUTURE: handle kBloch 
     int nc = 3;
     ParmColumns[1] = nc++; ParmNames[1]=strdup("x");
     ParmColumns[2] = nc++; ParmNames[2]=strdup("y");
     ParmColumns[3] = nc++; ParmNames[3]=strdup("z");
     ParmColumns[4] = nc++; ParmNames[4]=strdup("sourceSurface");
     SDColumn       = ParmColumns[4];
     NumParms=nParmNames=5;
     
     int ndc=0;
     char SC[3]="Sx", TC[4]="Txx";
     for(SC[1]='x'; SC[1]<='z'; SC[1]++)
      { DataNames[ndc]     = strdup(SC);
        DataColumns[ndc++] = nc++; 
      };
     for(TC[1]='x'; TC[1]<='z'; TC[1]++)
      for(TC[2]='x'; TC[2]<='z'; TC[2]++)
       { DataNames[ndc] = strdup(TC);
         DataColumns[ndc++]  = nc++; 
       };
     NumData=nDataNames=nc;
   }
  else if( strcasestr(DataFileName,".byXi") )
   {
     Log("Autodetecting .byXi data file");
     ErrExit("--byXi mode not yet implemented");
   };

  /***************************************************************/
  /* sanity checks ***********************************************/
  /***************************************************************/
  if (NumData==0)
   OSUsage(argv[0], VERSION, OSArray,"you must specify at least one --DataColumn");
  if (nDataNames!=0 && nDataNames!=NumData)
   OSUsage(argv[0], VERSION, OSArray,"numbers of --DataNames and --DataColumns do not agree");
  if (nParmNames!=0 && nParmNames!=NumParms)
   OSUsage(argv[0], VERSION, OSArray,"numbers of --ParmNames and --ParmColumns do not agree");

  if (nDataNames==0)
   for(int nd=0; nd<NumData; nd++)
    DataNames[nd]=vstrdup("Integral of data #%i\n",nd);
  if (nParmNames==0)
   for(int np=0; np<NumParms; np++)
    ParmNames[np]=vstrdup("Parameter #%i\n",np);
  

  /***************************************************************/
  /* read in the data and sanity check column specifications     */
  /***************************************************************/
  HMatrix *DataMatrix=new HMatrix(DataFileName);
  int NR = DataMatrix->NR;
  if ( ! (1<=FreqColumn && FreqColumn<=NR) )
   ErrExit("--FreqColumn not in range [1,%i]",NR);
  for(int np=0; np<NumParms; np++)
   if ( ! (1<=ParmColumns[np] && ParmColumns[np]<=NR) )
    ErrExit("--ParmColumn %i not in range [1,%i]",ParmColumns[np],NR);
  for(int nd=0; nd<NumData; nd++)
   if ( ! (1<=DataColumns[nd] && DataColumns[nd]<=NR) )
    ErrExit("--DataColumn %i not in range [1,%i]",DataColumns[nd],NR);

  /***************************************************************/
  /* sort by parameter values, then by frequency                 */
  /***************************************************************/
  ivec SortColumns(NumParms+1);
  for(int n=0; n<NumParms; n++)
   SortColumns[n]=ParmColumns[n]-1;
  SortColumns[NumParms]=FreqColumn-1; 
  DataMatrix->Sort(SortColumns);

  /***************************************************************/
  /* look at the source/dest index (if present) to figure out how*/
  /* many objects are in the geometry                            */
  /***************************************************************/
  int NSource=1, NDest=0;
  if (SDColumn)
   for(int nr=0; nr<DataMatrix->NR; nr++)
    { int SD = (int)DataMatrix->GetEntryD(nr, SDColumn-1);
      if (SD<11 || SD>99)
       ErrExit("%s:%i: invalid SDColumn %i",DataFileName,nr,SD);
      int SourceIndex = SD/10, DestIndex=SD%10;
      if (NSource < SourceIndex) NSource=SourceIndex;
      if (NDest   < DestIndex  ) NDest=DestIndex;
    };
  int MinDestIndex = (NDest==0 ? 0 : 1);
  if (NDest==0) NDest=1;

  /***************************************************************/
  /* process temperature specifications if any.                  */
  /* TemperatureMatrix column  0     = environment temperatures  */
  /* TemperatureMatrix columns 1..NO = object 1..NO temperatures */
  /***************************************************************/
  HMatrix *TemperatureMatrix=0;
  if ( TemperatureFile )
   { if (nTempStrings>0)
      ErrExit("--Temperature is incompatible with --TemperatureFile");
     TemperatureMatrix=new HMatrix(TemperatureFile);
     if ( TemperatureMatrix->NC != (NSource+1) )
      ErrExit("Data file %s describes %i objects, but temperature file %s describes environment plus %i objects",
               DataFileName, NSource, DataFileName, TemperatureMatrix->NC-1);
   }
  else if ( nTempStrings>0 )
   {
     TemperatureMatrix=new HMatrix(1,NSource+1);

     for(int nts=0; nts<nTempStrings; nts++)
      { int Index; 
        int Status=sscanf(TempStrings[2*nts+0],"%i",&Index);
        if (Status!=1 || Index<0 || Index>NSource)
         ErrExit("invalid object index %s",TempStrings[2*nts+0]);
        double T;
        Status=sscanf(TempStrings[2*nts+1],"%le",&T);
        if (Status!=1) 
         ErrExit("invalid temperature %s",TempStrings[2*nts+1]);
        TemperatureMatrix->SetEntry(0,Index,T);
      };
     printf("Temperatures: \n");
     printf(" Environment: %e\n",TemperatureMatrix->GetEntryD(0,0));
     for(int ns=0; ns<NSource; ns++)
      printf(" Object %2i  : %e\n",ns+1,TemperatureMatrix->GetEntryD(0,ns+1));
   };
  int NumTemperatures = TemperatureMatrix ? TemperatureMatrix->NR : 1;

  if (PreFactor==0.0)
   PreFactor = (TemperatureMatrix==0) ? 1.0 : HBAROMEGA02;

  /***************************************************************/
  /* NumParmValues is the total number of distinct sets of       */
  /* parameter values.                                           */
  /* The (x,y1,y2,...) data for e.g. the 3rd set of parameter    */
  /* values start on row # RowStarts[2]                          */
  /*      and end on row # RowStarts[3] - 1                      */
  /***************************************************************/
  int NumParmValues=1;
  ivec RowStarts(1,0);
  if (NumParms)
   for(int nr=1; nr<DataMatrix->NR; nr++)
    if ( !ParmsMatch(DataMatrix,nr,nr-1,ParmColumns,NumParms) )
     RowStarts[NumParmValues++]=nr;
  RowStarts.push_back(DataMatrix->NR);
  
  int MaxNumFreqs = RowStarts[1] - RowStarts[0];
  for(int npv=2; npv<NumParmValues; npv++)
   { int ThisNumFreqs=RowStarts[npv]-RowStarts[npv-1];
     if (ThisNumFreqs>MaxNumFreqs)
      MaxNumFreqs=ThisNumFreqs;
   };

  /***************************************************************/
  /* create output file and write preamble ***********************/
  /***************************************************************/
  if (!OutFileName)
   OutFileName=vstrdup("%s.Integrated",GetFileBase(DataFileName));
  FILE *OutFile=fopen(OutFileName,"a");
  if (!OutFile)
   ErrExit("could not open output file %s",OutFileName);
  fprintf(OutFile,"# columns: \n");
  int ncol=1;
  for(int np=0; np<NumParms; np++)
   fprintf(OutFile,"# %i   %s\n",ncol++,ParmNames[np]);
  if (TemperatureMatrix)
   { fprintf(OutFile,"# %i: environment temperature \n",ncol++);
     for(int ns=0; ns<NSource; ns++)
      fprintf(OutFile,"# %i: object %i temperature \n",ncol++,ns+1);
   };
  for(int nd=0; nd<NumData; nd++)
   fprintf(OutFile,"# %i   integral of %s\n",ncol++,DataNames[nd]);

  /***************************************************************/
  /* allocate storage buffers                                    */
  /***************************************************************/
  double *XValues = (double *)mallocEC(MaxNumFreqs*sizeof(double));
  double *YValues = (double *)mallocEC(MaxNumFreqs*NumData*sizeof(double));

  double *Result = (double *)mallocEC(NumData*sizeof(double));
  double *Error  = (double *)mallocEC(NumData*sizeof(double));

  /* array to accumulate contributions of all source objects to */
  /* quantities for individual destination objects              */
  double *TotalByDest=0;
  if (NSource>1)
   TotalByDest=(double *)mallocEC(NDest*NumTemperatures*NumData*sizeof(double));
     
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int npv=0; npv<NumParmValues; npv++)
   {
     int NR0=RowStarts[npv];

     int SourceObject=-1, DestObject=-1;
     if (NSource>1)
      { int SD = (int)DataMatrix->GetEntryD(NR0, SDColumn-1);
        SourceObject = SD/10;
        DestObject   = SD%10;
        if (SourceObject==0 && DestObject==MinDestIndex)
         memset(TotalByDest,0,NDest*NumTemperatures*NumData*sizeof(double));
      };

     // populate arrays of x,y points for this set of parameter vals
     int NSamples=0;
     for(int nr=RowStarts[npv]; nr<RowStarts[npv+1]; nr++)
      { 
        double ThisFreq=DataMatrix->GetEntryD(nr, FreqColumn-1);
        if (NSamples>1 && Equal(ThisFreq, XValues[NSamples-1], FreqRelTol))
         continue;
        XValues[NSamples]=ThisFreq;
        for(int nd=0; nd<NumData; nd++)
         YValues[NSamples*NumData + nd]=DataMatrix->GetEntryD(nr, DataColumns[nd]-1);
        NSamples++;
      };
     double MinFreq = XValues[0], MaxFreq = XValues[NSamples-1];
     if (NumParms>0)
      { Log("Integrating for parm values { ");
        for(int np=0; np<NumParms; np++)
         LogC("%g ",DataMatrix->GetEntryD(NR0,ParmColumns[np]-1));
        LogC("}");
      };
     Log("  %i data points in range (%e,%e)\n",NSamples,MinFreq,MaxFreq);
       
     /***************************************************************/
     /* initialize interpolator                                     */
     /***************************************************************/
     IntegrandData MyID;
     MyID.I1D=new Interp1D(XValues, YValues, NSamples, NumData);
     MyID.PreFactor    = PreFactor;
     MyID.LogFile      = 0;
     MyID.Temperature  = -1.0;
     MyID.TEnvironment = -1.0;
  
     /***************************************************************/
     /* evaluate integrals at all temperature combinations          */
     /***************************************************************/
     for(int nTemp=0; nTemp<NumTemperatures; nTemp++)
      { 
        if (TemperatureMatrix)
         { MyID.TEnvironment = TemperatureMatrix->GetEntryD(nTemp,0);
           MyID.Temperature  = TemperatureMatrix->GetEntryD(nTemp,SourceObject+1);
         };
           
        MyID.LogFile=0;
#if 0
        if (IntegrandFileName)
         { MyID.LogFile=fopen(IntegrandFileName,"a");
           fprintf(MyID.LogFile,"\n\n");
           for(int np=0; np<
             if (TemperatureMatrix)
               { fprintf(MyID.LogFile,"# environment / object temperatures: ");
                 fprintf(MyID.LogFile,"%e ",MyID.TEnvironment);
                 for(int no=0; no<NumObjects; no++)
                  fprintf(MyID.LogFile,"%e ",TemperatureMatrix->GetEntryD(nTemp,no+1));
                 fprintf(MyID.LogFile,"\n");
               };
            };
#endif

        /*--------------------------------------------------------------*/
        /*- evaluate the integral --------------------------------------*/
        /*--------------------------------------------------------------*/
        pcubature(NumData, Integrand, (void *)&MyID, 1,
 	          &MinFreq, &MaxFreq, 0, AbsTol, RelTol,
 	          ERROR_INDIVIDUAL, Result, Error);

        if (MyID.LogFile)
         { fclose(MyID.LogFile);
           MyID.LogFile=0;
         };

        /*--------------------------------------------------------------*/
        /*- write results to .Integrated file                          -*/
        /*--------------------------------------------------------------*/
        int EffectiveNumParms = NumParms - (SDColumn ? 1 : 0);
        for(int np=0; np<EffectiveNumParms; np++)
         fprintf(OutFile,"%e ",DataMatrix->GetEntryD(NR0,ParmColumns[np]-1));
        if (SDColumn) 
         fprintf(OutFile,"%i%i",SourceObject,DestObject);
        if (TemperatureMatrix)
         for(int nc=0; nc<TemperatureMatrix->NC; nc++)
          fprintf(OutFile,"%e ",TemperatureMatrix->GetEntryD(nTemp,nc));
        for(int nd=0; nd<NumData; nd++)
         fprintf(OutFile,"%e ",Result[nd]);
        fprintf(OutFile,"\n");
        fflush(OutFile);
  
        /*--------------------------------------------------------------*/
        /*- accumulate per-destination-object totals and print when we -*/
        /*- have added the contributions of all source objects         -*/
        /*--------------------------------------------------------------*/
        if (NSource>1)
         { 
           int EffDestObject = (DestObject==0 ? 0 : DestObject-1);
           int Offset=EffDestObject*NumTemperatures*NumData + nTemp*NumData;
           for(int nd=0; nd<NumData; nd++)
            TotalByDest[Offset + nd] += Result[nd];

           if (SourceObject==NSource-1)
            { for(int np=0; np<NumParms-1; np++)
               fprintf(OutFile,"%e ",DataMatrix->GetEntryD(NR0,ParmColumns[np]));
              fprintf(OutFile,"0%i",DestObject);
              if (TemperatureMatrix)
               for(int nc=0; nc<TemperatureMatrix->NC; nc++)
                fprintf(OutFile,"%e ",TemperatureMatrix->GetEntryD(nTemp,nc));
              for(int nd=0; nd<NumData; nd++)
               fprintf(OutFile,"%e ",TotalByDest[Offset+nd]);
              fprintf(OutFile,"\n");
              fflush(OutFile);
            }

         } // if (NSource>1)

      } // for(int nTemp=0; nTemp<NumTemperatures; nTemp++)

     delete MyID.I1D;

   } // for(int npv=0; npv<NumParmValues; npv++)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  fclose(OutFile);
  printf("Integrated data written to %s.\n",OutFileName);
  printf("Thank you for your support.\n");
 
}
