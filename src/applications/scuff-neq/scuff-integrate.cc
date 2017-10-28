#include <stdio.h>
#include <stdlib.h>

#include <string.h> 
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <libhrutil.h>
#include <libSGJC.h>
#include <libMDInterp.h>
#include <libscuff.h>
#include <fenv.h>

#define II cdouble(0.0,1.0)

using namespace scuff;

#define MAXDATA     20
#define MAXTAGS     100
#define MAXOBJ      9

/***************************************************************/
/***************************************************************/
/***************************************************************/
char *GetTagString(cdouble Tag)
{ 
  static char TagString[30];
  if ( imag(Tag)==0.0 )
   sprintf(TagString,"%5.3f ",real(Tag));
  else if (imag(Tag)<0.0)
   sprintf(TagString,"%5.3f-%5.3f ",real(Tag),-imag(Tag));
  else
   sprintf(TagString,"%5.3f+%5.3f ",real(Tag),imag(Tag));
  return TagString;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct FileData
 { 
   HMatrix *DataMatrix;
   int NumObjects;
   int NumTags;
   cdouble Tags[MAXTAGS];
   char *TagStrings[MAXTAGS];
  
 } FileData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
FileData *GetFileData(char *DataFileName, int TagColumn, int SDColumn)
{
  HMatrix *DataMatrix=new HMatrix(DataFileName);

  FileData *FD  = (FileData *)mallocEC(sizeof(FileData));
  cdouble *Tags = FD->Tags;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumTags=0;
  if (TagColumn)
   for(int nr=0; nr<DataMatrix->NR; nr++)
    { cdouble Tag=DataMatrix->GetEntry(nr, TagColumn-1);
      bool HaveTag=false;
      for(int nTag=0; nTag<NumTags && !HaveTag; nTag++)
       if (Tags[nTag]==Tag)
        HaveTag=true;
      if (HaveTag)
       continue;
      if (NumTags==MAXTAGS-1)
       Warn("Too many tags (skipping %s)",z2s(Tag));
      else
       Tags[NumTags++]=Tag;
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumObjects=0;
  if (SDColumn)
   for(int nr=0; nr<DataMatrix->NR; nr++)
    { int SD = (int)DataMatrix->GetEntryD(nr, SDColumn-1);
      if (SD<0 || SD>99) 
       ErrExit("%s:%i: invalid SDColumn %i",DataFileName,nr,SD);
      if (NumObjects < (SD%10)) 
       NumObjects=SD%10;
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FD->NumObjects = NumObjects;
  FD->NumTags    = NumTags;
  FD->DataMatrix = DataMatrix;

  return FD;
}

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
bool Equal(double X, double Y, double RelTol, double AbsTol)
{
  if ( fabs(X)<=AbsTol && fabs(Y)<=AbsTol )
   return true;
  double Delta = fabs(X-Y);
  if ( Delta<=AbsTol || Delta<=RelTol*fmax(fabs(X),fabs(Y)) )
   return true;
  return false;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RemoveDuplicates(char *FileName,
                      double *XValues, double *YValues, int NumData, int *NumPoints,
                      double RelTol=1.0e-6, double AbsTol=1.0e-8)
{
  int N = *NumPoints;
  int NumDuplicates = 0;
  for(int m=0; m<N; m++)
   for(int n=m+1; n<N; n++)
    if ( Equal(XValues[m], XValues[n], RelTol, AbsTol) )
     { 
       for(int nd=0; nd<NumData; nd++)
        if ( !Equal(YValues[ m*NumData + nd], YValues[n*NumData+nd], RelTol, AbsTol) )
         ErrExit("%s: found incompatible {frequency,data[%i]} pairs {%e,%e}, {%e,%e}",
                  FileName,nd,XValues[m],YValues[m*NumData+nd],XValues[n],YValues[n*NumData+nd]);
       
       // remove column #nd from X and Y arrays
       NumDuplicates++;
       for(int p=n; p<(N-1); p++)
        { XValues[p]=XValues[p+1];
          for(int nd=0; nd<NumData; nd++)
           YValues[p*NumData+nd]=YValues[(p+1)*NumData+nd];
        };
       XValues[N-1]=0.0;
       for(int nd=0; nd<NumData; nd++)
        YValues[(N-1)*NumData+nd]=0.0;
     };
  *NumPoints -= NumDuplicates;
  if (NumDuplicates>0)
   Warn("removed %i duplicate frequency points",NumDuplicates);
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
  /* process command-line arguments ******************************/
  /***************************************************************/
  char *DataFileName=0;
  int FreqColumn=1;
  int DataColumns[MAXDATA];         int nDataColumns=0;
  char *DataNames[MAXDATA];         int nDataNames=0;
  int TagColumn=0;
  int SDColumn=0;
  char *TemperatureFile=0;
  char *TempStrings[2*MAXOBJ];      int nTempStrings;
  double PreFactor=0.0;
  double AbsTol=0.0;
  double RelTol=1.0e-4;
  char *IntegrandFileName=0;
  char *OutFileName=0;
  double DupAbsTol=0.0;
  double DupRelTol=1.0e-6;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"DataFile",        PA_STRING,  1, 1,       (void *)&DataFileName,       0,   "data file"},
     {"FreqColumn",      PA_INT,     1, 1,       (void *)&FreqColumn,         0,   "frequency column index in data file"},
     {"TagColumn",       PA_INT,     1, 1,       (void *)&TagColumn,          0,   "tag column index in data file"},
     {"SDColumn",        PA_INT,     1, 1,       (void *)&SDColumn,           0,   "source/dest column index in data file"},
     {"DataColumn",      PA_INT,     1, MAXDATA, (void *)DataColumns,   &nDataColumns,   "data column index in data file"},
     {"DataName",        PA_STRING,  1, MAXDATA, (void *)DataNames,     &nDataNames,     "name for data column"},
     {"TemperatureFile", PA_STRING,  1, 1,       (void *)&TemperatureFile,    0,   "list of object temperatures"},
     {"Temperature",     PA_STRING,  2, MAXOBJ,  (void *)TempStrings,   &nTempStrings,   "set object temperature"},
     {"PreFactor",       PA_DOUBLE,  1, 1,       (void *)&PreFactor,          0,   "overall multiplicative prefactor"},
     {"IntegrandFile",   PA_STRING,  1, 1,       (void *)&IntegrandFileName,  0,   "name of file for integrand data"},
     {"OutFile",         PA_STRING,  1, 1,       (void *)&OutFileName,        0,   "name of output file"},
     {"DupRelTol",       PA_DOUBLE,  1, 1,       (void *)&DupRelTol,          0,   "relative tolerance for identifying equivalent data points"},
     {"DupAbsTol",       PA_DOUBLE,  1, 1,       (void *)&DupAbsTol,          0,   "absolute tolerance for identifying equivalent data points"},
//
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (DataFileName==0)
   OSUsage(argv[0],OSArray,"--DataFile option is mandatory");

  /***************************************************************/
  /* auto-detect special known file types and autoset values of  */
  /* input parameters                                            */
  /***************************************************************/
  if( strcasestr(DataFileName,"SIFlux") && nDataNames==0 )
   { 
     HMatrix *DataMatrix=new HMatrix(DataFileName);
     int nc=1;
     if( DataMatrix->NC==11 )
      { TagColumn=1;
        nc=2;
      };
     delete DataMatrix;
     Log("Autodetecting spatially-integrated flux data file");
     FreqColumn    = nc++;
     SDColumn      = nc++;
     DataColumns[0] = nc++; DataNames[0]=strdup("PAbs");
     DataColumns[1] = nc++; DataNames[1]=strdup("PRad");
     DataColumns[2] = nc++; DataNames[2]=strdup("XForce");
     DataColumns[3] = nc++; DataNames[3]=strdup("YForce");
     DataColumns[4] = nc++; DataNames[4]=strdup("ZForce");
     DataColumns[5] = nc++; DataNames[5]=strdup("XTorque");
     DataColumns[6] = nc++; DataNames[6]=strdup("YTorque");
     DataColumns[7] = nc++; DataNames[7]=strdup("ZTorque");
     nDataNames=nDataColumns=8;
   };

  if (nDataColumns==0)
   OSUsage(argv[0],OSArray,"you must specify at least one --DataColumn");
  if (nDataNames!=0 && nDataNames!=nDataColumns)
   OSUsage(argv[0],OSArray,"incorrect number of --DataNames");
  if (FreqColumn<1)
   ErrExit("invalid --FreqColumn");
  if (TagColumn<0)
   ErrExit("invalid --TagColumn");
  if (SDColumn<0)
   ErrExit("invalid --SDColumn");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FileData *FD=GetFileData(DataFileName, TagColumn, SDColumn);
  HMatrix *DataMatrix   = FD->DataMatrix;
  int NumObjects        = FD->NumObjects;
  int NumTags           = FD->NumTags;
  cdouble *Tags         = FD->Tags;

  if (NumObjects==0)
   NumObjects=1;
  if (NumTags==0)
   NumTags=1;

  /***************************************************************/
  /* TemperatureMatrix column  0     = environment temperatures  */
  /* TemperatureMatrix columns 1..NO = object 1..NO temperatures */
  /***************************************************************/
  HMatrix *TemperatureMatrix=0;
  if ( TemperatureFile )
   { if (nTempStrings>0)
      ErrExit("--Temperature is incompatible with --TemperatureFile");
     TemperatureMatrix=new HMatrix(TemperatureFile);
     if ( TemperatureMatrix->NC != (NumObjects+1) )
      ErrExit("Data file %s describes %i objects, but temperature file %s describes environment plus %i objects",
               DataFileName, NumObjects, DataFileName, TemperatureMatrix->NC-1);
   }
  else if ( nTempStrings>0 )
   {
     TemperatureMatrix=new HMatrix(1,NumObjects+1);

     for(int nts=0; nts<nTempStrings; nts++)
      { int Index; 
        int Status=sscanf(TempStrings[2*nts+0],"%i",&Index);
        if (Status!=1 || Index<0 || Index>NumObjects )
         ErrExit("invalid object index %s",TempStrings[2*nts+0]);
        double T;
        Status=sscanf(TempStrings[2*nts+1],"%le",&T);
        if (Status!=1) 
         ErrExit("invalid temperature %s",TempStrings[2*nts+1]);
        TemperatureMatrix->SetEntry(0,Index,T);
      };
     printf("Temperatures: \n");
     printf(" Environment: %e\n",TemperatureMatrix->GetEntryD(0,0));
     for(int no=0; no<NumObjects; no++)
      printf(" Object %2i  : %e\n",no+1,TemperatureMatrix->GetEntryD(0,no+1));
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (PreFactor==0.0)
   PreFactor = (TemperatureMatrix==0) ? 1.0 : HBAROMEGA02;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int TotalFreqs  = DataMatrix->NR;
  int NumData     = nDataColumns;
  double *XValues = (double *)mallocEC(TotalFreqs * sizeof(double));
  double *YValues = (double *)mallocEC(TotalFreqs * NumData *sizeof(double));

  int NumTemperatures = TemperatureMatrix ? TemperatureMatrix->NR : 1;

  double *TotalByDest=0;
  if (NumObjects>1)
   TotalByDest=(double *)mallocEC(NumTemperatures*NumData*sizeof(double));

  double *Result = (double *)mallocEC(NumData*sizeof(double));
  double *Error  = (double *)mallocEC(NumData*sizeof(double));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (!OutFileName)
   OutFileName=vstrdup("%s.Integrated",GetFileBase(DataFileName));
  FILE *OutFile=fopen(OutFileName,"a");
  if (!OutFile)
   ErrExit("could not open output file %s",OutFileName);
  fprintf(OutFile,"# data columns: \n");
  int nColumn=0;
  if (TagColumn)
   fprintf(OutFile,"# %i: tag \n",++nColumn);
  if (SDColumn)
   fprintf(OutFile,"# %i: source/dest index \n",++nColumn);
  if (TemperatureMatrix)
   { fprintf(OutFile,"# %i: environment temperature \n",++nColumn);
     for(int no=0; no<NumObjects; no++)
      fprintf(OutFile,"# %i: object %i temperature \n",++nColumn,no+1);
   };
  if (nDataNames)
   for(int nd=0; nd<NumData; nd++)
    fprintf(OutFile,"# %i: %s \n",++nColumn,DataNames[nd]);
  else
   for(int nd=0; nd<NumData; nd++)
    fprintf(OutFile,"# %i: data %i \n",++nColumn,nd);
     
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nTag=0; nTag<NumTags; nTag++)
   for(int nDest=0; nDest<NumObjects; nDest++)
    for(int nSource=0; nSource<NumObjects; nSource++)
      {
        /***************************************************************/
        /* populate XValues and YValues arrays with data for this      */
        /* (tag, sourceObject, destObject) tuple                       */
        /***************************************************************/
        int N=0;
        int SDIndex=10*(nSource+1) + (nDest+1);
        for(int nr=0; nr<DataMatrix->NR; nr++)
         { if (TagColumn && Tags[nTag]!=DataMatrix->GetEntry(nr, TagColumn-1))
            continue;
           if (SDColumn && SDIndex!=((int)DataMatrix->GetEntryD(nr, SDColumn-1)))
            continue;
 
           XValues[N]=DataMatrix->GetEntryD(nr, FreqColumn-1);
           for(int nd=0; nd<NumData; nd++)
            YValues[N*NumData + nd]=DataMatrix->GetEntryD(nr, DataColumns[nd]-1);
           N++;
         };
        RemoveDuplicates(DataFileName,XValues,YValues,NumData,&N,DupRelTol,DupAbsTol);
        if (TagColumn)
         printf("tag %s | ",GetTagString(Tags[nTag]));
        if (SDColumn)
         printf("SD  %i | ",SDIndex);
        double MinFreq = XValues[0];
        double MaxFreq = XValues[N-1];
        printf("%i data points in range (%e,%e)\n",N,MinFreq,MaxFreq);
       
        /***************************************************************/
        /* initialize interpolator                                     */
        /***************************************************************/
        IntegrandData MyID;
        MyID.I1D=new Interp1D(XValues, YValues, N, NumData);
        MyID.PreFactor    = PreFactor;
        MyID.LogFile      = 0;
        MyID.Temperature  = -1.0;
        MyID.TEnvironment = -1.0;

        if (nSource==0 && NumObjects>1)
         memset(TotalByDest,0,NumTemperatures*NumData*sizeof(double));
  
        /***************************************************************/
        /* evaluate integrals at all temperature combinations          */
        /***************************************************************/
        for(int nTemp=0; nTemp<NumTemperatures; nTemp++)
         { 
           if (TemperatureMatrix)
            { MyID.TEnvironment = TemperatureMatrix->GetEntryD(nTemp,0);
              MyID.Temperature  = TemperatureMatrix->GetEntryD(nTemp,nSource+1);
            };
           
           MyID.LogFile=0;
           if (IntegrandFileName)
            { MyID.LogFile=fopen(IntegrandFileName,"a");
              fprintf(MyID.LogFile,"\n\n");
              if (TagColumn)
               fprintf(MyID.LogFile,"# tag %s \n",GetTagString(Tags[nTag]));
              if (SDColumn)
               fprintf(MyID.LogFile,"# SD %i \n",SDIndex);
              if (TemperatureMatrix)
               { fprintf(MyID.LogFile,"# environment / object temperatures: ");
                 fprintf(MyID.LogFile,"%e ",MyID.TEnvironment);
                 for(int no=0; no<NumObjects; no++)
                  fprintf(MyID.LogFile,"%e ",TemperatureMatrix->GetEntryD(nTemp,no+1));
                 fprintf(MyID.LogFile,"\n");
               };
            };

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
           if (TagColumn)
            fprintf(OutFile, "%s ",GetTagString(Tags[nTag]));
           if (SDColumn)
            fprintf(OutFile, "%i%i ",nSource+1,nDest+1);
           if (TemperatureMatrix)
            for(int nc=0; nc<TemperatureMatrix->NC; nc++)
             fprintf(OutFile,"%e ",TemperatureMatrix->GetEntryD(nTemp,nc));
           for(int nd=0; nd<NumData; nd++)
            fprintf(OutFile,"%e ",Result[nd]);
           fprintf(OutFile,"\n");
           fflush(OutFile);
  
           /*--------------------------------------------------------------*/
           /*- accumulate per-destination-object totals -------------------*/
           /*--------------------------------------------------------------*/
           if (NumObjects>1)
            { 
              int Offset=nTemp*NumData;
              for(int nd=0; nd<NumData; nd++)
               TotalByDest[Offset + nd] += Result[nd];

              if (nSource==(NumObjects-1))
               { if (TagColumn)
                  fprintf(OutFile, "%s ",GetTagString(Tags[nTag]));
                 fprintf(OutFile, "0%i ",nDest+1);
                 if (TemperatureMatrix)
                  for(int nc=0; nc<TemperatureMatrix->NC; nc++)
                   fprintf(OutFile,"%e ",TemperatureMatrix->GetEntryD(nTemp,nc));
                 for(int nd=0; nd<NumData; nd++)
                  fprintf(OutFile,"%e ",TotalByDest[Offset+nd]);
                 fprintf(OutFile,"\n");
                 fflush(OutFile);
               };
            };
 
         }; // for(int nTemp=0; nTemp<NumTemperatures; nTemp++)

        delete MyID.I1D;
 
      }; // for(int nTag... for(int nSource... for(int nDest...)
  fclose(OutFile);
 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Integrated data written to %s.\n",OutFileName);
  printf("Thank you for your support.\n");
 
 }
