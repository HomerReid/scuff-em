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
 * CheckData.cc
 *
 * homer reid        -- 3/2016
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include "libhrutil.h"

#define MAXSTR 1000
#define MAXTOK 50

#define LOGLEVEL_SILENT  0
#define LOGLEVEL_DEFAULT 1
#define LOGLEVEL_VERBOSE 2

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
typedef struct Checklist
 { 
   int NumKeys;
   int *KeyColumns;
   char **KeyNames;

   int NumData;
   int *DataColumns;
   char **DataNames;

   bool CaseSensitive;
   double AbsTol, RelTol;
   int LogLevel;
   int MaxColumn;
   
 } Checklist;

typedef struct DataRecord
 { char **KeyValues;
   char **DataValues;
   int FileLine;
 } DataRecord;

typedef struct DataSet
 { DataRecord **DataRecords;
   int NumRecords;
   int TotalDataItems;
   char *FileName;
 } DataSet;

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
Checklist *ReadChecklist(char *FileName, int LogLevel)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumKeys=0;
  int *KeyColumns=0;
  char **KeyNames=0;
  int NumData=0;
  int *DataColumns=0;
  char **DataNames=0;
  bool CaseSensitive=false;
  double AbsTol=0.0;
  double RelTol=1.0e-2;
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *f=fopen(FileName,"r");
  if (!f) ErrExit("could not open %s",FileName);
  char Line[MAXSTR];
  int LineNum=0;
  int MaxColumn=0;
  while (fgets(Line, MAXSTR, f))
   {
     LineNum++;

     char *Tokens[MAXTOK];
     int NumTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

     if (!strcasecmp(Tokens[0],"KEY") || !strcasecmp(Tokens[0],"DATA"))
      { 
        char *Name = Tokens[1];
        int Column;
        if ( NumTokens!=3 || sscanf(Tokens[2],"%i",&Column)!=1 )
         ErrExit("%s:%i: syntax error",FileName,LineNum);
        Column-=1;

        if (Column>MaxColumn) 
         MaxColumn=Column;

        switch(toupper(Tokens[0][0]))
         { case 'K': 
            KeyColumns = (int *)reallocEC((void *)KeyColumns, (NumKeys+1)*sizeof(int));
            KeyNames   = (char **)reallocEC((void *)KeyNames, (NumKeys+1)*sizeof(char *));
            KeyColumns[NumKeys] = Column;
	    KeyNames[NumKeys]   = strdup(Name);
            NumKeys++;
            break;

           case 'D':
            DataColumns = (int *)reallocEC((void *)DataColumns, (NumData+1)*sizeof(int));
            DataNames   = (char **)reallocEC((void *)DataNames, (NumData+1)*sizeof(char *));
            DataColumns[NumData] = Column;
            DataNames[NumData]   = strdup(Name);
            NumData++;
            break;
         };
      }
     else if (!strcasecmp(Tokens[0],"CASE_SENSITIVE"))
      { 
        CaseSensitive=true;
      }
     else if (!strcasecmp(Tokens[0],"RELTOL") || !strcasecmp(Tokens[0],"ABSTOL") )
      { 
        double *Dest = (toupper(Tokens[0][0])=='R') ? &RelTol : &AbsTol;
        if (NumTokens!=2 || (sscanf(Tokens[1],"%le",Dest)!=1) )
         ErrExit("%s:%i: syntax error",FileName,LineNum);
      }
     else
      ErrExit("%s:%i: syntax error",FileName,LineNum);
   };
  fclose(f);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Checklist *CL=(Checklist *)mallocEC(sizeof(Checklist));
  CL->NumKeys=NumKeys;
  CL->KeyColumns=KeyColumns;
  CL->KeyNames=KeyNames;
  CL->NumData=NumData;
  CL->DataColumns=DataColumns;
  CL->DataNames=DataNames;
  CL->CaseSensitive=CaseSensitive;
  CL->AbsTol=AbsTol;
  CL->RelTol=RelTol;
  CL->LogLevel=LogLevel;
  CL->MaxColumn=MaxColumn;

  if (LogLevel>=LOGLEVEL_VERBOSE)
   { printf("Checklist %s: %i keys, %i data items per line\n",
             FileName, NumKeys, NumData);
     printf(" RelTol=%e, AbsTol=%e, %scase sensitive\n",
             RelTol, AbsTol, CaseSensitive ? "" : "not ");
   };

  return CL;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
DataRecord *CreateDataRecord(Checklist *CL, int FileLine )
{
  DataRecord *DR = (DataRecord *)mallocEC(sizeof(DataRecord));
  DR->KeyValues  = (char **)mallocEC(CL->NumKeys);
  DR->DataValues = (char **)mallocEC(CL->NumData);
  DR->FileLine   = FileLine;
  return DR;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
DataSet *ReadDataSet(char *FileName, Checklist *CL)
{
  int NumKeys      = CL->NumKeys;
  int *KeyColumns  = CL->KeyColumns;
  char **KeyNames  = CL->KeyNames;
  int NumData      = CL->NumData;
  int *DataColumns = CL->DataColumns;
  char **DataNames = CL->DataNames;

  DataRecord **DataRecords=0;
  int NumRecords=0;
  int NumWildcards=0;

  int MinTokens = NumKeys + NumData;
  if (MinTokens < CL->MaxColumn)
   MinTokens = CL->MaxColumn;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *f=fopen(FileName,"r");
  if (!f) ErrExit("could not open %s",FileName);
  char Line[MAXSTR];
  int LineNum=0;
  while (fgets(Line, MAXSTR, f))
   {
     LineNum++;

     char *Tokens[MAXTOK];
     int NumTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue;

     if ( NumTokens < MinTokens )
      { Warn("%s:%i: not enough tokens on line (skipping)...");
        continue;
      };

     DataRecord *DR = CreateDataRecord(CL, LineNum);

     for(int nk=0; nk<NumKeys; nk++)
      DR->KeyValues[nk] = strdup(Tokens[ KeyColumns[nk] ]);

     for(int nd=0; nd<NumData; nd++)
      { DR->DataValues[nd] = strdup(Tokens[ DataColumns[nd] ]);
        if ( DR->DataValues[nd][0]=='*' )
         NumWildcards++;
      };

     DataRecords = (DataRecord **)reallocEC(DataRecords, (NumRecords+1)*sizeof(DataRecord *));
     DataRecords[NumRecords]=DR;
     NumRecords++;
   };
  fclose(f);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  DataSet *DS = (DataSet *)mallocEC(sizeof(DataSet));
  DS->DataRecords  = DataRecords;
  DS->NumRecords   = NumRecords;
  DS->TotalDataItems = NumRecords*NumData - NumWildcards;
  DS->FileName     = strdup(FileName); 

  if (CL->LogLevel>=LOGLEVEL_DEFAULT)
   printf("Data set %s: %i total data items (%i wildcards)\n",
           FileName,DS->TotalDataItems,NumWildcards);

  return DS;

}

/***************************************************************/
/* returns true if the data agree ******************************/
/***************************************************************/
bool Agree(char *TestDatum, char *RefDatum, Checklist *CL)
{
  if (RefDatum==0 || TestDatum==0)
   return false;

  if (TestDatum[0]=='*' || RefDatum[0]=='*')
   return true;

  double TestValue, RefValue;
  bool TestIsNumerical = (1==sscanf(TestDatum,"%le",&TestValue));
  bool  RefIsNumerical = (1==sscanf( RefDatum,"%le", &RefValue));

  if (RefIsNumerical!=TestIsNumerical)
   return false;

  if (RefIsNumerical==false)
   return CL->CaseSensitive ? !strcmp(TestDatum,RefDatum)
                            : !strcasecmp(TestDatum,RefDatum);

  if ( IsFinite(TestValue) != IsFinite(RefValue) )
   return false;

  if ( RefValue==0.0 )
   return TestValue==RefValue;

  double Diff=fabs(TestValue-RefValue);
  bool AbsMatch = (Diff < CL->AbsTol);
  bool RelMatch = (Diff < (CL->RelTol) * fabs(RefValue));

  return AbsMatch || RelMatch;

}

DataRecord *FindDataRecord(DataSet *DS, char **KeyValues, Checklist *CL)
{
  for(int ndr=0; ndr<DS->NumRecords; ndr++)
   { DataRecord *DR = DS->DataRecords[ndr];

     bool IsAMatch=true;
     for(int nd=0; IsAMatch && nd<CL->NumKeys; nd++)
      IsAMatch &= Agree(DR->KeyValues[nd], KeyValues[nd], CL);

     if (IsAMatch)
      return DR;    
   };
  return 0;
    
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int CompareDataSets(DataSet *TestDS, DataSet *RefDS, Checklist *CL, FILE *ErrorFile)
{
  int NumKeys      = CL->NumKeys;
  int *KeyColumns  = CL->KeyColumns;
  char **KeyNames  = CL->KeyNames;
  int NumData      = CL->NumData;
  int *DataColumns = CL->DataColumns;
  char **DataNames = CL->DataNames;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int DataMatches=0;
  for(int ndr=0; ndr<RefDS->NumRecords; ndr++)
   { 
     DataRecord *RefDR = RefDS->DataRecords[ndr];

     DataRecord *TestDR = FindDataRecord(TestDS, RefDR->KeyValues, CL);

     if (!TestDR)
      { if (ErrorFile)
         fprintf(ErrorFile,"Data file has no line matching ref file line %i\n",
                            RefDR->FileLine);
        continue;
      };

     for(int nd=0; nd<NumData; nd++)
      { 
        // skip wildcards in the reference data
        if (RefDR->DataValues[nd][0]=='*')
         continue;

        if ( Agree(TestDR->DataValues[nd], RefDR->DataValues[nd], CL) )
         { 
           DataMatches++;
         }
        else
         {
           if (ErrorFile)
            { fprintf(ErrorFile,"\n*********************************\n\n");
              fprintf(ErrorFile,"Data item %s for keys ",DataNames[nd]);
              for(int nk=0; nk<NumKeys; nk++)
               fprintf(ErrorFile,"%s=%s%c ",KeyNames[nk],RefDR->KeyValues[nk],(nk==NumKeys-1)?':':',');
              fprintf(ErrorFile,"\n");
              fprintf(ErrorFile," should be: %s \n",RefDR->DataValues[nd]);
              fprintf(ErrorFile," is:        %s \n",TestDR->DataValues[nd]);
            };
         }; // if(Match)...else

      }; // for(int nd=0 ...)

   }; // for(int ndr=0; nds<RefDS->NumRecords; nds++)

  return DataMatches;

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
  SetLogFileName("CheckSCUFFData.log");
//
  char *DataFile=0;
  char *RefFile=0;
  char *CLFile=0;
  char *ErrorFileName=0;
  bool Silent=false;
  bool Verbose=false;
//
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"Data",           PA_STRING,  1, 1,       (void *)&DataFile,   0,             "data file"},
     {"Reference",      PA_STRING,  1, 1,       (void *)&RefFile,    0,             "reference file"},
     {"Checklist",      PA_STRING,  1, 1,       (void *)&CLFile,     0,             "checklist file"},
     {"ErrorFile",      PA_STRING,  1, 1,       (void *)&ErrorFileName, 0,         "name of error output file"},
     {"silent",         PA_BOOL,    0, 1,       (void *)&Silent,     0,             "print no output"},
     {"verbose",        PA_BOOL,    0, 1,       (void *)&Verbose,    0,             "print copious output"},
/**/
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (DataFile==0)
   OSUsage(argv[0], OSArray, "--data option is mandatory");
  if (RefFile==0)
   OSUsage(argv[0], OSArray, "--reference option is mandatory");
  if (CLFile==0)
   OSUsage(argv[0], OSArray, "--Checklist option is mandatory");

  int LogLevel = LOGLEVEL_DEFAULT;
  if (Silent)
   LogLevel = LOGLEVEL_SILENT;
  if (Verbose)
   LogLevel = LOGLEVEL_VERBOSE;
  
  Checklist *CL   = ReadChecklist(CLFile, LogLevel);
  DataSet *TestDS = ReadDataSet(DataFile, CL);
  DataSet *RefDS  = ReadDataSet(RefFile, CL);
 
  if (ErrorFileName==0)
   ErrorFileName=vstrdup("%s.errors",GetFileBase(DataFile));
  FILE *ErrorFile = 0;
  if (ErrorFileName)
   { ErrorFile=fopen(ErrorFileName,"w");
     fprintf(ErrorFile,"CheckSCUFFData running on %s (%s)\n",GetHostName(),GetTimeString());
     fprintf(ErrorFile,"Checklist:      %s \n",CLFile);
     fprintf(ErrorFile,"Data      file: %s \n",DataFile);
     fprintf(ErrorFile,"Reference file: %s \n",RefFile);
   };

  int DataMatches = CompareDataSets(TestDS, RefDS, CL, ErrorFile);

  if (ErrorFile)
   fclose(ErrorFile);

  printf("%i/%i matches.\n",DataMatches,RefDS->TotalDataItems);


}
