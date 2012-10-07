/*
 * Casimir2D.cc -- compute Casimir energy and/or force for a given 2D 
 *              -- geometry under a sequence of geometrical transformations
 * 
 * homer reid -- 11/2008 -- 10/2010
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdarg.h>

#include <libhrutil.h>

#include "libTDRT.h"
#include "scuff-cas2D.h"

/***************************************************************/
/* print usage message *****************************************/
/***************************************************************/
void Usage(const char *ErrMsg, ...)
{ 
  va_list ap;
  char buffer[1000];

  if (ErrMsg)
   { va_start(ap,ErrMsg);
     vsnprintf(buffer,1000,ErrMsg,ap);
     fprintf(stderr,"\nerror: %s (aborting)\n\n",buffer);
     va_end(ap);
   };

  fprintf(stderr,"usage: Casimir2D [options] \n");
  fprintf(stderr," options:                 \n");
  fprintf(stderr,"                           \n");
  fprintf(stderr,"  --geometry  MyGeo.tdgeo  \n");
  fprintf(stderr,"  --translist TransList    \n");
  fprintf(stderr,"                           \n");
  fprintf(stderr,"  --energy                 \n");
  fprintf(stderr,"  --xforce                 \n");
  fprintf(stderr,"  --yforce                 \n");
  fprintf(stderr,"                           \n");
  fprintf(stderr,"  --NumThreads xx             \n");
  fprintf(stderr,"  --TETM                   \n");
  fprintf(stderr,"  --WriteHDF5              \n");
  fprintf(stderr,"                           \n");
  fprintf(stderr,"  --VisualizeOnly          \n");
  fprintf(stderr,"                           \n");
  fprintf(stderr,"  --Xi xx                  \n");
  fprintf(stderr,"  --Q xx                   \n");
  fprintf(stderr,"  --XQList xx              \n");
  fprintf(stderr,"                           \n");
  fprintf(stderr,"  --T xx                   \n");
  fprintf(stderr,"                           \n");
  fprintf(stderr,"  --GroundPlane            \n");
  fprintf(stderr,"                           \n");
  fprintf(stderr,"  --Rectangle x1 y1 x2 y2  \n");
  fprintf(stderr,"                           \n");
  fprintf(stderr,"  --LengthUnit xx          \n");
  fprintf(stderr,"                           \n");
  fprintf(stderr," Cubature options:         \n");
  fprintf(stderr,"  --abstol xx              \n");
  fprintf(stderr,"  --reltol xx              \n");
  fprintf(stderr,"                           \n");
  fprintf(stderr," Logging options:          \n");
  fprintf(stderr,"  --loglevel 0  (no logging)\n");
  fprintf(stderr,"  --loglevel 1  (terse logging)\n");
  fprintf(stderr,"  --loglevel 2  (full logging)\n");
  fprintf(stderr,"\n");
  exit(1);
}

/***************************************************************/
/* write some quick identifying data to the top of an output   */
/* file so we can know what command-line options were used to  */
/* generate the data.                                          */
/*                                                             */
/* FileType==0:  the file is the .log file                     */
/* FileType==1:  the file is the .byXQ output file             */
/* FileType==2:  the file is the .byXi output file             */
/* FileType==3:  the file is the .out output file              */
/***************************************************************/
void WriteFilePreamble(C2DWorkspace *W, FILE *f, int FileType, 
                       double Xi, double q, double T, char *XQListName)
{ 
  char TimeString[200];
  time_t MyTime;
  struct tm *MyTm;
  int no, nt, nc=0;

  MyTime=time(0);
  MyTm=localtime(&MyTime);
  strftime(TimeString,30,"%D::%T",MyTm);

  fprintf(f,"#\n");

  fprintf(f,"# Casimir2D running on %s (%s)\n",getenv("HOST"),TimeString);

  fprintf(f,"# NumThreads: %i \n",W->NumThreads);
  fprintf(f,"#\n");

  fprintf(f,"# reading geometry from file: %s\n",W->G->GeoFileName);
  fprintf(f,"# which describes the following geometry: \n");
  fprintf(f,"#  external medium: %s \n",W->G->MP->Name);
  for(no=0; no<W->G->NumObjects; no++)
   fprintf(f,"#  object %i: %s (file %s, material %s)\n",
                   no,W->G->Objects[no]->Label,
                      W->G->Objects[no]->MeshFileName,
                      W->G->Objects[no]->MP->Name);
  fprintf(f,"#\n");
  fprintf(f,"# total number of basis functions: %i\n",W->G->TotalBFs);
//  fprintf(f,"# average segment length : %.3e\n",W->G->AverageSegmentLength);
  fprintf(f,"#\n");

  if (W->GroundPlane)
   fprintf(f,"#\n# Infinite conducting ground plane at y=0.\n#\n");

  if (W->NumContributingIVs>0)
   { fprintf(f,"#\n");
     fprintf(f,"# Retaining only the contributions of basis functions \n");
     fprintf(f,"# associated with vertices contained in the rectangle \n");
     fprintf(f,"# (%g,%g) -- (%g,%g)\n",W->Rectangle[0],W->Rectangle[1],W->Rectangle[2],W->Rectangle[3]);
     fprintf(f,"# (%i contributing vertices found)\n",W->NumContributingIVs);
     fprintf(f,"#\n");
   };

  fprintf(f,"# reading transformations from file %s\n",W->TransListName);
  fprintf(f,"# which contains the following %i geometrical transformations:\n",W->NumTransforms);
  for(nt=0; nt<W->NumTransforms; nt++)
   fprintf(f,"# %s %s",W->Tags[nt],W->TransLines[nt]);
  fprintf(f,"# \n");

  fprintf(f,"# Quantities requested: ");
  if ( W->WhichQuantities & QUANTITY_ENERGY )
   fprintf(f," energy");
  if ( W->WhichQuantities & QUANTITY_XFORCE )
   fprintf(f," xforce");
  if ( W->WhichQuantities & QUANTITY_YFORCE )
   fprintf(f," yforce");
  fprintf(f,"\n");

  fprintf(f,"# Frequency behavior: ");
  if (Xi!=-1.0 && q!=-1.0)
   fprintf(f,"single-point calculation at (Xi,q)=(%.12g,%.12g)\n",Xi,q);
  else if (Xi!=-1.0)
   fprintf(f,"evaluating q integral at Xi=%.12g\n",Xi);
  else if (q!=-1.0)
   fprintf(f,"evaluating Xi integral at q=%.12g\n",q);
  else if (T!=-1.0)
   fprintf(f,"evaluating matsubara sum at T=%g kelvin\n",T);
  else if (XQListName)
   fprintf(f,"reading (Xi,q) points from file %s\n",XQListName);
  else
   fprintf(f,"evaluating full (Xi,q) integral \n");

  if ( XQListName==0 && !(Xi==-1.0 && q==-1.0) )
   fprintf(f,"# integration tolerances: %.1e (relative), %.1e (absolute)\n",
              W->RelTol, W->AbsTol);

  fprintf(f,"# Other options: ");
  if (W->WriteHDF5)
   fprintf(f," --WriteHDF5");
  if (W->TETM)
   fprintf(f," --TETM");
  fprintf(f,"\n#\n");
  
  if (FileType==1) 
   { 
     fprintf(f,"# data file columns: \n");
     fprintf(f,"# 1: transform tag \n");
     fprintf(f,"# 2: Xi \n");
     fprintf(f,"# 3: q \n");
     nc=4; 
   }
  else if (FileType==2) 
   { fprintf(f,"# data file columns: \n");
     fprintf(f,"# 1: Xi \n");
     fprintf(f,"# 2: temperature\n");
     nc=2;
   }
  else if (FileType==3) 
   { fprintf(f,"# data file columns: \n");
     fprintf(f,"# 1: transform tag \n");
     fprintf(f,"# 2: temperature\n");
     nc=3;
   };

  if (FileType>0)
   { if (W->TETM)
      { if ( W->WhichQuantities & QUANTITY_ENERGY )
         { fprintf(f,"# %i: energy (TE)\n",nc++);
           fprintf(f,"# %i: energy (TM)\n",nc++);
         }
        if ( W->WhichQuantities & QUANTITY_XFORCE )
         { fprintf(f,"# %i: xforce (TE)\n",nc++);
           fprintf(f,"# %i: xforce (TM)\n",nc++);
         }
        if ( W->WhichQuantities & QUANTITY_YFORCE )
         { fprintf(f,"# %i: yforce (TE)\n",nc++);
           fprintf(f,"# %i: yforce (TM)\n",nc++);
         };
      }
     else
      { if ( W->WhichQuantities & QUANTITY_ENERGY )
         fprintf(f,"# %i: energy \n",nc++);
        if ( W->WhichQuantities & QUANTITY_XFORCE )
         fprintf(f,"# %i: xforce \n",nc++);
        if ( W->WhichQuantities & QUANTITY_YFORCE )
         fprintf(f,"# %i: yforce \n",nc++);
      };
   };
  
}

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{

  int narg, NumThreads, WhichQuantities;
  int CubatureMethod, TETM, GroundPlane;
  int WriteHDF5;
  TDRTGeometry *G;
  double Xi, Q, T;
  double AbsTol, RelTol;
  double LengthUnit;
  double RectangleBuffer[4], *Rectangle;

  char *GeoFileBase, *TransListName, *XQListName;
  char OutFileName[200];
  int nt, nq, ntnq, NTNQ;
  int VisualizeOnly;

  C2DWorkspace *W;
  FILE *f;
  double *I, *E;
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  InstallHRSignalHandler();

  /***************************************************************/
  /***************************************************************/
  /* 1. process command-line arguments ***************************/
  /***************************************************************/
  /***************************************************************/
  CubatureMethod=CMETHOD_DC;
  G=0;
  TransListName=0;
  Xi=Q=T=-1.0;
  NumThreads=0;
  WhichQuantities=0;
  TETM=GroundPlane=WriteHDF5=0;
  XQListName=0;
  AbsTol=1.0e-8;
  RelTol=1.0e-3;
  VisualizeOnly=0;
  LengthUnit=0.0;
  Rectangle=0;
  for(narg=1; narg<argc; narg++)
   { 
     if ( !StrCaseCmp(argv[narg],"--geometry") )
      { if (narg+1>=argc)
         ErrExit("--geometry option requires one argument");
        if (G)
         ErrExit("only one --geometry may be specified");
        G=new TDRTGeometry(argv[++narg]);
      }
     else if ( !StrCaseCmp(argv[narg],"--translist") )
      { if (narg+1>=argc)
         ErrExit("--translist option requires one argument");
        if (TransListName)
         ErrExit("only one --translist may be specified");
        TransListName=argv[++narg];
      }
     else if ( !StrCaseCmp(argv[narg],"--energy") )
      { WhichQuantities|=QUANTITY_ENERGY;
        printf("Computing energy.\n");
      }
     else if ( !StrCaseCmp(argv[narg],"--xforce") )
      { WhichQuantities|=QUANTITY_XFORCE;
        printf("Computing X force.\n");
      }
     else if ( !StrCaseCmp(argv[narg],"--yforce") )
      {
        WhichQuantities|=QUANTITY_YFORCE;
        printf("Computing Y force.\n");
      }
     else if ( !StrCaseCmp(argv[narg],"--TETM") )
      { TETM=1;
        printf("Separating TE and TM contributions.\n");
      }
     else if ( !StrCaseCmp(argv[narg],"--WriteHDF5") )
      { WriteHDF5=1;
        printf("Exporting matrices to MATLAB files.\n");
      }
     else if ( !StrCaseCmp(argv[narg],"--Xi") )
      { if (narg+1>=argc)
         ErrExit("--Xi option requires one argument");
        if (    1!=sscanf(argv[++narg],"%le",&Xi) || Xi<0.0 )
         ErrExit("invalid Xi value %s",argv[narg]);
      }
     else if ( !StrCaseCmp(argv[narg],"--Q") )
      { if (narg+1>=argc)
         ErrExit("--Q option requires one argument");
        if (    1!=sscanf(argv[++narg],"%le",&Q) )
         ErrExit("invalid Q value %s",argv[narg]);
      }
     else if ( !StrCaseCmp(argv[narg],"--XQList") )
      { if (narg+1>=argc)
         ErrExit("--XQList option requires an argument");
        XQListName=argv[++narg];
        if ( !(f=fopen(XQListName,"r")) )
         ErrExit("could not open file %s",XQListName);
        fclose(f);
      } 
     else if ( !StrCaseCmp(argv[narg],"--T") )
      { if (narg+1>=argc)
         ErrExit("--T option requires one argument");
        if (    1!=sscanf(argv[++narg],"%le",&T) || T<0.0 )
         ErrExit("invalid temperature value %s",argv[narg]);
      }
     else if ( !StrCaseCmp(argv[narg],"--LengthUnit") )
      { if (narg+1>=argc)
         ErrExit("--LengthUnit option requires one argument");
        if (    1!=sscanf(argv[++narg],"%le",&LengthUnit) || LengthUnit<=0.0 )
         ErrExit("invalid length unit %s",argv[narg]);
        printf("Setting length unit to %e.\n",LengthUnit);
      }
     else if ( !StrCaseCmp(argv[narg],"--NumThreads") )
      { if (++narg==argc) 
         ErrExit("--NumThreads option requires an argument");
        NumThreads=-1;
        sscanf(argv[narg],"%i",&NumThreads);
        if ( NumThreads<0 || NumThreads>100 ) 
         ErrExit("invalid number of threads %s",argv[narg]);
        printf("Using %i threads.\n",NumThreads);
      }
     else if ( !StrCaseCmp(argv[narg],"--LogLevel") )
      { if (++narg==argc) 
         ErrExit("--LogLevel option requires an argument");
        sscanf(argv[narg],"%i",&(TDRTGeometry::LogLevel));
        if ( TDRTGeometry::LogLevel<0 || TDRTGeometry::LogLevel>2) 
         ErrExit("invalid log level %s",argv[narg]);
        printf("Setting log level to %i.\n",TDRTGeometry::LogLevel);
      }
     else if ( !StrCaseCmp(argv[narg],"--abstol") )
      { if (narg+1>=argc)
         ErrExit("--abstol option requires one argument");
        if (    1!=sscanf(argv[++narg],"%le",&AbsTol) )
         ErrExit("invalid AbsTol value %s",argv[narg]);
      }
     else if ( !StrCaseCmp(argv[narg],"--reltol") )
      { if (narg+1>=argc)
         ErrExit("--abstol option requires one argument");
        if (    1!=sscanf(argv[++narg],"%le",&RelTol) )
         ErrExit("invalid RelTol value %s",argv[narg]);
      }
     else if ( !StrCaseCmp(argv[narg],"--VisualizeOnly") )
      { 
        VisualizeOnly=1;
        fprintf(stderr,"Generating visualization files only.\n");
      }
     else if ( !StrCaseCmp(argv[narg],"--GroundPlane") )
      { 
        GroundPlane=1;
        fprintf(stderr,"Infinite conducting plane at y=0.\n");
      }
     else if ( !StrCaseCmp(argv[narg],"--Rectangle") )
      { 
        fprintf(stderr,"Infinite conducting plane at y=0.\n");
        if (narg+4>=argc)
         ErrExit("--rectangle option requires four arguments");
        if (    1!=sscanf(argv[++narg],"%le",RectangleBuffer) 
             || 1!=sscanf(argv[++narg],"%le",RectangleBuffer+1) 
             || 1!=sscanf(argv[++narg],"%le",RectangleBuffer+2) 
             || 1!=sscanf(argv[++narg],"%le",RectangleBuffer+3)
           ) ErrExit("invalid argument passed to --rectangle");
       Rectangle=RectangleBuffer;
       fprintf(stderr,"Retaining contributions of vertices within rectangle\n");
       fprintf(stderr," (%g,%g) -- ",Rectangle[0],Rectangle[1]);
       fprintf(stderr," (%g,%g).\n",Rectangle[2],Rectangle[3]);
      }
     else
      Usage("unknown command-line argument %s",argv[narg]);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (LengthUnit!=0.0)
   MatProp::SetLengthUnit(LengthUnit);

  /***************************************************************/
  /* some quick sanity checks on input arguments******************/
  /***************************************************************/
  if (G==0)
   Usage("--geometry argument is mandatory");
  if (TransListName==0)
   Usage("--translist argument is mandatory");
  if ( WhichQuantities==0 && VisualizeOnly!=1 )
   Usage("you must specify at least one quantity to compute");

  /***************************************************************/
  /* 1b. tell user what we understood her to be asking us to do **/
  /***************************************************************/
  if (T==0.0)     // in case the user said --T 0.0
   T=Xi=Q=-1.0;
  
  if ( XQListName )
   { if ( Xi!=-1.0 || Q!=-1.0 )
      ErrExit("--XQList option may not be used with --Xi or --Q");
     printf("Reading (Xi,Q) points from file %s.\n",XQListName);
   }
  else if ( Xi==-1.0 && Q!=-1.0 )
   printf("Evaluating Xi integral at Q=%g.\n",Q);
  else if ( Xi!=-1.0 && Q==-1.0 )
   printf("Evaluating Q integral at Xi=%g.\n",Xi);
  else if ( Xi!=-1.0 && Q!=-1.0 )
   printf("Single-point calculation at (Xi,Q)=(%g,%g).\n",Xi,Q);
  else if ( T!=-1.0 )
   printf("Evaluating Matsubara sum at T=%g kelvin.\n",T);
  else
   printf("Evaluating integral over (Xi,Q) quarter-plane.\n");

  /***************************************************************/
  /* 2. initialize workspace structure and other data structures */
  /*    and arrays                                               */
  /***************************************************************/
  GeoFileBase=strdup(GetFileBase(G->GeoFileName));
  W=CreateC2DWorkspace(G, TransListName, WhichQuantities, Rectangle,
                       NumThreads, TETM, GroundPlane, WriteHDF5, VisualizeOnly);

  if (VisualizeOnly)
   exit(1);

  W->AbsTol=AbsTol;
  W->RelTol=RelTol;
  
  NTNQ=W->NTNQ;

  if (TDRTGeometry::LogLevel>=1)
   { SetLogFileName("%s.log",GeoFileBase);
     f=vfopen("%s.log","a",GeoFileBase);
     WriteFilePreamble(W, f, 0, Xi, Q, T, XQListName);
     fclose(f);
   };

  f=fopen(W->ByXQFileName,"a");
  WriteFilePreamble(W, f, 1, Xi, Q, T, XQListName);
  fclose(f);

  if (T!=-1.0)
   { 
     f=fopen(W->ByXiFileName,"a");
     WriteFilePreamble(W, f, 2, Xi, Q, T, XQListName);
     fclose(f);
   };

  I=(double *)malloc(NTNQ*sizeof(double));
  E=(double *)malloc(NTNQ*sizeof(double));
  
  /***************************************************************/
  /***************************************************************/
  /* 3. now switch off to determine what to do.                  */
  /***************************************************************/
  /***************************************************************/

  /*-------------------------------------------------------------*/
  /* 3A. if a single-point calculation was requested, compute    */
  /* integrand at a single point.                                */
  /*-------------------------------------------------------------*/
  if ( (Xi!=-1.0) && (Q!=-1.0) )
   { 
     XQIntegrand(W, Xi, Q, I);
   }
  /*-------------------------------------------------------------*/
  /* 3B. if a Xi point was specified, evaluate the Q integral    */
  /*-------------------------------------------------------------*/
  else if ( Xi!=-1.0 )
   { 
     EvaluateQIntegral(W, Xi, I, E);
   }
  /*-------------------------------------------------------------*/
  /* 3C. if a Q point was specified, evaluate the Xi integral    */
  /*-------------------------------------------------------------*/
  else if ( Q !=-1.0 )
   {    
     EvaluateXiIntegral(W, Q, I, E);

     FILE *ByQFile=vfopen("%s.byQ","a",GeoFileBase);
     if (ByQFile)
      { 
        for (ntnq=nt=0; nt<W->NumTransforms; nt++) 
         { fprintf(ByQFile,"%s %e ",W->Tags[nt],Q);
           for (nq=0; nq<W->NumQuantities; nq++, ntnq++) 
            fprintf(ByQFile,"%.15e ",I[ntnq]);
           ntnq-=W->NumQuantities;
           for(nq=0; nq<W->NumQuantities; nq++, ntnq++)
            fprintf(ByQFile,"%.2e ",E[ntnq]);
           fprintf(ByQFile,"\n");
         };
        fclose(ByQFile);
      };

   }
  /*-------------------------------------------------------------*/
  /* 3D. if a temperature was specified, evaluate the matsubara  */
  /*     sum                                                     */
  /*-------------------------------------------------------------*/
  else if ( T !=-1.0 )
   {
     EvaluateMatsubaraSum(W, T, I, E);
   }
  /*-------------------------------------------------------------*/
  /* 3E. if an XQ file was specified, read (Xi,Q) points from    */
  /*     the XQ file and evaluate the (Xi,q) integrand at        */
  /*     each point.                                             */
  /*     Note in this case we exit without printing output to    */
  /*     the console or to the .out file (console output is      */
  /*     handled by ProcessXQList, and there is no .out file     */
  /*     generated in this case).                                */
  /*-------------------------------------------------------------*/
  else if ( XQListName )
   { ProcessXQList(W, XQListName);

     printf("Data written to file %s.\n",W->ByXQFileName);
     printf("Thank you for your support.\n");

     exit(1);
   }
  /*-------------------------------------------------------------*/
  /* 3F. otherwise evaluate the full Xi-Q integral using the     */
  /*     requested cubature scheme.                              */
  /*-------------------------------------------------------------*/
  else
   { 
     EvaluateXQIntegral(W,I,E);

   }; // if .. else 

  /***************************************************************/
  /* 4. write results to console     *****************************/
  /***************************************************************/
  printf("\n");
  printf("   TAG     ");
  if (W->TETM)
   { 
     if ( W->WhichQuantities & QUANTITY_ENERGY )
      { printf("     ENERGY (TE)    ");
        printf("     ENERGY (TM)    "); 
      }
     if ( W->WhichQuantities & QUANTITY_XFORCE )
      { printf("     XFORCE (TE)    ");
        printf("     XFORCE (TM)    "); 
      }
     if ( W->WhichQuantities & QUANTITY_YFORCE )
      { printf("     YFORCE (TE)    "); 
        printf("     YFORCE (TM)    "); 
      };
   }
  else
   { 
     if ( W->WhichQuantities & QUANTITY_ENERGY )
      printf("       ENERGY       ");
     if ( W->WhichQuantities & QUANTITY_XFORCE )
      printf("       XFORCE       ");
     if ( W->WhichQuantities & QUANTITY_YFORCE )
      printf("       YFORCE       ");
   };
  printf("\n");

  printf("---------  ");
  for(nt=0; nt<W->NumQuantities; nt++)
   printf("  ----------------  ");
  printf("\n");

  for(ntnq=nt=0; nt<W->NumTransforms; nt++)
   { printf("%-9s  ",W->Tags[nt]); 
     for(nq=0; nq<W->NumQuantities; nq++)
      printf("  %+15.9e  ",I[ntnq++]);
     printf("\n");
   };
    
  /***************************************************************/
  /* 5. write results to output file if we did an integral or sum*/
  /***************************************************************/
  if ( (Xi==-1.0) || (Q==-1.0) || (T!=-1.0) )
   { 
     snprintf(OutFileName,200,"%s.out",GeoFileBase);
     nt=0;
     while( (f=fopen(OutFileName,"r"))!=0 )
      { fclose(f);
        sprintf(OutFileName,"%s.out.%i",GeoFileBase,++nt);
      };
     if (nt!=0)
      { fprintf(stderr,"\n** WARNING: file %s.out exists; ",GeoFileBase);
        fprintf(stderr,"writing output to %s\n\n",OutFileName);
      }; 
     f=fopen(OutFileName,"a");
     WriteFilePreamble(W, f, 3, Xi, Q, T, XQListName);
     for(ntnq=nt=0; nt<W->NumTransforms; nt++)
      { fprintf(f,"%s ",W->Tags[nt]);
        if (T==-1.0)
         fprintf(f,"0.0 ");
        else
         fprintf(f,"%e ",T);
        for(nq=0; nq<W->NumQuantities; nq++, ntnq++)
         fprintf(f,"%.15e ",I[ntnq]);
        ntnq-=W->NumQuantities;
        for(nq=0; nq<W->NumQuantities; nq++, ntnq++)
         fprintf(f,"%.2e ",E[ntnq]);
        fprintf(f,"\n");
      };
     fclose(f);
   };

  printf("Thank you for your support.\n");
}  

/***************************************************************/
/* 4. write results to console     *****************************/
/***************************************************************/
void PrintConsoleOutput(C2DWorkspace *W, double *I)
{ 
  int nt, nq, ntnq;

  printf("\n");

  printf("   TAG     ");
  if (W->TETM)
   { 
     if ( W->WhichQuantities & QUANTITY_ENERGY )
      { printf("     ENERGY (TE)    ");
        printf("     ENERGY (TM)    "); 
      }
     if ( W->WhichQuantities & QUANTITY_XFORCE )
      { printf("     XFORCE (TE)    ");
        printf("     XFORCE (TM)    "); 
      }
     if ( W->WhichQuantities & QUANTITY_YFORCE )
      { printf("     YFORCE (TE)    "); 
        printf("     YFORCE (TM)    "); 
      };
   }
  else
   { 
     if ( W->WhichQuantities & QUANTITY_ENERGY )
      printf("       ENERGY       ");
     if ( W->WhichQuantities & QUANTITY_XFORCE )
      printf("       XFORCE       ");
     if ( W->WhichQuantities & QUANTITY_YFORCE )
      printf("       YFORCE       ");
   };
  printf("\n");

  printf("---------  ");
  for(nt=0; nt<W->NumQuantities; nt++)
   printf("  ----------------  ");
  printf("\n");

  for(ntnq=nt=0; nt<W->NumTransforms; nt++)
   { printf("%-9s  ",W->Tags[nt]); 
     for(nq=0; nq<W->NumQuantities; nq++)
      printf("  %+15.9e  ",I[ntnq++]);
     printf("\n");
   };
}
