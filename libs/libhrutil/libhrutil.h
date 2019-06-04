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
 * libhrutil.h  -- header for general-purpose utility library
 *
 * homer reid   -- 10/2008
 */
#ifndef LIBHRUTIL_H
#define LIBHRUTIL_H 

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <complex>
#include <cmath>
#include <vector>
using namespace std;

/***************************************************************/
/* various random constants    *********************************/
/***************************************************************/
#define PLUSMINUS 177   // ascii plus-minus symbol

#ifndef M_PI
#define M_PI		3.1415926535897932384
#endif

#define DEG2RAD M_PI/(180.0) 
#define RAD2DEG 180.0/(M_PI) 

/* isnan was standardized in C in 1999, but in C++ only in 2011,
   and is still not universally available ... punt, since
   checking for nan in IEEE arithmetic is trivial: */
#ifndef ISNAN
#define ISNAN(x) ((x) != (x))
#endif

/***************************************************************/
/* shorthand type definitions **********************************/
/***************************************************************/
//typedef _Complex double cdouble;
#ifndef cdouble
  typedef complex<double> cdouble;
#endif
#ifndef iVec
  typedef vector<int> iVec;
#endif
#ifndef dVec
  typedef vector<double> dVec;
#endif
#ifndef zVec
  typedef vector<cdouble> zVec;
#endif
#ifndef sVec
  typedef vector<char *> sVec;
#endif
#ifndef strVec
  typedef vector<string> strVec;
#endif
#ifndef bVec
  typedef vector<bool> bVec;
#endif

/***************************************************************/
/* useful macros ***********************************************/
/***************************************************************/
#define VARARG_BUFLEN 1000
#define COMPLETE_VARARGS2(format, lastarg, buffer) \
  va_list ap;                                     \
  char buffer[VARARG_BUFLEN];                     \
  va_start(ap,lastarg);                           \
  vsnprintfEC(buffer,VARARG_BUFLEN,format,ap);    \
  va_end(ap);                                     \

#define COMPLETE_VARARGS(format, buffer) COMPLETE_VARARGS2(format, format, buffer)

#define STATICBUFFER(Name,Num,Type)                       \
 static int Name ## Size=0;                               \
 static Type *Name=0;                                     \
 if ( Name ## Size   < (Num) )                            \
  { Name ## Size   = (Num);                               \
    if (Name) free(Name);                                 \
    Name=(Type *)mallocEC((Num)*sizeof(Type));            \
  }; 

/***************************************************************/
/* timing functions  *******************************************/
/***************************************************************/
double Secs();
void Tic(bool MeasureBytesAllocated=false);
double Toc(unsigned long *BytesAllocated=0);

/***************************************************************/
/* string functions  *******************************************/
/***************************************************************/
char *RemoveDirectories(char *s);
char *RemoveExtension(char *s);
char *GetFileExtension(char *s);
char *GetFileBase(const char *s);
int Tokenize(char *s, char **Tokens, int MaxTokens, const char *Separators=" \t\n");
sVec Tokenize(char *s, const char *Separators=" \t\n");
int StrCaseCmp(const char *s1, const char *s2);

FILE *fopenPath(const char *Path, const char *FileName, const char *Mode, char **WhichDir=0);

/***************************************************************/
/* Vararg versions of common functions *************************/
/***************************************************************/
int vsnprintfEC(char *str, size_t size, const char *format, va_list ap);
FILE *vfopen(const char *format, const char *mode, ...);
int vmkdir(const char *format, ...);
int vsystem(const char *format, ...);
void vsetenv(const char *VariableName, const char *format, ...);
char *vstrdup(const char *format, ...);
char *vstrappend(char *s, const char *format, ...);
void vstrncat(char *s, size_t n, const char *format, ...);
void ErrExit(const char *format, ...);
void Warn(const char *format, ...);

/***************************************************************/
/* General-purpose status logging ******************************/
/***************************************************************/
char *GetHostName(); 
char *GetTimeString(char *s=0);
void SetLogFileName(const char *format, ...);
void InitializeLog(char *argv0, const char *Path=0);
void SetConsoleLogging();
void Log(const char *format, ...);
void LogC(const char *format, ...);
void LogThread(const char *format, ...);
void MutexLog(const char *format, ...);
void LogPercent(int n, int N, int Gradations=10);

/***************************************************************/
/* POSIX-specific nonportable system information functions *****/
/***************************************************************/
int GetNumProcs();
void SetNumThreads(int pNumThreads);
int GetNumThreads();
int GetThreadNum();
unsigned long GetMemoryUsage(unsigned long *MemoryUsage=0);
void SetCPUAffinity(int WhichProcessor);

/***************************************************************/
/* complex arithmetic ******************************************/
/***************************************************************/
cdouble expi(double x);
cdouble csqrt2(cdouble z);
int S2CD(const char *str, cdouble *z);
char *CD2S(cdouble z, const char *format);
char *CD2S(cdouble z);
void SetDefaultCD2SFormat(const char *format);
void z2s(cdouble z, char *zStr);
char *z2s(cdouble z);

/***************************************************************/
/* real-valued vector arithmetic                               */
/***************************************************************/
void VecZero(double *v, int N=3);
double *VecCopy(const double *v1, double *v2, int N=3);
double *VecLinComb(double Alpha, const double *v1, double Beta, const double *v2, double *v3, int N=3);
double *VecScale(const double *v1, double Alpha, double *v2, int N=3);
double *VecScale(double *v, double Alpha, int N=3);
double *VecScaleAdd(const double *v1, double Alpha, const double *v2, double *v3, int N=3);
double *VecAdd(const double *v1, const double *v2, double *v3, int N=3);
double *VecSub(const double *v1, const double *v2, double *v3, int N=3);
double *VecPlusEquals(double *v1, double Alpha, const double *v2, int N=3);
double VecDot(const double *v1, const double *v2, int N=3);
double VecNorm(const double *v, int N=3);
double VecNorm2(const double *v, int N=3);
double VecNormalize(double *v, int N=3);
double VecDistance(const double *v1, const double *v2, int N=3);
double VecDistance2(const double *v1, const double *v2, int N=3);

// the following three routines are for the 3D case only
double *VecCross(const double *v1, const double *v2, double *v3);

/***************************************************************/
/* complex-valued vector arithmetic                            */
/***************************************************************/
cdouble *VecLinComb(cdouble Alpha, const cdouble *v1, cdouble Beta, const cdouble *v2, cdouble *v3, int N=3);
cdouble *VecScale(const cdouble *v1, cdouble Alpha, cdouble *v2, int N=3);
cdouble *VecScale(cdouble *v, cdouble Alpha, int N=3);
cdouble *VecScaleAdd(const cdouble *v1, cdouble Alpha, const cdouble *v2, cdouble *v3, int N=3);
cdouble *VecPlusEquals(const cdouble *v1, cdouble Beta, const cdouble *v2, int N=3);
cdouble *VecPlusEquals(cdouble *v1, cdouble Alpha, const cdouble *v2, int N=3);
cdouble *VecPlusEquals(cdouble *v1, cdouble Alpha, const double *v2, int N=3);
cdouble *VecSub(const cdouble *v1, const cdouble *v2, cdouble *v3, int N=3);
cdouble *VecAdd(const cdouble *v1, const cdouble *v2, cdouble *v3, int N=3);
cdouble VecHDot(const cdouble *v1, const cdouble *v2, int N=3);
double VecNorm2(const cdouble *v, int N);
double VecNorm(const cdouble *v, int N);
double VecNormalize(cdouble *v, int N);
double VecDistance2(const cdouble *v1, const cdouble *v2, int N);
double VecDistance(const cdouble *v1, const cdouble *v2, int N);

/***************************************************************/
/* random numbers                                              */
/***************************************************************/
int irand(int A, int B);                       // uniform integer-value 
double randU(double A=0.0, double B=1.0);      // uniform real-valued
double randN(double Sigma=1.0, double Mu=0.0); // normal real-valued

/***************************************************************/
/* finite-difference derivatives *******************************/
/***************************************************************/
typedef void (*VVFunction) (void *UserData, double X[3], cdouble V[3]);
cdouble GetDivCurl(VVFunction VVFun, void *UserData, int Order,
                   double X[3], double Delta, cdouble CurlF[3]);

/***************************************************************/
/* single-precision comparisons of double-precision numbers    */
/***************************************************************/
bool EqualFloat(const double a, const double b);
bool EqualFloat(const cdouble a, const cdouble b);
bool VecEqualFloat(const double *a, const double *b, int N=3);
bool VecEqualFloat(const cdouble *a, const cdouble *b, int N=3);
bool VecClose(const double *a, const double *b, double abstol);

// relative differences |x-y| / (avg( |x|, |y| ) )
double RD(double x, double y);
double RD(cdouble x, cdouble y);

bool IsFinite(double d);
bool IsFinite(cdouble z);

/***************************************************************/
/* matrix arithmetic *******************************************/
/***************************************************************/
bool Matrix2x2_Inverse(double *a[2],double ainv[2][2]);

/***************************************************************/
/* pretty-printing of vectors and comparisons ******************/
/***************************************************************/
void Compare(double *V1, double *V2, int N, const char *str1, const char *str2);
void Compare(cdouble *V1, cdouble *V2, int N, const char *str1, const char *str2);
void Compare(double *V1, double *V2, double *V3, int N, const char *str1, const char *str2, const char *str3);
void Compare(cdouble *V1, cdouble *V2, cdouble *V3, int N, const char *str1, const char *str2, const char *str3);

void fprintVec(FILE *f, double *v, int Length=3, const char *format="%+.8e");
void fprintVecCR(FILE *f, double *v, int Length=3, const char *format="%+.8e");
void fprintVec(FILE *f, cdouble *v, int Length=3, const char *format="%+.8e %+.8e");
void fprintVecCR(FILE *f, cdouble *v, int Length=3, const char *format="%+.8e %+.8e"); 

/***************************************************************/
/* Other functions  ********************************************/
/***************************************************************/
void *mallocEC(size_t size);
void *reallocEC(void *v, size_t size);
char *strdupEC(const char *s);
void *memdup(void *v, size_t size);
void KeyPause();

bool CheckEnv(const char *Name, const char *fmt, void *Destination, bool LogSuccess=true);
bool CheckEnv(const char *Name, double *Destination, bool LogSuccess=true);
bool CheckEnv(const char *Name, int *Destination, bool LogSuccess=true);
bool CheckEnv(const char *Name, char **Destination, bool LogSuccess=true);
bool CheckEnv(const char *Name, bool LogSuccess=true);
void AppendEnv(const char *Name, const char *Separator, const char *format, ...);

FILE *CreateUniqueFile(const char *Base, int ConsoleMessage, char *FileName);
FILE *CreateUniqueFile(const char *Base, int ConsoleMessage); 
FILE *CreateUniqueFile(const char *Base);

void SetCodeMarker(const char *Marker);
void InstallHRSignalHandler();


/***************************************************************/
/* command-line argument processing ****************************/
/***************************************************************/

/*--------------------------------------------------------------*/
/* values for the TYPE field in the ArgStruct structure        -*/
/*--------------------------------------------------------------*/
#define PA_DOUBLE  0
#define PA_INT     1
#define PA_STRING  2
#define PA_BOOL    3
#define PA_CDOUBLE 4

/*--------------------------------------------------------------*/
/* user creates an array of structures like this for passage   */
/* to ProcessArguments().                                      */
/*--------------------------------------------------------------*/
typedef struct ArgStruct
 { const char *Name;
   int Type;
   void *Storage;
   const char *Default;
   const char *Description;
 } ArgStruct;

void ASUsage(char *ProgName, ArgStruct *ASArray, const char *format, ...);
void ProcessArguments(int argc, char *argv[], ArgStruct *ASArray);

/***************************************************************/
/* 20120123: i am expanding 'argument' processing into 'option'*/
/*           processing, which is different in several ways:   */
/*                                                             */
/*            (1) in addition to command-line arguments,       */
/*                options and their arguments can be specified */
/*                in a little text file piped into stdin       */
/*                                                             */
/*            (2) options may take more than one argument.     */
/*                                                             */
/*            (3) options can be specified more than once with */
/*                different arguments                          */
/***************************************************************/
typedef struct OptStruct
 { const char *Name;
   int Type;
   int NumArgs;
   int MaxInstances;
   void *Storage;
   int *NumInstances;
   const char *Description;
 } OptStruct;

void AppendOSUsageMessage(const char *Message);
void OSUsage(const char *ProgName, OptStruct *OSArray, const char *format, ...);
void OSUsage(const char *ProgName, const char *version, OptStruct *OSArray, const char *format, ...);
void ProcessOptions(int argc, char *argv[], OptStruct *OSArray,
                    bool AbortOnUnknown=true, bool ZeroArgs=false);
void ProcessOptions(const char *ArgString, OptStruct *OSArray);

#endif // LIBHRUTIL_H
