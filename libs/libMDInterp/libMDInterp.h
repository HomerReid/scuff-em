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
 *  libMDInterp.h -- header file for multidimensional interpolation 
 *                -- library
 *
 *  homer reid    -- 3/2011
 */

#ifndef LIBMDINTERP_H
#define LIBMDINTERP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include <libhmat.h> // only needed for GetInterpolationError

#define LMDI_LOGLEVEL_NONE    0
#define LMDI_LOGLEVEL_TERSE   1
#define LMDI_LOGLEVEL_VERBOSE 2

/***************************************************************/
/* prototype for a user-supplied function of one variable.     */
/*                                                             */
/* X is the input to the function (the point at which the      */
/* function is to be evaluated.)                               */
/*                                                             */
/* PhiVD stands for 'phi value and derivative.' on return,     */
/* this array must be filled in as follows:                    */
/*                                                             */
/* PhiVD[0]  = Phi_1                                           */
/* PhiVD[1]  = dPhi_1 /dX                                      */
/*                                                             */
/* PhiVD[2]  = Phi_2                                           */
/* PhiVD[3]  = dPhi_2 /dX                                      */
/*                                                             */
/* ... and so on for all components of the Phi function vector.*/
/* (the number of components is specified by the nFun parameter*/
/* to the Interp1D class constructor.)                         */
/*                                                             */
/* note: the function should be thread-safe, as it may be      */
/* called simultaneously by multiple concurrent threads.       */
/***************************************************************/ 
typedef void (*Phi1D)(double X1, void *UserData, double *PhiVD);

/***************************************************************/ 
/* prototype for a user-supplied function of two variables.    */
/*                                                             */
/* X1, X2 are the inputs to the function (the coordinates      */
/* of the point at which the function is to be evaluated.)     */
/*                                                             */
/* PhiVD stands for 'phi values and derivatives.' on return,   */
/* this array must be filled in as follows:                    */
/*                                                             */
/* PhiVD[0]  = Phi_1                                           */
/* PhiVD[1]  = dPhi_1 /dX1                                     */
/* PhiVD[2]  = dPhi_1 /dX2                                     */
/* PhiVD[3]  = d^2 Phi_1 / dX1 dX2                             */
/*                                                             */
/* PhiVD[4]  = Phi_2                                           */
/* PhiVD[5]  = dPhi_2 /dX1                                     */
/* PhiVD[6]  = dPhi_2 /dX2                                     */
/* PhiVD[7]  = d^2 Phi_2 / dX1 dX2                             */
/*                                                             */
/* ... and so on for all components of the Phi function vector.*/
/* (the number of components is specified by the nFun parameter*/
/* to the Interp2D class constructor.)                         */
/*                                                             */
/* note: the function should be thread-safe, as it may be      */
/* called simultaneously by multiple concurrent threads.       */
/***************************************************************/ 
typedef void (*Phi2D)(double X1, double X2, void *UserData, double *PhiVD);

/***************************************************************/
/* prototype for a user-supplied function of three variables.  */
/*                                                             */
/* X1, X2, X3 are the inputs to the function (the coordinates  */
/* of the point at which the function is to be evaluated.)     */
/*                                                             */
/* PhiVD stands for 'phi values and derivatives.' on return,   */
/* this array must be filled in as follows:                    */
/*                                                             */
/* PhiVD[0]  = Phi_1                                           */
/* PhiVD[1]  = dPhi_1 /dX1                                     */
/* PhiVD[2]  = dPhi_1 /dX2                                     */
/* PhiVD[3]  = dPhi_1 /dX3                                     */
/* PhiVD[4]  = d^2 Phi_1 / dX1 dX2                             */
/* PhiVD[5]  = d^2 Phi_1 / dX1 dX3                             */
/* PhiVD[6]  = d^2 Phi_1 / dX2 dX3                             */
/* PhiVD[7]  = d^3 Phi_1 / dX1 dX2 dX3                         */
/*                                                             */
/* PhiVD[8]  = Phi_2                                           */
/* PhiVD[9]  = dPhi_2 / dX1                                    */
/* ...                                                         */
/* PhiVD[15] = d^3 Phi_2 / dX1 dX2 dX3                         */
/*                                                             */
/* ... and so on for all components of the Phi function vector.*/
/* (the number of components is specified by the nFun parameter*/
/* to the Interp3D class constructor.)                         */
/*                                                             */
/* note: the function should be thread-safe, as it may be      */
/* called simultaneously by multiple concurrent threads.       */
/***************************************************************/ 
typedef void (*Phi3D)(double X1, double X2, double X3, 
                      void *UserData, double *PhiVD);

/***************************************************************/
/* prototype for a user-supplied function of four variables.   */
/*                                                             */
/* X1, X2, X3, X4 are the inputs to the function (the coords   */
/* of the point at which the function is to be evaluated.)     */
/*                                                             */
/* PhiVD stands for 'phi values and derivatives.' on return,   */
/* this array must be filled in as follows:                    */
/*                                                             */
/* PhiVD[0]  = Phi_1                                           */
/* PhiVD[1]  = dPhi_1 /dX1                                     */
/* PhiVD[2]  = dPhi_1 /dX2                                     */
/* PhiVD[3]  = dPhi_1 /dX3                                     */
/* PhiVD[4]  = dPhi_1 /dX4                                     */
/* PhiVD[5]  = d^2 Phi_1 / dX1 dX2                             */
/* PhiVD[6]  = d^2 Phi_1 / dX1 dX3                             */
/* PhiVD[7]  = d^2 Phi_1 / dX1 dX4                             */
/* PhiVD[8]  = d^2 Phi_1 / dX2 dX3                             */
/* PhiVD[9]  = d^2 Phi_1 / dX2 dX4                             */
/* PhiVD[10] = d^2 Phi_1 / dX3 dX4                             */
/* PhiVD[11] = d^3 Phi_1 / dX1 dX2 dX3                         */
/* PhiVD[12] = d^3 Phi_1 / dX1 dX2 dX4                         */
/* PhiVD[13] = d^3 Phi_1 / dX1 dX3 dX4                         */
/* PhiVD[14] = d^3 Phi_1 / dX2 dX3 dX4                         */
/* PhiVD[15] = d^4 Phi_1 / dX1 dX2 dX3 dX4                     */
/*                                                             */
/* PhiVD[16] = Phi_2                                           */
/* PhiVD[17] = dPhi_2 / dX1                                    */
/* ...                                                         */
/* PhiVD[31] = d^4 Phi_2 / dX1 dX2 dX3 dX4                     */
/*                                                             */
/* ... and so on for all components of the Phi function vector.*/
/* (the number of components is specified by the nFun parameter*/
/* to the Interp4D class constructor.)                         */
/*                                                             */
/* note: the function should be thread-safe, as it may be      */
/* called simultaneously by multiple concurrent threads.       */
/***************************************************************/
typedef void (*Phi4D)(double X1, double X2, double X3, double X4,
                      void *UserData, double *PhiVD);

/***************************************************************/
/* Interp1D is a class which stores an internal grid of        */
/* tabulated function values and provides a class method for   */
/* approximating the value of the function at an arbitrary     */
/* point by interpolating.                                     */
/***************************************************************/
class Interp1D
 { 
   public:

    /*--------------------------------------------------------------*/
    /*- class constructor 1: construct from a user-supplied function*/
    /*- and user-supplied array of points                           */
    /*--------------------------------------------------------------*/
    Interp1D(double *XPoints, int N,
             int nFun, Phi1D PhiFunc, void *UserData,
             int LogLevel=LMDI_LOGLEVEL_TERSE);

    /*--------------------------------------------------------------*/
    /*- class constructor 2: construct from a user-supplied function*/
    /*- and user-supplied uniformly-spaced grid                     */
    /*--------------------------------------------------------------*/
    Interp1D(double XMin, double XMax, int N,
             int nFun, Phi1D PhiFunc, void *UserData,
             int LogLevel=LMDI_LOGLEVEL_TERSE);

    /*--------------------------------------------------------------*/
    /*- body of class constructor for the above two cases ----------*/
    /*--------------------------------------------------------------*/
    void InitInterp1D(Phi1D PhiFunc, void *UserData);

    /*--------------------------------------------------------------*/
    /*- class constructor 3: construct from a user-specified set of */
    /*- X and Y data values.                                        */
    /*-                                                             */
    /*-  YPoints[ nx*nFun + nf ] = f_{nf}( x_{nx} )                 */
    /*--------------------------------------------------------------*/
    Interp1D(double *XPoints, double *YPoints, int N, int nFun,
             int LogLevel=LMDI_LOGLEVEL_TERSE);

    /*--------------------------------------------------------------*/
    /*- class constructor 4: construct from a data file previously  */
    /*- generated by a call to WriteToFile()                        */
    /*--------------------------------------------------------------*/
    Interp1D(const char *FileName);

    /*--------------------------------------------------------------*/
    /*- copy constructor                                            */
    /*--------------------------------------------------------------*/
    Interp1D(Interp1D *Original);

    /*--------------------------------------------------------------*/
    /*- class destructor -------------------------------------------*/
    /*--------------------------------------------------------------*/
    ~Interp1D();

    /*--------------------------------------------------------------*/
    /*- class method that does the interpolation to return an      -*/
    /*- approximate value for the function at the given eval point -*/
    /*--------------------------------------------------------------*/
    bool Evaluate(double X, double *Phi);
    double Evaluate(double X);             // returns Phi[0]

    bool EvaluateVD(double *X0, double *PhiVD);  // get values and 1st derivatives
    bool EvaluateVDD(double *X0, double *PhiVD); // get values and 1st,2nd derivatives

    bool PointInGrid(double X);

    /*--------------------------------------------------------------*/
    /*- class method that writes all internal data to a binary file */
    /*- that may be subsequently used to reconstruct the class     -*/
    /*--------------------------------------------------------------*/
    void WriteToFile(const char *FileName);

    /*----------------------------------------------------------------*/
    /*- internal class data that should be private but i don't bother */
    /*----------------------------------------------------------------*/
    double *XPoints;
    double XMin, DX;    // these are only used if XPoints==0

    int N, nFun;
    int LogLevel;
    int Type;

    double *CTable;
};
    

/***************************************************************/
/* Interp2D is a class which stores an internal 2D grid of     */
/* tabulated function values and provides a class method for   */
/* approximating the value of the function at an arbitrary     */
/* point by interpolating.                                     */
/***************************************************************/
class Interp2D
 { 
   public:

    /*--------------------------------------------------------------*/
    /*- class constructor 1: construct from a user-supplied function*/
    /*- and nonuniform grid                                         */
    /*--------------------------------------------------------------*/
    Interp2D(double *X1Points, int N1, double *X2Points, int N2,
             int nFun, Phi2D PhiFunc=0, void *UserData=0,
             int LogLevel=LMDI_LOGLEVEL_TERSE);

    /*--------------------------------------------------------------*/
    /*- class constructor 2: construct from a user-supplied function*/
    /*- and uniform grid                                            */
    /*--------------------------------------------------------------*/
    Interp2D(double X1Min, double X1Max, int N1, 
             double X2Min, double X2Max, int N2, 
             int nFun, Phi2D PhiFunc=0, void *UserData=0,
             int LogLevel=LMDI_LOGLEVEL_TERSE);

    /*--------------------------------------------------------------*/
    /*- the body of the class constructor for the above two entry  -*/
    /*- points                                                     -*/
    /*--------------------------------------------------------------*/
    void InitInterp2D(Phi2D PhiFunc, void *UserData);

    /*--------------------------------------------------------------*/
    /*- reinitialize an Interp2D object using the same interpolation*/
    /*- grid but with a different function and/or different data    */
    /*--------------------------------------------------------------*/
    void ReInitialize(Phi2D PhiFunc, void *UserData);

    /*--------------------------------------------------------------*/
    /*- class constructor 3: construct from a data file previously  */
    /*- generated by a call to WriteToFile()                        */
    /*--------------------------------------------------------------*/
    Interp2D(const char *FileName);

    /*--------------------------------------------------------------*/
    /*- class destructor -------------------------------------------*/
    /*--------------------------------------------------------------*/
    ~Interp2D();

    /*--------------------------------------------------------------*/
    /*- class method that does the interpolation to return an      -*/
    /*- approximate value for the function at the given eval point -*/
    /*--------------------------------------------------------------*/
    void Evaluate(double X1, double X2, double *Phi);
    void EvaluatePlus(double X1, double X2, double *Phi);
    void EvaluatePlusPlus(double X1, double X2, double *Phi);

    /*--------------------------------------------------------------*/
    /*- class method that writes all internal data to a binary file */
    /*- that may be subsequently used to reconstruct the class     -*/
    /*--------------------------------------------------------------*/
    void WriteToFile(const char *FileName);

    /*----------------------------------------------------------------*/
    /*- internal class data that should be private but i don't bother */
    /*----------------------------------------------------------------*/
    double *X1Points, *X2Points;
    int N1, N2;
    double X1Min, X2Min;   // only used if X1Points==X2Points==0
    double DX1, DX2;

    int nFun; 
    int LogLevel;

    double *CTable;
    
 };

/***************************************************************/
/* Interp3D is a class which stores an internal 3D grid of     */
/* tabulated function values and provides a class method for   */
/* approximating the value of the function at an arbitrary     */
/* point by interpolating.                                     */
/***************************************************************/
class Interp3D
 { 
   public:

    /*--------------------------------------------------------------*/
    /*- user-supplied function, nonuniform grid                     */
    /*--------------------------------------------------------------*/
    Interp3D(double *X1Points, int N1,
             double *X2Points, int N2,
             double *X3Points, int N3,
             int nFun, Phi3D PhiFunc=0, void *UserData=0,
             int LogLevel=LMDI_LOGLEVEL_TERSE);

    /*--------------------------------------------------------------*/
    /*- user-supplied data table, nonuniform grid                   */
    /*--------------------------------------------------------------*/
    Interp3D(double *PhiVDTable,
             double *X1Points, int N1,
             double *X2Points, int N2,
             double *X3Points, int N3,
             int nFun, int LogLevel=LMDI_LOGLEVEL_TERSE);

    /*--------------------------------------------------------------*/
    /*- user-supplied function, uniform grid                        */
    /*--------------------------------------------------------------*/
    Interp3D(double X1Min, double X1Max, int N1,
             double X2Min, double X2Max, int N2,
             double X3Min, double X3Max, int N3,
             int nFun, Phi3D PhiFunc=0, void *UserData=0,
             int LogLevel=LMDI_LOGLEVEL_TERSE);

    /*--------------------------------------------------------------*/
    /*- user-supplied data table, uniform grid                      */
    /*--------------------------------------------------------------*/
    Interp3D(double *PhiVDTable, 
             double X1Min, double X1Max, int N1,
             double X2Min, double X2Max, int N2,
             double X3Min, double X3Max, int N3,
             int nFun, int LogLevel=LMDI_LOGLEVEL_TERSE);

    /*--------------------------------------------------------------*/
    /*- body of class constructor                                   */
    /*--------------------------------------------------------------*/
    void InitInterp3D(
      double pX1Min, double X1Max, double *pX1Points, int pN1,
      double pX2Min, double X2Max, double *pX2Points, int pN2,
      double pX3Min, double X3Max, double *pX3Points, int pN3,
      int pnFun, Phi3D PhiFunc, void *UserData, double *PhiVDTable,
      int pLogLevel);

    /*--------------------------------------------------------------*/
    /*- reinitialize an Interp3D object using the same interpolation*/
    /*- grid but with a different function and/or different data    */
    /*--------------------------------------------------------------*/
    void ReInitialize(Phi3D PhiFunc, void *UserData, double *PhiVDTable);

    /*--------------------------------------------------------------*/
    /*- class constructor 3: construct from a data file previously  */
    /*- generated by a call to WriteToFile()                        */
    /*--------------------------------------------------------------*/
    Interp3D(const char *FileName);

    /*--------------------------------------------------------------*/
    /*- class destructor -------------------------------------------*/
    /*--------------------------------------------------------------*/
    ~Interp3D();

    /*--------------------------------------------------------------*/
    /*- class method that does the interpolation to return an      -*/
    /*- approximate value for the function at the given eval point -*/
    /*--------------------------------------------------------------*/
    void Evaluate(double X1, double X2, double X3, double *Phi);
    void EvaluatePlus(double X1, double X2, double X3, double *PhiVD);
    void EvaluatePlusPlus(double X1, double X2, double X3, double *PhiVD);

    /*--------------------------------------------------------------*/
    /*- return true if point lies in the interior or on the boundary*/
    /*- of the interpolation grid; false otherwise.                 */
    /*- on return, pn (if non-null) has the indices of the grid cell*/
    /*- and pxBar (if non-null) has the normalized point coordinates*/
    /*- within that cell (0<=pXbar[i]<=1).                          */
    /*--------------------------------------------------------------*/
    bool PointInGrid(double X1, double X2, double X3, 
                     int *pn=0, double *pxBar=0);

    /*--------------------------------------------------------------*/
    /*- class method that writes all internal data to a binary file */
    /*- that may be subsequently used to reconstruct the class     -*/
    /*--------------------------------------------------------------*/
    void WriteToFile(const char *FileName);

    /*----------------------------------------------------------------*/
    /*- internal class data that should be private but i don't bother */
    /*----------------------------------------------------------------*/
    double *X1Points, *X2Points, *X3Points;
    int N1, N2, N3;
    double X1Min, X2Min, X3Min;
    double DX1, DX2, DX3;
    int nFun; 
    int LogLevel;

    double *CTable;
    
 };

/***************************************************************/
/* Interp4D is a class which stores an internal 4D grid of     */
/* tabulated function values and provides a class method for   */
/* approximating the value of the function at an arbitrary     */
/* point by interpolating.                                     */
/***************************************************************/
class Interp4D
 { 
   public:

    /*--------------------------------------------------------------*/
    /*- class constructor 1: construct from a user-supplied function*/
    /*- and nonuniform grid                                         */
    /*--------------------------------------------------------------*/
    Interp4D(double *X1Points, int N1,
             double *X2Points, int N2,
             double *X3Points, int N3,
             double *X4Points, int N4,
             int nFun, Phi4D PhiFunc, void *UserData);

    /*--------------------------------------------------------------*/
    /*- class constructor 2: construct from a user-supplied function*/
    /*- and uniform grid                                            */
    /*--------------------------------------------------------------*/
    Interp4D(double pX1Min, double X1Max, int pN1,
             double pX2Min, double X2Max, int pN2,
             double pX3Min, double X3Max, int pN3,
             double pX4Min, double X4Max, int pN4,
             int pnFun, Phi4D PhiFunc, void *UserData);

    /*--------------------------------------------------------------*/
    /*- body of class constructor for the previous two entry points-*/
    /*--------------------------------------------------------------*/
    void InitInterp4D(Phi4D PhiFunc, void *UserData);

    /*--------------------------------------------------------------*/
    /*- reinitialize an Interp4D object using the same interpolation*/
    /*- grid but with a different function and/or different data    */
    /*--------------------------------------------------------------*/
    void ReInitialize(Phi4D PhiFunc, void *UserData);

    /*--------------------------------------------------------------*/
    /*- class constructor 3: construct from a data file previously  */
    /*- generated by a call to WriteToFile()                        */
    /*--------------------------------------------------------------*/
    Interp4D(const char *FileName);

    /*--------------------------------------------------------------*/
    /*- class destructor -------------------------------------------*/
    /*--------------------------------------------------------------*/
    ~Interp4D();

    /*--------------------------------------------------------------*/
    /*- class method that does the interpolation to return an      -*/
    /*- approximate value for the function at the given eval point -*/
    /*--------------------------------------------------------------*/
    void Evaluate(double X1, double X2, double X3, double X4, double *Phi);

    /*--------------------------------------------------------------*/
    /*- class method that writes all internal data to a binary file */
    /*- that may be subsequently used to reconstruct the class     -*/
    /*--------------------------------------------------------------*/
    void WriteToFile(const char *FileName);

    /*----------------------------------------------------------------*/
    /*- internal class data that should be private but i don't bother */
    /*----------------------------------------------------------------*/
    double *X1Points, *X2Points, *X3Points, *X4Points;
    int N1, N2, N3, N4; 
    double X1Min, X2Min, X3Min, X4Min;
    double DX1, DX2, DX3, DX4;
    int nFun; 

    double *CTable;
    
 };

/***************************************************************/
/* 3/2018: InterpND is a single general-purpose interpolator   */
/*         that supports arbitrary dimensions and should       */
/*         replace all the above classes.                      */
/***************************************************************/
#ifndef iVec
  typedef std::vector<int>    iVec;
#endif
#ifndef dVec
  typedef std::vector<double> dVec;
#endif
bool Increment(iVec &n, iVec &N);
double Monomial(dVec xVec, iVec pVec);
#define LOOP_OVER_IVECS(np, nVec, NVec)     \
   iVec nVec(NVec.size());                  \
   int np=0;                                \
   for(bool np ## _Done=false; !np ## _Done; np++, np ## _Done=Increment(nVec,NVec))
#define MAXDIM 4

typedef void (*PhiVDFunc)(double *X, void *UserData, double *PhiVD);

class InterpND
 { 
   public:

    /*--------------------------------------------------------------*/
    /*- user-supplied function, uniform grid                        */
    /*--------------------------------------------------------------*/
    InterpND(PhiVDFunc Phi, void *UserData, int NF,
             dVec X0Min, dVec X0Max, iVec N0Vec, bool Verbose=false);

    /*--------------------------------------------------------------*/
    /*- user-supplied function, user-supplied grids in each dimension*/
    /*--------------------------------------------------------------*/
    InterpND(PhiVDFunc UserFunc, void *UserData, int NF, 
             vector<dVec> &X0Grids, bool Verbose=false);

    /*--------------------------------------------------------------*/
    /*- user-supplied function with requested error tolerance       */
    /*- (non-uniform grids determined automatically)                */
    /*--------------------------------------------------------------*/
    InterpND(PhiVDFunc UserFunc, void *UserData, int NF,
             dVec XMin, dVec XMax, double MaxRelError=1.0e-4, bool Verbose=false);

    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
    /*- class destructor -------------------------------------------*/
    /*--------------------------------------------------------------*/
    ~InterpND();

    /*--------------------------------------------------------------*/
    /*- return true/false if X0 does/does not lie within the grid   */
    /*--------------------------------------------------------------*/
    bool PointInGrid(double *X0, int *nVec=0, double *XBarVec=0);

    /*--------------------------------------------------------------*/
    /*- perform interpolation to return approximate values for the  */
    /*- interpolated functions at the given eval point             -*/
    /*--------------------------------------------------------------*/
    bool Evaluate(double *X0, double *Phi);
    bool EvaluateVD(double *X0, double *PhiVD);
    bool EvaluateVDD(double *X0, double *PhiVD);

    double PlotInterpolationError(PhiVDFunc UserFunc, void *UserData, char *OutFileName, bool CentersOnly=false);

    /*--------------------------------------------------------------*/
    /*- return true if point lies in the interior or on the boundary*/
    /*- of the interpolation grid; false otherwise.                 */
    /*- on return, pn (if non-null) has the indices of the grid cell*/
    /*- and pxBar (if non-null) has the normalized point coordinates*/
    /*- within that cell (0<=pXbar[i]<=1).                          */
    /*--------------------------------------------------------------*/
    //bool PointInGrid(double *X, int *pn=0, double *pxBar=0);

    /*--------------------------------------------------------------*/
    /*- class method that writes all internal data to a binary file */
    /*- that may be subsequently used to reconstruct the class     -*/
    /*--------------------------------------------------------------*/
    //void WriteToFile(const char *FileName);

    /*----------------------------------------------------------------*/
    /*- private  class methods ---------------------------------------*/
    /*----------------------------------------------------------------*/
//private:

    // convert back and forth between the index of a grid
    // point (0...NumPoints-1) and its integer coordinates
    iVec GetPoint(size_t nPoint);
    size_t GetPointIndex(iVec nVec);

    // convert back and forth between the index of a grid
    // cell (0...NumCell-1) and the integer coordinates
    // of its lower-left corner
    iVec GetCell(size_t nVec);
    size_t GetCellIndex(iVec nVec);

    // convert grid-point indices to cartesian coordinates
    dVec n2X(iVec nVec);
    dVec n2X0(iVec nVec);

    // convert between reduced and full coordinates of an evaluation point
    dVec X02X(dVec X0Vec);
    dVec X2X0(dVec XVec);

    // utility methods for handling internal storage tables
    size_t GetCTableOffset(int nCell, int nFun);
    size_t GetCTableOffset(iVec nVec, int nFun);
    size_t GetPhiVDTableOffset(int nPoint, int nFun);
    size_t GetPhiVDTableOffset(iVec nVec, iVec tauVec, int nFun);
    void FillPhiVDBuffer(double *PhiVD, double *PhiVD0);

    // constructor helper method
    void Initialize(PhiVDFunc UserFunc, void *UserData, bool Verbose=false);

    /*----------------------------------------------------------------*/
    /*- internal class data that should be private but i don't bother */
    /*----------------------------------------------------------------*/
    int D;                   // number of dimensions
    iVec NVec;               // # grid points in each dimension
    int NF;                  // number of functions
    dVec XMin, DX;           // grid points (uniform spacing)
    vector<dVec> XGrids;     // grid points (non-uniform spacing)

    iVec CellStride;         // strides for computing grid-cell index
    iVec PointStride;        // strides for computing grid-point index
    int NVD;                 // # function values, derivatives per grid point
    int NCoeff;              // # coefficients per grid cell

    dVec FixedCoordinates;
    int D0;
    int NVD0;

    double *CTable;          // polynomial coefficients
 };

/***************************************************************/
/* routines for automatically determining optimal grid spacing */
/* required for given error tolerances                         */
/***************************************************************/
double GetInterpolationError(PhiVDFunc UserFunc, void *UserData, int NF,
                             dVec X0Vec, dVec DeltaVec, double AbsTol=1.0e-8, char *LogFile=0);

double GetInterpolationError(PhiVDFunc UserFunc, void *UserData, int NF,
                             int d0, double X, double Delta,
                             dVec X0Vec, dVec dX0Vec, double AbsTol=1.0e-8, char *LogFile=0);

dVec GetXGrid(PhiVDFunc UserFunc, void *UserData, int NF, dVec XMin, dVec XMax, int d,
              double DesiredMaxRE);

/***************************************************************/
/* non-member function for binary searching                    */
/***************************************************************/
void BinSearch(double X, double *XPoints, int N);

bool FindInterval(double X, double *XPoints, int N, double XMin, double DX,
                  int *n, double *XBar);

/***************************************************************/
/* a version of fread() with simple error-checking *************/
/***************************************************************/
void freadEC(void *p, size_t size, size_t nmemb, FILE *f, const char *FileName); 

#endif // #ifndef LIBMDINTERP_H
