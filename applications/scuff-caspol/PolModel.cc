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
 * PolModel.cc   -- implementation of a simple class for describing
 *               -- the frequency-dependent polarizability of atoms
 *               -- and molecules
 *
 * homer reid    -- 1/2012
 *
 */

#include "scuff-caspol.h"
#include <libMatProp.h>
#include <libhmat.h>

/***************************************************************/
/* Data tables from Derevianko et. al., "Electric dipole       */
/* polarizabilities at imaginary frequencies for hydrogen, the */
/* alkali–metal, alkaline–earth, and noble gas atoms,"         */
/* Atomic Data and Nuclear Data Tables 96 (2010) 323–331.      */
/* Note: DPB stands for "Derevianko, Porsev, Babb."            */
/***************************************************************/
#define NUMPOINTS_DPB 51

// original table from DPB paper with frequencies in atomic units
#if 0
double XiPoints_DPB[NUMPOINTS_DPB] =
{ 0.0, 
  0.00178065, 0.00937463, 0.0230068, 0.042631, 0.0681817,
  0.0995817, 0.136748, 0.179602, 0.228071, 0.282107,
  0.341686, 0.406823, 0.477579, 0.554072, 0.636488,
  0.725091, 0.820235, 0.922382, 1.03212, 1.15017,
  1.27743, 1.41501, 1.56425, 1.72679, 1.90461,
  2.10016, 2.31644, 2.55713, 2.82683, 3.13128,
  3.47776, 3.87553, 4.3366, 4.87665, 5.51655,
  6.28448, 7.21927, 8.37558, 9.83228, 11.7066,
  14.179, 17.5384, 22.2715, 29.2508, 40.168,
  58.6667, 93.8285, 173.862, 426.684, 2246.37
};
#endif

// table from DPB paper with frequencies in SCUFF units.
// note: in SCUFF units, 'Xi==1'  <--> Xi=3e14 rad /sec
//       in atomic units, 'Xi==1' <--> Xi=4.13393964448119e+16 rad/sec
// the conversion factor is 137.79798814937303
#define ATOMIC_TO_SCUFF 137.79798814937303
// 20140731 note: many thanks to Tom Judd and Tino Sering for
// helping to correct a mistake in this calculation!

double XiPoints_DPB[NUMPOINTS_DPB] =
{ 0.00000000e+00,   2.45369988e-01,   1.29180515e+00,
  3.17029075e+00,   5.87446603e+00,   9.39530109e+00,
  1.37221579e+01,   1.88435993e+01,   2.47487943e+01,
  3.14277250e+01,   3.88737770e+01,   4.70836434e+01,
  5.60593909e+01,   6.58094254e+01,   7.63500069e+01,
  8.77067659e+01,   9.99160810e+01,   1.13026733e+02,
  1.27102384e+02,   1.42224060e+02,   1.58491112e+02,
  1.76027284e+02,   1.94985531e+02,   2.15550503e+02,
  2.37948188e+02,   2.62451426e+02,   2.89397823e+02,
  3.19200772e+02,   3.52367369e+02,   3.89531487e+02,
  4.31484084e+02,   4.79228331e+02,   5.34040237e+02,
  5.97574755e+02,   6.71992559e+02,   7.60169492e+02,
  8.65988701e+02,   9.94800882e+02,   1.15413807e+03,
  1.35486840e+03,   1.61314593e+03,   1.95383767e+03,
  2.41675624e+03,   3.06896789e+03,   4.03070139e+03,
  5.53506959e+03,   8.08415323e+03,   1.29293785e+04,
  2.39578338e+04,   5.87961968e+04,   3.09545267e+05
};

double wPoints_DPB[NUMPOINTS_DPB] =
{ 0.0, 
  0.00456886, 0.0106185, 0.0166378, 0.0225996, 0.0284889, 
  0.0342971, 0.0400227, 0.0456718, 0.0512592, 0.0568084, 
  0.0623519, 0.0679318, 0.0735998, 0.0794176, 0.085458, 
  0.0918056, 0.0985589, 0.105832, 0.113758, 0.122493, 
  0.132222, 0.143164, 0.155587, 0.169813, 0.18624, 
  0.205361, 0.227799, 0.254343, 0.286005, 0.324106, 
  0.370381, 0.427154, 0.497572, 0.585976, 0.698465, 
  0.843785, 1.03477, 1.29076, 1.64181, 2.13626, 
  2.85525, 3.94176, 5.66354, 8.56095, 13.8343, 
  24.5132, 49.7405, 125.732, 483.298, 5763.83
};

double PolPoints_DPB_Hydrogen[NUMPOINTS_DPB] = 
{ 4.5, 
  4.4997,    4.4974,    4.4857,    4.452,     4.3798,
  4.252,     4.0561,    3.7895,    3.4624,    3.0959,
  2.7157,    2.3453,    2.0018,    1.6948,    1.4274,
  1.1985,    1.0049,    8.4223e-1, 7.06e-1,   5.9205e-1, 
  4.9668e-1, 4.1676e-1, 3.4963e-1, 2.9313e-1, 2.4546e-1, 
  2.0515e-1, 1.7102e-1, 1.4207e-1, 1.1751e-1, 9.6682e-2, 
  7.9029e-2, 6.4103e-2, 5.1524e-2, 4.0971e-2, 3.2172e-2, 
  2.4893e-2, 1.8931e-2, 1.4108e-2, 1.0263e-2, 7.2545e-3, 
  4.9534e-3, 3.2418e-3, 2.0123e-3, 1.1674e-3, 6.1937e-4, 
  2.9045e-4, 1.1357e-4, 3.3079e-5, 5.4925e-6, 1.9816e-7 
 };
double LargeXiCoefficient_DPB_Hydrogen = 1.33367e-06;

double PolPoints_DPB_Lithium[NUMPOINTS_DPB] = 
{ 164.0, 
  1.6386e+2, 1.6094e+2, 1.4727e+2, 1.181e+2,  8.2483e+1, 
  5.3094e+1, 3.355e+1,  2.1584e+1, 1.4337e+1, 9.8574, 
  7.0022,    5.1216,    3.844,     2.9512,    2.3111, 
  1.8417,    1.4902,    1.2221,    1.0139,    8.4974e-1, 
  7.1822e-1, 6.1136e-1, 5.2334e-1, 4.4989e-1, 3.8785e-1, 
  3.3485e-1, 2.891e-1,  2.4925e-1, 2.1428e-1, 1.834e-1, 
  1.5602e-1, 1.3171e-1, 1.1013e-1, 9.1031e-2, 7.4209e-2, 
  5.9517e-2, 4.6826e-2, 3.6019e-2, 2.6979e-2, 1.9581e-2, 
  1.3686e-2, 9.1389e-3, 5.7691e-3, 3.3923e-3, 1.8183e-3, 
  8.5884e-4, 3.3735e-4, 9.8504e-5, 1.6372e-5, 5.9081e-7 
};
double LargeXiCoefficient_DPB_Lithium = 3.96234e-06;

double PolPoints_DPB_Sodium[NUMPOINTS_DPB] = 
{ 162.6, 
  1.625e+2, 1.6025e+2, 1.4948e+2, 1.2503e+2, 9.2082e+1, 
  6.1986e+1, 4.0336e+1, 2.6414e+1, 1.7756e+1, 1.2345e+1, 
  8.8866, 6.6142, 5.0767, 4.0061, 3.2398, 
  2.6767, 2.2523, 1.9245, 1.6651, 1.4553, 
  1.2819, 1.1356, 1.0098, 8.9998e-1, 8.0261e-1, 
  7.1525e-1, 6.3608e-1, 5.6376e-1, 4.9735e-1, 4.3614e-1, 
  3.7966e-1, 3.2757e-1, 2.7968e-1, 2.3587e-1, 1.9607e-1, 
  1.6028e-1, 1.2848e-1, 1.0066e-1, 7.6781e-2, 5.6743e-2, 
  4.0389e-2, 2.7483e-2, 1.7703e-2, 1.0652e-2, 5.8682e-3, 
  2.8645e-3, 1.1666e-3, 3.5183e-4, 5.9695e-5, 2.1693e-6 
};
double LargeXiCoefficient_DPB_Sodium = 1.37391e-05;

double PolPoints_DPB_Potassium[NUMPOINTS_DPB] = 
{ 290.2, 
  2.8996e+2, 2.8329e+2, 2.5306e+2, 1.934e+2, 1.2835e+2, 
  8.0193e+1, 5.06e+1, 3.3376e+1, 2.326e+1, 1.7106e+1, 
  1.3188e+1, 1.0567e+1, 8.7219, 7.3564, 6.2984, 
  5.4454, 4.7358, 4.1313, 3.6081, 3.1504, 
  2.7476, 2.3919, 2.0774, 1.7993, 1.5535, 
  1.3367, 1.1456, 9.7766e-1, 8.3028e-1, 7.013e-1, 
  5.8879e-1, 4.9097e-1, 4.0628e-1, 3.333e-1, 2.7077e-1, 
  2.1753e-1, 1.7253e-1, 1.3479e-1, 1.0341e-1, 7.7582e-2, 
  5.6563e-2, 3.9715e-2, 2.6512e-2, 1.6525e-2, 9.3754e-3, 
  4.6648e-3, 1.919e-3, 5.8509e-4, 1.012e-4, 3.7292e-6
};
double LargeXiCoefficient_DPB_Potassium = 2.26217e-05;

double PolPoints_DPB_Rubidium[NUMPOINTS_DPB] = 
{ 318.6, 
  3.1832e+2, 3.1076e+2, 2.7667e+2, 2.1035e+2, 1.3939e+2, 
  8.7721e+1, 5.627e+1, 3.802e+1, 2.7268e+1, 2.066e+1, 
  1.6373e+1, 1.3424e+1, 1.1272e+1, 9.6174, 8.2871, 
  7.1818, 6.2421, 5.4316, 4.7265, 4.1101, 
  3.5702, 3.097, 2.6821, 2.3186, 2.0004, 
  1.722, 1.4788, 1.2664, 1.0812, 9.1991e-1, 
  7.7965e-1, 6.5784e-1, 5.522e-1, 4.6071e-1, 3.8156e-1, 
  3.1318e-1, 2.5416e-1, 2.0332e-1, 1.5966e-1, 1.2241e-1, 
  9.0954e-2, 6.4885e-2, 4.3891e-2, 2.7695e-2, 1.5939e-2, 
  8.0975e-3, 3.4329e-3, 1.0815e-3, 1.9126e-4, 7.2051e-6
};
double LargeXiCoefficient_DPB_Rubidium = 4.05829e-05;

double PolPoints_DPB_Cesium[NUMPOINTS_DPB] = 
{ 399.8, 
  3.9938e+2, 3.88e+2, 3.3819e+2, 2.4771e+2, 1.5921e+2, 
  9.9455e+1, 6.4749e+1, 4.5037e+1, 3.343e+1, 2.6166e+1, 
  2.1283e+1, 1.776e+1, 1.5058e+1, 1.2886e+1, 1.1086e+1, 
  9.5655, 8.2676, 7.1542, 6.1966, 5.372, 
  4.6611, 4.0476, 3.5171, 3.0576, 2.6583, 
  2.3104, 2.0061, 1.739, 1.5036, 1.2953, 
  1.1104, 9.4575e-1, 7.9902e-1, 6.683e-1, 5.5221e-1, 
  4.4973e-1, 3.6009e-1, 2.827e-1, 2.1698e-1, 1.6227e-1, 
  1.1779e-1, 8.2549e-2, 5.5437e-2, 3.5233e-2, 2.0733e-2, 
  1.0867e-2, 4.7421e-3, 1.5266e-3, 2.7665e-4, 1.0596e-5
};
double LargeXiCoefficient_DPB_Cesium = 5.61692e-05;

double PolPoints_DPB_Francium[NUMPOINTS_DPB] = 
{ 317.8, 
  3.1752e+2, 3.1082e+2, 2.8035e+2, 2.1973e+2, 1.5244e+2, 
  1.0137e+2, 6.9056e+1, 4.9579e+1, 3.7583e+1, 2.9784e+1, 
  2.4368e+1, 2.0363e+1, 1.7246e+1, 1.4729e+1, 1.2648e+1, 
  1.0902e+1, 9.4258, 8.17, 7.098, 6.18, 
  5.3912, 4.711, 4.1221, 3.61, 3.1623, 
  2.7689, 2.4215, 2.113, 1.8379, 1.5914, 
  1.3701, 1.1712, 9.9237e-1, 8.3219e-1, 6.8947e-1, 
  5.6332e-1, 4.5304e-1, 3.579e-1, 2.7711e-1, 2.0971e-1, 
  1.5452e-1, 1.1023e-1, 7.5407e-2, 4.8728e-2, 2.9029e-2, 
  1.535e-2, 6.7792e-3, 2.2355e-3, 4.1737e-4, 1.6446e-5
};
double LargeXiCoefficient_DPB_Francium =  8.04767e-5;

/***************************************************************/
/* PolModel class constructor entry point                      */
/***************************************************************/
PolModel::PolModel(char *String, int Which)
{
  if (Which==PM_BUILTIN)
   InitPolModel_BI(String);
  else
   InitPolModel_UD(String);
   
}

/***************************************************************/
/* PolModel class constructor for built-in polarizabilities.   */
/***************************************************************/
void PolModel::InitPolModel_BI(char *Atom)
{
  ErrMsg=0;

  Name  = strdup(Atom);

  NumPoints = NUMPOINTS_DPB;
  XiPoints = XiPoints_DPB;

  if ( !strcasecmp(Atom,"H")  || !strcasecmp(Atom,"Hydrogen") )
   { PolPoints = PolPoints_DPB_Hydrogen;
     LargeXiCoefficient = LargeXiCoefficient_DPB_Hydrogen;
   }
  else if ( !strcasecmp(Atom,"Li") || !strcasecmp(Atom,"Lithium") )
   {
     PolPoints = PolPoints_DPB_Lithium;
     LargeXiCoefficient = LargeXiCoefficient_DPB_Lithium;
   }
  else if ( !strcasecmp(Atom,"Na") || !strcasecmp(Atom,"Sodium") )
   {
     PolPoints = PolPoints_DPB_Sodium;
     LargeXiCoefficient = LargeXiCoefficient_DPB_Sodium;
   }
  else if ( !strcasecmp(Atom,"K")  || !strcasecmp(Atom,"Potassium") )
   {
     PolPoints = PolPoints_DPB_Potassium;
     LargeXiCoefficient = LargeXiCoefficient_DPB_Potassium;
   }
  else if ( !strcasecmp(Atom,"Rb") || !strcasecmp(Atom,"Rubidium") )
   {
     PolPoints = PolPoints_DPB_Rubidium;
     LargeXiCoefficient = LargeXiCoefficient_DPB_Rubidium;
   }
  else if ( !strcasecmp(Atom,"Cs") || !strcasecmp(Atom,"Cesium") )
   {
     PolPoints = PolPoints_DPB_Cesium;
     LargeXiCoefficient = LargeXiCoefficient_DPB_Cesium;
   }
  else if ( !strcasecmp(Atom,"Fr") || !strcasecmp(Atom,"Francium") )
   {
     PolPoints = PolPoints_DPB_Francium;
     LargeXiCoefficient = LargeXiCoefficient_DPB_Francium;
   }
  else
   { ErrMsg=vstrdup("unknown atom type %s",Atom);
     return;
   };

   // initialize MP, otherwise will be used in PolModel::GetPolarizability w/o assignment
   MP = 0;
   // initialize the interpolator
   PolInterp = new Interp1D(XiPoints, PolPoints, NumPoints, 1);

}

/***************************************************************/
/* PolModel class constructor for user-defined polarizabilities*/
/***************************************************************/
void PolModel::InitPolModel_UD(char *FileName)
{
  ErrMsg=0;
  LargeXiCoefficient = 0.0;
  Name = strdup(GetFileBase(FileName));

  MP=new MatProp(FileName);
  if (MP->ErrMsg==0)
   { PolInterp=0;
     MatProp::SetLengthUnit(1.0e-6);
     return;
   }
  else
   { delete MP;
     MP=0;
   };

  HMatrix *PolData = new HMatrix(FileName,LHM_TEXT,"--ncol 2 --strict");
  if (PolData->ErrMsg)
   ErrExit(PolData->ErrMsg);

  NumPoints = PolData->NR;
  XiPoints  = (double *)malloc(NumPoints * sizeof(double));
  PolPoints = (double *)malloc(NumPoints * sizeof(double));

  for(int np=0; np<NumPoints; np++)
   { XiPoints[np]  = ATOMIC_TO_SCUFF * PolData->GetEntryD(np,0);
     PolPoints[np] = PolData->GetEntryD(np,1);
   };

  PolInterp = new Interp1D(XiPoints, PolPoints, NumPoints, 1);
 
}

/***************************************************************/
/* routine to compute polarizability tensor at a given         */
/* frequency.                                                  */
/*                                                             */
/* on input, Xi is the imaginary frequency in units of         */
/* 3e14 rad/sec.                                               */
/*                                                             */
/* Alpha points to a preallocated 3x3 real-valued HMatrix.     */
/*                                                             */
/* On return, Alpha is filled in with the components of the    */
/* polarizability at this imaginary frequency.                 */
/***************************************************************/
void PolModel::GetPolarizability(double Xi, HMatrix *Alpha)
 {
   double AlphaDiag;

   if (MP)
    AlphaDiag = real(MP->GetEps(cdouble(0.0,Xi)));
   else
    PolInterp->Evaluate(Xi, &AlphaDiag);

   Alpha->Zero();
   Alpha->SetEntry(0,0,AlphaDiag);
   Alpha->SetEntry(1,1,AlphaDiag);
   Alpha->SetEntry(2,2,AlphaDiag);

 }
