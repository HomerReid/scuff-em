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
 * tlibMatProp.cc -- demonstration program for libMatProp that plots 
 *                   epsilon and mu as functions of real and imaginary
 *                   frequencies for a single material in a libMatProp 
 *                   database
 * 
 * homer reid     -- 12/2008
 * 
 * --------------------------------------------------------------
 * 
 * usage:
 *  
 *   tlibMatProp MaterialName [options]
 *  
 * options:
 *  
 *  --MinFreq    xx  (specify lower limit of frequency range to plot)
 *  --MaxFreq    xx  (specify upper limit of frequency range to plot)
 *  --LengthUnit xx  (specify libMatProp length unit)
 *  
 * --------------------------------------------------------------
 *  
 * if the program runs successfully, the result will be a data file
 * named 'MaterialName.out.' 
 * each line in this data file will have 6 columns: 
 *  
 * w  |  real(eps(w))  |  imag(eps(w))  |  mu(w)  |  eps(iw)  | mu(iw)
 *  
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include "libMatProp.h"

#define NUMPTS 100
#define REAL_FREQ 1
#define IMAG_FREQ 0

#define MINFREQ 1.0e6
#define MAXFREQ 1.0e16

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  MatProp *MP;
  cdouble Eps; 
  cdouble Mu;
  double Freq, MinFreq, MaxFreq, FMult, LengthUnit;
  int nArg, nConv;
  FILE *f;

  SetLogFileName("tlibMatProp.log");

  if (argc<2)
   { fprintf(stderr,"usage: %s MaterialName [options] \n",argv[0]);
     fprintf(stderr,"                                 \n");
     fprintf(stderr," options:                        \n");
     fprintf(stderr,"  --LengthUnit xx                \n");
     fprintf(stderr,"  --MinFreq xx                   \n");
     fprintf(stderr,"  --MaxFreq xx                   \n");
     fprintf(stderr,"                                 \n");
     exit(1);
   };

  /***************************************************************/
  /* try to create a MatProp for the specified material          */
  /***************************************************************/
  MP=new MatProp(argv[1]);
  if (MP->ErrMsg)
   ErrExit(MP->ErrMsg);

  /***************************************************************/
  /* process other command-line options **************************/
  /***************************************************************/
  MinFreq=MINFREQ;
  MaxFreq=MAXFREQ;
  for(nArg=2; nArg<argc; nArg++)
   { 
     if ( !StrCaseCmp(argv[nArg],"--MinFreq") )
      { nArg++;
        if (nArg==argc) 
         ErrExit("%s requires an option",argv[nArg-1]);
        nConv=sscanf(argv[nArg],"%le",&MinFreq);
        if ( nConv!=1 || MinFreq<0.0 )
         ErrExit("invalid value %s specified for %s",argv[nArg],argv[nArg-1]);
        printf("Using min freq=%e.\n",MinFreq);
      }
     else if ( !StrCaseCmp(argv[nArg],"--MaxFreq") )
      { nArg++;
        if (nArg==argc) 
         ErrExit("%s requires an option",argv[nArg-1]);
        nConv=sscanf(argv[nArg],"%le",&MaxFreq);
        if ( nConv!=1 || MaxFreq<0.0 )
         ErrExit("invalid value %s specified for %s",argv[nArg],argv[nArg-1]);
        printf("Using max freq=%e.\n",MaxFreq);
      }
     else if ( !StrCaseCmp(argv[nArg],"--LengthUnit") )
      { nArg++;
        if (nArg==argc) 
         ErrExit("%s requires an option",argv[nArg-1]);
        nConv=sscanf(argv[nArg],"%le",&LengthUnit);
        if ( nConv!=1 || LengthUnit <0.0 )
         ErrExit("invalid value %s specified for %s",argv[nArg],argv[nArg-1]);
        MatProp::SetLengthUnit(LengthUnit);
        printf("Using length unit %e.\n",LengthUnit);
      }
     else 
      ErrExit("invalid command-line option %s",argv[nArg]);
   };

  /***************************************************************/
  /* plot eps and mu over the specified frequency range **********/
  /***************************************************************/
  FMult=pow(MaxFreq/MinFreq,1.0/((double)NUMPTS ));  
  f=vfopen("%s.epsmu","w",MP->Name);
  for(Freq=MinFreq; Freq<=MaxFreq; Freq*=FMult)
   { 
      fprintf(f,"%e ",Freq);

      MP->GetEpsMu(Freq, REAL_FREQ, &Eps, &Mu);
      fprintf(f,"%+15.12e %+15.12e %+15.12e %+15.12e",real(Eps), imag(Eps), 
	      real(Mu), imag(Mu));

      MP->GetEpsMu(Freq, IMAG_FREQ, &Eps, &Mu);
      fprintf(f," %+15.12e %+15.12e \n",real(Eps), real(Mu));

   };
  fclose(f);

  /***************************************************************/
  /* generate gnuplot plots **************************************/
  /***************************************************************/
  f=popen("gnuplot -persist","w");
  if (f)
   { 
     fprintf(f,"set logscale x\n");

     fprintf(f,"set title 'Eps,Mu for %s along the real frequency axis\n",MP->Name);
     fprintf(f,"set xlabel 'Real(Omega)'\n");
     fprintf(f,"plot '%s.epsmu' u 1:2 t 'Re(Eps)' w lp,",MP->Name); 
     fprintf(f,"     '%s.epsmu' u 1:3 t 'Im(Eps)' w lp,",MP->Name); 
     fprintf(f,"     '%s.epsmu' u 1:4 t 'Re(Mu)' w lp,",MP->Name); 
     fprintf(f,"     '%s.epsmu' u 1:5 t 'Im(Mu)' w lp\n",MP->Name); 

     fprintf(f,"set title 'Eps,Mu for %s along the imaginary frequency axis\n",MP->Name);
     fprintf(f,"set terminal x11 2\n");
     fprintf(f,"set xlabel 'Imag(Omega)'\n");
     fprintf(f,"plot '%s.epsmu' u 1:6 t 'Eps' w lp,",MP->Name); 
     fprintf(f,"     '%s.epsmu' u 1:7 t 'Mu' w lp\n",MP->Name); 

     pclose(f);
   };

  /***************************************************************/
  /* that was fun, thanks   **************************************/
  /***************************************************************/
  printf("Epsilon-Mu vs. frequency written to file %s.epsmu.\n",MP->Name); 
  printf("Thank you for your support.\n");
   
}
