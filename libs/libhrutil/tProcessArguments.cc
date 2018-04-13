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
 * tProcessArguments.cc -- a demonstration of how to use the command-line- 
 *                      -- argument-processing functionality of libhrutil
 *
 * homer reid           -- 3/2008
 */

#include <stdio.h>
#include <stdlib.h>

#include "libhrutil.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  int lMax, RealFreq;
  double L, H, W, Frequency;
  char *Name, *Name2;

  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  ArgStruct ASArray[]=
   { {"L",         PA_DOUBLE, (void *)&L,         "0.1",  "L parameter"},
     {"H",         PA_DOUBLE, (void *)&H,         "0.1",  "H parameter"},
     {"W",         PA_DOUBLE, (void *)&W,         "0.05", "W parameter"}, 
     {"Frequency", PA_DOUBLE, (void *)&Frequency, "0.1",  "frequency"}, 
     {"RF",        PA_BOOL,   (void *)&RealFreq,  "1",    "1=real freq, 0=imag freq"},
     {"lMax",      PA_INT,    (void *)&lMax,      "11",   "1=real freq, 0=imag freq"},
     {"Name",      PA_STRING, (void *)&Name,      "Howd", "name of operation"},
     {"Name2",     PA_STRING, (void *)&Name2,     0,      "name of operation, 2"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);
  if (Name2==0)
   ASUsage(argv[0],ASArray,"--Name2 option is mandatory");


  printf("%-15s  =   %+12.9e\n","L",L);
  printf("%-15s  =   %+12.9e\n","H",H);
  printf("%-15s  =   %+12.9e\n","W",W);
  printf("%-15s  =   %+12.9e\n","Frequency",Frequency);
  printf("%-15s  =   %i     \n","RF",RealFreq);
  printf("%-15s  =   %i     \n","lMax",lMax);
  printf("%-15s  =   %s     \n","Name",Name);
  printf("%-15s  =   %s     \n","Name2",Name2==0 ? "null" : Name2);
 

} 
