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
 * tProcessOptions.cc -- a demonstration of how to use the options-      
 *                    -- processing functionality of libhrutil
 *
 * homer reid         -- 1/2012
 */

#include <stdio.h>
#include <stdlib.h>

#include "libhrutil.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  /*--------------------------------------------------------------*/
  /*- process options  -------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble pwPol[3*3]; int npwPol;
  double pwDir[3*3];  int npwDir;
  char *pwName[3];    int npwName;
  char *Other1=0;     int nOther1;
  char *Other2=0;
  OptStruct OSArray[]=
   { {"pwPol",   PA_CDOUBLE, 3, 3, (void *)pwPol,   &npwPol,  "plane wave polarization"},
     {"pwDir",   PA_DOUBLE,  3, 3, (void *)pwDir,   &npwDir,  "plane wave direction"},
     {"pwName",  PA_STRING,  1, 3, (void *)pwName,  &npwName, "plane wave name"},
     {"Other1",  PA_STRING,  1, 1, (void *)&Other1, &nOther1, "other parameter 1"},
     {"Other2",  PA_STRING,  1, 1, (void *)&Other2, 0,        "other parameter 2"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (npwName==0)
   OSUsage(argv[0],OSArray,"--pwName must be specified at least once");

  if ( npwPol != npwDir || npwDir!=npwName || npwName!=npwPol )
   ErrExit("numbers of --pwPol, --pwDir, and --pwName options must agree");

  if (Other1)
   printf("Value of Other1 is %s\n",Other1);
  if (Other2)
   printf("Value of Other2 is %s\n",Other2);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SetDefaultCD2SFormat("(%g,%g)");
  int npw;
  for(npw=0; npw<npwPol; npw++)
   { printf("Plane wave %i: (name %s): \n",npw, pwName[npw]);
     printf(" pol (%s, %s, %s)\n", CD2S(pwPol[3*npw+0]), CD2S(pwPol[3*npw+1]), CD2S(pwPol[3*npw+2]));
     printf(" dir (%g, %g, %g)\n", pwDir[3*npw + 0], pwDir[3*npw + 1], pwDir[3*npw + 2]);
     printf("\n");
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("Thank you for your support.\n");

}
