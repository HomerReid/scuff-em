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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "libhrutil.h"
#include "libhmat.h"

#define DIM 10
#define II cdouble (0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool ParseMatrixIndexString(const char *MIS, int End,
                            int *Start, int *Stop, int *Inc);

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  int Start, Stop, Inc;
  bool Status = ParseMatrixIndexString(argv[1], 10, &Start, &Stop, &Inc);
 
  if (Status==false)
   printf("Failed! \n");
  else
   { printf("Start=%i\n",Start);
     printf("Stop=%i\n",Stop);
     printf("Inc=%i\n",Inc);
   };
}
