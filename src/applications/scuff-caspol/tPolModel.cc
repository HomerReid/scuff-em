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

#include "libhrutil.h"
#include "libhmat.h"
#include "scuff-caspol.h"

#include <time.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXSTR  1000
#define MAXXI   10
#define MAXATOM 10

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *Atoms[MAXATOM];    int NumAtoms;
  /* name           type  #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { 
     {"Atom",        PA_STRING,  1, 10,      (void *)Atoms,        &NumAtoms,   "type of atom"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  PolModel **PolModels = (PolModel **)malloc(NumAtoms * sizeof(PolModel *));
  HMatrix **Alphas = (HMatrix **)malloc(NumAtoms * sizeof(HMatrix *));
  for(int na=0; na<NumAtoms; na++)
   { PolModels[na] = new PolModel(Atoms[na]);
     if (PolModels[na]->ErrMsg)
      ErrExit(PolModels[na]->ErrMsg);
     Alphas[na] = new HMatrix(3,3);
   };

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  FILE *f=fopen("tPolModel.out","w");
  for( double Xi = 1.0e-6; Xi<=10.0; Xi*=exp(0.1*log(10.0)) )
   { 
     fprintf(f,"%e ",Xi);
     for(int na=0; na<NumAtoms; na++)
      { PolModels[na] -> GetPolarizability(Xi, Alphas[na]);
        fprintf(f,"%e ",Alphas[na]->GetEntryD(0,0) );
      };
     fprintf(f,"\n");
   };
  fclose(f);

}
