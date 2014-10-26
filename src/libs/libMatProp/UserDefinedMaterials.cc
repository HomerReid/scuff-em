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
 * UserDefinedMaterials.cc -- libMatProp module for handling materials defined 
 *                         -- by user-specified functions of frequency
 *
 * homer reid  3/2009 -- 12/2011
 * SGJ 3/2012 - change to use bundled cmatheval instead of muParser
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <libhrutil.h>
#include "libMatProp.h"

#include "cmatheval.h"

#define MAXCONSTANTS 25
#define MAXSTR       1000

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddMPToMatPropDataBase(MatProp *MP);

/***************************************************************/
/* body of the class constructor for the case of user-defined  */
/* materials.                                                  */
/***************************************************************/
void MatProp::CreateUserDefinedMaterial(const char *MatPropFileName, 
                                        const char *MaterialName)
{
  int status=0;
  char FileName[MAXSTR];

  EpsExpression = MuExpression = NULL;

  if (MatPropFileName)
   { snprintf(FileName,MAXSTR,"%s",MatPropFileName);
     status=ReadMaterialFromFile(FileName, MaterialName);
   };
  if( status==0 && getenv("SCUFF_MATPROPFILE") )
   { snprintf(FileName,MAXSTR,"%s",getenv("SCUFF_MATPROPFILE"));
     status=ReadMaterialFromFile(FileName, MaterialName);
   };
  if (status==0)
   { snprintf(FileName,MAXSTR,"%s","matprop.dat");
     status=ReadMaterialFromFile(FileName, MaterialName);
   };
  if (status==0)
   { snprintf(FileName,MAXSTR,"%s/.matprop.dat",getenv("HOME"));
     status=ReadMaterialFromFile(FileName, MaterialName);
   };
  
  if (status==1)
   { Log("Found material %s in file %s.",MaterialName,FileName);
     return;
   };

  if (status==-1)
   ErrExit(ErrMsg);

  if (status==0)
   ErrExit("could not find material %s in any matprop file",MaterialName);
}

/***************************************************************/
/* constructor helper function that tries to find a specified  */
/* material name in a big data file.                           */
/*                                                             */
/* returns:                                                    */
/*  1: success                                                 */
/*  0: file not found or material not found in file            */
/* -1: material found in file but syntax is incorrect          */
/*                                                             */
/* if the return code is -1, ErrMsg contains an error message. */
/***************************************************************/
int MatProp::ReadMaterialFromFile(const char *FileName, const char *MaterialName)
{ 
  FILE *f;
  char Line[MAXSTR];
  char *Tokens[3];
  int LineNum, NumTokens;
   
  if ( FileName==0 ) return 0;
  if ( !(f=fopen(FileName,"r")) ) return 0;

  /*--------------------------------------------------------------*/
  /*- read lines until we encounter one like 'MATERIAL ETHANOL'---*/
  /*--------------------------------------------------------------*/
  LineNum=0;
  while ( fgets(Line,MAXSTR,f) )
   { 
     LineNum++;
     NumTokens=Tokenize(Line,Tokens,3);

     if ( NumTokens==2 
          && !StrCaseCmp(Tokens[0],"MATERIAL") 
          && !StrCaseCmp(Tokens[1],MaterialName) 
        ) 
      { int Status=ParseMaterialSectionInFile(f, FileName, &LineNum);
        fclose(f);
        return Status;
      };
   };

  /*--------------------------------------------------------------*/
  /* didn't find it                                               */
  /*--------------------------------------------------------------*/
  fclose(f);
  return 0;
  
} 

/***************************************************************/
/* on entry, FileName is open with read pointer pointing to    */
/*  the line immediately following 'MATERIAL ETHANOL'          */
/*                                                             */
/* subsequent lines must have the structure                    */
/*                                                             */
/*  wp    = 1e-3 ;            (actually semicolon is optional) */
/*  gamma = 1e-4 ;                                             */
/*    ...    ...                                               */
/*  Eps(w) =  (w-wp) / (1+gamma)                               */
/*  Mu(w)  =  1.0 + wp*wp/(w*(w+i*Gamma)                       */
/*                                                             */
/*  ENDMATERIAL                                                */
/*                                                             */
/* if the material is successfully parsed, we return 1.        */
/* on failure, we return -1 and ErrMsg contains an error msg.  */
/***************************************************************/
int MatProp::ParseMaterialSectionInFile(FILE *f, const char *FileName, int *LineNum)
{ 
  char Line[MAXSTR];
  char *p, *pp, *ppp;

  int NumConstants;
  char ConstantNames[MAXCONSTANTS][MAXSTR]; 
  double ConstantValues[MAXCONSTANTS];
  int n, nc, nConv;
  double CV;
   
  /*--------------------------------------------------------------*/
  /*- read lines until we encounter 'ENDMATERIAL'               --*/
  /*--------------------------------------------------------------*/
  NumConstants=0;
  ErrMsg=0;
  while ( fgets(Line,MAXSTR,f) )
   { 
     (*LineNum)++;

     /*--------------------------------------------------------------*/
     /*- skip blank lines and comments  -----------------------------*/
     /*--------------------------------------------------------------*/
     p=Line; 
     while ( *p && isspace(*p) )
      p++;
     if ( *p==0 || *p=='#' )
      continue;

     /*--------------------------------------------------------------*/
     /*- look for end-of-material keyword   -------------------------*/
     /*--------------------------------------------------------------*/
     if ( !strncasecmp(p,"ENDMATERIAL",11) )
      { 
        return 1;
      };

     /*--------------------------------------------------------------*/
     /*- process lines that define one of the two material properties*/
     /*--------------------------------------------------------------*/
     if ( !strncasecmp(p,"Eps(w)",6) )
      { 
        pp=strchr(p,'=');
        if ( !pp )
         { ErrMsg=vstrdup("%s:%i: syntax error",FileName,*LineNum);
           return -1;
         };

        /* strip off the trailing carriage-return and semicolon */
        ppp=strrchr(pp,';');
        if (ppp)
         *ppp=0;
        else
         { ppp=strrchr(pp,'\n');
           if (ppp) 
            *ppp=0;
         };

        /* attempt to create a new cevaluator for the Eps(w) function */
	EpsExpression = cevaluator_create(pp+1);
	if (!EpsExpression) {
	  ErrMsg=vstrdup("%s:%i: invalid expression",FileName,*LineNum);
	  return -1;
	}
	cevaluator_set_var_index(EpsExpression, "w", 0); // w is thread-safe
        for(nc=0; nc<NumConstants; nc++)
         cevaluator_set_var(EpsExpression, ConstantNames[nc],ConstantValues[nc]);

      }
     else if ( !strncasecmp(p,"Mu(w)",5) )
      { 
        pp=strchr(p,'=');
        if ( !pp )
         { ErrMsg=vstrdup("%s:%i: syntax error",FileName,*LineNum);
           return -1;
         };

        /* strip off the trailing carriage-return and semicolon */
        ppp=strrchr(pp,';');
        if (ppp)
         *ppp=0;
        else
         { ppp=strrchr(pp,'\n');
           if (ppp) 
            *ppp=0;
         };

        /* attempt to create a new cevaluator for the Mu(w) function */
	MuExpression = cevaluator_create(pp+1);
	if (!MuExpression) {
	  ErrMsg=vstrdup("%s:%i: invalid expression",FileName,*LineNum);
	  return -1;
	}
	cevaluator_set_var_index(MuExpression, "w", 0); // w is thread-safe
        for(nc=0; nc<NumConstants; nc++)
	  cevaluator_set_var(MuExpression, ConstantNames[nc],ConstantValues[nc]);

      }
     else
      { 
        /*--------------------------------------------------------------*/
        /*- otherwise try to interpret the line as a constant definition*/
        /*-  which is to say as a line of the form 'name = value;'      */
        /*--------------------------------------------------------------*/

        /* look for an equals sign */
        pp=strchr(p,'=');
        if ( !pp )
         { ErrMsg=vstrdup("%s:%i: unrecognized input",FileName,*LineNum);
           return -1;
         };

        /* try to interpret the characters after the = sign as a number */
        nConv=sscanf(pp+1,"%le%n",&CV,&n);
        if (nConv<1)
         { ErrMsg=vstrdup("%s:%i: invalid constant value",FileName,*LineNum);
           return -1;
         };

        /* make sure there isn't any junk (other than a possible '; or # comment') after the number*/ 
        ppp=pp+1+n;
        while( *ppp!=0 && *ppp!='#' && (isspace(*ppp) || *ppp==';') ) 
         ppp++;
        if ( *ppp!=0 && *ppp!='#' )
         { ErrMsg=vstrdup("%s:%i: junk at end of line",FileName,*LineNum);
           return -1;
         };

        /* conversion succeeded; note the name and value of the new constant */ 
        while ( isspace (*(pp-1) ) ) 
         pp--;
        *pp=0; 
        if ( NumConstants==MAXCONSTANTS )
         { ErrMsg=vstrdup("%s:%i: too many constants (max allowed is %i)",FileName,*LineNum,MAXCONSTANTS);
           return -1;
         };
        if ( strlen(p) > (MAXSTR-1) )
         { ErrMsg=vstrdup("%s:%i: name of constant is too long", FileName,*LineNum);
           return -1;
         };
        strcpy(ConstantNames[NumConstants],p);
        ConstantValues[NumConstants]=CV;
        NumConstants++;
      };
   };

  ErrMsg=vstrdup("%s: unexpected end of file",FileName);
  return -1;

}

/***************************************************************/
/* this routine allows user-defined materials to be defined    */
/* on-the-fly in places other than the standard database       */
/* files.                                                      */
/* how it works: on entry, the file read pointer is assumed to */
/* point immediately after a line reading                      */
/*                                                             */
/*  MATERIAL MaterialName                                      */
/*                                                             */
/* this routine parses the ensuing lines as if it were reading */
/* a MATERIAL...ENDMATERIAL section in an ordinary MatProp     */
/* database file.                                              */
/*                                                             */
/* if an error is encountered, the routine stops and returns   */
/* an error message.                                           */
/*                                                             */
/* otherwise, a new MatProp is created for the user-defined    */
/* material and added to the internal database of MatProp      */
/* structures maintained within libMatProp. then, when a user  */
/* subsequently calls the MatProp() constructor to create an   */
/* instance of the material in question, it is available and   */
/* the MatProp() constructor succeeds, which it otherwise      */
/* wouldn't since the material cannot be found in the usual    */
/* database files.                                             */
/*                                                             */
/* in this latter (successful) case, on return from this       */
/* routine the return value is 0 and the file read pointer     */
/* points to the line immediately following ENDMATERIAL.       */
/***************************************************************/
char *AddMaterialToMatPropDataBase(FILE *f, char *FileName,
                                   char *MaterialName, int *LineNum)
{ 
  MatProp *MP=new MatProp(); // call the default constructor

  MP->Type=MP_PARSED;
  MP->Zeroed=0;
  free(MP->Name);
  MP->Name=strdupEC(MaterialName);

  MP->ErrMsg=0;
  MP->EpsExpression=MP->MuExpression=NULL;
  MP->InterpReal=MP->InterpImag=NULL;
  MP->ParseMaterialSectionInFile(f, FileName, LineNum);

  if (MP->ErrMsg)
   return MP->ErrMsg;

  AddMPToMatPropDataBase(MP);

  return 0;
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void MatProp::GetEpsMu_Parsed(cdouble Omega, cdouble *pEps, cdouble *pMu)
{
  static char *OmegaVar= "w";
  Omega *= FreqUnit;

  // Note that this is thread-safe since "w" is an indexed variable,
  // which means that the internal symbol table is not modified, but
  // we MUST pass its value as the 0th array entry in cevaluator_evaluate.

  if (pEps) {
    if (EpsExpression) {
      *pEps = cevaluator_evaluate(EpsExpression, 1,&OmegaVar,&Omega);
    }
    else
      *pEps = 1.0;
  }

  if (pMu) {
    if (MuExpression) {
      *pMu = cevaluator_evaluate(MuExpression, 1,&OmegaVar,&Omega);
    }
    else
      *pMu = 1.0;
  }
}
