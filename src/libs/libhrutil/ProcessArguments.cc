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
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "libhrutil.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ASUsage(char *ProgName, ArgStruct *ASArray, const char *format, ...)
{
  va_list ap;
  char buffer[1000];
  ArgStruct *AS;

  if (format)
   { va_start(ap,format);
     vsnprintfEC(buffer,1000,format,ap);
     fprintf(stderr,"\nerror: %s (aborting)\n\n",buffer);
     va_end(ap);
   };

  fprintf(stderr,"usage: %s [options]\n\n",ProgName);
  fprintf(stderr," options: \n\n");

  for(AS=ASArray; AS->Name!=0; AS++)
   { 
     if (AS->Type==PA_BOOL) 
      sprintf(buffer,"--%s ",AS->Name);
     else
      sprintf(buffer,"--%s xx ",AS->Name);
     fprintf(stderr,"  %-15s",buffer);
     if ( AS->Description )
      fprintf(stderr,"(%s)",AS->Description);
     fprintf(stderr,"\n");
   };

  fprintf(stderr,"\n");
  
  exit(1);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ProcessArguments(int argc, char *argv[], ArgStruct *ASArray)
{
  int narg, nConv=0, ArgHandled;
  int Verbose=1;
  ArgStruct *AS;

  /*--------------------------------------------------------------*/
  /*- set default option values.  --------------------------------*/
  /*--------------------------------------------------------------*/
  for(AS=ASArray; AS->Name!=0; AS++)
   { 
     switch(AS->Type)
      { 
        case PA_DOUBLE: 
          if (AS->Default)
           sscanf(AS->Default,"%le",(double *)AS->Storage); 
          break;

        case PA_INT:
        case PA_BOOL:
          if (AS->Default)
           sscanf(AS->Default,"%i",(int *)AS->Storage); 
          break;

        case PA_STRING: 
          *(char **)AS->Storage=AS->Default ? strdupEC(AS->Default) : 0;
          break;

        default: 
          ErrExit("%s:%i: argument %s: unknown type %i",__FILE__,__LINE__,AS->Name,AS->Type);

      };
   };

  /*--------------------------------------------------------------*/
  /*- now process arguments, one at a time -----------------------*/
  /*--------------------------------------------------------------*/
  for (narg=1; narg<argc; narg++)
   { 
     /*--------------------------------------------------------------*/
     /*- make sure the argument starts with '--'                     */
     /*--------------------------------------------------------------*/
     if ( strncmp(argv[narg],"--",2) )
      ASUsage(argv[0],ASArray,"unknown option %s",argv[narg]);

     /*--------------------------------------------------------------*/
     /*- look for an option in the ASArray whose name matches the argument */
     /*--------------------------------------------------------------*/
     ArgHandled=0;
     for(AS=ASArray; !ArgHandled && AS->Name!=0; AS++)
      if ( !StrCaseCmp(argv[narg]+2,AS->Name) )
       { 
         ArgHandled=1;
         if ( AS->Type==PA_BOOL )
          { 
            *(int *)(AS->Storage)=1;

            if (Verbose)
             printf("--%s argument specified.\n",AS->Name);
          }
         else
          { 
            if ( narg+1 == argc )
             ASUsage(argv[0],ASArray,"%s requires an argument",argv[narg]);

            switch( AS->Type )
             { case PA_DOUBLE: 
                 nConv=sscanf(argv[narg+1],"%le",(double *)AS->Storage);
                 break;
               case PA_INT:
                 nConv=sscanf(argv[narg+1],"%i",(int *)AS->Storage);
                 break;
               case PA_STRING:
                 *((char **)AS->Storage)=strdupEC(argv[narg+1]);
                 nConv=1;
                 break;
             };
            if ( nConv != 1)
             ASUsage(argv[0],ASArray,"bad value %s specified for %s",argv[narg+1],argv[narg]);

            if (Verbose)
             printf("Setting --%s to %s.\n",AS->Name,argv[narg+1]);

           narg++;

          }; // if ( AS->Type==PA_BOOL ... else ... )
       }; // for(AS=ASArray...)

     /*--------------------------------------------------------------*/
     /*- handle unmatched arguments ---------------------------------*/
     /*--------------------------------------------------------------*/
     if ( AS->Name==0 && !ArgHandled )
      ASUsage(argv[0],ASArray,"unknown option %s",argv[narg]);

   };
}
