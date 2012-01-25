/*
 * changes required to make this work: 
 *  1) add a new class constructor to RWGObject that 
 *     takes an open file as input and ...
 *
 *  2) add a new class constructor to MatProp that 
 *     takes an open file as input and ...
 *
 *  3) finish implementing GTransformation
 *
 *  
 */

/*
 *
 */

#define MAXSTR 1000
#define MAXTOK 50  

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
RWGGeometry::RWGGeometry(const char *pGeoFileName)
{ 
  /***************************************************************/
  /* NOTE: i am not sure where to put this. put it here for now. */
  /***************************************************************/
  MatProp::SetLengthUnit(1.0e-6);
   
  /***************************************************************/
  /* storage for material properties defined on-the-fly in the   */
  /* .scuffgeo file.                                             */
  /* minor garbage-collection issue: MPs allocated in this       */
  /* routine are never freed.                                    */
  /***************************************************************/
  MatProp *MP, **MPs=0;
  int NumMPs;

  /***************************************************************/
  /* initialize simple fields ************************************/
  /***************************************************************/
  NumObjects=TotalBFs=TotalPanels=0;
  GeoFileName=strdup(pGeoFileName);
  ExternalMP=0;
  Objects=0;

  /***************************************************************/
  /* try to open input file **************************************/
  /***************************************************************/
  FILE *f=fopen(GeoFileName,"r");
  if (!f)
   RWGErrExit("could not open %s",GeoFileName);

  /***************************************************************/
  /* read and process lines from input file one at a time        */
  /***************************************************************/
  RWGObject *O;
  int no, nop;
  char Line[MAXSTR], Label[MAXSTR];
  int LineNum=0; 
  int nt, nTokens;
  char *p, *Tokens[MAXTOK];
  while( fgets(Line,MAXSTR,f) )
   { 
     LineNum++;

     /*--------------------------------------------------------------*/
     /*- break up line into tokens; skip blank lines and comments ---*/
     /*--------------------------------------------------------------*/
     nTokens=Tokenize(Line, Tokens, MAXTOK);
     if ( nTokens==0 || Tokens[0][0]=='#' )
      continue; 
    
     /*--------------------------------------------------------------*/
     /*- switch off based on first token on the line ----------------*/
     /*--------------------------------------------------------------*/
     if ( !strcasecmp(Tokens[0],"MEDIUM") )
      { 
       if ( nTokens!=3 || strcasecmp(Tokens[1],"MATERIAL") )
         RWGErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
        ExteriorMP=new MatProp(Tokens[2]);
        if (ExteriorMP->ErrMsg)
         RWGErrExit("file %s:%i: %s",GeoFileName,LineNum,ExteriorMP->ErrMsg);
      }
     else if ( !strcasecmp(Tokens[0],"MATERIAL") )
      {
        /*--------------------------------------------------------------*/
        /* hand off to MatProp class constructor to parse this section  */
        /*--------------------------------------------------------------*/
        if ( nTokens==1 )
         ErrExit("file %s:%i: no name given for MATERIAL ",GeoFileName,LineNum);
        else if ( nTokens>2 )
         ErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
         
        MP = new MatProp(f,Label);
        if (MP->ErrMsg)
         ErrExit("file %s:%i: %s",GeoFileName,LineNum,MP->ErrMsg); 

        NumMPs++;
        MPs=(MatProp *)realloc(MPs, NumMPs*sizeof(MatProp *));
        MPs[NumMPs-1]=MP;
           
      };
     else if ( !strcasecmp(Tokens[0],"OBJECT") )
      { 
        /*--------------------------------------------------------------*/
        /* hand off to RWGObject class constructor to parse this section*/
        /*--------------------------------------------------------------*/
        if ( nTokens>2 )
         ErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
        else if ( nTokens==2 )
         O=new RWGObject(f,Tokens[1],&LineNum);
        else if ( nTokens==1 )
         { snprintf(Label,MAXSTR,"Object_%i",NumObjects+1);
           O=new RWGObject(f,Label,&LineNum);
         };

        if (O->ErrMsg)
         ErrExit("file %s:%i: %s",GeoFileName,LineNum,O->ErrMsg); 

        NumObjects++;
        Objects=(RWGObject **)realloc(Objects, NumObjects*sizeof(RWGObject *) );
        Objects[NumObjects-1]=O;

        TotalBFs+=O->NumBFs;
        TotalPanels+=O->NumPanels;
        O->Index=NumObjects-1;
      }
     else 
       { 
         /*--------------------------------------------------------------*/
         /* unknown keyword                                              */
         /*--------------------------------------------------------------*/
         ErrExit("file %s:%i: syntax error",GeoFileName,LineNum);
       };

   }; // while( fgets(Line,MAXSTR,f) )

  /***************************************************************/
  /* process material properties of exterior medium             */
  /***************************************************************/
  int nmp;
  if (ExteriorMPName )
   {
     for(nmp=0; ExteriorMP==0 && nmp<NumMPs; nmp++)
      if ( !strcasecmp(ExteriorMPName, MPs[nmp]->Name) )
       ExteriorMP=MPs[nmp];

     if (ExteriorMP==0)
      { ExteriorMP = new MatProp(O->MPName);
        if (ExteriorMP->ErrMsg)
        ErrExit("%s: medium: error in MATERIAL value: %s",GeoFileName,O->MP->ErrMsg);
      };
   }
  else
   ExteriorMP=new MatProp(MP_VACUUM);

  /***************************************************************/
  /* process material properties of objects                      */
  /***************************************************************/
  for(no=0; no<NumObjects; no++)
   { 
     O=Objects[no];

     if ( O->MPName )
      { 
        // deallocate the default MP created by the RWGObject constructor
        delete O->MP;

        /* look first to see if the requested material was one of */
        /* the materials defined on-the-fly in the .scuffgeo file */
        for(nmp=0; O->MP==0 && nmp<NumMPs; nmp++)
         if ( !strcasecmp(O->MPName, MPs[nmp]->Name) )
          O->MP=MPs[nmp];

        /* otherwise ... */ 
        if (O->MP==0)
         { O->MP = new MatProp(O->MPName);
           if (O->MP->ErrMsg)
            ErrExit("%s: object %s (%s): error in MATERIAL value: %s",
                     GeoFileName,O->Label,O->MeshFileName,O->MP->ErrMsg); 
         };

        free(O->MPName);
      }
     else
      O->MP = new MatProp(MP_PEC);
   };

  /***************************************************************/
  /* process object nesting relationships                        */
  /***************************************************************/
  for(no=0; no<NumObjects; no++)
   { 
     O=Objects[no];

     if ( O->ContainingObjectName )
      { 
        /* look for an object appearing earlier in the .scuffgeo file */
        /* whose label matches the requested object label             */
        for(nop=0; O->ContainingObject==0 && nop<no; nop++)
         if ( !strcasecmp(O->ContainingObjectName, Objects[nop]->Label ) )
          O->ContainingObject=Objects[nop];

        /* if no matching object was found, it's an error */
        if (O->ContainingObject==0)
         { 
           for(nop=no+1; nop<NumObjects; nop++)
            if ( !strcasecmp(O->ContainingObjectName, Objects[nop]->Label ) )
             ErrExit("%s: object %s: containing object %s must appear earlier in file",
                     GeoFileName, O->Label, O->ContainingObjectName);

           ErrExit("%s: object %s: containing object %s not found",
                    GeoFileName, O->Label, O->ContainingObjectName);
         };

        free( O->ContainingObjectName );

      };
   };
 
  /*******************************************************************/
  /* compute average panel area for statistical bookkeeping purposes */
  /*******************************************************************/
  AveragePanelArea=0.0; 
  for(no=0; no<NumObjects; no++)
   for(np=0; np<Objects[no]->NumPanels; np++)
    AveragePanelArea+=Objects[no]->Panels[np]->Area;
  AveragePanelArea/=((double) TotalPanels);

  /*******************************************************************/
  /* set the AllPEC flag based on whether or not all material objects*/
  /* are PEC bodies                                                  */
  /*******************************************************************/
  AllPEC=1;
  for(no=0; no<NumObjects && AllPEC; no++)
   if ( !(Objects[no]->MP->IsPEC()) )
    AllPEC=0;

  /***************************************************************/
  /* initialize arrays of basis-function and panel index offsets */
  /***************************************************************/
  BFIndexOffset=(int *)RWGMalloc(NumObjects*sizeof(int) );
  PanelIndexOffset=(int *)RWGMalloc(NumObjects*sizeof(int) );
  BFIndexOffset[0]=PanelIndexOffset[0]=0;
  for(no=1; no<NumObjects; no++)
   { BFIndexOffset[no]=BFIndexOffset[no-1] + Objects[no-1]->NumBFs;
     PanelIndexOffset[no]=PanelIndexOffset[no-1] + Objects[no-1]->NumPanels;
   };

  /***************************************************************/
  /* initialize Identical[][] and Mate[] arrays.                 */
  /*                                                             */
  /* how it works:                                               */
  /*                                                             */
  /* (1) two objects are considered identical if                 */
  /*     (a) they have the same mesh file, and                   */
  /*     (b) they have the same material properties (i.e. they   */
  /*         were given identical values for the MATERIAL        */
  /*         keyword in the .rwggeo file.)                       */
  /*                                                             */
  /* (2) Identical[][] array: We set Identical[i][j] = 1 if      */
  /*                          objects i and j are identical, =0  */
  /*                          otherwise.                         */
  /*                                                             */
  /* (3) Mate[] array: If objects i, j, k, ... are identical and */
  /*                   i<j<k<..., then we set                    */
  /*                   Mate[i] = -1                              */
  /*                   Mate[j] = i                               */
  /*                   Mate[k] = i                               */
  /***************************************************************/
  Mate=(int *)malloc(NumObjects*sizeof(int));
  Mate[0]=-1;
  for(i=1; i<NumObjects; i++)
   { Mate[i]=-1;
     for(j=0; j<i && Mate[i]==-1; j++)
      if (    !strcmp(Objects[i]->MeshFileName, Objects[j]->MeshFileName)
           && !strcmp(Objects[i]->MP->Name    , Objects[j]->MP->Name)
         ) 
       Mate[i]=j;
   };

  /***************************************************************/
  /* initialize ObjectMoved[] array                              */
  /***************************************************************/
  ObjectMoved=(int *)RWGMalloc(NumObjects*sizeof(int));

}
