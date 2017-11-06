##################################################
##################################################
##################################################
function LMP2Alpha(L, M, P)

  int( 2*( L*(L+1) + M - 1 ) + (P-1) + 1 );

end

##################################################
# This function reads a .TMatrix file produced by
# scuff-tmatrix and returns a list of N angular
# frequencies and N matrices. Each matrix is
# dimension DxD and corresponds to the T-matrix
# for a single angular frequency.
#
# inputs:
#
#  FileName: a string like "SiO2Sphere.TMatrix" 
#
# outputs:
#
#  Omega:    a vector whose entries are the
#            angular frequencies at which the 
#            input file contained T-matrix data
#
#  TMatrix:  array of T-matrix data
#
#  TMatrix[n, Alpha, Beta] = 
#    (Alpha, Beta) element of T-matrix 
#    corresponding to angular frequency 
#    Omega[n]                             
#
##################################################
function ReadTMatrix(FileName)

  ##################################################
  # try to open the file ###########################
  ##################################################
  RawData=0;
  try
    RawData = readdlm(FileName)
  catch
    @printf("error: could not read file %s (aborting)\n",FileName);
    return
  end

  ##################################################
  # try to convert the first column of the file to   
  # floating-point numbers and extract the unique 
  # frequencies
  ##################################################
  Omega=0;
  try
    Omega = unique(float(RawData[:,1]));
  catch
     @printf("error: file %s contained bad data (aborting)\n",FileName);
     return
  end
  nOmega = length(Omega);

  ##################################################
  # try to determine the value of LMax               
  ##################################################
  TMatrixRemainder = length(RawData[:,1]) % nOmega;

  if (TMatrixRemainder!=0)
     @printf("error: not all frequencies have the same number of entries (aborting)\n")
     return
  end

  TMatrixDimension= int( sqrt( length(RawData[:,1]) / nOmega ) );
  LMax = sqrt( TMatrixDimension/2 + 1 ) - 1;
  if (floor(LMax) != LMax )
     @printf("error: incorrect matrix size (aborting)\n")
     return
  end

  @printf("%i frequencies, LMax=%i, TMatrix dimension=%i\n",nOmega,LMax,TMatrixDimension);

  ##################################################
  # Now go through the table and put each entry in its 
  # proper placeatrid
  ##################################################
  nRow=1
  TMatrix=(0.0 + 0.0im)*zeros(nOmega, TMatrixDimension, TMatrixDimension);
  for n=1:nOmega

    for L=1:LMax, M=-L:L, P=1:2, LP=1:LMax, MP=-LP:LP, PP=1:2

        ##################################################
        ##################################################
        ##################################################
        if (    float(RawData[nRow,1])!=float(Omega[n])
             ||   int(RawData[nRow,3])!=L
             ||   int(RawData[nRow,4])!=M
             ||   int(RawData[nRow,6])!=LP
             ||   int(RawData[nRow,7])!=MP
           ) 
          @printf("Line %i: should be %e (%i %i %+i) (%i %i %+i), ",
                   nRow,Omega[n],P,L,M,PP,LP,MP);
          @printf(" is %s (%s %s %s) (%s %s %s)\n",
                   RawData[nRow,1],
                   RawData[nRow,2], RawData[nRow,3], RawData[nRow,4],
                   RawData[nRow,5], RawData[nRow,6], RawData[nRow,7]);
                   
          return
        end
      
        ##################################################
        ##################################################
        ##################################################
        @printf("%i (%i %i %i) (%i %i %+i) (%i,%i) \n",nRow,
                 L,M,P,LP,MP,PP,LMP2Mu(L,M,P),LMP2Mu(LP,MP,PP));
        TMatrix[n, LMP2Mu(L, M, P), LMP2Mu(LP, MP, PP) ] = float(RawData[nRow,8]) + 1.0im*float(RawData[nRow,9]);

        nRow+=1
    end
  end

  #################################################
  ##################################################
  ##################################################
  Omega, TMatrix

end
