file ../../../applications/scuff-caspol/scuff-caspol
set environment SCUFF_MATRIX_2018=1
set args --geometry PECSlab_2640.scuffgeo --Atom Rb --EPFile EPFile4 --XiFile ../XiFile2
break GetDyadicGFs.cc:127
break ErrExit
run
