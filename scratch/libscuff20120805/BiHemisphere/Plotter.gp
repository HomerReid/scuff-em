FILE1='WithoutHalfEdges/BiHemisphere_3348.X0P1Y0P1.out'
FILE2='New/BiHemisphere_3348.X0P1Y0P.WithoutHalfEdges.out'
FILE3='New/BiHemisphere_3348.X0P1Y0P1.WithHalfEdgesWithoutLineCharges.out'
FILE4='New/BiHemisphere_3348.X0P1Y0P.WithHalfEdgesWithLineCharges.out'
FILE9='X0P1Y0P1.20.2'		

plot   FILE1   u (-$3):8 t '(Old) Without half-RWG functions' w p pt 7 ps 1.5 \
     , FILE2   u 3:8 t '(New) Without half-RWG functions' w p pt 7 ps 1.5 \
     , FILE3   u 3:8 t '(New) With half-RWG functions, without line charges' w p pt 7 ps 1.5 \
     , FILE4   u 3:8 t '(New) With half-RWG functions, with line charges' w p pt 7 ps 1.5  \
     , FILE9   u 3:6 t 'Exact calculation' w l

set terminal x11 2

FILE10 = 'New/E10Sphere_3348.X0P1Y0P1.WithHalfBFsWithLineCharges.out'
FILE11 = 'New/E10Sphere_3348.X0P1Y0P1.WithHalfBFsWithoutLineCharges.out'
FILE12 = 'New/E10Sphere_3348.X0P1Y0P1.WithoutHalfBFs.out'
FILE19 = 'X0P1Y0P1.10.10'

plot   FILE10  u 3:8 t 'With half-RWG functions, with line charges'    w p pt 7 ps 1.5 \
     , FILE11  u 3:8 t 'With half-RWG functions, without line charges' w p pt 7 ps 1.5 \
     , FILE12  u 3:8 t 'WithoutHalfBFs.out'                            w p pt 7 ps 1.5 \
     , FILE19  u 3:6 t 'Exact calculation' w l

