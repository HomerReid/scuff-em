FILE1='OutFiles/Sphere_E1P1E1P2_762.YN.M11_M11'
FILE2='OutFiles/Sphere_E1P1E1P2_3348.YN.M11_M11'
FILE3='OutFiles/Sphere_E1P1E1P1_762.YN.M11_M11'
FILE4='OutFiles/Sphere_E1P1E1P1_3348.YN.M11_M11'
FILE5='../CONST_EPS_1.1.SphereTMatrix'

set logscale xy

plot  FILE1 u 1:(D2($2,$3)) t 'EpsUpper=1.1, EpsLower=1.2, N=762'  w lp pt 7 ps 1.5  \
    , FILE2 u 1:(D2($2,$3)) t 'EpsUpper=1.1, EpsLower=1.2, N=3348' w lp pt 7 ps 1.5  \
    , FILE3 u 1:(D2($2,$3)) t 'EpsUpper=1.1, EpsLower=1.1, N=762'  w lp pt 7 ps 1.5  \
    , FILE4 u 1:(D2($2,$3)) t 'EpsUpper=1.1, EpsLower=1.1, N=3348' w lp pt 7 ps 1.5  \
    , FILE5 u 1:(D2($2,$3)) t 'EpsUpper=1.1, EpsLower=1.1, Exact'  w l
