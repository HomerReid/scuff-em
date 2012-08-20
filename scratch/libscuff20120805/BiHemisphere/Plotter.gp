FILE00='X0P1Y0P1.20.2'
FILE01='X0P5Y0P1.20.2'
FILE10='BiHemisphere_762.X0P1Y0P1.out'
FILE11='BiHemisphere_762.X0P5Y0P1.out'
FILE20='BiHemisphere_3348.X0P1Y0P1.out'
FILE21='BiHemisphere_3348.X0P5Y0P1.out'

set terminal x11 1
set xlabel 'Z'
set ylabel 'E_x/E_0'
call 'vline' '-0.99'
call 'vline' '0.99'
set title 'X component of E-field along (X,Y,Z)=(0.1,0.1,Z)'
plot   FILE10  u 3:4 t 'scuff, N=762' w p pt 7 ps 1  \
     , FILE20  u 3:4 t 'scuff, N=3348' w p pt 7 ps 1 \
     , FILE00  u 3:4 t 'Exact' w l
unset arrow

set terminal x11 2
set xlabel 'Z'
set ylabel 'E_x/E_0'
call 'vline' '-0.8602'
call 'vline' '+0.8602'
set title 'X component of E-field along (X,Y,Z)=(0.5,0.1,Z)'
plot   FILE11  u 3:4 t 'scuff, N=762' w p pt 7 ps 1  \
     , FILE21  u 3:4 t 'scuff, N=3348' w p pt 7 ps 1 \
     , FILE01  u 3:4 t 'Exact' w l
unset arrow

set terminal x11 3
set xlabel 'Z'
set ylabel 'E_z/E_0'
call 'vline' '-0.99'
call 'vline' '+0.99'
set title 'Z component of E-field along (X,Y,Z)=(0.1,0.1,Z)'
plot   FILE10  u 3:8 t 'scuff, N=762' w p pt 7 ps 1  \
     , FILE20  u 3:8 t 'scuff, N=3348' w p pt 7 ps 1 \
     , FILE00  u 3:6 t 'Exact' w l
unset arrow

set terminal x11 4
set xlabel 'Z'
set ylabel 'E_z/E_0'
call 'vline' '-0.8602'
call 'vline' '+0.8602'
set title 'Z component of E-field along (X,Y,Z)=(0.5,0.1,Z)'
plot   FILE11  u 3:8 t 'scuff, N=762' w p pt 7 ps 1  \
     , FILE21  u 3:8 t 'scuff, N=3348' w p pt 7 ps 1 \
     , FILE01  u 3:6 t 'Exact' w l
unset arrow
