set xlabel '$z$'
set ylabel '$\phi(z)$'
set key at 0.15,0.5 spacing 4

plot \
  'None.EPFile.out'         u 3:4 t 'No substrate'      w lp pt 7 ps 1 \
 ,'E10HalfSpace.EPFile.out' u 3:4 t 'Eps=10 half space' w lp pt 7 ps 1 \
 ,'E10Slab.EPFile.out'      u 3:4 t 'Eps=10 slab'       w lp pt 7 ps 1 \
 ,'E10SlabGP.EPFile.out'    u 3:4 t 'Eps=10 slab w GP' w lp pt 7 ps 1 \

call 'latex' 'E10Substrate'
