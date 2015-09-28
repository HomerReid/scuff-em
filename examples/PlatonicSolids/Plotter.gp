FILE1='Sphere_498.AlphaVsTau'
FILE2='Sphere_924.AlphaVsTau'
FILE3='Octahedron_456.AlphaVsTau'
FILE4='Octahedron_975.AlphaVsTau'
FILE5='Tetrahedron_615.AlphaVsTau'
FILE6='Tetrahedron_1179.AlphaVsTau'

V=4*pi/3;

set terminal x11 1
set logscale x
unset logscale y
set xrange [1e-3:1e3]
set yrange [-2:5]
set key at 0.5,4
plot FILE1 u ($1/$2):($12/V) t 'Sphere, N=498'      w lp pt 7 ps 1.0  \
    ,FILE2 u ($1/$2):($12/V) t 'Sphere, N=924'      w lp pt 6 ps 1.25 \
    ,FILE3 u ($1/$2):($12/V) t 'Octahedron, N=456'  w lp pt 7 ps 1.0  \
    ,FILE4 u ($1/$2):($12/V) t 'Octahedron, N=975'  w lp pt 6 ps 1.25 \
    ,FILE5 u ($1/$2):($12/V) t 'Tetrahedron, N=456' w lp pt 7 ps 1.0  \
    ,FILE6 u ($1/$2):($12/V) t 'Tetrahedron, N=975' w lp pt 6 ps 1.25 
