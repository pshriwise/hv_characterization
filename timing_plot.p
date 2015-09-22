set ylabel 'Area Fraction (A_f)' 
set xlabel 'Valence (n)'
set format x "%2.0tE%L"
set format cb "%2.0tE%L"
set palette rgb 33,13,10 negative
#set log x
set zlabel 'Avg Ray Fire Time (ms)' 
set pm3d map
set dgrid3d 
set cbrange [1e-06:1e-4]
set cblabel 'time (s)'
splot 'data.dat' using 1:2:4
pause -1
#
#unset pm3d
#set xlabel 'Triangle Density (HV Area / Valence)'
#set ylabel 'Avg Ray Fire Time (ms)'
#set log xy
#plot 'data.dat' matrix using (($1*2)/($2*60000)):3
#
#pause -1
