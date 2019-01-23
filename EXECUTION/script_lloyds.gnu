 xmax=   1.0000000000000000     
 xmin=   0.0000000000000000     
 ymax=   1.0000000000000000     
 ymin=   0.0000000000000000     
 dx=  0.12000000000000001     
 dy=  0.1200000000000001     

set term  png enhanced

set output 'vor0.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part0.10" t "" w p pt 7 ps 1 lc 1,"mesh0.10" t "Generator#=          4      Cycle#           1 " w l lt 1 lc 'black', "cent0.10" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor1.png'
  p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part1.10" t "" w p pt 7 ps 1 lc 1,"mesh1.10" t "Generator#=          5      Cycle#           1 " w l lt 1 lc 'black', "cent1.10" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1   
set output 'vor2.png'
  p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part2.10" t "" w p pt 7 ps 1 lc 1,"mesh2.10" t "Generator#=          6      Cycle#           1 " w l lt 1 lc 'black', "cent2.10" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor3.png'
  p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part3.10" t "" w p pt 7 ps 1 lc 1,"mesh3.10" t "Generator#=          7      Cycle#           1 " w l lt 1 lc 'black', "cent3.10" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1  
set output 'vor4.png'
  p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part4.10" t "" w p pt 7 ps 1 lc 1,"mesh4.10" t "Generator#=          8      Cycle#           1 " w l lt 1 lc 'black', "cent4.10" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor5.png'
  p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part5.10" t "" w p pt 7 ps 1 lc 1,"mesh5.10" t "Generator#=          9      Cycle#           1 " w l lt 1 lc 'black', "cent5.10" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1  
set output 'vor6.png' 
  p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part6.10" t "" w p pt 7 ps 1 lc 1,"mesh6.10" t "Generator#=         10      Cycle#           1 " w l lt 1 lc 'black', "cent6.10" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor7.png'
  p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part7.10" t "" w p pt 7 ps 1 lc 1,"mesh7.10" t "Generator#=         11      Cycle#           1 " w l lt 1 lc 'black', "cent7.10" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor8.png'
  p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part8.10" t "" w p pt 7 ps 1 lc 1,"mesh8.10" t "Generator#=         12      Cycle#           1 " w l lt 1 lc 'black', "cent8.10" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor9.png'
  p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part9.10" t "" w p pt 7 ps 1 lc 1,"mesh9.10" t "Generator#=         13      Cycle#           1 " w l lt 1 lc 'black', "cent9.10" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1



set output 'vor_c01.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.11" t "" w p pt 7 ps 1 lc 1,"mesh.11" t "Generator#=          13      Cycle#           1 " w l lt 1 lc 'black', "cent.11" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c02.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.12" t "" w p pt 7 ps 1 lc 1,"mesh.12" t "Generator#=          13      Cycle#           2 " w l lt 1 lc 'black', "cent.12" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c03.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.13" t "" w p pt 7 ps 1 lc 1,"mesh.13" t "Generator#=          13      Cycle#           3 " w l lt 1 lc 'black', "cent.13" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c04.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.14" t "" w p pt 7 ps 1 lc 1,"mesh.14" t "Generator#=          13      Cycle#           4 " w l lt 1 lc 'black', "cent.14" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c05.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.15" t "" w p pt 7 ps 1 lc 1,"mesh.15" t "Generator#=          13      Cycle#           5 " w l lt 1 lc 'black', "cent.15" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c06.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.16" t "" w p pt 7 ps 1 lc 1,"mesh.16" t "Generator#=          13      Cycle#           6 " w l lt 1 lc 'black', "cent.16" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c07.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.17" t "" w p pt 7 ps 1 lc 1,"mesh.17" t "Generator#=          13      Cycle#           7 " w l lt 1 lc 'black', "cent.17" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c08.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.18" t "" w p pt 7 ps 1 lc 1,"mesh.18" t "Generator#=          13      Cycle#           8 " w l lt 1 lc 'black', "cent.18" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c09.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.19" t "" w p pt 7 ps 1 lc 1,"mesh.19" t "Generator#=          13      Cycle#           9 " w l lt 1 lc 'black', "cent.19" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c10.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.20" t "" w p pt 7 ps 1 lc 1,"mesh.20" t "Generator#=          13      Cycle#          10 " w l lt 1 lc 'black', "cent.20" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c11.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.21" t "" w p pt 7 ps 1 lc 1,"mesh.21" t "Generator#=          13      Cycle#          11 " w l lt 1 lc 'black', "cent.21" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c12.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.22" t "" w p pt 7 ps 1 lc 1,"mesh.22" t "Generator#=          13      Cycle#          12 " w l lt 1 lc 'black', "cent.22" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c13.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.23" t "" w p pt 7 ps 1 lc 1,"mesh.23" t "Generator#=          13      Cycle#          13 " w l lt 1 lc 'black', "cent.23" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c14.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.24" t "" w p pt 7 ps 1 lc 1,"mesh.24" t "Generator#=          13      Cycle#          14 " w l lt 1 lc 'black', "cent.24" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c15.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.25" t "" w p pt 7 ps 1 lc 1,"mesh.25" t "Generator#=          13      Cycle#          15 " w l lt 1 lc 'black', "cent.25" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c16.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.26" t "" w p pt 7 ps 1 lc 1,"mesh.26" t "Generator#=          13      Cycle#          16 " w l lt 1 lc 'black', "cent.26" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c17.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.27" t "" w p pt 7 ps 1 lc 1,"mesh.27" t "Generator#=          13      Cycle#          17 " w l lt 1 lc 'black', "cent.27" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c18.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.28" t "" w p pt 7 ps 1 lc 1,"mesh.28" t "Generator#=          13      Cycle#          18 " w l lt 1 lc 'black', "cent.28" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c19.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.29" t "" w p pt 7 ps 1 lc 1,"mesh.29" t "Generator#=          13      Cycle#          19 " w l lt 1 lc 'black', "cent.29" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c20.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.30" t "" w p pt 7 ps 1 lc 1,"mesh.30" t "Generator#=          13      Cycle#          20 " w l lt 1 lc 'black', "cent.30" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c21.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.31" t "" w p pt 7 ps 1 lc 1,"mesh.31" t "Generator#=          13      Cycle#          21 " w l lt 1 lc 'black', "cent.31" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c22.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.32" t "" w p pt 7 ps 1 lc 1,"mesh.32" t "Generator#=          13      Cycle#          22 " w l lt 1 lc 'black', "cent.32" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c23.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.33" t "" w p pt 7 ps 1 lc 1,"mesh.33" t "Generator#=          13      Cycle#          23 " w l lt 1 lc 'black', "cent.33" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c24.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.34" t "" w p pt 7 ps 1 lc 1,"mesh.34" t "Generator#=          13      Cycle#          24 " w l lt 1 lc 'black', "cent.34" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c25.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.35" t "" w p pt 7 ps 1 lc 1,"mesh.35" t "Generator#=          13      Cycle#          25 " w l lt 1 lc 'black', "cent.35" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c26.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.36" t "" w p pt 7 ps 1 lc 1,"mesh.36" t "Generator#=          13      Cycle#          26 " w l lt 1 lc 'black', "cent.36" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c27.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.37" t "" w p pt 7 ps 1 lc 1,"mesh.37" t "Generator#=          13      Cycle#          27 " w l lt 1 lc 'black', "cent.37" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c28.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.38" t "" w p pt 7 ps 1 lc 1,"mesh.38" t "Generator#=          13      Cycle#          28 " w l lt 1 lc 'black', "cent.38" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c29.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.39" t "" w p pt 7 ps 1 lc 1,"mesh.39" t "Generator#=          13      Cycle#          29 " w l lt 1 lc 'black', "cent.39" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c30.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.40" t "" w p pt 7 ps 1 lc 1,"mesh.40" t "Generator#=          13      Cycle#          30 " w l lt 1 lc 'black', "cent.40" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c31.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.41" t "" w p pt 7 ps 1 lc 1,"mesh.41" t "Generator#=          13      Cycle#          31 " w l lt 1 lc 'black', "cent.41" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c32.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.42" t "" w p pt 7 ps 1 lc 1,"mesh.42" t "Generator#=          13      Cycle#          32 " w l lt 1 lc 'black', "cent.42" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c33.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.43" t "" w p pt 7 ps 1 lc 1,"mesh.43" t "Generator#=          13      Cycle#          33 " w l lt 1 lc 'black', "cent.43" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c34.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.44" t "" w p pt 7 ps 1 lc 1,"mesh.44" t "Generator#=          13      Cycle#          34 " w l lt 1 lc 'black', "cent.44" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c35.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.45" t "" w p pt 7 ps 1 lc 1,"mesh.45" t "Generator#=          13      Cycle#          35 " w l lt 1 lc 'black', "cent.45" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c36.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.46" t "" w p pt 7 ps 1 lc 1,"mesh.46" t "Generator#=          13      Cycle#          36 " w l lt 1 lc 'black', "cent.46" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c37.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.47" t "" w p pt 7 ps 1 lc 1,"mesh.47" t "Generator#=          13      Cycle#          37 " w l lt 1 lc 'black', "cent.47" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c38.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.48" t "" w p pt 7 ps 1 lc 1,"mesh.48" t "Generator#=          13      Cycle#          38 " w l lt 1 lc 'black', "cent.48" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c39.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.49" t "" w p pt 7 ps 1 lc 1,"mesh.49" t "Generator#=          13      Cycle#          39 " w l lt 1 lc 'black', "cent.49" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c40.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.50" t "" w p pt 7 ps 1 lc 1,"mesh.50" t "Generator#=          13      Cycle#          40 " w l lt 1 lc 'black', "cent.50" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c41.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.51" t "" w p pt 7 ps 1 lc 1,"mesh.51" t "Generator#=          13      Cycle#          41 " w l lt 1 lc 'black', "cent.51" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1
set output 'vor_c42.png'
   p [xmin-dx:xmax+dx][ymin-dy:ymax+dy] "part.52" t "" w p pt 7 ps 1 lc 1,"mesh.52" t "Generator#=          13      Cycle#          42 " w l lt 1 lc 'black', "cent.52" t "Centroids" w p pt 2 ps 1 lc 3
  pause -1

set term x11
 p 1
