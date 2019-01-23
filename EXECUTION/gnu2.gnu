reset

xmin = -0.05
xmax = 1.05
ymin = -0.05
ymax = 1.1

set term  png truecolor
set output 'vor0.png'
p [xmin:xmax][ymin:ymax] './mesh.0' t '' w l lw 2 lt 1,'./cent.0' t '' w p ps 2 lt 3,'./part.0' t 'Generator' w p ps 2 lt 1 pt 7
pause -1

set output 'vor1.png'
p [xmin:xmax][ymin:ymax] './mesh.1' t '' w l lw 2 lt 1,'./cent.1' t '' w p ps 2 lt 3,'./part.1' t 'Generator' w p ps 2 lt 1 pt 7
pause -1

set output 'vor2.png'
p [xmin:xmax][ymin:ymax] './mesh.2' t '' w l lw 2 lt 1,'./cent.2' t '' w p ps 2 lt 3,'./part.2' t 'Generator' w p ps 2 lt 1 pt 7
pause -1


set output 'vor3.png'
p [xmin:xmax][ymin:ymax] './mesh.3' t '' w l lw 2 lt 1,'./cent.3' t '' w p ps 2 lt 3,'./part.3' t 'Generator' w p ps 2 lt 1 pt 7
pause -1


set output 'vor4.png'
p [xmin:xmax][ymin:ymax] './mesh.4' t '' w l lw 2 lt 1,'./cent.4' t '' w p ps 2 lt 3,'./part.4' t 'Generator' w p ps 2 lt 1 pt 7
pause -1


set output 'vor5.png'
p [xmin:xmax][ymin:ymax] './mesh.5' t '' w l lw 2 lt 1,'./cent.5' t '' w p ps 2 lt 3,'./part.5' t 'Generator' w p ps 2 lt 1 pt 7
pause -1

set output 'vor6.png'
p [xmin:xmax][ymin:ymax] './mesh.6' t '' w l lw 2 lt 1,'./cent.6' t '' w p ps 2 lt 3,'./part.6' t 'Generator' w p ps 2 lt 1 pt 7
pause -1

set output 'vor7.png'
p [xmin:xmax][ymin:ymax] './mesh.7' t '' w l lw 2 lt 1,'./cent.7' t '' w p ps 2 lt 3,'./part.7' t 'Generator' w p ps 2 lt 1 pt 7
pause -1

set output 'vor8.png'
p [xmin:xmax][ymin:ymax] './mesh.8' t '' w l lw 2 lt 1,'./cent.8' t '' w p ps 2 lt 3,'./part.8' t 'Generator' w p ps 2 lt 1 pt 7
pause -1

set output 'vor9.png'
p [xmin:xmax][ymin:ymax] './mesh.9' t '' w l lw 2 lt 1,'./cent.9' t '' w p ps 2 lt 3,'./part.9' t 'Generator' w p ps 2 lt 1 pt 7
pause -1

p 1
set term x11
