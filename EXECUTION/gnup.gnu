 p './accuracy.plt' u 1:2 w lp,'./accuracy.plt' u 1:3 w l lt 3,'' u 1:4 w l lt 3


p './bounds.plt' u 1:2 w l,'./nodal_resolution.plt' u 4:3 w p lt 3


p  './mesh.plt' w l,'./shock.plt' w p ps 0.2 pt 7 lt 3
