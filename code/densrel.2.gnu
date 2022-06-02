
   
    set terminal pngcairo enhanced font "arial, 20.0" fontscale 1 size 1000, 900 
    show terminal
    set output 'densrel.2.png'
    
    unset key
    set grid
    set xlabel "Shoot density C. nodosa (m^{-2})"
    set ylabel "Shoot density C. prolifera (m^{-2})"


 
    m = "./hyst.txt"
    g(x) = b2 + a2*x
    fit g(x) m using 2:3 via a2,b2
    
    plot m using 2:3 w p pt 7 ps 2
 
