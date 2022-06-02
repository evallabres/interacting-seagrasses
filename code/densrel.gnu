
   
    set terminal pngcairo enhanced font "arial, 20.0" fontscale 1 size 1000, 900 
    show terminal
    set output 'densrel.png'
    
    unset key
    set grid
    set ylabel "Shoot density C. nodosa (m^{-2})"
    set xlabel "Shoot density C. prolifera (m^{-2})"


 
    m = "./hyst.txt"
    g(x) = b + a*x
    fit g(x) m using 3:2 via a,b
    
    plot m using 2:3 w p pt 7 ps 2
 
