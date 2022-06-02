
   
    set terminal pngcairo enhanced font "arial, 20.0" fontscale 1 size 1000, 900 
    show terminal
    set output 'hyst.png'
    
    unset key
    set grid
    set xlabel "Mortality rate C. prolifera (yr^{-1}) "
    set ylabel "Shoot density (m^{-2})"
    set ytics 500
 
    m = "./hyst.txt"
    
    plot m using 1:2 w p pt 7 ps 2, m u 1:3 w p pt 7 ps 2
    
   


    
 
