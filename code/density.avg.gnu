
   
    set terminal pngcairo enhanced font "arial, 20.0" fontscale 1 size 1000, 900 
    show terminal
    set output 'density.avg.png'
    
    unset key
    set grid
    set xlabel "Age (years) " 
    set ylabel "Shoot density (m^{-2})"
    set ytics 500
 
    m = "./density.avg.txt"
    
    plot m using 2:3 w l lw 2.5, m u 2:4 w l lw 2.5
    

    
 
