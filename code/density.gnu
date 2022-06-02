
   
    set terminal pngcairo enhanced font "arial, 20.0" fontscale 1 size 1000, 900 
    show terminal
    set output 'density.png'
    
    unset key
    set grid
    set xlabel "Age (years) " 
    set ylabel "Shoot density (m^{-2})" 
    set ytics 500

 
    m = "./density.txt" 
    n = "./pics.txt"

    plot m using 1:2 w l lw 2.5, m u 1:3 w l lw 2.5, n u 1:2 w p pt 7
    

    
 
