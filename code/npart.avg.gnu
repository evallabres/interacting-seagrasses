
   
    set terminal pngcairo enhanced font "arial, 20.0" fontscale 1 size 1000, 900 
    show terminal
    set output 'npart.avg.png'
    
    unset key
    set grid
    set xlabel "Age (years) " 
    set ylabel "Number of shoots" 
    set ytics 200000 

 
    m = "./npart.avg.txt"

    plot m using 2:3 w l lw 2.5, m u 2:4 w l lw 2.5
    

    
 
