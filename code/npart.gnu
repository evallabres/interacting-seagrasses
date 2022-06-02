
   
    set terminal pngcairo enhanced font "arial, 20.0" fontscale 1 size 1000, 900 
    show terminal
    set output 'npart.png'
    
    unset key
    set grid
    set xlabel "Age (years) " 
    set ylabel "Number of shoots"
    set ytics 200000 

 
    m = "./npart.txt"

    plot m using 1:2 w l lw 2.5, m u 1:3 w l lw 2.5
    

    
 
