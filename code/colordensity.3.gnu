
    set terminal pngcairo enhanced font "arial, 20.0"  fontscale 1.25 size 1500, 1500  
    show terminal
    set output 'colordensity.3.png'
    
    set xlabel 'x (m)'
    set ylabel 'y (m)'
    set xtics 5
    set ytics 5

    m = "./colordensity.3.txt"
    
    set autoscale xfix
    set autoscale yfix
    set autoscale cbfix
    set pm3d map interpolate 0,0
    !set isosample 25, 25
    unset key
    set palette rgbformulae 33,13,10
    splot m u 1:2:3 with pm3d
    
    reset
    
