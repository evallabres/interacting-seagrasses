
    set terminal pngcairo enhanced font "arial,12.0" fontscale 0.5
    show terminal
    set output 'colordensity.0.png'
    
    set title 'Shoot Density 0'
    set xlabel 'meters'

    m = "./densitycolor.txt"
    
    set autoscale xfix
    set autoscale yfix
    set autoscale cbfix
    set pm3d map interpolate 0,0
    !set isosample 25, 25
    unset key
    set palette rgbformulae 33,13,10
    splot m u 1:2:3 with pm3d
    
    reset
    
