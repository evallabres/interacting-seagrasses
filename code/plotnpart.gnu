
    set terminal png nocrop enhanced size 800,600 font "arial,12.0"
    show terminal
    set output 'npart.png'

    
    
    set xlabel "time (years)"
    set ylabel "Number of shoot (log scales)"
    set title 'Number of shoots'
    !set logscale x
    set logscale y
    set grid
    
    m = "./npart.txt"
    
    f(x) = b*exp(a*x)
    fit [0:8] f(x) m using 1:2 via a,b
    

    
       
    plot m using 1:2 w p pt 7 ps 0.25, [0:8] f(x)
    !, [7:16] g(x)
    reset
    
