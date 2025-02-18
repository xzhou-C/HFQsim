#!/bin/bash


set -e

#---------------------------------------------------------------------
function run_prog ()
{
    prog=$1

    echo "Compiling..."
    make $prog
    
    last_line=$(tail -n 1 input_parameters.txt)
    echo "$last_line"

    echo "Running program..."
    ./${prog} > ${last_line}_eqcat.out

    echo "Plotting..."
    gnuplot <<EOF
    set terminal pngcairo enhanced font 'Verdana,8'
    set output '${prog}_timemag.png'
    set xrange [10:20]
    set xlabel "Time [yr]"
    set ylabel "M_L"
    plot '${prog}_eqcat.out' u 1:6 w points
EOF
    

    echo "Displaying"
    if [ -x /usr/bin/eog ]; then
        eog ${prog}_timemag.png
    elif [ -x /usr/bin/ristretto ]; then
        ristretto ${prog}_timemag.png
    else
        echo Output check plot written to $prog.png
    fi
}
#---------------------------------------------------------------------

run_prog test_model
