#!/bin/bash
mkdir experiments
mkdir experiments/u0.01
mkdir experiments/u0.01875
mkdir experiments/u0.01transition

for mutation_rate in "u0.01" "u0.01875" "u0.01transition"; do
    for filename in $(echo "./parameter_files/$mutation_rate/*.h" | cut -d ' ' -f 1); do
        echo "compiling $filename";
        # mutarate=$(echo $filename | cut -d '/' -f 3);
        # echo $substring;
        fileroot=$(basename $filename .h);
        # echo $fileroot;
        cp $filename ./params.h
        g++ simulation.cpp main.cpp functions.cpp -w -O3 -o ./$fileroot.exe
        mv $fileroot.exe "./experiments/$mutation_rate/$fileroot.exe"
    done
done