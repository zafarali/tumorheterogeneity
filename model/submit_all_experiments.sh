#!/bin/bash
# argument1 == mutation rate
# argument2 == seed
# argument3 == dryrun
# usage bash submit_all_experiments.sh u0.01 1 DRY
if [ $# <= 2 ]
then
    echo "not enough arguments"
    exit 1
fi
mutation_rate=${1}
seed=${2}

for executable in $(find "./experiments/${mutation_rate}/" -name "*.exe"); do
#    echo "executable name: $(echo $executable | cut -d '/' -f 4)";
    executable=$(echo $executable | cut -d '/' -f 4)
    executable=$(basename $executable .exe)
    if [ $# -eq 3 ]
    then
        echo "Command prepared: Qmsub -h 48 -n 4 -q sw run.sh ${executable} ${mutation_rate} ${seed}"
    else
        echo "Submitting: Qmsub -h 48 -n 4 -q sw run.sh ${executable} ${mutation_rate} ${seed}"
        Qmsub -h 48 -n 4 -q sw run.sh ${executable} ${mutation_rate} ${seed}
    fi
done

