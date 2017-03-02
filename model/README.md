# TumorSimulator

Using [TumorSimulator](http://www2.ph.ed.ac.uk/~bwaclaw/cancer-code/).

_A spatial model predicts that dispersal and cell turnover limit intratumour heterogeneity_, Bartlomiej Waclaw, Ivana Bozic, Meredith E. Pittman, Ralph H. Hruban, Bert Vogelstein, and Martin A. Nowak, Nature 525, no. 7568 (September 10, 2015): 261-64. doi:10.1038/nature14971

# Reproducing Results

To compile these simulations, one must change the `params.h` file and then run `./compile.sh SIMULATION_NAME`. You can then use `/run.sh SIMULATION_NAME SEED` to run a simulation with the `SEED`.

Our parameters are shown in `./params.h` on lines `56-57`. To use one of these parameters, modify line `63` to the death rate you want and uncomment `68` to use surface turnover mode.

To modify the mutation rate to reproduce the data from the HCC paper, change `gama` on line `76` to what is written on `78`

We used the following seeds:

1. `0`
1. `1`
1. `2`
1. `3`
1. `4`
1. `5`
1. `6`
1. `101`
1. `102`
1. `12151`
1. `12152`