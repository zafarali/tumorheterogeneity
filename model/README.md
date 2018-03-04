# TumorSimulator

Using [TumorSimulator](http://www2.ph.ed.ac.uk/~bwaclaw/cancer-code/).

_A spatial model predicts that dispersal and cell turnover limit intratumour heterogeneity_, Bartlomiej Waclaw, Ivana Bozic, Meredith E. Pittman, Ralph H. Hruban, Bert Vogelstein, and Martin A. Nowak, Nature 525, no. 7568 (September 10, 2015): 261-64. doi:10.1038/nature14971

# Reproducing Results

The first step is to compile the executable files:

```
bash compile_all_experiments.sh
```

You will now see `./experiments/` folder with all the compiled executables.


To run experiment you will need to do the following:

```
bash run.sh 1_0_0 u0.01 1
```

where `1_0_0` is the name of the executable, `u0.01` is the mutation rate and `1` is the seed.
This will run the file in `./experiments/u0.01/1_0_0.exe` and save data into `./experiments/u0.01/1_0_0_outs_1/`


To run all of these on a cluster in a batch use:

```
bash submit_all_experiments.sh u0.01 1 DRY
```

which will run all the experiments with mutation rate 0.01 and seed 1.

Note this uses an internal Qmsub (quick msub) submission script so you will need to modify this based on your cluster.

