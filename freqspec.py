#!/usr/bin/env python

import numpy as np
import pandas as pd
import glob
import sys

folder = sys.argv[1]
bins = np.linspace(0,0.5,num=500)

to_be_meaned = []
for path in glob.glob(folder+'*all_PMs*'):
    data = pd.read_csv(path, sep=' ', names=['Num', 'SNP', 'abundancy'])
    y,_ = np.histogram(data.abundancy, bins=bins)
    to_be_meaned.append(y)

np.save(folder+'PMfreqs_averaged.npy',np.mean(to_be_meaned, axis=0))

np.save(folder+'bins.npy',bins)


to_be_meaned = []
for path in glob.glob(folder+'*drv_PMs*'):
    data = pd.read_csv(path, sep=' ', names=['Num', 'SNP', 'abundancy'])
    y,_ = np.histogram(data.abundancy, bins=bins)
    to_be_meaned.append(y)

np.save(folder+'drvfreqs_averaged.npy',np.mean(to_be_meaned, axis=0))
