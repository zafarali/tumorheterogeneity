# modules for pipeline
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import csv
import os
import itertools
from BaseObjects import Tumor
from Samplers import SphericalSampler, KDTSphericalSampler
from PipelineModules import sample_coordinate
import Statistics
import random
import time
import json
import multiprocessing.dummy as multiprocessing
from collections import namedtuple
import logging
import sys

CPU_COUNT = 12
EPS = np.finfo(float).eps
MAX_CELLS_RANGE = 50
INTERVAL = 5
LARGE_CLUSTERS = [100, 1000, 10000, 20000]
SMALL_CLUSTERS = [2, 6, 27, 50]
CLUSTER_SIZES = SMALL_CLUSTERS + LARGE_CLUSTERS

def load_snps(datafile):
    mutations = pd.read_csv(datafile, sep=' ', names=['Num', 'SNP', 'abundancy'])
    mutations['abundancy'] = mutations['abundancy'] + EPS 
    mutations['sampling_prob'] = mutations['abundancy']/mutations['abundancy'].sum()
    return mutations 

# convinience functions to store information and save it in np format
cluster_information = namedtuple('cluster_information','distance_from_tumor_COM cell_distance_from_cluster_COM p_shared_gas n_shared_gas p_private_gas n_private_gas cumulative_S')
def save_cluster_information(ci):
    to_save = []
    to_save.append([ci.distance_from_tumor_COM] * len(ci.cell_distance_from_cluster_COM))
    to_save.append(ci.cell_distance_from_cluster_COM)
    to_save.append(ci.p_shared_gas)
    to_save.append(ci.n_shared_gas)
    to_save.append(ci.p_private_gas)
    to_save.append(ci.n_private_gas)
    to_save.append(ci.cumulative_S)
    return np.array(to_save)


def perform_mixing_analysis(pipeline):
    results = []
    pipeline.print2('Performing mixing analysis:')
    rand_coordinate = sample_coordinate(pipeline.sampler.cell_positions)
    for i in range(100):
        try:
            outs = pipeline.sampler.sample_fixed_points(n=100,
                                                        centre=rand_coordinate.next(),
                                                        with_genotypes=True,
                                                        with_distances=True)

            cluster_COM, _, genotypes, cell_distances = outs
            distance_to_tumor_COM = np.sqrt(np.sum((np.array(cluster_COM) - np.array(pipeline.tumor.COM))**2))
            psgas, nsgas, ppgas, npgas = Statistics.diff_GAs(genotypes[0], genotypes[1:])

            cumulative_S = [0]+[len(Statistics.SNP_count(genotypes[:k]).keys()) for k in range(1,len(genotypes))]
            ci = cluster_information(distance_to_tumor_COM, cell_distances, psgas, nsgas, ppgas, npgas, cumulative_S)
            results.append(save_cluster_information(ci))
        except Exception as e:
            pipeline.print2('1/2 Exception occured: '+str(sys.exc_info()))
            pipeline.print2('2/2 Exception occured:'+str(e))
    
    results = np.stack(results)

    save_loc = pipeline.FILES['out_directory'] +'/mixing_analysis.npy'
    pipeline.print2('Saving mixing analysis...')
    np.save(save_loc, results)
    pipeline.print2('Saved mixing analysis!')

