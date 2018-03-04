import sys
import numpy as np


sys.path.append('..')


import matplotlib
matplotlib.use('Agg') # prevent use of a display

import matplotlib.pyplot as plt

from analysis import Statistics
from analysis.Pipeline import Pipeline
from analysis.PipelineModules import load_tumor, create_kdsampler, sample_coordinate

if len(sys.argv) < 4:

	sys.exit('not enough arguments: (need: \n sys.argv[1] = turnover/noturnover \n sys.argv[2] = close/notclose, sys.argv[3] = path to folder containing all simulations')

print 'mode=',sys.argv[1]
print 'distance=',sys.argv[2]

if sys.argv[1] == 'turnover':
	filename = sys.argv[3] + '/1_0_02_outs_10'
else:
	filename = sys.argv[3] + '/1_0_0_outs_10'

# surface_turnover = sys.argv[1]
p = Pipeline(filename,append_name='AFS_big_test',modules=[load_tumor,create_kdsampler])
p.execute()

def neighbour_iterator(arr):
    index = 0
    while index < len(arr)-1:
        yield (arr[index], arr[index+1])
        index += 1

BINS = np.linspace(0, 0.3, num=100)

class AFS(object):
    """
        Holds an Allele-Frequency Spectra for sample
    """
    
    def __init__(self, genotypes, bins=BINS):
        """
            @params
                genotypes: list of all genotypes in the sample
                n: size of the sample
        """
        self.n = len(genotypes)
        self.counts = Statistics.SNP_count(genotypes) # returns { SNPID : SNPCOUNT, ...}
        self.freqs = np.array(self.counts.values(), dtype=np.float) / np.float(self.n)  # returns [ SNPCOUNT / n, ... ]
        self.AFS = np.histogram(self.freqs, bins=bins)
        
    def plot(self, ax=None, rebin=False, loglog=False, color=(0,0,0), stop_plot=False):
        y1, x1 = self.AFS if not rebin else np.histogram(self.freqs, bins=rebin)
        
        x1_meaned = np.mean(np.array(list(neighbour_iterator(x1))), axis=1)
        y1_keep = y1!=0

        x1_plt = np.log10(x1_meaned[y1_keep]) if loglog else x1_meaned[y1_keep]
        y1_plt = np.log(y1[y1_keep]) if loglog else y1[y1_keep]
        
        if not ax:
            f = plt.figure()
            ax = f.add_subplot(111)
        
        s001 = np.sum(self.freqs > 0.01)
        s01 = np.sum(self.freqs > 0.1)
        s = self.freqs.shape[0]
        
        if not stop_plot:
            ax.scatter(x1_plt, y1_plt, c=color,alpha=0.5)
            ax.set_title('sample size = '+str(self.n)+'\n{S(>0) = '+str(s) + ', S(>0.01) = '+str(s001)+' , S(>0.10) = '+str(s01)+'}')
            ax.set_xlabel('Frequency')
            ax.set_ylabel('Count')
        
        
        return ax, s, s001, s01
        

if sys.argv[2] == 'close':
    coord = (0,0,0)
else:
    coord = (-145.35, 201.64, 163.64)


sampler = p.sampler
rand_coordinate = sample_coordinate(sampler.cell_positions, deviate=True)
# coord = (100,100,0)
sample_COM, cell_positions, genotypes_sample_sorted = sampler.sample_fixed_points(100000, centre=coord)
s = np.zeros(100000)
s001 = np.zeros(100000)
s01 = np.zeros(100000)

print 'distance found='+str(np.sqrt( np.sum( ( np.array(sample_COM) - np.array(p.tumor.COM) )**2) ))
# for i in [4, 6, 10]:
for i in range(1, 1000, 1) + range(1000, 100000, 100):
    
    genotypes_selected = genotypes_sample_sorted[:i]
    afs = AFS(genotypes_selected)
    _, _s, _s001, _s01 = afs.plot(stop_plot=True)
    s[i] = _s
    s001[i] = _s001
    s01[i] = _s01


np.save(sys.argv[1]+'_'+sys.argv[2]+'_0.npy', s)
np.save(sys.argv[1]+'_'+sys.argv[2]+'_001.npy', s001)
np.save(sys.argv[1]+'_'+sys.argv[2]+'_01.npy', s01)

print 'DONE'

