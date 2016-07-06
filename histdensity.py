#!/usr/bin/env python
import sys

import matplotlib
matplotlib.use('Agg') # prevent use of a display
import matplotlib.pyplot as plt

from analysis.BaseObjects import Tumor
from analysis.Samplers import KDTSphericalSampler
import numpy as np

from glob import glob
print('Loading...')
folder = sys.argv[1]
cell_file = glob(folder+'/cells*')[0]
genotype_file = glob(folder+'/genotypes*')[0]


print(cell_file)
print(genotype_file)

t = Tumor.from_files(cell_file, genotype_file)
ss = KDTSphericalSampler(t)
print('Loaded.')

print('Plotting...')
# plt.figure(1,figsize=(10,10))
# plt.hist(np.sqrt(np.sum((ss.cell_positions-ss.tumor.COM)**2,axis=1)), bins=100)
# plt.xlabel('distance from COM r=sqrt((x-x_COM)^2+(y-y_COM)^2+(z-z_COM)^2)')
# plt.ylabel('counts')
# plt.savefig(folder+'_'+'r.pdf')


# plt.figure(2,figsize=(7.5,10))
# # plt.suptitle('Distribution of x,y,z coordinates and distance from COM')
# plt.subplot(3,1,1)
# plt.title('Distribution of cell x coordinates')
# plt.hist(ss.cell_positions[:,0], bins=100)
# plt.xlabel('x')
# plt.subplot(3,1,2)
# plt.title('Distribution of cell y coordinates')
# plt.hist(ss.cell_positions[:,1],bins=100)
# plt.xlabel('y')
# plt.subplot(3,1,3)
# plt.title('Distribution of cell z coordinates')
# plt.hist(ss.cell_positions[:,2],bins=100)
# plt.xlabel('z')
# plt.tight_layout()
# plt.savefig(folder+'_'+'xyz.pdf')


rho, r = np.histogram(np.sqrt(np.sum((ss.cell_positions-ss.tumor.COM)**2,axis=1)), bins=np.linspace(0,1000,1000))

# np.save(folder+'_rho', rho)
# np.save(folder+'_r', r)

# plt.figure(3,figsize=(7.5,10))
# rho_r = rho/ (4*np.pi*r[1:]**2)
# np.save(folder+'_rho_r1', rho_r)
# plt.plot(r[1:],rho_r)
# plt.ylim([0,1])
# plt.xlabel('r')
# plt.ylabel('density')
# plt.title('using r[1:]')
# plt.savefig(folder+'_rho_r1.pdf')

# plt.figure(4,figsize=(7.5,10))
# rho_r = rho/ (4*np.pi*r[:-1]**2)
# np.save(folder+'_rho_r-1', rho_r)
# plt.plot(r[:-1],rho_r)
# plt.ylim([0,1])
# plt.xlabel('r')
# plt.ylabel('density')
# plt.title('using r[:-1]')
# plt.savefig(folder+'_rho_r0.pdf')
# print('Plotted')


def neighbour_iterator(arr):
    index = 0
    while index < len(arr)-1:
        yield (arr[index], arr[index+1])
        index += 1


r_meaned = np.mean(np.array(list(neighbour_iterator(r))), axis=1)
rho2 = rho/ (4*np.pi*r_meaned**2)        
np.save('./'+folder+'_outs_r_meaned', r_meaned)
np.save('./'+folder+'_outs_rho_corrected', rho2)

plt.figure(3,figsize=(10,7.5))
plt.plot(r_meaned,rho2)
plt.xlabel('r')
plt.ylabel('density')
plt.title('Tumor Density')
plt.savefig(folder+'_rho_corrected.pdf')