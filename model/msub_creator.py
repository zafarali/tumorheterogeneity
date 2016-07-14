#!/usr/bin/env python

import glob
import pandas as pd
import numpy as np
import random
import os
import stat
import sys

folder = sys.argv[1]
model = folder[:5] if len(sys.argv) < 2 else sys.argv[2]
if model[:2] == './':
    raise Exception("Make sure to include absolute path only...")
    
bins = [(2,7), (8,12), (13,17), (18,23), (23,30)]
all_samples = glob.glob('./'+folder+'/pipe*/samples/*')
print all_samples

all_data = []
for pathname in all_samples:
    details = pathname.split('/')[-1]
    detail_split =  details.split('size')[-1].split('dist')
    sz = int(detail_split[0])
    detail_split = detail_split[-1].split('haplo')
    dist = float(detail_split[0])
    nhaplo = int(detail_split[-1].split('i')[0])
    all_data.append([sz,dist,nhaplo,pathname])
    
all_data = pd.DataFrame(all_data, columns=['samplesize', 'distance', 'nhaplotypes', 'pathname'])
all_data = all_data[all_data['samplesize'] < 30]



commands = []
for binnum,limit in enumerate(bins):
    subsample = all_data[(all_data['samplesize'] <= limit[1]) & (all_data['samplesize'] >= limit[0])]
    selected = random.sample(subsample.values, 3)
    local_commands = []
    for select in selected:
        cells_file = glob.glob(select[3]+'/cells*')[0]
        gen_file = glob.glob(select[3]+'/gen*')[0]
        det_file = glob.glob(select[3]+'/det*')[0]
        # local_commands.append('cp '+cells_file+' ./')
        # local_commands.append('cp '+gen_file+' ./')
        local_commands.append('./'+model+'_reseed.exe reseeded_'+model+'-'+str(binnum)+' 1 1 '+cells_file+' '+gen_file+' \n')
        local_commands.append('cp '+det_file+' ./reseeded_'+model+'-'+str(binnum)+'/ \n')
    local_commands.append('../freqspec.py reseeded_'+model+'-'+str(binnum)+'/ \n')
    commands.append(local_commands)

    
for i,commands in enumerate(commands):
    filename = './'+model+'_executable'+str(i)+'.sh'
    with open(filename, 'w') as f:
        f.writelines(commands)
    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IEXEC)
