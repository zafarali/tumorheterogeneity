import sys

import matplotlib
matplotlib.use('Agg') # prevent use of a display

from analysis.Pipeline import Pipeline
from analysis.PipelineModules import load_tumor, create_kdsampler, marginal_counts_ordered, big_samples
from analysis.MixingModules import perform_mixing_analysis

print('arguments:'+str(sys.argv))

p = Pipeline(sys.argv[1],append_name='MixingAnalysisOnly',modules=[
        load_tumor,
        create_kdsampler,
        perform_mixing_analysis
    ])
p.execute()
