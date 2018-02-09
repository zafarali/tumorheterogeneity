import sys

import matplotlib
matplotlib.use('Agg') # prevent use of a display

from analysis.Pipeline import Pipeline
from analysis.PipelineModules import load_tumor, create_kdsampler, marginal_counts_ordered, big_samples
from analysis.MixingModules import prepare_tumor_for_mixing_analysis, perform_mixing_analysis

print('arguments:'+str(sys.argv))

p = Pipeline(sys.argv[1],append_name='Jan18_MixingAnalysis',modules=[
        load_tumor,
        create_kdsampler,
        # marginal_counts_ordered,
        # big_samples,
        prepare_tumor_for_mixing_analysis,
        perform_mixing_analysis
    ])
p.execute()
