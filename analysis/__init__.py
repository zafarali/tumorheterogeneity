from MixingModules import perform_mixing_analysis
from PipelineModules import (
    load_tumor,
    create_kdsampler,
    marginal_counts_unordered,
    marginal_counts_ordered,
    density_plot,
    big_samples)

PAPER_PIPELINE = [ 
    load_tumor,
    create_kdsampler,
    marginal_counts_unordered,
    marginal_counts_ordered, 
    density_plot,
    big_samples,
    perform_mixing_analysis
]
