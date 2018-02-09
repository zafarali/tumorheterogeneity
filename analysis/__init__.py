from MixingModules import prepare_tumor_for_mixing_analysis, perform_mixing_analysis
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
    prepare_tumor_for_mixing_analysis,
    perform_mixing_analysis
]
