import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import numpy as np
import sys
import seaborn as sns
from analysis.BaseObjects import Tumor
from analysis.MixingModules import load_snps
path_to_folder = sys.argv[1]

SNPS = load_snps(path_to_folder+'/all_PMs_1_0.dat')
tumor = Tumor.from_files(path_to_folder+'/cells_1000000.dat', path_to_folder+'/genotypes_1000000.dat')

def cmap(tumor, all_snps, top_n, palette='colorblind'):
    """
    Creates a color map
    """
    color_palette = sns.color_palette(palette, top_n) + [(0, 0, 0)]
    most_abundant_snps = np.array(all_snps).astype(int)[:top_n].tolist()
    genotype_to_color = {}
    for i, snp in enumerate(most_abundant_snps):
        try:
            genotype_idx = tumor.genotype_idx_with_snp(snp)
        except KeyError as e:
            pass
        for genotype in genotype_idx:
            genotype_to_color[genotype] = i
    def cmap_(i_):
        try:
            print(genotype_to_color[i_])
            return color_palette[genotype_to_color[i_]]
        except Exception as e:
            print('Skipped key:'+str(e))
            return color_palette[-1]
    return cmap_

cmap_vectorized = np.vectorize(cmap(tumor, SNPS.SNP, 50))
cells_to_plot = tumor.cells[np.where(tumor.cells[:, 2]==0)]
cells_to_plot_x_y = cells_to_plot[:, [0, 1]]
colors = cells_to_plot[:, [3]]
sns.set_style('white')
plt.scatter(*zip(*cells_to_plot_x_y), c=zip(*cmap_vectorized(colors)))
plt.xlabel('x')
plt.ylabel('y')
plt.savefig(path_to_folder+'/tumor_image.pdf')
