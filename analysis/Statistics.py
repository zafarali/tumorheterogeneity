import numpy as np
from collections import Counter

GAMMA = 0.5772156649015328606065120900824024310421

def SNP_count(genotypes_list):
	"""
		returns the counts for SNPs in the sample
		@params:
			genotypes_list : list of Genotypes in the sample
		@returns the SNP count
	"""

	counts = Counter()

	map(lambda genotype: counts.update(genotype.snps), genotypes_list)

	return counts

def proportion_of_pairwise_differences(SNP_counts, sample_size):
	"""
		Calculates k or E(pi) the proportion of pairwise differences 
		according to the following equation:
		x_l = proportion of chromosomes where locus l is mutated

		pi = sum_(all loci l) 2*x_l*(1-x_l)
		sample_size = the number of cells in the sample

		@params:
			SNP_counts returned from SNP_count function
	"""
	# x holds the proportion of chromosomes at which a locus is mutated
	x = np.array(SNP_counts.values(), dtype=np.float)/sample_size
	return np.sum( 2 * x * (1-x) )

def number_of_segregating_sites(SNP_counts):
	"""
		Calculates S the number of segregating sites
		(basically the number of mutations in the sample)
		@params:
			SNP_counts returned from SNP_count function
	"""
	return len(SNP_counts.keys())

def number_of_singletons(SNP_counts):
	"""
		Calculates the number of mutations appearing in only one sampled cell
		@params:
			SNP_counts returned from SNP_count function
	"""
	return np.sum(np.array(SNP_counts.values()) == 1)

def centre_of_mass(sample):
	"""
		Calculates the center of mass of a sample of cells
		@params:
			sample containing a list of Cell objects or raw list of integers corresponding
			to the same
	"""
	sample = np.array(sample, dtype=np.float) if type(sample) is not np.ndarray else sample

	R_CM = np.sum( sample[:,:3] , axis = 0 ) / float( len( sample ) )

	return R_CM

def harmonic_number(n):
	"""
		Returns the nth-harmonic number according to:
		sum_{i=1, n-1} 1/i
		implemented as https://en.wikipedia.org/wiki/Harmonic_number#Calculation

	"""
	return GAMMA + np.log(n) + 0.5/n - 1./(12*n**2) + 1./(120*n**4)

def normalized_segregating_sites(SNP_counts, sample_size):
	"""
		Calculates the normalized number of segreating sites: S / H(n-1)
		@params:
			SNP_counts: returned from the SNP_count function
			sample_size: the number of samples
	"""
	return number_of_segregating_sites(SNP_counts) / float(harmonic_number(sample_size-1))
