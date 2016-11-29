import numpy as np
from collections import Counter

GAMMA = 0.5772156649015328606065120900824024310421
L_PM = (1<<30) -1 
D_PM = 1<<30
R_PM = 1<<31

def SNP_count(genotypes_list, min_freq=0):
	"""
		returns the counts for SNPs in the sample
		@params:
			genotypes_list : list of Genotypes in the sample
		@returns the SNP count
	"""

	counts = Counter()

	number_of_cells = float(len(genotypes_list))

	map(lambda genotype: counts.update(genotype.snps), genotypes_list)
	
	if min_freq >0:
		uc = []
		# c[0] = the SNP id
		# c[1] = the count of that SNP
		for f in filter( lambda c: c[1]/number_of_cells > min_freq, counts.items() ):
			# f[0] = the SNP id
			# f[1] = the count of that SNP
			uc.extend( [ f[0] for i in range(f[1]) ] )
		counts = Counter(uc)

	return counts

def drivers_count(genotypes_list, min_freq=0):
	"""
		returns the counts for driver SNPs in the sample
		@params:	
			genotypes_list : list of Genotypes in the sample
		@returns the driverSNP count
	"""
	counts = Counter()
	# map over all genotypes
	map(\
		lambda genotype: counts.update( \
			filter( lambda snp: (snp & D_PM) > 0, genotype.snps ) ), \
		genotypes_list ) 		# ^ checks if a PM is a driver, filter snps accordingly.
		# only counts the SNPs that pass the check.
	
	number_of_cells = float(len(genotypes_list))

	if min_freq >0:
		uc = []
		for f in filter( lambda c: c[1]/number_of_cells > min_freq, counts.items() ):
			uc.extend( [ f[0] for i in range(f[1]) ] )
		counts = Counter(uc)

	return counts

def unique_driver_combinations(genotypes_list, min_freq=0):

	# here we hash a tuple of the driver snps
	# this hash is unique based on the combination of driver snps
	# 
	driver_combos = set( map(\
		lambda genotype: \
		tuple( sorted( filter( lambda snp: (snp & D_PM ) > 0 , genotype.snps ) ) ) , \
		genotypes_list ) )

	return len(driver_combos)

def unique_snp_combinations(genotypes_list, min_freq=0):

	# here we hash a tuple of the snps
	# this hash is unique based on the combination of snps
	# 
	snp_combos = set( map(\
		lambda genotype: \
		tuple( sorted(  genotype.snps  ) ) , \
		genotypes_list ) )

	return len(snp_combos)

def unique_res_combinations(genotypes_list):

	# here we hash a tuple of the snps
	# this hash is unique based on the combination of snps
	# 
	snp_combos = set( map(\
		lambda genotype: \
		tuple( sorted( filter( lambda snp: (snp & R_PM ) > 0 , genotype.snps ) ) ) , \
		genotypes_list ) )
	return len(snp_combos)

def unique_driver_proportion(driver_counts):
	"""
		Obtains the proportion of driver SNPs that are unique
	"""
	counts = np.array( driver_counts.values() )
	return np.sum( counts == 1 ) / float(len( driver_counts.keys() ))

def proportion_of_pairwise_differences(SNP_counts, sample_size):
	"""
		Calculates k or E(pi) the proportion of pairwise differences 
		according to the following equation:
		x_l = proportion of chromosomes where locus l is mutated

		using exact formula:
			2*(m/n)*(n-m)/(n-1)
		sample_size = the number of cells in the sample

		@params:
			SNP_counts returned from SNP_count function
	"""
	# x holds the proportion of chromosomes at which a locus is mutated
	m = np.array(SNP_counts.values(), dtype=np.float)
	x = m/float(sample_size)
	y = (sample_size - m)/float(sample_size-1)
	return np.sum( 2 * x * y )

def number_of_segregating_sites(SNP_counts, sample_size):
	"""
		Calculates S the number of segregating sites
		(basically the number of mutations not shared by all in the sample)
		@params:
			SNP_counts returned from SNP_count function
			sample_size: the size of the sample
	"""
	S = 0
	for SNP, count in SNP_counts.items():
		if count != sample_size:
			S += 1
		# only add to the number of segregating sites if it is not
		# shared by the whole sample
	#endfor
	return S


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
		sum_{i=1, @n} 1/i
		implemented as https://en.wikipedia.org/wiki/Harmonic_number#Calculation

	"""
	assert n>0, 'n must be bigger than 0'
	return GAMMA + np.log(n) + 0.5/n - 1./(12*n**2) + 1./(120*n**4)

def generalized_harmonic_number(n, power):
	"""
		Returns the nth-generalized 2nd order harmonic number according to:
		sum_{i=1, @n} 1/i^@power
		@TODO: is there a better way to do this?
	"""
	assert n>0, 'n must be bigger than 0'
	# use n+1 because arange doesn't include the n+1st integer.
	if power == 1 and n > 100: 
		return harmonic_number(n)
	return np.sum(1./(np.arange(1,n+1)**power))

def normalized_segregating_sites(SNP_counts, sample_size):
	"""
		Calculates the normalized number of segreating sites: S / H(n-1)
		@params:
			SNP_counts: returned from the SNP_count function
			sample_size: the number of samples
	"""
	return number_of_segregating_sites(SNP_counts) / float(harmonic_number(sample_size-1))

def tajimas_D(SNP_counts, n, return_parts=False):
	"""
		Calculates Tajimas D according to reference below.
		@params:
			SNP_counts : SNP counts
			n : size of the sample
			return_parts: returns intermediate calculations (S, S/H and pi)
		Simonsen, K. L., Churchill, G. A., & Aquadro, C. F. (1995). 
		Properties of statistical tests of neutrality for DNA polymorphism data. 
		Genetics, 141(1), 413-429.
	"""
	# the harmonic numbers
	# an = harmonic_number(n-1)
	an = generalized_harmonic_number(n-1, 1) # 1st harmonic number
	bn = generalized_harmonic_number(n-1, 2) # 2nd harmonic number
	
	vT = (float(2*(n**2 + n + 3))/(9*n*(n-1)) - float(n+2)/(n*an) + bn/an)/(an**2 + bn)
	uT = float(float(n+1)/(3*(n-1)) - 1./an)/an - vT

	S = number_of_segregating_sites(SNP_counts, n)
	pi = proportion_of_pairwise_differences(SNP_counts, n)

	denom = np.sqrt(uT*S + vT*(S**2))
	D = float(pi - float(S)/an)/denom if denom != 0 else 0
	if np.isnan(D):
		print 'vt:',vT
		print 'ut:',uT
		print 'n:',n
		print 'an:',an
		print 'bn:',bn
		print 'S:',S
		print 'pi:',pi
		raise Exception('NAN found during D calculation')

	return (pi, S, float(S)/an, D) if return_parts else D

