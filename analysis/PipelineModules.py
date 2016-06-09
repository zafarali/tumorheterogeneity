# modules for pipeline
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import csv
import itertools
from BaseObjects import Tumor
from Samplers import SphericalSampler, KDTSphericalSampler
import Statistics


# generates a coordinate
def generate_coordinate(max_value, min_value):
	while True:
		yield (max_value-min_value)*np.random.random() + max_value

def sample_coordinate(coordinates, deviate=False):
	"""
		Returns a random coordinate from a sample of coordintes
		@params:
			coordinates : the array of shape (num_coordinates, dimension) to sample from 
			deviate[=False] : perturbates the sampled coordinate by N(0,1)
		@returns:
			coordinate of shape (1, dimension)
	"""
	max_size = coordinates.shape[0]
	while True:
		sampled_coordinate= coordinates[np.random.randint(max_size),:]
		yield ( sampled_coordinate + np.random.rand() ) if deviate else sampled_coordinate

def generate_angle():
	while True:
		yield np.random.random() * 2* np.pi


def load_tumor(pipeline):
	"""
		Loads Tumor Into the Pipeline
	"""
	pipeline.tumor = Tumor.from_files(pipeline.FILES['cells'], \
		pipeline.FILES['genotypes'] )
	pipeline.print2('Tumor was loaded')

def random_spherical_samples(pipeline):
	"""
		generates samples from the tumor according to the specs
	"""
	pipeline.tumor.cells
	samples = []
	if not pipeline.tumor:
		raise Exception('Tumor doesnt exist')
	sampler = SphericalSampler(pipeline.tumor)
	# rand_coordinate = generate_coordinate( np.max( sampler.cell_positions ), np.min( sampler.cell_positions ) )
	rand_coordinate = sample_coordinate(sampler.cell_positions, deviate=True)

	pipeline.print2('Begining Tumor Sampling')
	for radius in pipeline.specs['RADII']:
		# generate 3 random coodinate
		pipeline.print2('Sampling radius:'+str(radius))
		for i in xrange( pipeline.specs.get('repeats', 25) ):
			# generate a new coordinate
			# centre = ( rand_coordinate.next() , rand_coordinate.next() , rand_coordinate.next() )
			centre = rand_coordinate.next()

			# conduct the sample
			sample = sampler.sample(radius=radius, centre=centre, with_genotypes=True)
			
			print Statistics.get_drivers_only(sample[1],sampler.tumor.drivers)

			# insert a tuple of (radius, center, sample)
			samples.append( ( radius, centre, sample ) )
			
		pipeline.print2(str(pipeline.specs.get('repeats', 25))+' samples conducted')
	pipeline.print2('Sampling completed')
	pipeline.samples = samples

def calculate_statistics(pipeline):
	"""
		Calculates Statistics
	"""
	if not pipeline.samples:
		raise Exception('Tumor must be sampled first!')

	stats = [['radius', 'distance_to_COM', 'sample_size', 'pi', 'S', 'S/H', 'Tajimas D', 'D','total_SNPs', 'unique_drivers', 'total_drivers']]

	pipeline.print2('Calculating Statistcs')
	for radius, centre, sample in pipeline.samples:
		# first index contains the genotypes
		distance_to_COM = np.sqrt(np.sum((np.array(centre) - np.array(pipeline.tumor.COM))**2))
		
		SNP_counts = Statistics.SNP_count(sample[1])
		n = len(sample[1])
		if n < 2:
			pipeline.print2('Skipped sample due to zero or one individuals')
			continue
		pi, S, SH, D = Statistics.tajimas_D(SNP_counts, n, return_parts=True)
		
		raise DeprecationWarning('This module is going to be deprecated due to outdated function signatures')

		unique_drivers, total_drivers = Statistics.unique_driver_proportion(SNP_counts, pipeline.tumor.drivers)

		stats.append( [radius, distance_to_COM, n, pi, S, SH, D, pi-SH, SNP_counts, unique_drivers, total_drivers ] )

	pipeline.stats = stats
	pipeline.print2('Statistics Calculated')


def inline_statistics(pipeline):

	stats = [['radius', 'distance_to_COM', 'sample_size', \
	'pi', 'S', 'S/H', 'Tajimas D', 'D', \
	'total_SNPs', 'proportion_driver', 'unique_driver_propotion',\
	'total_drivers', 'unique_combos']]

	pipeline.print2('Creating KDTSphericalSampler')
	sampler = KDTSphericalSampler(pipeline.tumor)
	pipeline.print2('KDTSphericalSampler created.')
	rand_coordinate = sample_coordinate(sampler.cell_positions, deviate=True)
	number_of_samples = pipeline.specs.get('repeats', 100)

	for radius in pipeline.specs['RADII']:
		pipeline.print2('Sampling radius:'+str(radius))
		i = 0
		# for i in xrange( number_of_samples ):
		while i != number_of_samples:
			if i > 1000:
				pipeline.print2('Did not find enough samples for radius='+str(radius)+'.')
				break
			# generate a new coordinate
			centre = rand_coordinate.next()

			# conduct the sample
			sample = sampler.sample(radius=radius, centre=centre, with_genotypes=True)
			
			# calculate statistics

			distance_to_COM = np.sqrt( np.sum( ( np.array(centre) - np.array(pipeline.tumor.COM ) )**2 ) )
			
			SNP_counts = Statistics.SNP_count(sample[1])

			n = len(sample[1])
			if n < 2:
				pipeline.print2('Skipped sample due to only ' + str(n)+ ' individuals')
				continue
			i+=1
			# base statistics pi, S, SH, D
			pi, S, SH, D = Statistics.tajimas_D(SNP_counts, n, return_parts=True)

			# number of SNPs in the sample
			total_SNPs = len( SNP_counts.keys() )

			# frequency of driver mutations
			drivers_count = Statistics.drivers_count(sample[1])

			# number of driver mutations in the sample
			total_drivers = len( drivers_count.keys() ) 

			unique_drivers = Statistics.unique_driver_proportion(drivers_count)

			unique_combos = Statistics.unique_driver_combinations(sample[1])
			
			driver_to_SNP_ratio = total_drivers/float(total_SNPs) if total_SNPs > 0 else 0

			stats.append( [radius, distance_to_COM, n, pi, S, SH, D, pi-SH, \
				total_SNPs, driver_to_SNP_ratio,\
				unique_drivers, total_drivers, unique_combos ] )
		
		#end sampling for loop

		pipeline.print2( str(number_of_samples) +' samples conducted and statistics too')
	# end

	pipeline.stats = stats
	pipeline.print2('Sampling/Statistics Completed')




def save_statistics(pipeline):
	"""
		Saves Statistics
	"""

	if not pipeline.stats:
		raise Exception('Tumor must first have statistics')

	pipeline.print2('Saving Statistics...')
	save_location = pipeline.FILES['out_directory']+'/statistics.csv'
	with open(save_location, 'w') as f:
		writer = csv.writer(f)
		for row in pipeline.stats:
			writer.writerow(row)
		pipeline.print2('Statistics saved')

def all_plot_combos(pipeline):

	indep = ['radius', 'distance_to_COM', 'sample_size']
	dep = ['pi', 'S', 'S/H', 'Tajimas D', 'D', \
		'total_SNPs', 'proportion_driver', 'unique_driver_propotion',\
		'total_drivers', 'unique_combos']
	
	# generates all pairs of indep and dep variables
	pairs = list(itertools.product(indep, dep))

	# creates a dataframe for easy plotting
	df = pd.DataFrame(pipeline.stats[1:], columns=pipeline.stats[0])
	radii = pipeline.specs['RADII']

	pipeline.print2('Plotting')
	for x_label,y_label in pairs:

		colors = iter( cm.rainbow( np.linspace( 0, 1, len(radii) ) ) ) 
		plt.figure(figsize=(10,10))
		plt.title(x_label+' vs '+y_label)
		for radius in radii:
			# subsample the data according to the radius so we can colour it
			subsample = df[df['radius'] == radius]
			plt.scatter( subsample[x_label], \
				subsample[y_label], color=next(colors), label='radius='+str(radius))
			plt.xlabel(x_label)
			plt.ylabel(y_label)
		plt.legend(bbox_to_anchor=(1.25,0.75))
		y_label_serialized = y_label.replace('/', '')
		plt.savefig( pipeline.FILES['out_directory']+'/'+x_label+'_vs_'+y_label_serialized+'.pdf',\
		bbox_inches='tight')
	pipeline.print2('Plotting done')


BASE = [ load_tumor, random_spherical_samples, calculate_statistics, save_statistics ]
KD_SAMPLING = [ load_tumor, inline_statistics, save_statistics, all_plot_combos ]