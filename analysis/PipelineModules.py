# modules for pipeline
import numpy as np
import csv
from BaseObjects import Tumor
from Samplers import SphericalSampler
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
		pipeline.FILES['genotypes'], drivers_file=pipeline.FILES['drivers'])
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

		unique_drivers, total_drivers = Statistics.driver_proportion(SNP_counts, pipeline.tumor.drivers)

		stats.append( [radius, distance_to_COM, n, pi, S, SH, D, pi-SH, SNP_counts, unique_drivers, total_drivers ] )

	pipeline.stats = stats
	pipeline.print2('Statistics Calculated')


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

def create_plots(pipeline):
	raise NotImplementedError('Not implemented yet!')

def create_plots(pipeline):
	raise NotImplementedError('Not implemented yet!')


BASE = [ load_tumor, random_spherical_samples, calculate_statistics, save_statistics ]
