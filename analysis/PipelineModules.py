# modules for pipeline

from BaseObjects import Tumor
from Samplers import SphericalSampler
import Statistics


# generates a coordinate
def generate_coordinate(max_value, min_value):
	while True:
		yield (max_value-min_value)*np.random.random() + max_value

def generate_angle():
	while True:
		yield np.random.random() * 2* np.pi


def load_tumor(pipeline):
	"""
		Loads Tumor Into the Pipeline
	"""
	pipeline.tumor = Tumor(pipeline.FILES['cells'], pipeline.FILES['genotypes'])
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
	rand_coordinate = generate_coordinate( np.max( sampler.cell_positions ), np.min( sampler.cell_positions ) )

	pipeline.print2('Begining Tumor Sampling')
	for radius in pipeline.specs['RADII']:
		# generate 3 random coodinate
		pipeline.print2('Sampling radius:',radius)
		for i in xrange( pipeline.specs['repeats'] or 10 ):
			# generate a new coordinate
			centre = ( rand_coordinate.next() , rand_coordinate.next() , rand_coordinate.next() )
			sample = sampler.sample(radius=radius, centre=centre, with_genotypes=True)
			# insert a tuple of (radius, center, sample)
			samples.append( ( radius, centre, sample ) )

	pipeline.samples = samples

def calculate_statistics(pipeline):
	"""
		Calculates Statistics
	"""
	if not pipeline.samples:
		raise Exception('Tumor must be sampled first!')

	stats = [['radius', 'distance_to_COM', 'sample_size', 'S', 'S/H', 'Tajimas D', 'D']]

	for radius, centre, sample in pipeline.samples:
		# first index contains the genotypes
		distance_to_COM = np.sqrt(np.sum((np.arrays(centre) - np.array(pipeline.COM))**2))
		
		SNP_counts = Statistics.SNP_count(sample[1])
		n = len(sample[1])
		pi, S, SH, D = Statistics.tajimas_D(SNP_counts, n, return_parts=True)

		stats.append( [radius, distance_to_COM, n, pi, S, SH, D, pi-SH ] )

	pipeline.stats = stats

