# modules for pipeline
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import csv
import os
import itertools
from BaseObjects import Tumor, CTC
from Samplers import SphericalSampler, KDTSphericalSampler
import Statistics
import random
import time

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

def create_sample_directory(pipeline):
	"""
		Creates a directory to store the samples
	"""
	time_info = '_'.join('_'.join(time.asctime().split(':')).split(' '))
	try:
		pipeline.FILES['sample_directory'] = pipeline.FILES['out_directory'] + '/samples'
		# create the output directory and save specs and sampling there
		if not os.path.exists( pipeline.FILES['sample_directory'] ):
			os.makedirs( pipeline.FILES['sample_directory'] )
			pipeline.print2('sample directory created')
	except Exception as e:
		raise Exception('Exception occured in create_sample_directory'+str(e))


def create_kdsampler(pipeline):
	pipeline.print2('Creating KDTSphericalSampler')
	pipeline.sampler = KDTSphericalSampler(pipeline.tumor)
	pipeline.print2('KDTSphericalSampler created.')

def inline_statistics(pipeline):

	stats = [['radius', 'distance_to_COM', 'sample_size', \
	'pi', 'S', 'S/H', 'Tajimas D', 'D', \
	'num_SNPs', 'proportion_driver', 'unique_driver_propotion',\
	'num_drivers', 'unique_combos','unique_combos_snps', 'uniqe_res_combos',\
	'mean_drivers', 'mean_SNPs', 'driver_enrichment']]

	rand_coordinate = sample_coordinate(pipeline.sampler.cell_positions, deviate=True)
	number_of_samples = pipeline.specs.get('repeats', 500)
	to_save = pipeline.specs.get('save_samples', False)

	for radius in pipeline.specs['RADII']:
		pipeline.print2('Sampling radius:'+str(radius))
		i = 0 # holds the number of samples so far
		# for i in xrange( number_of_samples ):
		i2=0 # holds the number of attempts so far
		while i != number_of_samples:
			try:
				if i2 > 2000:
					pipeline.print2('Did not find enough samples for radius='+str(radius)+'.')
					break
				# generate a new coordinate
				centre = rand_coordinate.next()

				# conduct the sample
				sample = pipeline.sampler.sample(radius=radius, centre=centre, with_genotypes=True)
				cell_data = pipeline.sampler.sample(radius=radius, centre=centre)
				genotypes = sample[1]
				# calculate statistics

				distance_to_COM = np.sqrt( np.sum( ( np.array(centre) - np.array(pipeline.tumor.COM ) )**2 ) )
				
				# frequency of snps
				SNP_counts = Statistics.SNP_count(sample[1])

				# frequency of driver mutations
				drivers_count = Statistics.drivers_count(sample[1])

				n = len(sample[1])
				
				i2 += 1

				if n < 2:

					if n > 0:
						# snp count
						num_SNPs = len(SNP_counts.keys())
						mean_SNPs = np.sum( SNP_counts.values() ) / float(n)

						num_drivers = len( drivers_count.keys() ) 
						mean_drivers = np.sum( drivers_count.values() ) / float(n)

						pipeline.print2('Only calculatd num snps/ and drivers due to only ' + str(n)+ ' individuals')
						stats.append( [radius, distance_to_COM, \
							n, 0, 0, 0, 0, 0, \
							num_SNPs, 0,\
							0, 0, \
							0, 0, 0,\
							mean_drivers, mean_SNPs, 0 ] )
					else:
						pipeline.print2('Skipped sampling due to only ' + str(n)+ ' individuals')
					continue

				i+=1
				# base statistics pi, S, SH, D
				pi, S, SH, D = Statistics.tajimas_D(SNP_counts, n, return_parts=True)

				# number of SNPs in the sample
				num_SNPs = len( SNP_counts.keys() )

				# total number of SNPs in the sample
				mean_SNPs = np.sum( SNP_counts.values() ) / float(n)

				# number of driver mutations in the sample
				num_drivers = len( drivers_count.keys() ) 

				# total number of driver mutations:
				mean_drivers = np.sum( drivers_count.values() ) / float(n)

				# the number of drivers that are unique / total number of drivers
				unique_drivers = Statistics.unique_driver_proportion(drivers_count)

				# the number of unique haplotypes (combinations of drivers)
				unique_combos = Statistics.unique_driver_combinations(sample[1])
				unique_combos_snps = Statistics.unique_snp_combinations(sample[1])
				unique_combos_res = Statistics.unique_res_combinations(sample[1])

				# the proportion of total SNPs in the data that are drivers
				driver_enrichment = mean_drivers/float(mean_SNPs) if mean_SNPs > 0 else 0
				
				# the number of SNPs that are driver
				driver_to_SNP_ratio = num_drivers/float(num_SNPs) if num_SNPs > 0 else 0

				stats.append( [radius, distance_to_COM, \
					n, pi, S, SH, D, pi-SH, \
					num_SNPs, driver_to_SNP_ratio,\
					unique_drivers, num_drivers, \
					unique_combos, unique_combos_snps, unique_combos_res,\
					mean_drivers, mean_SNPs, driver_enrichment ] )

				if to_save:
					save_loc = pipeline.FILES['sample_directory'] +'/size'+str(n)+'dist'+str(distance_to_COM)+'haplo'+str(unique_combos)+'i'+str(i)
					# create the output directory and save specs and sampling there
					if not os.path.exists( save_loc ):
						os.makedirs( save_loc )
						# pipeline.print2('sample dire directory created')


					ctc = CTC.prepare_CTC(cell_data, genotypes, dist=distance_to_COM,nhaplotypes=unique_combos,notes='snp_combos:'+str(unique_combos_snps))
					
					# using randint() to randomize the file names just in case
					identifier = str(np.random.randint(10000))
					with open(save_loc+'/gen'+identifier+'.dat', 'w') as f:
						for row in ctc.genomes_to_string():
							f.write(row+'\n')
					
					with open(save_loc+'/cells'+identifier+'.dat', 'w') as f:
						for row in ctc.cells_to_string():
							f.write(row+'\n')

					with open(save_loc+'/details'+identifier+'.dat', 'w') as f:
						f.write(ctc.details())

			except Exception as e:
				print("exception occured... graceful skip:"+str(e))
			#end try catch

			# end sampling for loop for certain radius
		
		#end radius for loop

		pipeline.print2( str(number_of_samples) +' samples conducted and statistics too')
	# end

	pipeline.stats = stats
	pipeline.print2('Sampling/Statistics Completed')

def density_plot(pipeline):
	rho, r = np.histogram(np.sqrt(np.sum((pipeline.sampler.cell_positions-pipeline.tumor.COM)**2,axis=1)), bins=np.linspace(0,1000,1000))
	def neighbour_iterator(arr):
		index = 0
		while index < len(arr)-1:
			yield (arr[index], arr[index+1])
			index += 1
	r_meaned = np.mean(np.array(list(neighbour_iterator(r))), axis=1)
	rho2 = rho/ (4*np.pi*r_meaned**2)  
	plt.figure(figsize=(10,7.5))
	plt.plot(r_meaned,rho2)
	plt.xlabel('r')
	plt.ylabel('density')
	plt.title('Tumor Density')
	np.save(pipeline.FILES['out_directory']+'/r_meaned.npy',r_meaned)
	np.save(pipeline.FILES['out_directory']+'/rho2.npy',rho2)
	plt.savefig(pipeline.FILES['out_directory']+'/density_corrected.pdf', bbox_inches='tight')

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
		'num_SNPs', 'proportion_driver', 'unique_driver_propotion',\
		'num_drivers', 'unique_combos', 'mean_drivers', \
		'mean_SNPs', 'driver_enrichment']
	
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

def mutation_count_plots(pipeline):
	pipeline.print2('Plotting mutation counts')
	_mutation_count_plots(pipeline, 'distance_to_COM')
	_mutation_count_plots(pipeline, 'sample_size')
	pipeline.print2('Plotting mutation counts done')

def driver_stats_plots(pipeline):
	pipeline.print2('Plotting driver stats')
	_driver_stats_plots(pipeline, 'distance_to_COM')
	_driver_stats_plots(pipeline, 'sample_size')
	pipeline.print2('Plotting driver stats done')

def pop_gen_plots(pipeline):
	pipeline.print2('Plotting pop_gen stats')
	_pop_gen_stats(pipeline, 'distance_to_COM')
	_pop_gen_stats(pipeline, 'sample_size')
	pipeline.print2('Plotting pop_gen stats done')

def _driver_stats_plots(pipeline, x_label):
	
	printlabel = 'Distance to Center of Mass' if x_label == 'distance_to_COM' else 'Sample Size'

	df = pd.DataFrame(pipeline.stats[1:], columns=pipeline.stats[0])
	plt.figure(figsize=(11,10))
	radii = np.unique(df['radius'].values)
	colors = itertools.cycle( cm.rainbow( np.linspace( 0, 1, len(radii) ) ) )

	plt.subplot(2,2,1)

	plt.title('Singleton Drivers vs Distance')
	for i in radii:
		subsample = df[df['radius'] == i]
		plt.scatter( subsample[x_label], subsample['unique_driver_propotion'], color=next(colors), label='Radius='+str(i))
		plt.xlabel(printlabel)
		plt.ylabel('Proportion of Drivers that are Singletons')

	plt.subplot(2,2,2)

	plt.title('Driver Enrichment vs Distance')
	for i in radii:
		subsample = df[df['radius'] == i]
		plt.scatter( subsample[x_label], subsample['driver_enrichment'], color=next(colors), label='Radius='+str(i))
		plt.xlabel(printlabel)
		plt.ylabel('Mean Drivers per Cell / Mean SNPs per Cell')
		
	plt.subplot(2,2,3)
	plt.title('Unique Haplotypes vs Distance')
	for i in radii:
		subsample = df[df['radius'] == i]
		plt.scatter( subsample[x_label], subsample['unique_combos'], color=next(colors), label='Radius='+str(i))
		plt.xlabel(printlabel)
		plt.ylabel('Number of Unique Haplotypes')

	plt.subplot(2,2,4)

	plt.title('Proportion of SNPs that are Driver')
	for i in radii:
		subsample = df[df['radius'] == i]
		plt.scatter( subsample[x_label], subsample['proportion_driver'], color=next(colors), label='Radius='+str(i))
		plt.xlabel(printlabel)
		plt.ylabel('Number of Drivers / Number of SNPs')

	plt.suptitle('Statistics about Driver SNPs in Samples', fontsize=15)
	plt.legend(ncol=5, bbox_to_anchor=(0.88,2.38))
	plt.savefig( pipeline.FILES['out_directory']+'/driver_stats_'+x_label+'.pdf',\
		bbox_inches='tight')

def _mutation_count_plots(pipeline, x_label):

	printlabel = 'Distance to Center of Mass' if x_label == 'distance_to_COM' else 'Sample Size'


	pipeline.print2('Special Plot being conducted.')

	df = pd.DataFrame(pipeline.stats[1:], columns=pipeline.stats[0])
	radii = np.unique(df['radius'].values)

	plt.figure(figsize=(11,10))
	colors = itertools.cycle( cm.rainbow( np.linspace( 0, 1, len(radii) ) ) )

	plt.subplot(2,2,1)
	plt.title('Mean number of SNPs per cell')

	for i in radii:
		subsample = df[df['radius'] == i]
		plt.scatter( subsample[x_label], subsample['mean_SNPs'], color=next(colors), label='Radius='+str(i))
		plt.xlabel(printlabel)
		plt.ylabel('Mean SNPs per Cell')

	plt.subplot(2,2,2)
	
	plt.title('Mean number of Drivers per cell')
	for i in radii:
		subsample = df[df['radius'] == i]
		plt.scatter( subsample[x_label], subsample['mean_drivers'], color=next(colors), label='Radius='+str(i))
		plt.xlabel(printlabel)
		plt.ylabel('Mean Drivers per Cell')

	plt.subplot(2,2,3)
	
	plt.title('Number of SNPs in the sample vs Distance')
	for i in radii:
		subsample = df[df['radius'] == i]
		plt.scatter( subsample[x_label], subsample['num_SNPs'], color=next(colors), label='Radius='+str(i))
		plt.xlabel(printlabel)
		plt.ylabel('Number of SNPs')

	plt.subplot(2,2,4)
	
	plt.title('Number of drivers in the sample vs Distance')
	for i in radii:
		subsample = df[df['radius'] == i]
		plt.scatter( subsample[x_label], subsample['num_drivers'], color=next(colors), label='Radius='+str(i))
		plt.xlabel(printlabel)
		plt.ylabel('Number of Drivers')
	
	plt.savefig( pipeline.FILES['out_directory']+'/mutation_count_plots_'+x_label+'.pdf',\
		bbox_inches='tight')
	
def _pop_gen_stats(pipeline,x_label):
	plt.figure(figsize=(15,10))
	printlabel = 'Distance to COM' if x_label == 'distance_to_COM' else 'Sample Size'
	df = pd.DataFrame(pipeline.stats[1:], columns=pipeline.stats[0])

	radii = np.unique(df['radius'].values)
	colors = itertools.cycle( cm.rainbow( np.linspace( 0, 1, len(radii) ) ) )

	plt.subplot(2,2,1)
	plt.title('D vs Radius')
	for i in radii:
		subsample = df[df['radius'] == i]
		plt.scatter( subsample[x_label], subsample['D'], color=next(colors))
		plt.xlabel(printlabel)
		plt.ylabel('D')
		
	plt.subplot(2,2,2)
	plt.title('Tajimas D vs Radius')
	for i in radii:
		subsample = df[df['radius'] == i]
		plt.scatter( subsample[x_label], subsample['Tajimas D'], color=next(colors))
		plt.xlabel(printlabel)
		plt.ylabel('Tajimas D')
		
	plt.subplot(2,2,3)
	plt.title('S/H vs Radius')
	for i in radii:
		subsample = df[df['radius'] == i]
		plt.scatter( subsample[x_label], subsample['S/H'], color=next(colors))
		plt.xlabel(printlabel)
		plt.ylabel('Normalized Number of Segregating Sites (S/H)')
		
	plt.subplot(2,2,4)
	plt.title('Pi vs Radius')
	for i in radii:
		subsample = df[df['radius'] == i]
		plt.scatter( subsample[x_label], subsample['pi'], color=next(colors))
		plt.xlabel(printlabel)
		plt.ylabel('Proportion of Pairwise Differences (Pi)')
		
	plt.suptitle('Statistics for Neutrality of Samples', fontsize=15)
	plt.savefig( pipeline.FILES['out_directory']+'/pop_gen_stats'+x_label+'.pdf',\
		bbox_inches='tight')

def marginal_counts_ordered(pipeline):
	pipeline.print2('Calculating ordered marginals...')
	sampler = pipeline.sampler
	rand_coordinate = sample_coordinate(sampler.cell_positions, deviate=True)
	all_deltas = []
	all_s_list = []
	all_d_list = []
	all_csnp_list = []
	all_cdrv_list = []

	for k in range(1000):
		coord = rand_coordinate.next()

		# here we update the "coord" based on the actual COM of the sample
		# get the 40 nearest cells so that we can sequentially increase the number of cells we use in the next part of the loop
		coord, cell_positions, cell_genotypes = sampler.sample_fixed_points(40, centre=coord, with_genotypes=True)
		samples = [ cell_positions, cell_genotypes ]

		# the distance from each cell to the COM of the sample
		distance_to_COM = np.sqrt( np.sum( ( np.array(samples[0]) - np.array( coord ) )**2, axis=1 ) )

		# use this to sort the cells in increasing distance
		# TODO: check if this is already done for us in the KdTree sampler?
		samples[0] = distance_to_COM.tolist()
		distances_sorted, samples_sorted = zip(*sorted(zip(*samples),key=lambda x:x[0])) # sorts the cells according to distance from COM

		# first entry is the distance between the tumor COM and sample COM
		delta = [ np.sqrt( np.sum( ( np.array(coord) - np.array( pipeline.tumor.COM ) )**2) ) ] 


		S_list = [ ] # holds the original S values, S=number of somatic mutations
		D_list = [ ] # holds the original D values, D=number of driver mutations
		combosnp_list = [ ] # holds the number of genotyps
		combodrv_list = [ ] # holds the number of driver genotypes
		try:
			for i in range(1,33): # loop through all sample sizes
				samples_selected = samples_sorted[:i]

				assert len(samples_selected) >= i, 'Did not have enough cells to calculate statistics. Wanted: '+str(i)+', Got:'+str(len(samples_selected))

				# calculate statistics
				S = len(Statistics.SNP_count(samples_selected).keys())
				D = len(Statistics.drivers_count(samples_selected).keys())
				S_list.append(S)
				D_list.append(D)

				if i == 1:
					delta.append(S)
				else:
					delta.append(float(S - S_list[-2])/float(S_list[-1]))

				combosnp = Statistics.unique_snp_combinations(samples_sorted[:i])
				combodrv = Statistics.unique_driver_combinations(samples_sorted[:i])

				combosnp_list.append(combosnp)
				combodrv_list.append(combodrv)

		except Exception as e:
			pipeline.print2('Exception occured in marginal_counts: '+str(e))
			continue
	#     print delta
		all_deltas.append(delta)
		all_s_list.append(S_list)
		all_d_list.append(D_list)
		all_csnp_list.append(combosnp_list)
		all_cdrv_list.append(combodrv_list)


	pipeline.print2('1000 ordered marginal trajectories calculated.')

	all_deltas = np.array(all_deltas)
	all_s_list = np.array(all_s_list)
	all_d_list = np.array(all_d_list)
	all_csnp_list = np.array(all_csnp_list)
	all_cdrv_list = np.array(all_cdrv_list)


	np.save(pipeline.FILES['out_directory']+'/deltas_ordered.npy',all_deltas)
	np.save(pipeline.FILES['out_directory']+'/S_list_ordered.npy',all_s_list)
	np.save(pipeline.FILES['out_directory']+'/D_list_ordered.npy',all_d_list)
	np.save(pipeline.FILES['out_directory']+'/csnp_list_ordered.npy',all_csnp_list)
	np.save(pipeline.FILES['out_directory']+'/cdrv_list_ordered.npy',all_cdrv_list)
	
def big_samples(pipeline):
	pipeline.print2('Calculating big samples...')
	sampler = pipeline.sampler
	rand_coordinate = sample_coordinate(sampler.cell_positions, deviate=True)
	
	# each of these are kxi matrix with the values of S
	# the index of the columns correspond to the size of the sample
	S_000_to_save = []
	S_001_to_save = []
	S_010_to_save = []

	# each of these are kxi matrices with the distance to com of the sample
	# the index of the columns correspond to the size of the sample in S
	d_000_to_save = []
	d_001_to_save = []
	d_010_to_save = []


	for k in range(1000):

		# global coordinate to use for sampling
		coord = rand_coordinate.next()

		S_000 = []
		S_001 = []
		S_010 = []

		d_000 = []
		d_001 = []
		d_010 = []

		try:
			for i in [100, 1000, 10000, 20000]: 

				# get i cells at this position
				sample_COM, cell_positions, samples_sorted = sampler.sample_fixed_points(i, centre=coord)

				assert len(samples_sorted) >= i, 'Did not have enough cells to calculate statistics. Wanted: ' + str(i) + ', Got:' + str( len(samples_sorted) )

				# distance from the COM of the sample to the COM of the tumor
				d_cm = np.sqrt( np.sum( ( np.array(sample_COM) - np.array(pipeline.tumor.COM) )**2) )

				try:
					# get the statistics necessary
					S00 = len( Statistics.SNP_count(samples_sorted, min_freq=0).keys() )
					S001 = len( Statistics.SNP_count(samples_sorted, min_freq=0.01).keys() )
					S01 = len( Statistics.SNP_count(samples_sorted, min_freq=0.1).keys() )

					# save into our larger array
					S_000.append( S00 )
					S_001.append( S001 )
					S_010.append( S01 )

					d_000.append( d_cm )
					d_001.append( d_cm )
					d_010.append( d_cm )

				except Exception as e2:
					pipeline.print2('failed to sample n='+str(i)+' exception:'+str(e2))

			# save relevant information
			S_000_to_save.append( S_000 )
			S_001_to_save.append( S_001 )
			S_010_to_save.append( S_010 )

			d_000_to_save.append( d_000 )
			d_001_to_save.append( d_001 )
			d_010_to_save.append( d_010 )

		except Exception as e:
			pipeline.print2('Exception occured in big_samples: '+str(e))
			continue

	pipeline.print2('100 big samples calculated!')

	np.save(pipeline.FILES['out_directory']+'/S_list00_big.npy', np.array(S_000_to_save))
	np.save(pipeline.FILES['out_directory']+'/S_list001_big.npy', np.array(S_001_to_save))
	np.save(pipeline.FILES['out_directory']+'/S_list01_big.npy', np.array(S_010_to_save))


	np.save(pipeline.FILES['out_directory']+'/dist_00_big.npy', np.array(d_000_to_save))
	np.save(pipeline.FILES['out_directory']+'/dist_001_big.npy', np.array(d_001_to_save))
	np.save(pipeline.FILES['out_directory']+'/dist_01_big.npy', np.array(d_010_to_save))
	


def marginal_counts_unordered(pipeline):
	pipeline.print2('Calculating unordered marginals...')
	sampler = pipeline.sampler
	rand_coordinate = sample_coordinate(sampler.cell_positions, deviate=True)
	all_deltas = []
	all_s_list = []
	all_d_list = []
	all_epsilons = []
	all_csnp_list = []
	all_cdrv_list = []


	for k in range(1000):
		
		## holds the marginal increases
		# first entry is the distance from COM of the sample in question
		# here there is no concept of COM since all the cells are randomly drawn.
		delta = [ 0 ] 
		S_list = [ ] # holds the original S values
		D_list = [ ] # holds the original D values
		combosnp_list = [ ] # holds the number of genotyps
		combodrv_list = [ ] # holds the number of driver genotypes

		try:
			for i in range(1,33): # loop through all sample sizes

				# randomly i cells sample from the tumor
				sample_indicies = random.sample(np.arange( 0, sampler.cell_data.shape[0] ), i)
				
				cell_locations, genotypes = sampler.cell_data[ sample_indicies, :3 ], \
				sampler.tumor.get_genotypes(sampler.cell_genotypes[sample_indicies])

				S = len(Statistics.SNP_count(genotypes).keys())
				D = len(Statistics.drivers_count(genotypes).keys())
				S_list.append(S)
				D_list.append(D)
				if i == 1:
					delta.append(S)
				else:
					delta.append(float(S - S_list[-2])/float(S_list[-1]))


				combosnp = Statistics.unique_snp_combinations(genotypes)
				combodrv = Statistics.unique_driver_combinations(genotypes)

				combosnp_list.append(combosnp)
				combodrv_list.append(combodrv)


		except Exception as e:
			pipeline.print2('Exception occured in marginal_counts: '+str(e))
			continue

		all_deltas.append(delta)
		all_s_list.append(S_list)
		all_d_list.append(D_list)
		all_csnp_list.append(combosnp_list)
		all_cdrv_list.append(combodrv_list)

	
	all_deltas = np.array(all_deltas)
	all_s_list = np.array(all_s_list)
	all_d_list = np.array(all_d_list)
	all_csnp_list = np.array(all_csnp_list)
	all_cdrv_list = np.array(all_cdrv_list)


	np.save(pipeline.FILES['out_directory']+'/deltas_unordered.npy',all_deltas)
	np.save(pipeline.FILES['out_directory']+'/S_list_unordered.npy',all_s_list)
	np.save(pipeline.FILES['out_directory']+'/D_list_unordered.npy',all_d_list)
	np.save(pipeline.FILES['out_directory']+'/csnp_list_unordered.npy',all_csnp_list)
	np.save(pipeline.FILES['out_directory']+'/cdrv_list_unordered.npy',all_cdrv_list)


	pipeline.print2('1000 unordered marginal trajectories calculated.')




BASE = [ load_tumor, random_spherical_samples, calculate_statistics, save_statistics ]
KD_SAMPLING = [ load_tumor, create_sample_directory, create_kdsampler, inline_statistics, save_statistics,\
	# all_plot_combos, mutation_count_plots, driver_stats_plots, pop_gen_plots, density_plot, \
	density_plot, \
	marginal_counts_ordered, marginal_counts_unordered, big_samples ]

ONLY_MARGINALS = [ load_tumor, create_kdsampler, marginal_counts_unordered, marginal_counts_ordered, density_plot, big_samples]

ONLY_BIG = [ load_tumor, create_kdsampler, big_samples, density_plot ]

