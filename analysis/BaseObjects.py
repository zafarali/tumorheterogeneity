from collections import namedtuple, Counter
from Statistics import centre_of_mass
import csv
import pandas as pd
import numpy as np
# make this really fast:
import multiprocessing.dummy as multiprocessing

D_PM = 1<<30 # convert to driver
R_PM = 1<<31 # convert to resistant
L_PM = (1<<30) - 1 # return from driver
L_PM2 = (1<<31) - 1 # return from resistant


Cell_ = namedtuple('Cell', 'x y z genotype')
Cell_.__new__.__defaults__ = ( 0, 0, 0, -1 )

class Cell(Cell_):
  """
	Generic holder for cell details.
	@params:
	  id = cell ID
	  x,y,z = coordinates of the cell in the tumor
	  genotype = the genotype ID of this cell
  """
  @staticmethod
  def load_cells(file_name):
	"""
	Returns a list of Cells
	@params:
	  file_name
	"""

	cell_info = pd.read_csv(file_name,sep=' ').values
	return map(lambda cell:Cell(*cell),cell_info)
	

Genotype_ = namedtuple('Genotype', 'original_id parent_genotype n_resistant n_driver frequency snps')
Genotype_.__new__.__defaults__ = (-1, -1, 0, 0, 1, [])

class Genotype(Genotype_):
  """
	Generic holder for Genotypes
	@params
	  orignal_id = the original ID of this genome during the simulation
	  parent_genotype = the parent_ID of this genome
	  n_resistant = number of resistant mutations
	  n_driver = number of driver mutations
	  frequency = the number of cells with this genotype
	  snps = the sequence of SNPs in this genotype.
  """
  @staticmethod
  def load_genotypes(file_name):
	"""
	Returns a list of genotypes
	@params:
	  file_name
	"""
	genotypes = []
	with open( file_name, 'r' ) as genotype_file:
	  for row in csv.reader( genotype_file ):
		details, sequence = row[0].split('\t') # split on the tab
		original_id, parent_genotype, numbers, frequency = details.split('  ')
		n_resistant, n_driver = numbers.split(' ')
		sequence = filter( lambda SNP: len(SNP), sequence.split(' ') ) # remove blank SNPs
		sequence = map( int , sequence ) # convert everything into ints
		genotypes.append(Genotype( int(original_id), int(parent_genotype), int(n_resistant), \
					   int(n_driver), int(frequency), sequence ) )
	return genotypes


class Tumor(object):
	def __init__( self, cells, genotypes ):
		"""
			Generic tumor object containing cells and genotypes.
			@params:
				cells: a list of cells
				genotypes: a list of genotypes
		"""
		self.cells = np.copy(cells)
		self.genotypes = list(genotypes)
		self.number_of_cells = len(cells)
		self.number_of_genotypes = len(genotypes)
		self.COM = centre_of_mass(self.cells)
		self._snp_lookup_efficient = False
		# print self.drivers
		
	def get_genotypes( self, genotype_indicies ):
		"""
			Returns a list of Genotype Objects based on the genotype_indicies
			@params:
				genotype_indicies: a list of indicies of genotypes required
		"""	
		return map(self.get_genotype, genotype_indicies)

	def get_genotype( self, genotype_index ):
		"""
			Returns a Genotype object based on the genotype_index
			@params:
				genotype_index: the index of the genotype wanted
		"""
		return self.genotypes[genotype_index]

	def _make_snp_lookup_efficient(self):
		"""
		Makes looking up snps efficient.
		"""
		print('Making snp lookup efficient.')
		snp_to_genotypes = {}

		# loop through genomes
		for gid, genotype in enumerate(self.genotypes):
			# loop through snps
			for snp_id in genotype.snps:
				# store in a list of genomes corresponding to that snp
				# the id of the genome is found at gid
				genotypes_with_snp = snp_to_genotypes.get(snp_id, [])
				genotypes_with_snp.append(gid)
				snp_to_genotypes[snp_id] = genotypes_with_snp

		self._snp_lookup_efficient = True
		self.snp_to_genotypes = snp_to_genotypes
		

	def genotype_idx_with_snp(self, snp_id ):
		"""
		Looks up the genotypes with the snp_id
		@params:
			snp_id: the snp_id to lookup
		@returns:
			list of tuples containing (index_in_internal_list, original_id)
		"""
		if not self._snp_lookup_efficient:
			self._make_snp_lookup_efficient()
		
		return self.snp_to_genotypes[snp_id]

	def _cells_with_genotype(self, genotype_index):
		"""
		Internal lookup function
		"""
		return set(np.where(self.cells[:, 3]==genotype_index)[0].tolist())

	def cells_with_snp(self, snp_id):
		"""
		Looks up the cells that contain snp_id
		@params:
			snp_id: the snp_id to lookup
		@returns:
			cell_idx a list containing ids
			np.array containing cell information
		"""
		if not self._snp_lookup_efficient:
			self._make_snp_lookup_efficient()

		indices = self.genotype_idx_with_snp(snp_id)
		cell_ids = set()
		
		pool = multiprocessing.Pool() 
		[cell_ids.update(cells_with_genotype) for cells_with_genotype in pool.map(self._cells_with_genotype, indices)]
		
		cell_ids = list(cell_ids)

		return indices, cell_ids, self.cells[cell_ids]


	@staticmethod
	def from_files( cell_file, genome_file ):
		"""
			Returns a tumor from cell and genotype data files
			@params:
				cell_file : file containing cell positions and genotypes
				genome_file : file containing genome data
		"""

		return Tumor( Cell.load_cells( cell_file ), \
					 Genotype.load_genotypes( genome_file ) )

class CTC(Tumor):
	@staticmethod
	def prepare_CTC(sample, genomes, dist=False,nhaplotypes=False,notes=''):
		
		# get important data
		cell_positions = np.copy(sample[:,:3])
		genome_ids = np.copy(sample[:,3].astype(int))
		
		# create a CTC
		tm = CTC(sample, genomes)

		# correct for the position relative to COM
		cell_positions_corrected = cell_positions-tm.COM
		tm.cells[:,:3] = cell_positions_corrected

		# calculate the new frequencies
		genome_counts = Counter(genome_ids.astype(int))

		# reprocess the way genomes are stored for reseeding
		genome_dict = {} #temporary holding for genomes.
		for gid, genome in zip(genome_ids, genomes):
			genome_dict[int(gid)] = genome

		# obtain the new ids for genomes
		new_mapping = { int(gid):new_id for new_id, gid in enumerate(genome_dict.keys()) }
		
		# create an empty array of the same size.
		genomes_modified = [0]*len(new_mapping.keys())


		# now rehash the snps
		snp_mapping = set()
		snp_types = {}
		new_snp_sequences = {}

		for gid, genotype in genome_dict.items():
		    for snp in genotype.snps:
		#         print 'SNP:',snp
		        if snp & D_PM:
		#             print('driver', snp & L_PM)
		            snp_types[snp] = 'driver'
		            if not snp & L_PM in snp_mapping:
		                snp_mapping.update([snp & L_PM])
		        elif snp & R_PM:
		#             print('resistant', snp & L_PM2)
		            snp_types[snp] = 'resistant'
		            if not snp & L_PM2 in snp_mapping:
		                snp_mapping.update([snp & L_PM2])
		        else:
		            if not snp in snp_mapping:
		                snp_mapping.update([snp])
		#         print '---'

		    
		snp_mapping = { mut: newid for newid, mut in enumerate(snp_mapping)}

		for gid, genotype in genome_dict.items():
		    new_snp_sequence = []
		    
		    for snp in genotype.snps:
		        snp_type = snp_types.get(snp, False)
		        if snp_type == 'resistant':
		            mapped = snp_mapping.get(snp & L_PM2)
		            new_snp_sequence.append(mapped | R_PM)
		        elif snp_type == 'driver':
		            mapped = snp_mapping.get( snp & L_PM )
		            new_snp_sequence.append(mapped | D_PM)
		        else:
		            mapped = snp_mapping.get(snp)
		            new_snp_sequence.append(mapped)
		            
		    new_snp_sequences[gid] = new_snp_sequence

		# now update what will be in the new CTC object.
		for original_id, new_id in new_mapping.items():

			genome_ids[np.where(genome_ids==original_id)] = new_id

			original_genome = genome_dict[original_id]
			new_frequency = genome_counts[original_id]
			genomes_modified[new_id] = Genotype(
				original_genome.original_id,
				-1,
				original_genome.n_resistant,
				original_genome.n_driver,
				new_frequency,
				new_snp_sequences[original_id]
				)

		tm.snp_mapping = snp_mapping
		tm.new_snp_sequences = new_snp_sequences
		tm.genotypes = genomes_modified
		tm.cells[:,3] = genome_ids
		tm.new_mapping = new_mapping
		tm.genome_counts = genome_counts
		tm.genome_dict = genome_dict
		tm.dist = dist
		tm.notes = notes
		tm.nhaplotypes = nhaplotypes
		return tm

	def cells_to_string(self):
		# cell positions and genome ids
		to_return = []
		for gid, cell in zip(self.cells[:,3], self.cells[:,:3]):
			x,y,z = cell
			to_return.append( str(x)+' '+str(y)+' '+str(z)+' '+str(int(gid)))

		return to_return

	def genomes_to_string(self):

		to_return = []
		for gid,genome in self.genome_dict.items():
			new_mapping = self.new_mapping[int(gid)]
			string_builder = str(self.new_mapping[int(gid)])+' '+str(self.genome_counts[int(gid)])+' '+str(genome.n_resistant)+' '+str(genome.n_driver)+' '
			for snp in self.genotypes[new_mapping].snps:
				string_builder+=str(snp)+','
			if int(new_mapping)+1 == len(self.genome_dict.keys()):
				string_builder+='-2'
			else:
				string_builder += '-1'
			to_return.append(string_builder)

		largest_snp = max(self.snp_mapping.values()) if len(self.snp_mapping.values()) > 0 else 0
		to_return.append(str(largest_snp+1))

		return to_return

	def details(self):
		string_builder = 'notes:\nNumber of Cells:'+str(self.number_of_cells)+'\nDistance from COM:'+str(self.dist)

		string_builder += '\nNumber of Unique Haplotypes:'+str(self.nhaplotypes) + '\nothernotes:'+str(self.notes)

		return string_builder

