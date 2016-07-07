from collections import namedtuple
from Statistics import centre_of_mass
import csv
import pandas as pd

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
		genotypes.append( Genotype( int(original_id), int(parent_genotype), int(n_resistant), \
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
		self.cells = cells
		self.genotypes = genotypes
		self.number_of_cells = len(cells)
		self.number_of_genotypes = len(genotypes)
		self.COM = centre_of_mass(self.cells)
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
	def prepare_CTC(sample, genomes):
		
		cell_positions = sample[:,:3]
		genome_ids = sample[:,3]
		tm = Tumor(sample, genomes)
		# need to calculate the COM of the sample and correct for that
		# COM = centre_of_mass(cell_positions)    
		
		cell_positions_corrected = cell_positions-COM

		tm.cells[:,:3] = cell_positions_corrected

		genome_counts = Counter(genome_ids.astype(int))
		# reprocess the way genomes are stored for reseeding
		genome_dict = {}
		for gid, genome in zip(genome_ids, genomes):
			genome_dict[int(gid)] = genome

		new_mapping = { int(gid):new_id for new_id, gid in enumerate(genome_dict.keys()) }
		
		genomes_modified = [0]*len(new_mapping.keys())
		for original_id, new_id in new_mapping.items():

			original_genome = genomes[original_id]
			new_frequency = genome_counts[original_id]
			original_genome.frequency = new_frequency
			genomes_modified[new_id] = original_genome

		tm.genotypes = genomes_modified

	def cells_to_save(self):
		# cell positions and genome ids
		to_return = []
		for gid, cell in zip(genome_ids,cell_positions_corrected):
			x,y,z = cell
			to_return.append([x,y,z, new_mapping[int(gid)]])

		return to_return

	def genomes_to_string(self):

		to_return = []
		for gid,genome in genome_dict.items():
			string_builder = str(new_mapping[int(gid)])+' '+str(genome_counts[int(gid)])+' '+str(genome.n_resistant)+' '+str(genome.n_driver)+' '
			for snp in genome.snps:
				string_builder+=str(snp)+','
			string_builder+='-1'
			
			to_return.append(string_builder)

		return to_return
		  