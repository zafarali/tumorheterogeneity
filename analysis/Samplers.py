import numpy as np
from BaseObjects import Tumor

class Sampler(object):
	def __init__(self, tumor):
		"""
			An object to sample a tumor
			@params:
				tumor: a Tumor object.
		"""
		self.tumor = tumor
		self.cell_data = np.array(self.tumor.cells, dtype=np.float)
		self.cell_positions = self.cell_data[:,0:3].astype(int)
		self.cell_genotypes = self.cell_data[:,3].astype(int)
	def all_cells(self, as_cells=False):
		"""
			Returns all the cells in the tumor
			@params:
				as_cells [=False] 
					by default it will return a numpy array of data
					if this flag is set to True then returns Cell objects
		"""
		return self.tumor.cells if as_cells else self.cell_data
	# @staticmethod
	# def from_file(cell_file, genotype_file):
	# 	return Sampler(Tumor.from_files(cell_file, genotype_file))
	
class EllipsoidSampler(Sampler):
	def __init__(self, tumor):
		"""
			Samples Ellipsoids from a Tumor
			@params:
				tumor: a Tumor object
		"""
		Sampler.__init__(self, tumor)
		
	def sample(self, radii=(1,1,1), centre=(0,0,0), as_cells = False, with_genotypes=False):
		"""
			Samples an ellipsoid in the tumor and return the cells
			@params:
				radii = a tuple of length 3 of the three axes length
				centre = a tuple of length 3 of the centre position
				as_cells [= False] returns Cell objects if set to true.
				with_genotypes [=False] also returns genotype indices
			@returns:
				a sample of the tumor that is within the @params
		"""
		
		# using the equation of an ellipsoid:
		# (x-x0)^2/a^2 + (y-y0)^2/b^2 + (z-z0)^2/c^2 = 1 
		sample_indicies = np.sum(((self.cell_positions - np.array(centre))**2)/np.array(radii)**2, axis=1) <= 1
		
		if not as_cells:
			if with_genotypes:
				return self.cell_data[sample_indicies,:3], self.cell_genotypes[sample_indicies] 
			else:
				return self.cell_data[sample_indicies,:]
		else:
			raise NotImplementedError('This hasn\'t been implemented yet and is not planned.')

class SphericalSampler(EllipsoidSampler):
	def __init__(self, tumor):
		"""
			Samples spheres from a Tumor
			@params:
				tumor: a Tumor object
		"""
		Sampler.__init__(self, tumor)

	def sample(self, radius=1, centre=(0,0,0), as_cells=False, with_genotypes=False):
		"""
			Samples a Sphere from the tumor
			@params:
				radius: the radius of the sphere
				centre: the centre of the sphere
				as_cells[=False]: returns cells rather than numpy array
				with_genotypes [=False]: also returns genotype indices
		"""
		return super(SphericalSampler, self).sample(radii=(radius, )*3, centre=centre, as_cells=as_cells, with_genotypes=with_genotypes)

