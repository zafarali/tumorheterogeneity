import numpy as np
from BaseObjects import Tumor
from Statistics import centre_of_mass
from scipy.spatial import cKDTree

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
				return self.cell_data[sample_indicies,:3], self.tumor.get_genotypes(self.cell_genotypes[sample_indicies])
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

class KDTSphericalSampler(Sampler):
	def __init__(self, tumor):
		"""
			Samples a Sphere from the tumor using a fast KDTree implementation
			@params:
				tumor: a Tumor Object
		"""
		Sampler.__init__(self, tumor)
		self.kdt = cKDTree(self.cell_positions, leafsize=1000)

	def sample(self, radius=1, centre=(0,0,0), with_genotypes=False):
		"""
			Samples a Sphere from the Tumor
			@params:
				radius: radius of the sphere
				centre: centre of the sphere
				with_genotypes [=False]: returns genotpe indicies as well
		"""
		sample_indicies = self.kdt.query_ball_point(centre, r=radius)
		if with_genotypes:
			return self.cell_data[sample_indicies, :3	], self.tumor.get_genotypes(self.cell_genotypes[sample_indicies])
		else:
			return self.cell_data[sample_indicies, :]

	def sample_fixed_points(self, n=1, centre=(0,0,0), with_genotypes=True):
		"""
			samples a fixed number of points from the tumor
			@params:
				n [=1]: the number of points you want
				centre [=(0,0,0)]: the location around which to sample
				with_genotypes [=True]: returns the genotypes as well
			@returns:
				if with_genotypes:
					COM[array], cell_positions[array] , genotypes[list[genotypes]]
		"""
		_, sample_indicies = self.kdt.query(centre, n) # returns distances, idx

		# make sure we got the right number of cells.
		assert sample_indicies.shape[0] == n, 'got ' + str(sample_indicies.shape[0]) + ' cells when asked for ' + str(n)

		cell_positions = self.cell_data[sample_indicies, :3] # cell positions.

		# calculate Center Of Mass from this
		_COM = centre_of_mass(cell_positions)

		if with_genotypes:
			return _COM, cell_positions, self.tumor.get_genotypes(self.cell_genotypes[sample_indicies])
		else:
			raise NotImplementedError('This method is not available for sample_fixed_points.')

