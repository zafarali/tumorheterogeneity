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
	
class EllipsoidSampler(Sampler):
	def __init__(self, tumor):
		"""
			Samples Ellipsoids from a Tumor
			@params:
				tumor: a Tumor object
		"""
		Sampler.__init__(self, tumor)
		
	def sample(self, radii=(1,1,1), centre=(0,0,0), as_cells = False):
		"""
			Samples an ellipsoid in the tumor and return the cells
			@params:
				radii = a tuple of length 3 of the three axes length
				centre = a tuple of length 3 of the centre position
				as_cells [= False] returns Cell objects if set to true.
			@returns:
				a sample of the tumor that is within the @params
		"""
		
		# using the equation of an ellipsoid:
		# (x-x0)^2/a^2 + (y-y0)^2/b^2 + (z-z0)^2/c^2 = 1 
		sample_indicies = np.sum(((self.cell_positions - np.array(centre))**2)/np.array(radii)**2, axis=1) <= 1
		
		if not as_cells:
			return self.cell_data[sample_indicies,:]
		else:
			raise NotImplementedError('Not yet implemented')
			
	