import numpy as np
import Statistics
from glob import glob
import sys
import time 
import os
import json


## set of tools to analyze a sample

class Pipeline(object):
	def __init__(self, folder, modules=[], **kwargs):
		self.modules = [ open_files, create_output_directory, spec_IO ]
		self.modules.extend(modules) # add new modules
		self.folder = folder # folder where the analysis is taking place
		self.kwargs = kwargs # extra arguments
		self.event_sequence = 0 # for printing
		self.FILES = {'specs':glob(self.folder+'/specs.json')[0]} # stores files in this pipeline
		self.specs = {} # stores specifications of this pipeline

	def execute(self):
		for module in self.modules:
			try:
				module(self)
			except Exception as e:
				self.print2('An exception occured in module:'+str(module)+'\n'+str(e))
				return False
				break
		return self

	def print2(self, string):
		string = str(string)
		print ' Out [' + str(self.event_sequence) + ']: ' + string
		sys.stdout.flush()
		self.event_sequence += 1

def open_files(pipeline):
	"""
		Opens files.
	"""
	try:
		pipeline.FILES.update({
			'cells' : glob(pipeline.folder+'/cells*')[0],
			'genotypes' : glob(pipeline.folder+'/genotypes*')[0],
		})
	except IndexError as e:
		sys.exit('(!) Certain files were missing')

	pipeline.print2('All files present')

def create_output_directory(pipeline):
	"""
		Creates an output director
	"""
	time_info = '_'.join('_'.join(time.asctime().split(':')).split(' '))
	try:
		pipeline.FILES['out_directory'] = pipeline.folder+'/pipe_out_'+time_info
		# create the output directory and save specs and sampling there
		if not os.path.exists( pipeline.FILES['out_directory'] ):
			os.makedirs( pipeline.FILES['out_directory'] )
			pipeline.print2('output directory created')
	except Exception as e:
		raise Exception('Exception occured in create_output_directory'+str(e))


def spec_IO(pipeline):
	"""
		Does IO for specs
	"""
	with open( pipeline.FILES['specs'] ,'r' ) as f:
		pipeline.specs = json.load( f )
		pipeline.print2('Sampling specifications loaded')

	with open(pipeline.FILES['out_directory']+'/specs.json', 'w') as f:
		json.dump(pipeline.specs, f)
		pipeline.print2('sampling specs saved in folder.')
