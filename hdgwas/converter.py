import os
import tables
import h5py
from numpy import genfromtxt
import bitarray as ba
import numpy as np
import gc
import subprocess
from tools import Timer
import pandas as pd
from data import MINIMACHDF5Folder
import shutil

class Genotype(object):
	def __init__(self):
		self.file_name = None
		self.reader=None
		self.probes = None
		self.individuals = None
		self.genotype = None
		self.out=None


class GenotypeHDF5(Genotype):
	def __init__(self, name,force=True):
		super(GenotypeHDF5, self).__init__()
		self.h5_name = '%s.h5' % name
		self.file_name = name
		self.pytable_filters = tables.Filters(complevel=9, complib='zlib')
		self.h5_gen_file = None
		self.h5_ind_file = None
		self.h5_pr_file = None
		self.gen_iter=0


	def write_data(self, type, overwrite=True):

		type_dic={'gen':['genotype',self.h5_gen_file],
				  'ind':['individuals',self.h5_ind_file],
				  'pr':['probes',self.h5_pr_file]}

		if (not overwrite) and os.path.isfile( os.path.join(self.out,type_dic[type][0],self.h5_name) ):
			print('File %s found. Please remove manually.' % self.h5_name)
			return
		else:
			if type=='pr':
				self.h5_pr_file = tables.openFile(os.path.join(self.out,type_dic[type][0],self.h5_name), 'w',
												  title=self.file_name)
				self.h5_pr_file.close() #need to close file before join data
			elif type=='ind':
				self.h5_ind_file = tables.openFile(os.path.join(self.out,type_dic[type][0],self.h5_name), 'w',
												   title=self.file_name)
			elif type=='gen':
				self.h5_gen_file = tables.openFile(os.path.join(self.out,type_dic[type][0],str(self.gen_iter)+'_'+self.h5_name),
												   'w', title=self.file_name)
				self.gen_iter+=1

	def close(self):
		self.h5_gen_file.close()
		self.h5_ind_file.close()
		self.h5_pr_file.close()

	def summary(self):
		raise (NotImplementedError)


class GenotypePLINK(GenotypeHDF5):

	def __init__(self, name, reader=None):
		super(GenotypePLINK, self).__init__(name)
		self.reader=reader
		self.split_size=None


	def convert_individuals(self):

		individuals=self.reader.folder.get_fam()

		self.h5_ind_file.createTable(self.h5_ind_file.root, 'individuals', individuals,
									 title='Individuals', filters=self.pytable_filters)
		self.h5_ind_file.root.individuals[:] = individuals
		self.individuals = self.h5_ind_file.root.individuals[:]
		self.n_ind=len(individuals)


	#@profile
	def convert_probes(self, chunk_size=100000):

		os.remove(os.path.join(self.out,'probes',self.h5_name))

		i=0
		chunk=np.array([])
		while True:
			chunk=self.reader.folder.get_bim(chunk_size)
			if isinstance(chunk,type(None)):
				break
			print i
			# WARNING!!! doesn't work on windows
			chunk.to_hdf(os.path.join(self.out,'probes',self.h5_name), key='probes',format='table', append=True,
						 names=['CHR', 'ID', 'distance', 'bp', 'allele1', 'allele2'],
						 complib='zlib',complevel=9, min_itemsize = 25)
			gc.collect()
			i+=1
		print('Number of Probes {} converted'.format(self.reader.folder.N_probes))

	#@profile
	def convert_genotypes(self):

		chunk_size=self.split_size
		if chunk_size is None:
			raise ValueError('CONVERTER_SPLIT_SIZE does not define in config file!')
		G=np.array([])
		#self.reader.folder.processed=0
		while True:
			with Timer() as t:
				G=self.reader.folder.get_bed(chunk_size)
				if isinstance(G,type(None)):
					break

			print ('Time to read {} SNPs is {} s'.format(G.shape[0], t.secs))

			self.write_data('gen')
			atom = tables.Int8Atom()
			self.genotype = self.h5_gen_file.createCArray(self.h5_gen_file.root, 'genotype', atom,
														  (G.shape),
														  title='Genotype', filters=self.pytable_filters)
			with Timer() as t:
				self.genotype[:] = G

			print ('Time to write {} SNPs is {} s'.format(G.shape[0], t.secs))

			self.h5_gen_file.close()
			G=None
			gc.collect()



	def plink2hdf5(self,out, force=True):

		if force:
			try:
				os.mkdir(os.path.join(out,'genotype') )
				os.mkdir(os.path.join(out,'individuals') )
				os.mkdir(os.path.join(out,'probes') )
			except:
				print('Directories "genotype","probes","individuals" are already exist in {}...'.format(self.out))

		self.out=out
		self.write_data('ind')
		self.convert_individuals()
		self.h5_ind_file.close()

		self.reader.folder.processed=0
		self.write_data('pr')
		self.convert_probes()

		self.reader.folder.processed=0
		self.convert_genotypes()


	def _summary(self, head=10):
		pass #TODO (low) rewrite to get statistic not load to memory

		print('Number of Probes: %d' % None)
		print('Number of Individuals: %d' % None)
		print('The genotype matrix is of size %d by %d' % None)

class GenotypeMINIMAC(object):

	def __init__(self, name, reader=None):
		self.reader=reader
		self.study_name=name
		self.split_size=None
		self.hdf5_iter=0
		self.pytable_filter=tables.Filters(complevel=9, complib='zlib')
		self.cluster=False

	def save_hdf5_chunk(self,data,out,name):
		print 'Saving chunk...{}'.format(os.path.join(out,'genotype',str(self.hdf5_iter)+'_'+name+'.h5'))
		h5_gen_file = tables.openFile(
			os.path.join(out,'genotype',str(self.hdf5_iter)+'_'+name+'.h5'), 'w', title=name)

		atom = tables.Float16Atom()  # TODO (low) check data format
		genotype = h5_gen_file.createCArray(h5_gen_file.root, 'genotype', atom,
											(data.shape),
											title='Genotype',
											filters=self.pytable_filter)
		genotype[:] = data
		h5_gen_file.close()
		genotype=None
		data=None
		gc.collect()
		self.hdf5_iter+=1

	def MACH2hdf5(self, out, id=False,remove_id=False):
		FNULL = open(os.devnull, 'w')
		subprocess.call(['bash',os.path.join(os.environ['HASEDIR'],'tools','minimac2hdf5.sh'),
						 self.reader.folder.path, out , os.environ['HASEDIR'], self.study_name ], shell=False,stderr=FNULL)

		if id:
			if self.cluster:
				ind=pd.read_hdf(os.path.join(out,'individuals',self.study_name+'.h5'),'individuals').individual
				N=ind.shape[0]
				print 'Submit to cluster!'
				cmd="qsub -sync y -t 1-{} {} {}".format(N,os.path.join(os.environ['HASEDIR'],'tools','qsub_helper.sh'),os.path.join( out,'id_convert.sh' ))
				print cmd
				proc=subprocess.Popen(cmd, shell=True,stderr=FNULL,stdout=subprocess.PIPE).communicate()
			else:
				proc=subprocess.Popen(['bash',os.path.join( out,'id_convert.sh' ) ], shell=False,stderr=FNULL)
				print proc.communicate()
			self.folder=MINIMACHDF5Folder(out,self.study_name)
			self.folder.pool.split_size=self.split_size
			self.folder.pool.chunk_size=self.split_size
			print('Start to convert id files to chunk files...')
			while True:
				data=self.folder.get_next()
				if data is None:
					break
				self.save_hdf5_chunk(data,out,self.study_name)
				gc.collect()

			print('Finished')
			if remove_id:
				self.folder.pool.remove(type='all')
			else:
				try:
					os.mkdir(os.path.join(out,'id_genotype') )
				except:
					print 'id_genotype folder already exist'

				self.folder.pool.move(os.path.join(out,'id_genotype'),type='all')

			os.remove(os.path.join(out,'id_convert.sh'))
			shutil.move(os.path.join(out,'SUB_FAM.txt'), os.path.join(out,'id_genotype','SUB_FAM.txt') )
			shutil.move(os.path.join(out,'SUB_ID.txt'), os.path.join(out,'id_genotype','SUB_ID.txt') )

		else:
			f=open(os.path.join( out,'minimac_convert.sh' ), 'w')
			N_probes=pd.read_hdf(os.path.join(out,'probes', self.study_name +'.h5'),'probes', where='columns=[ID]').ID.shape[0]
			chunk=np.vstack(((np.arange(0,N_probes,self.split_size)+1)[:-1],np.arange(0,N_probes,self.split_size)[1:]))
			for i_ch in range(chunk.shape[1]):
				ch=chunk[:,i_ch]
				print ch
				l='bash {} {} {} {} {} {} {} {} \n'.format(
					os.path.join(os.environ['HASEDIR'],'tools','minimacGenotype2hdf5.sh'),
					self.reader.folder.path,
					out,
					os.environ['HASEDIR'],
					self.study_name,
					ch[0],
					ch[1],
					i_ch
						)
				f.write(l)
			if ch[1]!=N_probes:
				l='bash {} {} {} {} {} {} {} {} \n'.format(
					os.path.join(os.environ['HASEDIR'],'tools','minimacGenotype2hdf5.sh'),
					self.reader.folder.path,
					out,
					os.environ['HASEDIR'],
					self.study_name,
					ch[1]+1,
					N_probes,
					i_ch+1
						)
				f.write(l)
			f.close()
			if self.cluster:
				print 'Submit to cluster!'
				cmd="qsub -sync y -t 1-{} {} {}".format(len(chunk),os.path.join(os.environ['HASEDIR'],'tools','qsub_helper.sh'),os.path.join( out,'minimac_convert.sh' ))
				print cmd
				proc=subprocess.Popen(cmd, shell=True,stderr=FNULL,stdout=subprocess.PIPE).communicate()
			else:
				proc=subprocess.Popen(['bash',os.path.join( out,'minimac_convert.sh' ) ], shell=False,stderr=FNULL)
				print proc.communicate()

			os.remove(os.path.join(out,'minimac_convert.sh'))

		os.remove(os.path.join(out,'files_order.txt'))


	def summary(self):
		pass

