
import time
import numpy as np
import pandas as pd
import tables
import h5py
import bitarray as ba
import gc
from collections import OrderedDict
import os
from scipy import stats
from hash import *
import glob

def timer(func):
	def f(*args,**kwargs):
		with Timer(verbose=True) as t:
			return func(*args,**kwargs)
	return f


class Timer(object):
	def __init__(self, verbose=False):
		self.verbose = verbose

	def __enter__(self):
		self.start = time.time()
		return self

	def __exit__(self, *args):
		self.end = time.time()
		self.secs = self.end - self.start
		self.msecs = self.secs * 1000  # millisecs
		if self.verbose:
			print 'elapsed time: %f ms' % self.msecs


########################################################################
########################################################################
########################################################################

class Analyser(object):

	def __init__(self,name):
		self.name=name


class HaseAnalyser(Analyser):

	def __init__(self):
		self.name=None
		self.t_stat=None
		self.p_value=None
		self.betas=None
		self.SE=None
		self.DF=None
		self.n_studies=None
		self.probes_path={}
		self.results={}
		self.split_size=None
		self.threshold=None
		self.result_path=None
		self.result_index=1
		self.out=None
		self.MAF=None
		self.rsid=None
		self.reference=None
		self.cluster=False
		self.permutation=False
		self.rsid_dic={}
		self.result_folder=None
		self.result_dump_size=3


	def summary(self):

		if isinstance(self.result_path,type(None)):
			raise ValueError('Please set result pathway to Analyser!')

		if isinstance(self.DF,type(None)):
			raise ValueError('Please set DF to Analyser!')

		if self.result_folder is None:
			self.result_folder=glob.glob(self.result_path + '*.npy')

		self.results['RSID']=np.array([])
		self.results['p_value']=np.array([])
		self.results['t-stat']=np.array([])
		self.results['phenotype']=np.array([])
		self.results['SE']=np.array([])
		self.results['MAF']=np.array([])

		files=[]
		for i in range(self.result_dump_size):
			try:
				files.append(self.result_folder.pop())
			except:
				break
		if len(files)!=0:
			for i in files:
				print i
				d=np.load(i).item()
				p_value=stats.t.sf(np.abs(d['t-stat']),self.DF)*2
				self.results['t-stat']=np.append(self.results['t-stat'],d['t-stat'].flatten())
				self.results['SE']=np.append(self.results['SE'],d['SE'].flatten())
				self.results['RSID']=np.append(self.results['RSID'],d['index'])
				self.results['phenotype']=np.append(self.results['phenotype'],d['phenotype'])
				self.results['MAF']=np.append(self.results['MAF'],d['MAF'])
				self.results['p_value']=np.append(self.results['p_value'],p_value.flatten())
		else:
			self.results=None

	def save_result(self , phen_names):

		t_threshold=self.threshold
		save_path=self.out

		mask=([],[],[])

		if t_threshold:
			mask = np.where(np.abs(self.t_stat) > t_threshold)
		else:
			mask = np.where(np.abs(self.t_stat) > 0)

		if (len(mask[0]) != 0):
			print ('Saving results to {}'.format(save_path))
			t_save = self.t_stat[mask[0],mask[1],mask[2]]
			se=self.SE[mask[0],mask[1],mask[2]]
			result = {'phenotype': phen_names[mask[2]], 't-stat': t_save,'index':self.rsid[mask[0]],'SE':se, 'MAF':self.MAF[mask[0]]}

			if not self.cluster:
				np.save(os.path.join(save_path, str(self.result_index) + 'result.npy'), result)
			else:
				if self.permutation:
					result=self.t_stat
				self.rsid_dic[str(self.chunk[0])+'_' +str(self.chunk[1])]=self.rsid
				np.save(os.path.join(save_path,'node'+str(self.node)+'_'+ str(self.chunk[0])+'_' +str(self.chunk[1])+'_'+ str(self.result_index) + 'result.npy'),result)
		self.result_index+=1



########################################################################
########################################################################
########################################################################
class Log (object): #TODO (low) write this class

	def __init__(self):
		pass


########################################################################
########################################################################
########################################################################
class Checker(object): #TODO (mid) finish or remove this class

	def __init__(self):
		self.mode=None
		self.name=None
		self.args=None

	def check(self,args, mode=None):

		if mode=='converting':
			return self.converting(args)
		elif mode=='encoding':
			self.encoding(args)
		elif mode=='single-meta':
			self.single_meta(args)
		elif mode=='meta-stage':
			self.meta_stage(args)
		elif mode=='regression':
			self.regression(args)


	def system_check(self,args):
		pass

	def converting(self,args):


		if len(args.study_name)>1:
			raise ValueError('Should be only one study name!')

		if len(args.genotype)>1:
			raise ValueError('Should be only one directory for genotype!')

		if isinstance(args.study_name, type(None)):
			raise ValueError('study name is not defined')

		args.genotype=args.genotype[0]
		if not os.path.isdir(args.genotype):
			raise ValueError('{} is not a directory!'.format(args.genotype))

	def encoding(self,args):

		if len(args.study_name)>1:
			raise ValueError('Should be only one study name!')

		if len(args.genotype)>1:
			raise ValueError('Should be only one directory for genotype!')

		if isinstance(args.study_name, type(None)):
			raise ValueError('study name is not defined')

		args.genotype=args.genotype[0]
		args.study_name=args.study_name[0]

		if not os.path.isdir(args.genotype) or not os.path.isdir(args.phenotype) or not os.path.isdir(args.out):
			raise ValueError('{},{} or {} is not a directory'.format(args.genotype, args.phenotype, args.out))

		if  not os.path.isdir(os.path.join(args.out, 'encode_genotype' )) or not os.path.isdir(os.path.join(args.out, 'encode_phenotype' )):
			raise ValueError('{} or {} not exist it is up to you to create these directory'.format(
				(os.path.join(args.out, 'encode_genotype' )),
				(os.path.join(args.out, 'encode_phenotype' ))
			)
			)


		id=np.array(pd.read_hdf(os.path.join(args.genotype,'individuals',args.study_name+'.h5'),'individuals').individual.tolist(), dtype=np.int)



	def single_meta(self,args):


		if len(args.study_name)>1:
			raise ValueError('Should be only one study name!')

		if len(args.genotype)>1:
			raise ValueError('Should be only one directory for genotype!')

		if isinstance(args.study_name, type(None)):
			raise ValueError('study name is not defined')

		args.genotype=args.genotype[0]
		args.study_name=args.study_name[0]

		if not os.path.isdir(args.genotype) or not os.path.isdir(args.phenotype) or not os.path.isdir(args.covariates) :
			raise ValueError('{},{} or {} is nor a directory'.format(os.path.isdir(args.genotype),os.path.isdir(args.phenotype), os.path.isdir(args.covariates )) )



	def meta_stage(self,args):

		if isinstance(args.mapper, type(None)):
			raise ValueError('You should define mapper folder --mapper, if it does not exist you should first run Mapper script')

		if isinstance(args.mapper_name, type(None)):
			raise ValueError('You should define mapper name --mapper_name, if it does not exist you should first run Mapper script')

		if len(args.genotype)<=1:
			raise ValueError('Should be only one directory for genotype!')

		g=[os.path.isdir(i) for i in args.genotype ]
		d=[os.path.isdir(i) for i in args.derivatives ]

		if np.sum(g)!=len(args.genotype) or not os.path.isdir(args.phenotype) or  np.sum(d)!=len(args.derivatives) :
			raise ValueError('{},{} or {} is nor a directory'.format(os.path.isdir(args.genotype),os.path.isdir(args.phenotype), os.path.isdir(args.derivatives )) )


	def regression(self,args):

		if isinstance(args.mapper, type(None)):
			raise ValueError('You should define mapper folder --mapper, if it does not exist you should first run Mapper script')

		if isinstance(args.mapper_name, type(None)):
			raise ValueError('You should define mapper name --mapper_name, if it does not exist you should first run Mapper script')

		g=[os.path.isdir(i) for i in args.genotype ]

		if np.sum(g)!=len(args.genotype) or not os.path.isdir(args.phenotype) or not os.path.isdir(args.covariates) :
			raise ValueError('{},{} or {} is nor a directory'.format(os.path.isdir(args.genotype),os.path.isdir(args.phenotype), os.path.isdir(args.derivatives )) )


########################################################################
########################################################################
########################################################################
class Mapper(object):

	def __init__(self):

		self.name=None
		self.genotype_names=[]
		self.dic=OrderedDict()
		self.n_study=0
		self.values=[]
		self.keys=None
		self.n_keys=None
		self.processed=0
		self._limit=None
		self.reference=False
		self.include=None
		self.exclude=None
		self.include_ind=np.array([])
		self.exclude_ind=np.array([])
		#self.hash=HashTable()
		self.cluster=None
		self.node=None
		self.chunk_pool=None
		self.chunk_size=None
		self.column_names=[]
		self.flip=[]
		self.probes=None


	def chunk_pop(self):
		if self.chunk_pool is None and self.cluster is None:
			raise ValueError('Cluster settings are not defined!')

		if self.chunk_pool is None:
			self.chunk_pool=[]
			if self.node is None:
				raise ValueError('cluster setting not defined! (number of nodes / this node number) ')
			self.node=[int(i) for i in self.node]
			N=self.n_keys/self.node[0]
			n=N/self.chunk_size
			if N<1:
				raise ValueError('Too many nodes! Change chunk size in mapper')
			elif N==1:
				self.chunk_pool=[[self.node[1]-1,self.node[1]]]
			else:
				if n==0:
					self.chunk_pool=[[N*(self.node[1]-1),N*(self.node[1]-1)+N]]
				else:
					self.chunk_pool=[]
					for i in range(n):
						self.chunk_pool.append(   [N*(self.node[1]-1)+ self.chunk_size*i, N*(self.node[1]-1) + self.chunk_size*(i+1)   ] )
					if (N-n*self.chunk_size)!=0:
						self.chunk_pool.append( [N*(self.node[1]-1) + self.chunk_size*(n), N*(self.node[1]-1) + self.chunk_size*(n)+ (N-n*self.chunk_size)  ]  )


			if self.node[1]==self.node[0]:
				if (self.n_keys-N*self.node[0])!=0:
					self.chunk_pool.append([self.chunk_pool[-1][1],self.n_keys])
			self.chunk_pool=self.chunk_pool[::-1]

		if len(self.chunk_pool)!=0:
			return self.chunk_pool.pop()
		else:
			return None


	def fill(self, keys, name ,repeats=False, reference=False):#TODO (middle) remove
		self.reference=reference
		self.reference_name=name
		self.column_names.append(name)
		l=len(keys)
		uniq=np.unique(keys)
		u=len(uniq)
		if not repeats:
			if l!=u:
				raise ValueError('length of keys {} does not =  {} number of uniq values!!!'.format(l,u))
		else:
			l=u
			keys=uniq
		self.n_keys=l
		if isinstance(self._limit, type(None)):
			for i,j in enumerate(keys):
				self.dic[j]=[i]
		else:
			for i,j in enumerate(keys):
				if i==self._limit:
					break
				self.dic[j]=[i]

	def push(self, new_keys,name=None, new_id=True):#TODO (middle) remove
		if not self.reference and len(self.dic)==0:
			raise ValueError('You should fill mapper first with ref panel or your own rsids!')
		self.n_study+=1
		if name is not None:
			try:
				self.genotype_names.index(name)
				raise ValueError('Trying to push the same study to mapper')
			except:
				self.genotype_names.append(name)
				self.column_names.append(name)

		if isinstance(self._limit, type(None)):
			for i,j in enumerate(new_keys):
				if self.dic.get(j):
					#self.dic[j][n]=i
					self.dic[j]=self.dic[j] + [i]
				else:
					if new_id:
						print ('WARNING! You are pushing ID {} which is not present in reference panel!'.format(j))
						self.dic[j]=[-1]*(self.n_study+1)
						self.dic[j][self.n_study]=i
					else:
						print ('WARNING! You are pushing ID {} which is not present in reference panel!'.format(j))
		else:
			for i,j in enumerate(new_keys):
				if i==self._limit:
					break
				if self.dic.get(j):
					self.dic[j]=self.dic[j] + [i]
				else:
					self.dic[j]=[-1]*(self.n_study+1)
					self.dic[j][self.n_study]=i

	def load_flip(self,folder):

		if folder is None:
			raise ValueError('Mapper is not defined!')
		if not os.path.isdir(folder):
			raise ValueError('{} is not a folder'.format(folder))

		if self.reference_name is None:
			raise ValueError('Reference name for mapper is not defined!')

		if len(glob.glob(folder + 'flip_*'))==0:
			raise ValueError('There is no flip mapper data in folder {}'.format(folder))

		for j,i in enumerate(self.genotype_names):
			self.flip.append(np.load(os.path.join(folder, 'flip_'+self.reference_name+'_'+i+'.npy')))
			if j==0:
				len_flip=self.flip[0].shape[0]
			else:
				if len_flip!=self.flip[j].shape[0]:
					raise ValueError('Different length of flip array between studies; used different ref panel!')
		self.flip=np.array(self.flip).T

	def load (self, folder):
		if folder is None:
			raise ValueError('Mapper is not defined!')
		if not os.path.isdir(folder):
			raise ValueError('{} is not a folder'.format(folder))

		if self.reference_name is None:
			raise ValueError('Reference name for mapper is not defined!')

		if len(glob.glob(folder + 'keys_*'))==0:
			raise ValueError('There is no mapper data in folder {}'.format(folder))

		keys=glob.glob(folder+ 'keys_*')
		if len(keys)>1:
			raise ValueError('There are more than one reference keys in folder {}'.format(folder))

		self.keys=np.load(os.path.join(keys[0]) )#TODO (middle) not safety to load only one file
		self.n_keys=self.keys.shape[0]

		values=glob.glob(folder+ 'values_*')
		if len(values)==0:
			raise ValueError('There is no mapper data in folder {}'.format(folder))

		if self.genotype_names is None:
			raise ValueError('Genotype names are not defined!')

		if len(values)!=len(self.genotype_names):
			raise ValueError('Number of files in mapper folder {} is not equal to number of studies {}'.format(len(values),len(self.genotype_names)))

		for j,i in enumerate(self.genotype_names):
			self.values.append(np.load(os.path.join(folder, 'values_'+self.reference_name+'_'+i+'.npy')))
			if j==0:
				len_value=self.values[0].shape[0]
			else:
				if len_value!=self.values[j].shape[0]:
					raise ValueError('Different length of values array between studies; used different ref panel!')

		self.values=np.array(self.values).T

		if self.n_keys<self.n_study:
			print ('WARNING!!! Normally should be more test values that studies! Check your saved data for Mapper')


		if self.n_keys!=self.values.shape[0]:
			raise ValueError('Number of indexes {} in mapper values is different from keys length {}!'.format(self.values.shape[0],self.n_keys ))

		self.n_study=self.values.shape[1]

		print ('You loaded values for {} studies and {} test values'.format(self.n_study, self.n_keys))

		if self.include is not None or self.exclude is not None:
			self.reference=Reference()
			self.reference.name=self.reference_name
			self.reference.load_index()

	def get(self,chunk_number=None):

		if isinstance(self.keys, type(None)) and isinstance(self.values, type(None)):
			raise ValueError('mapper data is not loaded!')

		if isinstance(self.chunk_size,type(None)):
			raise ValueError('chunk_size and should be defined')

		if self.processed==self.n_keys:
			return None, None

		if self.include is not None:
			if len(self.include_ind)==0:
				if 'ID' in self.include.columns:
					self.query_include=self.reference.index.select('reference',where='ID=self.include.ID')
					self.include_ind=self.query_include.index
				elif 'CHR' in self.include.columns and 'bp' in self.include.columns:
					self.query_include=self.reference.index.select('reference',where='CHR=self.include.CHR & bp=self.include.bp')
					self.include_ind=self.query_include.index

		if self.exclude is not None:
			if len(self.exclude_ind)==0:
					if 'ID' in self.exclude.columns:
						self.query_exclude=self.reference.index.select('reference',where='ID=self.exclude.ID')
						self.exclude_ind=self.query_exclude.index
					elif 'CHR' in self.exclude.columns and 'bp' in self.exclude.columns:
						self.query_exclude=self.reference.index.select('reference',where='CHR=self.exclude.CHR & bp=self.exclude.bp')
						self.exclude_ind=self.query_exclude.index

		if chunk_number is not None:
			start=chunk_number[0]
			finish=chunk_number[1]
		else:
			start=self.processed
			finish=self.processed+self.chunk_size if (self.processed+self.chunk_size)<self.n_keys else self.n_keys
			self.processed=finish
		self.include_ind=np.setxor1d(self.include_ind,self.exclude_ind)
		if len(self.include_ind)>0:
			ind=np.intersect1d(np.arange(start,finish),self.include_ind)
		else:
			ind=np.arange(start,finish)

		if chunk_number is not None and len(ind)==0:
			return None,None
		else:
			while len(ind)==0:
				start=self.processed
				finish=self.processed+self.chunk_size if (self.processed+self.chunk_size)<self.n_keys else self.n_keys
				self.processed=finish
				ind=np.intersect1d(np.arange(start,finish),self.include_ind)
				if self.processed==self.n_keys:
					return None, None

		ind=ind.astype('int')
		indexes=self.values[ind,:]
		keys=self.keys[ind]
		r=(indexes==-1).any(axis=1)
		indexes=indexes[~r]
		keys=keys[~r]
		return [indexes[:,i] for i in range(self.n_study)], keys



	def get_all(self, name, nonempty=True):
		ind=self.genotype_names.index(name)
		indexes=self.values[:,ind]
		if nonempty:
			r=(self.values==-1).any(axis=1)
			if len(r)==self.values.shape[0]:
				raise ValueError('There is no common names between studies')
			indexes=indexes[~r]
		return indexes


########################################################################
########################################################################
########################################################################

class Reference(object):

	def __init__(self):
		self.name=None
		self.path_default=os.path.join(os.environ['HASEDIR'], 'data')
		self.path={
			'1000Gp1v3_ref':{'table':os.path.join(os.environ['HASEDIR'], 'data', '1000Gp1v3.ref.gz'),
							 'index':os.path.join(os.environ['HASEDIR'], 'data', '1000Gp1v3.ref_info.h5') }
				  }
		self.dataframe=None
		self.loaded=False
		self.chunk=10000
		self.read=0
		self.columns=['ID','allele1','allele2','CHR','bp']
		self.index=None



	def load(self):
		if self.name is not None:
			if self.path.get(self.name) is not None:
				try:
					self.dataframe=pd.read_csv(self.path[self.name]['table'], compression='gzip', sep=' ',chunksize=self.chunk)
				except:
					self.dataframe=pd.read_csv(self.path[self.name]['table'], sep=' ',chunksize=self.chunk)
				self.loaded=True
			else:
				if os.path.isfile(os.path.join(self.path_default,self.name)):
					try:
						self.dataframe=pd.read_csv(os.path.join(self.path_default,self.name), compression='gzip', sep=' ',chunksize=self.chunk)
					except:
						self.dataframe=pd.read_csv(os.path.join(self.path_default,self.name), sep=' ',chunksize=self.chunk)
					self.loaded=True
					self.name=os.path.basename(self.name)
				else:
					raise ValueError('Unknown reference {}!'.format((self.name)))
		else:
			raise ValueError('Reference name is not define!')


	def load_index(self):
		try:
			self.index=pd.HDFStore( os.path.join( self.path[self.name]['index'] ), 'r' )
		except:
			if os.path.isfile(os.path.join(self.path_default,self.name)):
				self.index=pd.HDFStore(os.path.join(self.path_default,self.name), 'r' )
			else:
				raise ValueError('There is {} no index file {}'.format(self.path_default, self.name))


	def next(self):
		df=self.dataframe.get_chunk()
		self.read+=df.shape[0]
		return df














########################################################################
########################################################################
########################################################################

def study_indexes( args=None, genotype=None,phenotype=None,covariates=None):

	def _get_id(notype):
		if isinstance(notype,type(None)):
			id_p=np.array([])
		else:
			if isinstance(notype,tuple):
				id_p=np.array([])
				for i in notype:
					if isinstance(i,dict):
						id_p=np.append(id_p,i['id'])
					else:
						id_p=np.append(id_p,i.get_id())
			else:
				id_p=notype.get_id()
		return id_p.astype(np.str)

	id_c=_get_id(covariates)
	id_p=_get_id(phenotype)
	id_g=_get_id(genotype)

	index_g=np.array([])
	index_p=np.array([])
	index_c=np.array([])

	if args is not None:
		if not isinstance(args.ind_id_inc,type(None)) or not isinstance(args.ind_id_exc,type(None)):
			pass
		if not isinstance(args.cov_name_inc,type(None)) or not isinstance(args.cov_name_exc,type(None)):
			pass

	if len(id_g)!=0 and len(id_p)!=0 and len(id_c)!=0:
		common_id=(np.intersect1d(np.intersect1d(id_g,id_p),id_c))

		index_g=np.array([np.where(id_g==i)[0][0] for i in common_id])
		index_p=np.array([np.where(id_p==i)[0][0] for i in common_id])
		index_c=np.array([np.where(id_c==i)[0][0] for i in common_id])

	elif len(id_g)==0:
		common_id=(np.intersect1d(id_p,id_c))
		index_p=np.array([np.where(id_p==i)[0][0] for i in common_id])
		index_c=np.array([np.where(id_c==i)[0][0] for i in common_id])

	elif len(id_c)==0:
		common_id=(np.intersect1d(id_p,id_g))
		index_g=np.array([np.where(id_g==i)[0][0] for i in common_id])
		index_p=np.array([np.where(id_p==i)[0][0] for i in common_id])

	else:
		common_id=(np.intersect1d(id_g,id_c))
		index_g=np.array([np.where(id_g==i)[0][0] for i in common_id])
		index_c=np.array([np.where(id_c==i)[0][0] for i in common_id])

	print ('There are {} common ids'.format(len(common_id)))
	np.savetxt(os.path.join(os.environ['HASEOUT'],'study_common_id.txt'),common_id, fmt='%s')
	if len(common_id)==0:
		exit(0)

	return [index_g,index_p,index_c], np.array(common_id)


def merge_pard(pard, SNPs_index,B4=False):

	if len(pard)!=len(SNPs_index):
		raise ValueError('There are not equal number od PD and SNPs indexes {}!={}'.format(len(pard), len(SNPs_index)))

	a_test, b_cov, C, a_cov = pard[0].get(index=SNPs_index[0])

	if B4:
		b4=pard[0].folder.data.b4

	for i in range(1, len(pard)):

		a,b,c,a_c=pard[i].get(index=SNPs_index[i])

		a_test=a_test + a
		b_cov=b_cov + b
		C=C+c
		a_cov=a_cov + a_c
		if B4:
			b4=b4+pard[i].folder.data.b4

	if B4:
		return a_test, b_cov, C, a_cov,b4
	else:
		return a_test, b_cov, C, a_cov


def maf_pard(pard,SNPs_index): #TODO (middle) delete function

	samples=0
	maf=np.zeros( (SNPs_index[0]) )
	for j,i in enumerate(pard):
		n=len(pard.folder._data.metadata['id'])
		samples+=n
		maf=maf+n*pard.folder._data.metadata['maf'][SNPs_index[j]]

	maf=maf/np.float(samples)
	return maf



def merge_genotype(genotype, SNPs_index , mapper, flip_flag=True):

	if SNPs_index is None:
		gen=genotype[0].get_next()
		if gen is not None:
			for i in range(1, len(genotype)):
				gen=np.hstack( (gen,genotype[i].get_next() ) )
		return gen
	else:
		if len(genotype)!=len(SNPs_index):
			raise ValueError('There are not equal number of genotypes and SNPs indexes {}!={}'.format(len(genotype), len(SNPs_index)))
		gen=genotype[0].get(SNPs_index[0])
		if flip_flag:
			flip=mapper.flip[SNPs_index[0],0]
			flip_index=(flip==-1)
			gen=np.apply_along_axis(lambda x: flip*(x-2*flip_index) ,0,gen)
		for i in range(1, len(genotype)):
			g=genotype[i].get(SNPs_index[i])
			if flip_flag:
				flip=mapper.flip[SNPs_index[i],i]
				flip_index=(flip==-1)
				g=np.apply_along_axis(lambda x: flip*(x-2*flip_index) ,0,g)
			gen=np.hstack( (gen,g ) )
		return gen


if __name__=='__main__':
	print 'tools'


