
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
		self.result_dump_size=100


	def summary(self):

		if isinstance(self.result_path,type(None)):
			raise ValueError('Please set result pathway to Analyser!')

		if isinstance(self.DF,type(None)):
			raise ValueError('Please set DF to Analyser!')

		if self.result_folder is None:
			self.result_folder=os.listdir(self.result_path)

		self.results['RSID']=np.array([])
		self.results['p_value']=np.array([])
		self.results['t-stat']=np.array([])
		self.results['phenotype']=np.array([])
		self.results['SE']=np.array([])

		files=[]
		for i in range(self.result_dump_size):
			try:
				files.append(self.result_folder.pop())
			except:
				break
		if len(files)!=0:
			for i in files:
				print i
				d=np.load(os.path.join(self.result_path, i)).item()
				p_value=stats.t.sf(np.abs(d['t-stat']),self.DF)*2
				self.results['t-stat']=np.append(self.results['t-stat'],d['t-stat'].flatten())
				self.results['SE']=np.append(self.results['SE'],d['SE'].flatten())
				self.results['RSID']=np.append(self.results['RSID'],d['index'])
				self.results['phenotype']=np.append(self.results['phenotype'],d['phenotype'])
				self.results['p_value']=np.append(self.results['p_value'],p_value.flatten())
		else:
			self.results=None

	def save_result(self ,phen_names):

		t_threshold=self.threshold
		save_path=self.out

		mask=([],[],[])

		if t_threshold:
			mask = np.where(np.abs(self.t_stat) > t_threshold)
		else:
			mask = np.where(np.abs(self.t_stat) > 0)

		if (len(mask[0]) != 0):
			t_save = self.t_stat[mask[0],mask[1],mask[2]]
			se=self.SE[mask[0],mask[1],mask[2]]
			result = {'phenotype': phen_names[mask[2]], 't-stat': t_save,'index':self.rsid[mask[0]],'SE':se}

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

	def __init__(self,name):

		self.name=name
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
		#self.hash=HashTable()
		self.cluster=None
		self.node=None
		self.chunk_pool=None
		self.chunk_size=None
		self.column_names=[]


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
				raise ValueError('Too many nodes! Change chunk size in mapper.')
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
					self.chunk_pool.append([self.n_keys-N*self.node[0],self.n_keys])
			self.chunk_pool=self.chunk_pool[::-1]

		if len(self.chunk_pool)!=0:
			return self.chunk_pool.pop()
		else:
			return None


	def fill(self, keys, name ,repeats=False, reference=False):
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
				#self.dic[j]=[i]+[-1]*self.n_study
				self.dic[j]=[i]
		else:
			for i,j in enumerate(keys):
				if i==self._limit:
					break
				#self.dic[j]=[i]+[-1]*self.n_study
				self.dic[j]=[i]

	def push(self, new_keys,name=None, new_id=True):
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

	def save(self,path):

		np.save(os.path.join(path, 'keys_'+self.reference_name+'.npy'), self.dic.keys())
		values=np.array(self.dic.values() )
		n=1
		for i in range(values.shape[1]):
			if (self.genotype_names)!=0:
				name='values_'+self.genotype_names[i-1] +'.npy'
			else:
				name='values_'+str(n)+'_'+self.name+'.npy'
			if i==0:
				np.save(os.path.join(path, 'ref_'+self.reference_name+'.npy'), values[:,i] )
			else:
				np.save(os.path.join(path, name), values[:,i] )
			n+=1

		#TODO (low) to do through Log class
		f=open(os.path.join(path,'mapper_log.txt'))
		f.write('##############################\n')
		f.write(self.reference.reference_name + '\n')
		for i in self.genotype_names:
			f.write(i + '\n')
		f.write('##############################\n')
		f.close()



	def load (self, folder):

		if len(glob.glob(folder + 'key_*'))==0:
			raise ValueError('There is no mapper data in folder {}'.format(folder))

		keys=glob.glob(folder+ 'key_*')
		if len(keys)>1:
			raise ValueError('There are more than one reference keys in folder {}'.format(folder))

		self.keys=np.load(os.path.join(keys[0]) )
		self.n_keys=self.keys.shape[0]

		values=glob.glob(folder+ 'values_*')
		if len(values)==0:
			raise ValueError('There is no mapper data in folder {}'.format(folder))

		if self.genotype_names is None:
			raise ValueError('Genotype names are not defined!')

		if len(values)!=len(self.genotype_names):
			raise ValueError('Number of files in mapper folder {} is not equal to number of studies {}'.format(len(values),len(self.genotype_names)))

		for i in self.genotype_names:
			self.values.append(np.load(os.path.join(folder, 'values_'+i+'.npy')))

		self.values=np.array(self.values).T

		if self.n_keys<self.n_study:
			print ('WARNING!!! Normally should be more test values that studies! Check your saved data for Mapper')


		if self.n_keys!=self.values.shape[0]:
			raise ValueError('Number of indexes {} in mapper values is different from keys length {}!'.format(self.values.shape[0],self.n_keys ))

		self.n_study=self.values.shape[1]

		print ('You loaded values for {} studies and {} test values'.format(self.n_study, self.n_keys))

		if self.include is not None or self.exclude is not None:
			self.hash=HashTable()
			self.hash.fill(self.keys)
			#self.hash.save('pathway', 'rsid_hash')#TODO (low) check
			print ('Hash Table ready!')


	def get_chunk(self,chunk_number ):
		if isinstance(self.keys, type(None)) and isinstance(self.values, type(None)):
			raise ValueError('mapper data is not loaded!')

		if isinstance(self.chunk_size,type(None)):
			raise ValueError('chunk_size and should be defined')

		if self.include is None and self.exclude is None:
			start=chunk_number[0]
			finish=chunk_number[1]
			indexes=self.values[start:finish,:]
			keys=self.keys[start:finish]
			r=(indexes==-1).any(axis=1)
			indexes=indexes[~r]
			keys=keys[~r]
			#return [link(indexes[:,i], split_size) for i in range(self.n_study)], keys
			return [indexes[:,i] for i in range(self.n_study)], keys

		elif self.include is not None: #TODO (mid) check for millions snps
			ind=np.array([self.hash.get_index(i) for i in self.include])
			if self.exclude is not None:
				ind_exc=np.array([self.hash.get_index(i) for i in self.exclude])
				ind=np.setxor1d(ind,ind_exc)
			if ind.shape[0]==0:
				raise ValueError('No rsid for test!')
			indexes=self.values[ind,:]
			keys=self.include
			r=(indexes==-1).any(axis=1)
			indexes=indexes[~r]
			keys=keys[~r]
			self.processed=self.n_keys
			#return [link(indexes[:,i], split_size) for i in range(self.n_study)], keys
			return [indexes[:,i] for i in range(self.n_study)], keys
		elif self.exclude is not None:
			ind_exc=np.array([self.hash.get_index(i) for i in self.include])
			start=chunk_number[0]
			finish=chunk_number[1]
			ind=np.setxor1d(np.arange(start,finish),ind_exc)
			indexes=self.values[ind,:]
			keys=self.keys[ind]
			r=(indexes==-1).any(axis=1)
			indexes=indexes[~r]
			keys=keys[~r]
			#return [link(indexes[:,i], split_size) for i in range(self.n_study)], keys
			return [indexes[:,i] for i in range(self.n_study)], keys



	def get_next(self):

		if isinstance(self.keys, type(None)) and isinstance(self.values, type(None)):
			raise ValueError('mapper data is not loaded!')


		if isinstance(self.chunk_size,type(None)):
			raise ValueError('chunk_size   should be defined')

		if self.processed==self.n_keys:
			return None, None
		else:
			if self.include is None and self.exclude is None:
				start=self.processed
				finish=self.processed+self.chunk_size if (self.processed+self.chunk_size)<self.n_keys else self.n_keys
				self.processed=finish
				indexes=self.values[start:finish,:]
				keys=self.keys[start:finish]
				r=(indexes==-1).any(axis=1)
				indexes=indexes[~r]
				keys=keys[~r]
				#return [link(indexes[:,i], split_size) for i in range(self.n_study)], keys
				return [indexes[:,i] for i in range(self.n_study)], keys
			elif self.include is not None:
				ind=np.array([self.hash.get_index(i) for i in self.include])
				if self.exclude is not None:
					ind_exc=np.array([self.hash.get_index(i) for i in self.exclude])
					ind=np.setxor1d(ind,ind_exc)
				if ind.shape[0]==0:
					raise ValueError('No rsid for test!')
				indexes=self.values[ind,:]
				keys=self.include
				r=(indexes==-1).any(axis=1)
				indexes=indexes[~r]
				keys=keys[~r]
				self.processed=self.n_keys
				#return [link(indexes[:,i], split_size) for i in range(self.n_study)], keys
				return [indexes[:,i] for i in range(self.n_study)], keys

			elif self.exclude is not None:
				ind_exc=np.array([self.hash.get_index(i) for i in self.include])
				start=self.processed
				finish=self.processed+self.chunk_size if (self.processed+self.chunk_size)<self.n_keys else self.n_keys
				self.processed=finish
				ind=np.setxor1d(np.arange(start,finish),ind_exc)
				indexes=self.values[ind,:]
				keys=self.keys[ind]
				r=(indexes==-1).any(axis=1)
				indexes=indexes[~r]
				keys=keys[~r]
				#return [link(indexes[:,i], split_size) for i in range(self.n_study)], keys
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

	def __init__(self, name):
		self.name=name
		self.reference_name=None
		self.path={
			'1000Gp1v3_20101123_ALL_impute_legend':os.path.join(os.environ['HASEDIR'], 'data', '1000Gp1v3.20101123.ALL.impute.legend.final.gz')
				  }
		self.dataframe=None
		self.loaded=False
		self.chunk=100000
		self.read=0


	def load(self):
		if self.reference_name is not None:
			if self.path.get(self.reference_name) is not None:
				try:
					self.dataframe=pd.read_csv(self.path[self.reference_name], compression='gzip', sep=' ')
				except:
					self.dataframe=pd.read_csv(self.path[self.reference_name], sep=' ')
				self.loaded=True
			else:
				raise ValueError('Unknown reference!')
		else:
			raise ValueError('Reference name is not define!')

	def open(self):
		if self.reference_name is not None:
			if self.path.get(self.reference_name) is not None:
				self.dataframe=pd.read_csv(self.path[self.reference_name], compression='gzip', sep=' ', iterator =True)
			else:
				raise ValueError('Unknown reference!')
		else:
			raise ValueError('Reference name is not define!')

	def next(self):
		df=self.dataframe.get_chunk(self.chunk)
		self.read+=df.shape[0]
		return df














########################################################################
########################################################################
########################################################################

def study_indexes( args=None, genotype=None,phenotype=None,covariates=None):

	if isinstance(genotype,type(None)):
		id_g=np.array([])
	else:
		if isinstance(genotype,tuple):
			id_g=np.array([])
			for i in genotype:
				id_g=np.append(id_g,i.get_id())
		else:
			id_g=genotype.get_id()

	if isinstance(phenotype,type(None)):
		id_p=np.array([])
	else:
		id_p=phenotype.get_id()

	if isinstance(covariates,type(None)):
		id_c=np.array([])
	else:
		if isinstance(covariates,tuple):
			id_c=np.array([])
			for i in covariates:
				id_c=np.append(id_c,i['id'])
		else:
			id_c=covariates.get_id()

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


def maf_pard(pard,SNPs_index):

	samples=0
	maf=np.zeros( (SNPs_index[0]) )
	for j,i in enumerate(pard):
		n=len(pard.folder._data.metadata['id'])
		samples+=n
		maf=maf+n*pard.folder._data.metadata['maf'][SNPs_index[j]]

	maf=maf/np.float(samples)
	return maf



def merge_genotype(genotype, SNPs_index):

	if SNPs_index is None:
		gen=genotype[0].get_next()
		if gen is not None:
			for i in range(1, len(genotype)):
				gen=np.hstack( (gen,genotype[i].get_next() ) )
		return gen
	else:
		if len(genotype)!=len(SNPs_index):
			raise ValueError('There are not equal number od genotypes and SNPs indexes {}!={}'.format(len(genotype), len(SNPs_index)))
		gen=genotype[0].get(SNPs_index[0])
		for i in range(1, len(genotype)):
			gen=np.hstack( (gen,genotype[i].get(SNPs_index[i]) ) )
		return gen


if __name__=='__main__':
	print 'tools'


