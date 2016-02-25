import h5py
import tables
from  hdgwas.tools import Timer, Checker, study_indexes, Mapper,HaseAnalyser, merge_pard, merge_genotype, Reference, maf_pard
from hdgwas.converter import  GenotypePLINK, GenotypeMINIMAC
from hdgwas.data import Reader, MetaParData
from hdgwas.fake import Encoder
from hdgwas.hdregression import HASE, A_covariates, A_tests, B_covariates, C_matrix, A_inverse,B4
import argparse
import numpy as np
import gc
import sys
import os
from hdgwas.pard import partial_derivatives
from hdgwas.regression import haseregression
import pandas as pd
import time
from hdgwas.protocol import Protocol
from config import MAPPER_CHUNK_SIZE, basedir,CONVERTER_SPLIT_SIZE
os.environ['HASEDIR']=basedir
__version__ = '1.0.0'

HEAD = "*********************************************************************\n"
HEAD += "* HASE: Framework for efficient high-dimensional association analyses \n"
HEAD += "* Version {V}\n".format(V=__version__)
HEAD += "* (C) 2015-2016 Gennady Roshchupkin and Hieab Adams\n"
HEAD += "* Erasmus MC, Rotterdam /  Department of Medical Informatics, Radiology and Epidemiology \n"
HEAD += "* GNU General Public License v3\n"
HEAD += "*********************************************************************\n"




if __name__=='__main__':

	start = time.time()

	parser = argparse.ArgumentParser(description='Script to use HASE in command line')

	#MAIN SETTINGS
	parser.add_argument("-thr",type=float,
						help="predefined threshold for t-value, if not defined save everything")

	parser.add_argument("-o", "--out", type=str, required=True, help="path to save result folder")


	parser.add_argument("-mode", required=True,type=str,
						choices=['regression','converting','single-meta', 'meta-stage','encoding'],
						help='Main argument, specify type of analysis'
							 '*****************'
							 'converting - before start to use any kind of analysis you need to convert data, it could take some time, '
							 'you need to do it ones and then everything will be super fast'
							 '*****************'
							 'single-meta - if you are the part of meta-analysis and will send data to main center, use this mode'
							 '*****************'
							 'encoding - if you are doing HD association, before send data to main center create fake data'
							 '*****************'
							 'meta-stage - if you are the main center of meta-analysis and already got precomputed data from all sites, use it'
							 '*****************'
							 'regression - if you want run everything yourself and do not use precomputed data use it'

							 'If you want to run single site full association analysis run (1) single-meta mode and then (2) meta-stage'

							 '*****************'
						)

	parser.add_argument("-g", "--genotype",nargs='+', type=str, help="path/paths to genotype data folder")
	parser.add_argument("-ph", "--phenotype", type=str, help="path to phenotype data folder")
	parser.add_argument("-cov", "--covariates", type=str, help="path to covariates data folder")


	parser.add_argument( "--derivatives", nargs='+', type=str, help="path to derivatives data folder")

	parser.add_argument("-protocol", type=str, help="path to study protocol file")

	parser.add_argument('-study_name', type=str, required=True,nargs='+', help=' Name for saved genotype data, without ext')

	parser.add_argument('-mapper', type=str, help='Mapper data folder')
	parser.add_argument('-mapper_name', type=str, help='Mapper name')

	#ADDITIONAL SETTINGS
	parser.add_argument("-snp_id_inc", type=str, help="path to file with SNPs id to include to analysis")
	parser.add_argument("-snp_id_exc", type=str, help="path to file with SNPs id to exclude from analysis")


	#parser.add_argument("-ind_id_inc", type=str, help="path to file with individuals id to include to analysis") #TODO (low)
	#parser.add_argument("-ind_id_exc", type=str, help="path to file with individuals id to exclude from analysis")#TODO (low)
	#parser.add_argument("-cov_name_inc", type=str, help="path to file with covariates names to include to analysis")#TODO (low)
	#parser.add_argument("-cov_name_exc", type=str, help="path to file with covariates names to exclude from analysis")#TODO (low)

	parser.add_argument('-intercept', type=str, default='y', choices=['y','n'], help='include intercept to regression, default yes')

	parser.add_argument('-maf', type=float, default=0.05, help='MAF for genetics data')
	###

	#FLAGS
	parser.add_argument('-hdf5', type=bool,default=True, help='flag for genotype data format')
	parser.add_argument('-pd_full', type=bool, default=False, help='For not HD association study')
	parser.add_argument('-permute_ph', type=bool, default=False, help='Flag for phenotype permutation')
	#TODO (low) save genotype after MAF
	###

	###CLUSTER SETTING
	parser.add_argument('-cluster', type=str, default='n', choices=['y','n'], help=' Is it parallel cluster job, default no')
	parser.add_argument('-node', nargs='+',help='number of nodes / this node number, example: 10 2 ')
	###

	#ADVANCED SETTINGS
	parser.add_argument('-interaction', type=str, help='path to file with data for genotype interaction test') #TODO (low)
	###
	args = parser.parse_args()
	if not args.thr:
		print ('WARNING!!! You did not set threshold for t-value, all results will be saved')
	ARG_CHECKER=Checker()
	print args


	################################### CONVERTING ##############################

	if args.mode=='converting':

		#ARG_CHECKER.check(args,mode='converting')

		R=Reader('genotype')
		R.start(args.genotype[0])

		with Timer() as t:
			if R.format=='PLINK':
				G=GenotypePLINK(args.study_name[0],reader=R)
				G.split_size=CONVERTER_SPLIT_SIZE
				G.plink2hdf5(out=args.out)

			elif R.format=='MINIMAC':
				G=GenotypeMINIMAC(args.study_name,reader=R)
				G.split_size=CONVERTER_SPLIT_SIZE
				G.MACH2hdf5(args.out)
			else:
				raise ValueError('Genotype data should be in PLINK or MINIMAC format and alone in folder')
		print ('Time to convert all data: {} sec'.format(t.secs))

	################################### ENCODING ##############################

	elif args.mode=='encoding':

		#ARG_CHECKER.check(args,mode='encoding')

		phen=Reader('phenotype')
		phen.start(args.phenotype)

		gen=Reader('genotype')
		gen.start(args.genotype, hdf5=args.hdf5, study_name=args.study_name, ID=False)

		e=Encoder()
		e.out=args.out
		e.study_name=args.study_name

		row_index, ids =  study_indexes(phenotype=phen.folder._data,genotype=gen.folder._data)
		with Timer() as t:

			e.matrix(len(ids),save=True)

			while True:
				with Timer() as t_gen:
					genotype = gen.get_next()
					genotype=genotype[:,row_index[0]]
					if isinstance(genotype, type(None)):
						break
					encode_genotype=e.encode(genotype, data_type='genotype')
					e.save_hdf5(encode_genotype,os.path.join(args.out, 'encode_genotype' ) )
					encode_genotype=None
					gc.collect()

				print ('Time to create fake genotype {} is {}sec'.format(genotype.shape,t_gen.secs))

			while True:
				with Timer() as t_phen:
					phenotype=phen.get_next(index=row_index[1])
					if isinstance(phenotype, type(None)):
						break
					encode_phenotype=e.encode(phenotype, data_type='phenotype')
					if phen.folder.format=='.npy':
						e.save_npy(encode_phenotype,save_path=os.path.join(args.out, 'encode_phenotype'),info=phen.folder)
					if phen.folder.format in ['.csv', '.txt']:
						e.save_csv(encode_phenotype,save_path=os.path.join(args.out, 'encode_phenotype'),info=phen.folder)
					encode_phenotype=None
					gc.collect()
			if phen.folder.format=='.npy':
				np.save(os.path.join(os.path.join(args.out, 'encode_phenotype', 'info_dic.npy')),e.phen_info_dic )
			print ('Time to create fake genotype {} is {}sec'.format(phenotype.shape,t_phen.secs))

		print ('Time to encode all data: {} sec'.format(t.secs))

	################################### SINGLE META STAGE ##############################

	elif args.mode=='single-meta':

		#ARG_CHECKER.check(args,mode='single-meta')

		phen=Reader('phenotype')
		phen.start(args.phenotype)

		cov=Reader('covariates')
		cov.start(args.covariates)

		if cov.folder.n_files>1:
			raise ValueError('In covariates folder should be only one file!')

		gen=Reader('genotype')
		gen.start(args.genotype, hdf5=args.hdf5, study_name=args.study_name, ID=False)

		with Timer() as t:
			partial_derivatives(save_path=args.out,COV=cov,PHEN=phen, GEN=gen,
								MAP=None, MAF=args.maf, R2=None, B4_flag=args.pd_full,
								study_name=args.study_name, intercept=args.intercept)
		print ('Time to compute partial derivatives : {} sec'.format(t.secs))

	################################### MULTI META STAGE ##############################

	elif args.mode=='meta-stage':

		#ARG_CHECKER.check(args,mode='meta-stage')

		##### Init data readers #####

		mapper=Mapper(args.mapper_name)
		mapper.chunk_size=MAPPER_CHUNK_SIZE
		mapper.genotype_names=args.genotype
		mapper.load(args.mapper)
		mapper.cluster=args.cluster
		mapper.node=args.node

		Analyser=HaseAnalyser()

		pard=[]


		for i,j in enumerate(args.derivatives):
			pard.append(Reader('partial') )
			pard[i].start(j, study_name=args.study_name[i])
			pard[i].folder.load()


		PD=[False if isinstance(i.folder.data.b4,None ) else True for i in pard]

		if np.sum(PD)!=len(pard) and np.sum(PD)!=0:
			raise ValueError('All study should have b4 data for partial derivatives!')

		if args.protocol is not None:
			protocol=Protocol(args.protocol)
		else:
			protocol=None

		meta_pard=MetaParData(pard,protocol=protocol)

		if np.sum(PD)==0:
			phen=Reader('phenotype')

			phen.start(args.phenotype)
			N_studies=len(args.genotype)

			gen=[]
			for i,j in enumerate(args.genotype):
				gen.append(Reader('genotype'))
				gen[i].start(j,hdf5=args.hdf5, study_name=args.study_name[i], ID=False)

			for i in gen:
				i._data.link()
			row_index, ids =  study_indexes(phenotype=phen.folder._data,genotype=(i.folder._data for i in gen),covariates=(i.metadata for i in pard))

		while True:
			if mapper.cluster=='n':
				SNPs_index, keys=mapper.get_next()
			else:
				ch=mapper.chunk_pop()
				if ch is None:
					SNPs_index=None
					break
				print ch
				SNPs_index, keys=mapper.get_chunk(ch)

			if isinstance(SNPs_index, type(None)):
				break

			Analyser.rsid=keys
			if np.sum(PD)==0:
				genotype=np.array([])
				genotype=merge_genotype(genotype, SNPs_index)
				genotype=genotype[:,row_index[0]]

			#TODO (low) add interaction

			a_test=np.array([])
			b_cov=np.array([])
			C=np.array([])
			a_cov=np.array([])
			b4=np.array([])

			if protocol.enable:
				regression_model=protocol.regression_model()
			else:
				regression_model=None

			if np.sum(PD)==0:
				a_test, b_cov, C, a_cov = meta_pard.get( SNPs_index=SNPs_index, regression_model=regression_model)
			else:
				a_test, b_cov, C, a_cov, b4 = meta_pard.get( SNPs_index=SNPs_index, B4=True, regression_model=regression_model)

			MAF=meta_pard.maf_pard(SNPs_index=SNPs_index)

			if args.maf!=0:
				filter=(MAF>args.maf) & (MAF<1-args.maf) & (MAF!=0.5)
				a_test=a_test[filter,:]
				Analyser.MAF=MAF[filter]
				if np.sum(PD)==0:
					genotype=genotype[filter,:]
				else:
					b4=b4[filter,:]
				Analyser.rsid=Analyser.rsid[filter]
				if a_test.shape[0]==0:
					print 'NO SNPs > MAF'
					continue
			else:
				Analyser.MAF=MAF

			a_inv=A_inverse(a_cov, a_test)
			N_con=a_inv.shape[1] - 1
			DF=( meta_pard.get_n_id()- a_inv.shape[1]  )

			if np.sum(PD)==0:

				while True:
					phenotype=phen.get_next(index=row_index[1])
					if isinstance(phenotype, type(None)):
						phen.folder.processed=0
						break

					b4 = B4(phenotype,genotype)

					t_stat=HASE(b4, a_inv, b_cov, C, N_con, DF)

					if mapper.cluster=='y':
						Analyser.cluster=True
						Analyser.chunk=ch
						Analyser.node=mapper.node[1]
					Analyser.t_stat=t_stat
					Analyser.threshold=args.thr
					Analyser.out=args.out
					Analyser.save_result(  meta_pard.phen_names )

					t_stat=None
					Analyser.t_stat=None
					gc.collect()

			else:
				t_stat,SE=HASE(b4, a_inv, b_cov, C, N_con, DF)

				Analyser.t_stat=t_stat
				Analyser.SE=SE
				Analyser.threshold=args.thr
				Analyser.out=args.out
				Analyser.save_result(  meta_pard.phen_names )

				t_stat=None
				Analyser.t_stat=None
				gc.collect()



################################### TO DO EVERYTHING IN ONE GO ##############################

	elif args.mode=='regression':

		#ARG_CHECKER.check(args,mode='regression')

		print ('START regression mode...')

		phen=Reader('phenotype')
		phen.start(args.phenotype)
		phen.permutation=args.permute_ph

		cov=Reader('covariates')
		cov.start(args.covariates)

		if cov.folder.n_files>1:
			raise ValueError('In covariates folder should be only one file!')

		gen=[]
		for i,j in enumerate(args.genotype):
			gen.append(Reader('genotype'))
			gen[i].start(j,hdf5=args.hdf5, study_name=args.study_name[i], ID=False)


		if args.mapper_name is not None:
			mapper=Mapper(args.mapper_name)
			mapper.chunk_size=MAPPER_CHUNK_SIZE
			mapper.genotype_names=args.genotype
			if args.snp_id_inc is not None:
				mapper.include=np.array(pd.DataFrame.from_csv(args.snp_id_inc)['rsid'].tolist())
			if args.snp_id_exc is not None:
				mapper.exclude=np.array(pd.DataFrame.from_csv(args.snp_id_exc)['rsid'].tolist())
			mapper.load(args.mapper)
			mapper.cluster=args.cluster
			mapper.node=args.node
		else:
			mapper=None

		Analyser=HaseAnalyser()
		Analyser.threshold=args.thr
		Analyser.out=args.out
		haseregression(phen,gen,cov, mapper, Analyser,args.maf,intercept=args.intercept)

end = time.time()

print ('experiment finished in {} s'.format((end-start)))