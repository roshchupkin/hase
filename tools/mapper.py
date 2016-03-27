import pandas as pd
import numpy as np
import argparse
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from hdgwas.tools import Reference, Mapper, Timer
from hdgwas.hash import *
import gc

if __name__=='__main__':

	os.environ['HASEDIR']=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

	parser = argparse.ArgumentParser(description='Script to map studies for meta-stage')
	parser.add_argument("-g", "--genotype",nargs='+',required=True, type=str, help="path/paths to genotype data folder")
	parser.add_argument('-study_name', nargs='+',required=True, default=None, help=' Study names')
	parser.add_argument("-o", "--out", type=str,required=True, help="path to save result folder")
	parser.add_argument('-ref_name', type=str, default='1000Gp1v3_ref', help='Reference panel name')

	args = parser.parse_args()

	reference=Reference()
	print('Use {} as reference panel'.format(args.ref_name))
	reference.reference_name=args.ref_name
	reference.load()
	print reference.dataframe.head()
	print reference.dataframe.info()
	mapper=Mapper()

	for j,i in enumerate(args.genotype):
		mapper.hash=None
		index=[]
		flip=[]
		gc.collect()
		mapper.hash=HashTableRSID()
		with Timer() as t:
			df=pd.read_hdf(os.path.join(i,'probes', args.study_name[j] +'.h5'),'probes', where='columns=[ID,allele1,allele2]')
			df=df[['ID','allele1','allele2']]
			with Timer() as t1:
				mapper.hash.fill(df)
			print ('time to fill hash {}s'.format(t1.secs))
			for k,l in enumerate(reference.dataframe.iterrows()):
				ind,fl=mapper.hash.get_map( l[1].tolist() )
				index.append(ind)
				flip.append(fl)
		print ('time to add ID {}s'.format(t.secs))
		print ('Saving results for {} to {} ...'.format(i,args.out))
		np.save(os.path.join(args.out,'values_'+args.ref_name+'_'+args.study_name[j]+'.npy'),index)
		np.save(os.path.join(args.out,'flip_'+args.ref_name+'_'+args.study_name[j]+'.npy'),flip)

	np.save(os.path.join(args.out,'keys_'+args.ref_name+'.npy'),reference.dataframe['ID'].tolist())
	print ('Data successfully saved')


	overlap=np.where(index==1)
	not_overlap=np.where(index==-1)
	flip_index=np.where(flip==-1)
	print 'There are {} common variances with reference panel, which would be included in study'.format(len(overlap[0]))
	print 'There are {} variances excluded from study (not found in reference panel)'.format(df['ID'].shape[0] -len(overlap[0])   )
	if len(not_overlap[0])!=0:
		print 'Examples not found in reference panel:', df['ID'][not_overlap[0][:10]]
	print 'There are {} flipped variances',format(len(flip_index[0]))
	if len(flip_index[0])!=0:
		print 'Examples flipped:', df['ID'][flip_index[0][:10]]






