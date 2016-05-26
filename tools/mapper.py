import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import *
if PYTHON_PATH is not None:
	for i in PYTHON_PATH: sys.path.insert(0,i)
import h5py
import pandas as pd
import numpy as np
import argparse
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
		gc.collect()
		with Timer() as t:
			df=pd.read_hdf(os.path.join(i,'probes', args.study_name[j]+'.h5'),'probes')
			print df.head()
			print df.info()
			if "CHR" in df.columns and 'bp' in df.columns:
				mapper.hash=HashTablePOS()
			else:
				if ':' in df.ID.iloc[0] and ':' in df.ID.iloc[1]:
					CHR=[]
					bp=[]
					for i in df.ID:
						s=i.split(":")
						CHR.append(s[0])
						bp.append(s[1])
					CHR=np.array(CHR,dtype=np.int8)
					bp=np.array(bp)
					if np.max(CHR)<23 and np.min(CHR)>0:
						df['CHR']=CHR
						df['bp']=bp
						mapper.hash=HashTablePOS()
					else:
						print 'No CHR and pb info...'
						mapper.hash=HashTableRSID()
				else:
					print 'No CHR and pb info...'
					mapper.hash=HashTableRSID()

			with Timer() as t1:
				mapper.hash.fill(df)
			print ('time to fill hash {}s'.format(t1.secs))
			index=[]
			flip=[]
			for l in reference.dataframe.iterrows():
				ind,fl=mapper.hash.get_map( l[1])
				index.append(ind)
				flip.append(fl)
		print ('time to add ID {}s'.format(t.secs))
		print ('Saving results for {} to {} ...'.format(i,args.out))
		np.save(os.path.join(args.out,'values_'+args.ref_name+'_'+args.study_name[j]+'.npy'),index)
		np.save(os.path.join(args.out,'flip_'+args.ref_name+'_'+args.study_name[j]+'.npy'),flip)

	np.save(os.path.join(args.out,'keys_'+args.ref_name+'.npy'),reference.dataframe['ID'].tolist())
	print ('Data successfully saved')

	index=np.array(index)
	flip=np.array(flip)
	overlap=np.where(index!=-1)
	not_overlap=np.where(index==-1)
	flip_index=np.where(flip==-1)
	print 'There are {} common variances with reference panel, which would be included in study'.format(overlap[0].shape[0] )
	print 'There are {} variances excluded from study (not found in reference panel)'.format( df.ID.shape[0]-overlap[0].shape[0]  )
	print 'There are {} flipped variances'.format(len(flip_index[0]))







