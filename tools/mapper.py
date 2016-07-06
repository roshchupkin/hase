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
	parser.add_argument("-g",required=True, type=str, help="path/paths to genotype data folder")
	parser.add_argument('-study_name',type=str,required=True, default=None, help=' Study names')
	parser.add_argument("-o", "--out", type=str,required=True, help="path to save result folder")
	parser.add_argument('-ref_name', type=str, default='1000Gp1v3_ref', help='Reference panel name')
	parser.add_argument('-mismatch_table',action='store_true',default=False, help='Save table with mismatch IDs')
	parser.add_argument('-flipped_table',action='store_true',default=False, help='Save table with mismatch IDs')
	parser.add_argument('-probe_chunk',type=int,default=10000, help='Probes chunk')
	parser.add_argument('-ref_chunk',type=int,default=10000, help='Probes chunk')

	args = parser.parse_args()
	print args

	reference=Reference()
	print('Use {} as a reference panel'.format(args.ref_name))
	reference.reference_name=args.ref_name
	reference.chunk=args.ref_chunk
	reference.load()
	mapper=Mapper()


	mapper.hash=None
	gc.collect()
	df_ch=pd.read_hdf(os.path.join(args.g,'probes', args.study_name+'.h5'),'probes',chunksize=args.probe_chunk)
	ID=[]
	with Timer() as t1:
		for i,df in enumerate(df_ch):
			if i==0:
				print df.head()
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
							print 'No CHR and bp info...'
							mapper.hash=HashTableRSID()
					else:
						print 'No CHR and bp info...'
						mapper.hash=HashTableRSID()

			df.apply(mapper.hash.fill_id, axis=1)
	print ('time to fill hash by {} keys {}s'.format(mapper.hash.n_keys,t1.secs))
	ind_fl=np.array([])
	with Timer() as t:
		for k,ref in enumerate(reference.dataframe):
			ind_fl=np.append( ind_fl, ref.apply(mapper.hash.get_map, axis=1))
			ID=ID+ref.ID.tolist()

	print ('time to add {} IDs {}s'.format(len(ID),t.secs))
	index,flip=zip(*ind_fl)

	print ('Saving results for {} to {} ...'.format(i,args.out))
	np.save(os.path.join(args.out,'values_'+args.ref_name+'_'+args.study_name+'.npy'),index)
	np.save(os.path.join(args.out,'flip_'+args.ref_name+'_'+args.study_name+'.npy'),flip)

	np.save(os.path.join(args.out,'keys_'+args.ref_name+'.npy'),ID)

	full_index_list=np.arange(mapper.hash.n_keys)
	mismatch_index=np.setdiff1d(full_index_list,index )
	np.save(os.path.join(args.out,'excluded_'+args.ref_name+'_'+args.study_name+'.npy'),mismatch_index)
	print ('Data successfully saved')

	df_ch.close()
	store=pd.HDFStore(os.path.join(args.g,'probes', args.study_name+'.h5'),'r')


	index=np.array(index)
	flip=np.array(flip)
	overlap=np.where(index!=-1)
	not_overlap=np.where(index==-1)
	flip_index=np.where(flip==-1)
	print 'There are {} common variances with reference panel, which would be included in study'.format(overlap[0].shape[0] )
	print 'There are {} variances excluded from study (not found in reference panel)'.format( mapper.hash.n_keys-overlap[0].shape[0]  )
	if args.mismatch_table:
		df_mismatch=store.select('probes',where=mismatch_index)
		df_mismatch.to_csv(os.path.join(args.out,'excluded_ID_info.csv'))
	else:
		print store.select('probes',where=mismatch_index[:10])

	print 'There are {} flipped variances'.format(len(flip_index[0]))
	if args.flipped_table and flip_index[0].shape[0]!=0:
		df_flipped=store.select('probes',where=flip_index[0])
		df_flipped.to_csv(os.path.join(args.out,'flipped_ID_info.csv'))
	elif flip_index[0].shape[0]!=0:
		print store.select('probes',where=flip_index[:10])







