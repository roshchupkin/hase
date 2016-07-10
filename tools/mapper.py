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
	parser.add_argument('-collision',action='store_true',default=False, help='Print collision warning')
	parser.add_argument('-chunk',type=int,default=2000000, help='Chunk size')
	args = parser.parse_args()
	print args

	probes=pd.HDFStore(os.path.join(args.g,'probes', args.study_name+'.h5'),'r')
	probes_n_rows=probes.get_storer('probes').nrows
	chunk_size = np.min([args.chunk,probes_n_rows])

	print ('Merge chunk size {}'.format(chunk_size))
	match_key=np.array([])
	match_index=np.array([])
	flip_key=np.array([])
	flip_index=np.array([])
	ID=np.array([])

	del_counter_ref={}
	ID2CHR=False
	IDconv=False
	merge_on={
			'ID':{
				'straight':["ID",'allele1','allele2'],
				  'reverse':["ID",'allele2','allele1']
				 },

			'CHR':{
				 'straight':["CHR",'bp','allele1','allele2'],
				 'reverse':["CHR",'bp','allele2','allele1']

			      }
		  	}

	for p in xrange(int(np.ceil(probes_n_rows / chunk_size))):
		print 'p',p

		p_start_i = p * chunk_size
		p_stop_i  = min((p + 1) * chunk_size, probes_n_rows)

		a = probes.select('probes', start = p_start_i, stop = p_stop_i)

		if p==0:
			if "CHR" in a.columns and 'bp' in a.columns:
				ID2CHR=True
				merge=merge_on['CHR']
			else:
				if ':' in a.ID.iloc[0] and ':' in a.ID.iloc[1]:
					CHR=[]
					bp=[]
					for i in a.ID:
						s=i.split(":")
						CHR.append(s[0])
						bp.append(s[1])
					CHR=np.array(CHR,dtype=np.int8)
					bp=np.array(bp)
					if np.max(CHR)<23 and np.min(CHR)>0:
						a['CHR']=CHR
						a['bp']=bp
						ID2CHR=True
						IDconv=True
						merge=merge_on['CHR']
					else:
						print 'No CHR and bp info...'
						merge=merge_on['ID']
				else:
					print 'No CHR and bp info...'
					merge=merge_on['ID']

		elif IDconv:
			def f(x):
				s=x.ID.split(':')
				return s[0],s[1]
			CHR_bp=a.apply(f, s[0],s[1], axis=1 )
			a['CHR'],a['bp']=zip(*CHR_bp)

		a['counter_prob']=np.arange(p_start_i,p_stop_i,dtype='int32')

		reference=Reference()
		reference.reference_name=args.ref_name
		reference.chunk=chunk_size
		reference.load()
		counter_ref=0
		for r,b in enumerate(reference.dataframe):
			if r==0:
				if np.sum(np.array([ 1 if i in reference.columns else 0 for i in b.columns.tolist()  ]))!=len(reference.columns):
					raise ValueError('Reference table should have {} columns'.format(reference.columns))


			print 'r',r
			if p==0:
				ID=np.append(ID,b.ID)

			b['counter_ref']=np.arange(counter_ref,counter_ref+b.shape[0],dtype='int32')
			counter_ref+=b.shape[0]

			if len(match_index) or len(flip_index):
				print 'matched {}'.format(match_index.shape[0])
				print 'flipped {}'.format(flip_index.shape[0])
				if del_counter_ref.get(r) is not None:
					with Timer() as t:
						b=b[~b.counter_ref.isin(del_counter_ref[r])]
					print 'time {}'.format(t.secs)

			match_df = pd.merge(b,a, left_on=merge['straight'], right_on=merge['straight'])
			flip_df=pd.merge(b[~b.counter_ref.isin(match_df.counter_ref)],a, left_on=merge['reverse'], right_on=merge['straight'])

			if len(match_df):
				print match_df.head()
				match_key=np.append(match_key,match_df.counter_ref)
				match_index=np.append(match_index,match_df.counter_prob)
				if del_counter_ref.get(r) is None:
					del_counter_ref[r]=match_key
				else:
					del_counter_ref[r]=np.append(del_counter_ref[r], match_key)
			if len(flip_df):
				print flip_df.head()
				flip_key=np.append(flip_key,flip_df.counter_ref)
				flip_index=np.append(flip_index,flip_df.counter_prob)
				if del_counter_ref.get(r) is None:
					del_counter_ref[r]=flip_key
				else:
					del_counter_ref[r]=np.append(del_counter_ref[r], flip_key)
			gc.collect()
		reference=None

	match_key = match_key.astype(np.int32, copy=False)
	match_index = match_index.astype(np.int32, copy=False)
	flip_key = flip_key.astype(np.int32, copy=False)
	flip_index = flip_index.astype(np.int32, copy=False)

	index=np.ones(ID.shape[0])*-1
	flip=np.ones(ID.shape[0])
	index[match_key]=match_index
	index[flip_key]=flip_index
	flip[flip_key]=-1
	print ('Saving results for {} to {} ...'.format(args.study_name,args.out))
	np.save(os.path.join(args.out,'values_'+args.ref_name+'_'+args.study_name+'.npy'),index)
	np.save(os.path.join(args.out,'flip_'+args.ref_name+'_'+args.study_name+'.npy'),flip)
	np.save(os.path.join(args.out,'keys_'+args.ref_name+'.npy'),ID)
	print ('Data successfully saved')

	mismatch_index=np.setdiff1d(np.arange(probes_n_rows),np.append(match_index,flip_index) )


	print 'There are {} common variances with reference panel, which will be included in study'.format(np.where(index!=-1)[0].shape[0] )
	print 'There are {} variances from reference panel, which were not found in probes'.format(np.where(index==-1)[0].shape[0] )
	print 'There are {} variances excluded from study (not found in reference panel)'.format( probes_n_rows-np.where(index!=-1)[0].shape[0]  )
	if args.mismatch_table and mismatch_index.shape[0]!=0:
		df_mismatch=probes.select('probes',where=mismatch_index)
		df_mismatch.to_csv(os.path.join(args.out,'Mismatch_ID_info.csv'))
		print 'Mismatch ID info saved to {}'.format(os.path.join(args.out,'Mismatch_ID_info.csv'))
	elif mismatch_index.shape[0]!=0:
		print ('Mismatch examples:')
		print probes.select('probes',where=mismatch_index[:10])

	print 'There are {} flipped variances'.format(len(flip_index))
	if args.flipped_table and flip_index.shape[0]!=0:
		df_flipped=probes.select('probes',where=flip_index)
		df_flipped.to_csv(os.path.join(args.out,'flipped_ID_info.csv'))
		print 'Flipped ID info saved to {}'.format(os.path.join(args.out,'flipped_ID_info.csv'))
	elif flip_index.shape[0]!=0:
		print ('Flipped examples:')
		print probes.select('probes',where=flip_index[:10])



	exit(0)

	#N_ref=reference.get_storer('reference').nrows
	#index=np.ones(N_ref)*-1
	#flip=np.ones(N_ref)

	df_ch=pd.read_hdf(os.path.join(args.g,'probes', args.study_name+'.h5'),'probes',chunksize=args.probe_chunk)

	id_ref=reference.select('reference',where='columns=ID')
	id_ref=map(hash,id_ref)

	id_probes=[]
	with Timer() as t:
		for df in df_ch:
			id_probes=id_probes + map(hash,df.ID)
	print 'time {}'.format(t.secs)
	with Timer() as t:
		indeces=np.in1d(id_ref,id_probes)
	# for i, df in enumerate(df_ch):
	# 	ind_fl=df.apply(get_info, axis=1)
	# 	if i==50:
	# 		break
	print 'time {}'.format(t.secs)
	exit(0)

	reference=Reference()
	print ('********************************')
	print('Use {} as a reference panel'.format(args.ref_name))
	reference.reference_name=args.ref_name
	reference.chunk=args.ref_chunk
	reference.load()
	print reference._df.head()
	print ('********************************')
	mapper=Mapper()

	super_set=set()
	mapper.hash=None
	gc.collect()
	df_ch=pd.read_hdf(os.path.join(args.g,'probes', args.study_name+'.h5'),'probes',chunksize=args.probe_chunk)
	ID=[]
	with Timer() as t1:
		df_ch=reference.dataframe
		for i,df in enumerate(df_ch):
			#df=reference._df
			if i==0:
				print 'Probes table head:'
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
			mapper.hash._print=args.collision
			df.apply(mapper.hash.fill_id, axis=1)
			if i==2: break
	print ('time to fill hash by {} keys {}s'.format(mapper.hash.n_keys,t1.secs))
	ind_fl=np.array([])
	with Timer() as t:
		ind_fl=np.append( ind_fl, reference._df.apply(mapper.hash.get_map, axis=1))
		ID=ID+reference._df.ID.tolist()
		for k,ref in enumerate(reference.dataframe):
			ind_fl=np.append( ind_fl, df.apply(mapper.hash.get_map, axis=1))
			ID=ID+ref.ID.tolist()
			break

	print ('time to add {} IDs {}s'.format(len(ID),t.secs))
	index,flip=zip(*ind_fl)

	print ('Saving results for {} to {} ...'.format(args.study_name,args.out))
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

	print 'There are {} common variances with reference panel, which will be included in study'.format(overlap[0].shape[0] )
	print 'There are {} variances from reference panel, which were not found in probes'.format(not_overlap[0].shape[0] )
	print 'There are {} variances excluded from study (not found in reference panel)'.format( mapper.hash.n_keys-overlap[0].shape[0]  )
	if args.mismatch_table and mismatch_index.shape[0]!=0:
		df_mismatch=store.select('probes',where=mismatch_index)
		df_mismatch.to_csv(os.path.join(args.out,'Mismatch_ID_info.csv'))
		print 'Mismatch ID info saved to {}'.format(os.path.join(args.out,'Mismatch_ID_info.csv'))
	elif mismatch_index.shape[0]!=0:
		print ('Mismatch examples:')
		print store.select('probes',where=mismatch_index[:10])

	print 'There are {} flipped variances'.format(len(flip_index[0]))
	if args.flipped_table and flip_index[0].shape[0]!=0:
		df_flipped=store.select('probes',where=flip_index[0])
		df_flipped.to_csv(os.path.join(args.out,'flipped_ID_info.csv'))
		print 'Flipped ID info saved to {}'.format(os.path.join(args.out,'flipped_ID_info.csv'))
	elif flip_index[0].shape[0]!=0:
		print ('Flipped examples:')
		print store.select('probes',where=flip_index[0][:10])







