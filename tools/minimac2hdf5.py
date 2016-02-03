
import numpy as np

import os
import sys
import argparse
import h5py

def probes_minimac2hdf5(data_path, save_path,study_name):
	n=[]
	f=open(data_path,'r')
	for i,j in enumerate(f):
		n.append((j[:-1]))
	f.close()
	n=np.array(n)
	store=h5py.File(os.path.join(save_path,'probes','ID.h5'), 'w')
	store.create_dataset('RSID',data=n,compression='gzip',compression_opts=9 )
	store.close()


def ind_minimac2hdf5(data_path, save_path,study_name):
	n=[]
	f=open(data_path,'r')
	for i,j in enumerate(f):
		n.append((j[:-1]))
	f.close()
	n=np.array(n)
	store=h5py.File(os.path.join(save_path,'individuals',study_name+'.h5'), 'w')
	store.create_dataset('individuals',data=n,compression='gzip',compression_opts=9 )
	store.close()

def id_minimac2hdf5(data_path,id, save_path):

	n=[]
	f=open(data_path,'r')
	for i,j in enumerate(f):
		try:
			n.append(np.float(j))
		except:
			n.append(np.float(-1))
	n=np.array(n)
	f.close()
	store=h5py.File(os.path.join(save_path,'genotype',id+'.h5'), 'w')
	store.create_dataset(id,data=n,compression='gzip',compression_opts=9 )
	store.close()


if __name__=="__main__":

	parser = argparse.ArgumentParser(description='Script to convert Minimac data')
	parser.add_argument("-study_name", required=True, type=str, help="Study specific name")
	parser.add_argument("-id", type=str, help="subject id")
	parser.add_argument("-data",required=True, type=str, help="path to file")
	parser.add_argument("-out",required=True, type=str, help="path to results save folder")
	parser.add_argument("-flag",required=True,type=str,choices=['genotype','individuals','probes'], help="path to file with SNPs info")


	args = parser.parse_args()

	print args
	try:
		print ('Creating directories...')
		os.mkdir(os.path.join(args.out,'genotype') )
		os.mkdir(os.path.join(args.out,'individuals') )
		os.mkdir(os.path.join(args.out,'probes') )
	except:
		print('Directories "genotype","probes","individuals" are already exist in {}...'.format(args.out))

	if args.id is not None and args.flag=='genotype':
		id_minimac2hdf5(args.data, args.id, args.out)
	elif args.flag=='probes':
		probes_minimac2hdf5(args.data, args.out, args.study_name)
	elif args.flag=='individuals':
		ind_minimac2hdf5(args.data, args.out,args.study_name)

