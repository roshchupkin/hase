import pandas as pd
import numpy as np
import argparse
import os
from hdgwas.tools import Reference, Mapper, Timer
import string


if __name__=='__main__':

	os.environ['HASEDIR']=os.path.dirname(os.getcwd())

	parser = argparse.ArgumentParser(description='Script to map studies for meta-stage')
	parser.add_argument("-g", "--genotype",nargs='+',required=True, type=str, help="path/paths to genotype data folder")
	parser.add_argument('-study_name', type=str,required=True, default=None, help=' Study specific name')
	parser.add_argument("-o", "--out", type=str,required=True, help="path to save result folder")
	parser.add_argument("-dic",default=False ,type=bool, help="flag, if you already saved reference keys from mappper before")
	parser.add_argument('-ref_name', type=str, default='1000Gp1v3_20101123_ALL_impute_legend', help='Reference panel name')

	args = parser.parse_args()


	if not args.dic:
		reference=Reference(args.study_name)
		print('Use {} as reference panel'.format(args.ref_name))
		reference.reference_name=args.ref_name
		reference.load()
		keys=reference.dataframe['ID']
	elif args.dic:
		keys=np.load(os.path.join(args.out,'keys_'+args.study_name+'.npy'))


	mapper=Mapper(args.study_name)
	mapper.fill(keys, args.ref_name ,reference=True)


	for i in args.genotype:
		with Timer() as t:
			ID=pd.read_hdf(os.path.join(i,'probes', 'ID.h5'),'RSID')
			mapper.push(ID, genotype_name=i.split('/')[-1])
		print ('time to add ID {}s'.format(t.secs))

	print ('Saving results to {} ...'.format(args.out))
	mapper.save(args.out)
	print ('Data successfully saved')





