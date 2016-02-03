import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from hdgwas.tools import Timer,HaseAnalyser, Reference
import argparse
import pandas as pd
import numpy as np


if __name__=="__main__":

	os.environ['HASEDIR']=os.path.dirname(os.getcwd())
	parser = argparse.ArgumentParser(description='Script analyse results of HASE')
	parser.add_argument("-r", required=True,help="path to hase results")
	parser.add_argument("-o", "--out", type=str, required=True,help="path to save result folder")
	parser.add_argument("-df", type=str,required=True, help="degree of freedom = ( #subjects in study  - #covariates - 1 )")
	#TODO (low) add reference panel
	args = parser.parse_args()
	Analyser=HaseAnalyser()


	Analyser.DF=np.float(args.df)
	Analyser.result_path=args.r

	while True:
		Analyser.summary()
		if Analyser.results is None:
			break
		print('Saving data...')
		df=pd.DataFrame.from_dict(Analyser.results)
		df.to_csv( os.path.join(args.out,'results'+'.csv'), sep=" "  )

