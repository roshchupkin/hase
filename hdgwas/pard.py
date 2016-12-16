

from hdgwas.hdregression import HASE, A_covariates, A_tests, B_covariates, C_matrix, A_inverse,B4
from hdgwas.tools import study_indexes, Timer
import numpy as np
import os

def partial_derivatives(save_path=None,COV=None,PHEN=None, GEN=None,
                        MAP=None, MAF=None, R2=None, B4_flag=False, study_name=None,intercept=True):

    row_index, ids =  study_indexes(phenotype=PHEN.folder._data,genotype=GEN.folder._data,covariates=COV.folder._data)

    metadata={}



    #TODO (mid) add parameter to compute PD only for new phenotypes or cov
    metadata['id']=ids
    metadata['MAF']=[]
    metadata['filter']=[]
    metadata['names']=[] #TODO (low) change to cov_names
    metadata['phenotype']=[]
    b_cov=[]
    C=[]
    a_test=[]
    b4=[]

    covariates=COV.get_next(index=row_index[2])
    if intercept:
        metadata['names'].append(study_name+ '_intercept')
    metadata['names']=metadata['names']+[ study_name+ '_' + i for i in COV.folder._data.get_names() ]

    a_cov=A_covariates(covariates,intercept=intercept)
    np.save(os.path.join(save_path,study_name+'_a_cov.npy'),a_cov)

    i=1

    with Timer() as t_phen:

        while True:

            phenotype=PHEN.get_next(index=row_index[1])
            if isinstance(phenotype, type(None)):
                b_cov=np.concatenate(b_cov, axis=1)
                C=np.concatenate(C, axis=0)
                np.save(os.path.join(save_path,study_name+'_b_cov.npy'),b_cov)
                np.save(os.path.join(save_path,study_name+'_C.npy'),C)
                break

            metadata['phenotype']=metadata['phenotype']+ list(PHEN.folder._data.get_names())
            b_cov.append(B_covariates(covariates,phenotype,intercept=intercept))
            C.append(C_matrix(phenotype))

    print ('Time to PD phenotype {} is {} s'.format(np.array(C).shape, t_phen.secs))

    N_snps_read=0
    while True:
        with Timer() as t_gen:
            genotype=GEN.get_next()
            if isinstance(genotype, type(None)):
                np.save(os.path.join(save_path,study_name+'_a_test.npy'), np.concatenate(a_test) )
                np.save(os.path.join(save_path,study_name+'_metadata.npy'),metadata)
                if B4_flag:
                    b4=np.concatenate(b4, axis=0)
                    np.save(os.path.join(save_path,study_name+'_b4.npy'),b4)
                break
            flip = MAP.flip[N_snps_read:N_snps_read + genotype.shape[0], 0]
            N_snps_read += genotype.shape[0]
            flip_index=(flip==-1)
            genotype=np.apply_along_axis(lambda x: flip*(x-2*flip_index) ,0,genotype)
            genotype=genotype[:,row_index[0]]
            maf=np.mean(genotype, axis=1)/2
            metadata['MAF']=metadata['MAF']+list(maf)

            #TODO (low) add interaction
            a_test.append(A_tests(covariates,genotype,intercept=intercept))

            if B4_flag:
                #works only when all phenotypes in one chunk, if not, do not use this option!
                #it would use to much disk space anyway
                if len([f for f in PHEN.folder.files if f!='info_dic.npy' ])>1:
                    print 'pd_full flag disabled!'
                    B4_flag=False
                    continue
                PHEN.folder.processed=0
                if isinstance(phenotype, type(None)):
                    phenotype=PHEN.get_next(index=row_index[1])
                b4.append(B4(phenotype,genotype))

        print ('Time to PD genotype {} is {} s'.format(genotype.shape, t_gen.secs))







