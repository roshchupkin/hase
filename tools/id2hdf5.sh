#!/bin/sh

CWD=`pwd`

id=$1
GENOTYPE_DIR=$2
SAVE_DIR=$3
HASE_DIR=$4

cd $GENOTYPE_DIR

for i in $(cat ${SAVE_DIR}/files_order.txt ); do
zcat $i.dose.gz | awk "{if( match($1,${id}) ){\$1=\"\"; \$2=\"\"; printf \"%s \",\$0; exit}}" >>${SAVE_DIR}/${id}_row.txt
done

cat ${SAVE_DIR}/${id}_row.txt | gawk 'BEGIN{ind=1}{for(i=1;i<=NF;i+=1){SNPs[ind]=$i; ind+=1}}END{for (j=1;j<ind;j+=1){print SNPs[j]}}'>>${SAVE_DIR}/${id}.txt

rm ${SAVE_DIR}/${id}_row.txt

python ${HASE_DIR}/tools/minimac2hdf5.py -flag genotype -id ${id} -data ${SAVE_DIR}/${id}.txt -out ${SAVE_DIR}

rm ${SAVE_DIR}/${id}.txt

cd $CWD
