import numpy as np
import cPickle as pickle
import os


class HashTable(object):

    def __init__(self):
        self.hash_table=np.zeros([250*10**6])
        self.hash_dic={}
        self.protected=False

    def fill(self,df):
        if not self.protected:
            for j,i in enumerate(df['ID']):
                try:
                    if i[:2]=='rs':
                        m=int(i[2:])
                        if self.hash_table[m]!=0 or m>len(self.hash_table):
                            print 'collision {}!'.format(i)
                        else:
                            self.hash_table[m]=[j]
                    else:
                        self.hash_dic[i]=[j]
                except:
                    self.hash_dic[i]=[j]

            self.protected=True
        else:
            raise ValueError('Hash table already filled!')

    def get_index(self,value):
        try:
            if value[:2]=='rs':
                v=int(value[2:])
                r=self.hash_table[v]
                if len(r)==1: r=[-1,'R','R']
                return self.hash_table[v]
            else:
                return self.hash_dic[value]
        except:
            return -1,'R','R'

class HashTablePOS(HashTable):
    pass

class HashTableRSID(HashTable):

    def fill(self,df):
        if not self.protected:
            for j,i in enumerate(df['ID']):
                try:
                    if i[:2]=='rs':
                        m=int(i[2:])
                        if self.hash_table[m]!=0 or m>len(self.hash_table):
                            print 'collision {}!'.format(i)
                        else:
                            self.hash_table[m]=[j,df['allele1'][j],df['allele2'][j]]
                    else:
                        self.hash_dic[i]=[j,df['allele1'][j],df['allele2'][j]]
                except:
                    self.hash_dic[i]=[j,df['allele1'][j],df['allele2'][j]]

            self.protected=True
        else:
            raise ValueError('Hash table already filled!')


    def get_index(self,value):
        try:
            if value[:2]=='rs':
                v=int(value[2:])
                r=self.hash_table[v]
                if len(r)==1: r=[-1,'R','R']
                return self.hash_table[v]
            else:
                return self.hash_dic[value]
        except:
            return -1,'R','R'

    def get_map(self,info):

        ind,A1,A2=self.get_index(info[0])

        if ind==-1:
            return -1,0
        if A1==info[1] and A2==info[2]:
            return ind, 1
        if A1==info[2] and A2==info[1]:
            return ind, -1
        else:
            return -1,0

    def save(self,path, filename):
        with open(os.path.join(path,filename), 'wb') as output:
            pickle.dump(self, output, -1)

if __name__=="__main__":
    print 'Hash'


