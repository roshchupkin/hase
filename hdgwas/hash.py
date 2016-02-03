import numpy as np
import cPickle as pickle
import os
class HashTable(object):

    def __init__(self):
        self.hash_table=np.zeros([10**9])
        self.hash_dic={}
        self.protected=False

    def rsid_hash(self,value):
        pass

    def fill(self,ordered_values):
        if not self.protected:
            for j,i in enumerate(ordered_values):
                try:
                    if i[:2]=='rs':
                        m=int(i[2:])
                        if self.hash_table[m]!=0:
                            print 'collision!'
                        else:
                            self.hash_table[m]=j
                    else:
                        self.hash_dic[i]=j
                except:
                    self.hash_dic[i]=j

            self.protected=True
        else:
            raise ValueError('Hash table already filled!')


    def get_index(self,value):
            try:
                if value[:2]=='rs':
                    return int(self.hash_table[int(value[2:])])
                else:
                    return int(self.hash_dic[value])
            except:
                raise ValueError('Value {} does not include to table!'.format(value))


    def save(self,path, filename):
        with open(os.path.join(path,filename), 'wb') as output:
            pickle.dump(self, output, -1)



if __name__=="__main__":
    print 'Hash'


