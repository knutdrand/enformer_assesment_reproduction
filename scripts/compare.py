import bionumpy as bnp
import numpy as np
import sparse
# gene_name = 'ENSG00000228327.6'
gene_name = 'ENSG00000169885.10'
# gene_name = 'ENSG00000011009.11'
my_pathname= '../enformer_assesment_reproduction/%s.npz' % gene_name

def get_dense(filename):
    data = np.load(filename)
    return sparse.COO(data['coords'], data['data'], shape=data['shape']).todense()



my_dense  = get_dense(my_pathname)
my_seq = bnp.EncodedArray(np.argmax(my_dense[4:, :, 0], axis=0).astype(np.uint8), bnp.DNAEncoding)[:-1]

other_pathname = '../data/debug/sequence10K/chr1/%s.npz' % gene_name
other_dense = get_dense(other_pathname)
other_seq = bnp.EncodedArray(np.argmax(other_dense[[0, 3, 2, 1], :, 0], axis=0).astype(np.uint8), bnp.DNAEncoding)[1:]
print(my_seq)
print(other_seq)
idxs = np.flatnonzero(my_seq != other_seq)
print(my_seq[idxs])
print(other_seq[idxs])
