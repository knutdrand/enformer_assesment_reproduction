import os

import pandas as pd
import numpy as np
import gzip
import sparse
import sys

nSym = 8  # ATGCx2 for paternal and maternal
# nSub = 1161
nSub = 3

wpl = 50  # Number of words per line in hg38
import bionumpy as bnp


def extract_seq(gene_table, ch, refpath, filepath):
    win = 10000
    winSiz = '10K'
    genome = bnp.Genome.from_file(f'{refpath}/genome.fa')
    sequence = genome.read_sequence()
    gene_table.columns = ['ensg', 'chr', 'start', 'stop']
    gene_table = gene_table.loc[gene_table.chr == ch]
    #locations = bnp.LocationEntry([str(c) for c in gene_table['chr']], gene_table['start'].astype(int))
    #locations = genome.get_locations(locations, has_numeric_chromosomes=True)
    intervals = bnp.Interval([f'chr{i}' for i in gene_table['chr']], gene_table['start'].astype(int)-1, gene_table['start'].astype(int)+win*2)
    #windo = locations.get_windows(flank=int(win))
    gene_sequences = sequence[intervals] #.raw().view('S1'
    nGene = len(gene_table.index)
    for gene_id in np.arange(nGene):
        seqR = gene_sequences[gene_id]
        print(os.getcwd())
        specific_filepath = filepath + 'variantNucleotide' + winSiz + '/' + gene_table.iloc[gene_id, 0] + '.csv'
        f = open(specific_filepath)
        line = f.readline()
        f.close()
        seqP = np.concatenate([seqR] * nSub).reshape(nSub, -1)
        seqM = np.concatenate([seqR] * nSub).reshape(nSub, -1)
        if line != '':
            snp = pd.read_csv(specific_filepath, header=None,
                              encoding='latin1')
            positions = snp.iloc[:, 2]-int(intervals[gene_id].start)
            P_seqs = [''.join(str(snp)[2] for snp in column) for column in snp.iloc[:, 3:3+nSub].to_numpy().T]
            M_seqs = [''.join(str(snp)[2] for snp in column) for column in snp.iloc[:, 3:3+nSub].to_numpy().T]
            for k, (P_seq, M_seq) in enumerate(zip(P_seqs, M_seqs)):
                seqP[k, positions] = P_seq
                seqM[k, positions] = M_seq
        onehotM = seqM.T == bnp.as_encoded_array('ATGC').reshape(4, 1, 1)
        onehotP = seqP.T == bnp.as_encoded_array('ATGC').reshape(4, 1, 1)
        onehot = np.concatenate([onehotP, onehotM])
        print(os.getcwd())
        np.save(f'{gene_table.iloc[gene_id, 0]}_new.npy', onehot.astype('i8'))

        onehot = sparse.COO.from_numpy(onehot)
        folder = filepath + '/sequence' + winSiz + '/chr' + str(ch) + '/'
        os.makedirs(folder, exist_ok=True)
        sparse.save_npz(folder + gene_table.iloc[gene_id, 0], onehot)
