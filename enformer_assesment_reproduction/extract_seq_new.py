import os

import pandas as pd
import numpy as np
import gzip
import sparse
import sys

nSym = 8  # ATGCx2 for paternal and maternal
# nSub = 1161
nSub = 3
win = 1e5
winSiz='100K'
if win == 1e4:
    winSiz = '10K'
elif win == 1e5:
    winSiz = '100K'
wpl = 50  # Number of words per line in hg38
import bionumpy as bnp



def extract_seq(gene_table, ch, refpath, filepath):
    genome = bnp.Genome.from_file(f'{refpath}/genome.fa')
    sequence = genome.read_sequence()
    gene_table.columns = ['ensg', 'chr', 'start', 'stop']
    gene_table = gene_table.loc[gene_table.chr == ch]
    locations = bnp.LocationEntry([str(c) for c in gene_table['chr']], gene_table['start'].astype(int))
    locations = genome.get_locations(locations, has_numeric_chromosomes=True)
    #.columns = ['ensg', 'chr', 'winS', 'winE']

    #genes = bnp.Interval(['chr' + str(ch)] * len(gene_table.index),
    #gene_table['start'],
    #                    gene_table['stop'])# , gene_table['ensg'])

    windo = locations.get_windows(flank=int(win))
    gene_sequences =  sequence[windo] #.raw().view('S1'
    print(gene_sequences)
    #    .from_data_frame(gene_table)
    # data = data.sort_values(by='chr',axis=0)
    nGene = len(gene_table.index)
    # Loop through genes
    for gene_id in np.arange(nGene):
        print(gene_id)
        seqR = gene_sequences[gene_id]
        # Load reference sequence
        gene_table_i = gene_table.iloc[gene_id, 2]
        # Check if variants exist
        specific_filepath = filepath + 'variantNucleotide' + winSiz + '/' + gene_table.iloc[gene_id, 0] + '.csv'
        if os.path.exists(specific_filepath):
            f = open(specific_filepath)
            line = f.readline()
            f.close()
        else:
            line = ''
        onehot = np.zeros((nSym, len(seqR), nSub), dtype='i8')
        if line != '':
            # Load variant nucleotide
            snp = pd.read_csv(specific_filepath, header=None,
                              encoding='latin1')
            positions = snp.iloc[:, 2]-int(windo[gene_id].start)
            P_seqs = [''.join(str(snp)[2] for snp in column) for column in snp.iloc[:, 3:].to_numpy().T]
            M_seqs = [''.join(str(snp)[2] for snp in column) for column in snp.iloc[:, 3:].to_numpy().T]
            seq_size = len(seqR)
            # int(2 * win + 1)
            onehot = np.zeros((nSym, seq_size, len(P_seqs)), dtype='i8')

            for k, (P_seq, M_seq) in enumerate(zip(P_seqs, M_seqs)):
                seqP = seqR.copy()
                seqM = seqR.copy()
                seqP[positions] = P_seq
                seqM[positions] = M_seq

                #seqP[positions] = ''.join(snpP.iloc[:, k + 3].to_list())
                #seqM[positions] = ''.join(snpM.iloc[:, k + 3].to_list())

                # Convert to one-hot encoding
                onehot[:, :, k] = [seqP == 'A', seqP == 'T', seqP == 'G', seqP == 'C',
                                   seqM == 'A', seqM == 'T', seqM == 'G', seqM == 'C']

        else:
            print('no variants')

            # Loop over subjects

            for k in np.arange(nSub):
                # Convert to one-hot encoding
                onehot[:, :, k] = [seqR == 'A', seqR == 'T', seqR == 'G', seqR == 'C', seqR == 'A', seqR == 'T',
                                   seqR == 'G', seqR == 'C']

        # Convert to sparse and save
        print(os.getcwd())
        np.save(f'{gene_table.iloc[gene_id, 0]}_new.npy', onehot)

        onehot = sparse.COO.from_numpy(onehot)
        folder = filepath + '/sequence' + winSiz + '/chr' + str(ch) + '/'
        os.makedirs(folder, exist_ok=True)
        sparse.save_npz(folder + gene_table.iloc[gene_id, 0], onehot)


# Load gene names and windows
# data = pd.read_csv(filepath+'geneWin'+winSiz+'.txt',sep='\t',header=None)
if __name__ == '__main__':
    pittPath = '/bgfs/mchikina/byng/'
    filepath = pittPath + 'rosmapAD/projects/insilicoMutagenesis/extractSequence/results/'
    refpath = pittPath + 'humanRefGenome/data/hg38/'

    # Parameters
    # ch = 22 # Select chromosome
    ch = int(sys.argv[1])
    gList = sys.argv[2]


    gene_file_name = filepath + 'geneWin' + winSiz + gList + '.txt'
    data = pd.read_csv(gene_file_name, sep='\s+', header=None)
    extract_seq(data)
