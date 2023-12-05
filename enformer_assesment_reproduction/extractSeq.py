import os

import pandas as pd
import numpy as np
import gzip
import sparse
import sys

nSym = 8  # ATGCx2 for paternal and maternal
# nSub = 1161
nSub = 3
win = 1e4
winSiz='100K'
if win == 1e4:
    winSiz = '10K'
elif win == 1e5:
    winSiz = '100K'
wpl = 50  # Number of words per line in hg38


def extract_seq(gene_table, ch, refpath, filepath):
    print('running')
    gene_table.columns = ['ensg', 'chr', 'winS', 'winE']
    # data = data.sort_values(by='chr',axis=0)
    gene_table = gene_table.loc[gene_table.chr == ch]
    nGene = len(gene_table.index)
    # Loop through genes
    for gene_id in np.arange(nGene):
        print(gene_id)

        # Load reference sequence
        gene_table_i = gene_table.iloc[gene_id, 2]
        with gzip.open(refpath + 'chr' + str(ch) + '.fa.gz', 'rt') as f:
            f.readline()  # Remove non-sequence line
            # ref = f.read().replace('\n','') # Read file as a single string with \n removed but uses too much RAM
            start = np.mod(gene_table_i, wpl) - 1
            nLine = np.floor_divide(gene_table_i, wpl)
            if start == -1:  # gene j at the end of a line
                nLine -= 1
            for l in np.arange(nLine):
                f.readline()  # Get to line where gene j is located
            if start == -1:
                f.read(wpl - 1)  # Get to end of the current line
                ref = f.read(int(2 * win + 1 + np.floor_divide(2 * win + 1, wpl) + 1))  # Read 2*win+1 bases + \n's + 1
            else:
                f.read(start)  # Get to start of gene j location
                ref = f.read(int(2 * win + 1 + np.floor_divide(2 * win + 1, wpl)))  # Read 2*win+1 bases + \n's

        # Extract reference sequence within window
        # seqR = np.array(list(ref[data.iloc[j,2]-1:data.iloc[j,3]].upper()))
        seqR = np.array(list(ref.replace('\n', '').upper()))

        # Check if variants exist
        #specific_filepath = filepath + 'variantNucleotide' + winSiz + '/' + gene_table.iloc[gene_id, 0] + '.csv'
        specific_filepath = filepath + '/' + gene_table.iloc[gene_id, 0] + '.csv'
        if os.path.exists(specific_filepath):
            f = open(specific_filepath)
            line = f.readline()
            f.close()
        else:
            line = ''
        if line != '':
            # Load variant nucleotide
            snp = pd.read_csv(specific_filepath, header=None,
                              encoding='latin1')

            # Extract paternal variants
            snpP = snp.replace(['A|A', 'A|T', 'A|G', 'A|C'], 'A')
            snpP = snpP.replace(['T|A', 'T|T', 'T|G', 'T|C'], 'T')
            snpP = snpP.replace(['G|A', 'G|T', 'G|G', 'G|C'], 'G')
            snpP = snpP.replace(['C|A', 'C|T', 'C|G', 'C|C'], 'C')

            # Extract maternal variants
            snpM = snp.replace(['A|A', 'T|A', 'G|A', 'C|A'], 'A')
            snpM = snpM.replace(['A|T', 'T|T', 'G|T', 'C|T'], 'T')
            snpM = snpM.replace(['A|G', 'T|G', 'G|G', 'C|G'], 'G')
            snpM = snpM.replace(['A|C', 'T|C', 'G|C', 'C|C'], 'C')

            # Insert variants to reference sequence
            nVar = len(snp.index)

            # Loop over subjects
            onehot = np.zeros((nSym, int(2 * win + 1), nSub), dtype='i8')

            for k in np.arange(nSub):
                seqP = seqR.copy()
                seqM = seqR.copy()

                # Loop over variants
                for i in np.arange(nVar):
                    print(k, i)
                    snp_i = snpP.iloc[i, 2]
                    seqP[snp_i - gene_table_i] = snpP.iloc[i, k + 3]
                    seqM[snpM.iloc[i, 2] - gene_table_i] = snpM.iloc[i, k + 3]

                # Convert to one-hot encoding
                onehot[:, :, k] = [seqP == 'A', seqP == 'T', seqP == 'G', seqP == 'C', seqM == 'A', seqM == 'T',
                                   seqM == 'G', seqM == 'C']
        else:
            print('no variants')

            # Loop over subjects
            onehot = np.zeros((nSym, int(2 * win + 1), nSub), dtype='i8')
            for k in np.arange(nSub):
                # Convert to one-hot encoding
                onehot[:, :, k] = [seqR == 'A', seqR == 'T', seqR == 'G', seqR == 'C', seqR == 'A', seqR == 'T',
                                   seqR == 'G', seqR == 'C']

        # Convert to sparse and save
        print(os.getcwd())
        np.save(f'{gene_table.iloc[gene_id, 0]}.npy', onehot)

        onehot = sparse.COO.from_numpy(onehot)
        folder = filepath + '/sequence' + winSiz + '/chr' + str(ch) + '/'
        os.makedirs(folder, exist_ok=True)
        sparse.save_npz(folder + gene_table.iloc[gene_id, 0], onehot)
    print('finished')

# Load gene names and windows
# data = pd.read_csv(filepath+'geneWin'+winSiz+'.txt',sep='\t',header=None)
if __name__ == '__main__':
    # pittPath = '/bgfs/mchikina/byng/'
    # filepath = pittPath + 'rosmapAD/projects/insilicoMutagenesis/extractSequence/results/'
    # refpath = pittPath + 'humanRefGenome/data/hg38/'

    # Parameters
    # ch = 22 # Select chromosome
    ch = 1#  = int(sys.argv[1])
    #gList = sys.argv[2]
    #gene_file_name = filepath + 'geneWin' + winSiz + gList + '.txt'
    gene_file_name = 'enformer_assesment_reproduction/tmp.csv'
    data = pd.read_csv(gene_file_name, sep='\s+', header=None)
    extract_seq(data, ch,  refpath='/home/knut/Data/', filepath='data/debug/')
