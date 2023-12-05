import bionumpy as bnp
import pandas as pd
import sparse
from bionumpy import SequenceEntry
from bionumpy.variants import apply_variants
import numpy as np
from bionumpy.variants.consensus import apply_variants_to_sequence


# genome: bnp.Genome
    # returns: np.ndarray
def get_one_hot(genome_filename, variant_filename, annotation_filename):
    genes_iter = (chunk.get_genes() for chunk in bnp.open(annotation_filename, lazy=False).read_chunks())
    intervals = np.concatenate([bnp.datatypes.Bed6(genes.chromosome, genes.start, genes.stop, genes.gene_id, [0]*len(genes), ['+']*len(genes)) for genes in genes_iter])[::50]
    print('Reading intervals')
    # intervals = bnp.open('tmp.bed', buffer_type=bnp.io.Bed6Buffer, lazy=False).read()[::10]
    intervals.start = intervals.start - 10000
    intervals.stop = intervals.start + 10000+1
    bnp.open('tmp.bed', 'w').write(intervals)
    df = intervals.topandas()
    new_df = pd.DataFrame({'ensg': df['name'], 'chr': [s[3:] for s in df['chromosome']], 'winS': df['start'], 'winE': df['stop']})
    new_df.to_csv('tmp.csv', index=False, header=False, sep='\t')
    print('Reading genome')
    genome = bnp.Genome.from_file(genome_filename)
    # genomic_intervals= genome.get_intervals(intervals)
    # with bnp.open('tmp.vcf', 'w') as out:
    #     for i, variants in enumerate(bnp.open(variant_filename, buffer_type=bnp.io.vcf_buffers.VCFBuffer2).read_chunks(min_chunk_size=1000000)):
    #         print(f'Reading chunk {i}')
    #         locations = genomic_intervals.map_locations(bnp.replace(variants, chromosome=bnp.as_encoded_array(
    #             ['chr' + c for c in variants.chromosome.tolist()])))
    #         out.write(locations)
    print('Reading variants')
    variants = bnp.open('tmp.vcf', buffer_type=bnp.io.vcf_buffers.VCFBuffer2, lazy=False).read_chunk()
    print('Filtering variants')
    variants = variants[(variants.ref_seq.lengths == 1) & (variants.alt_seq.lengths == 1)]

#    print(variants)
#    print(intervals)
    print('Getting sequences')
    write_one_hot(genome, intervals, variants)
    # write_new_sequences(genome, intervals, variants)


def sparse_one_hot(paternal_sequences, maternal_sequences):
    if (np.any(paternal_sequences=='n') or np.any(maternal_sequences=='n')):
        print(np.sum(paternal_sequences=='n'), np.sum(maternal_sequences=='n'))
    n_seq = len(paternal_sequences)
    seq_l = len(paternal_sequences[0])
    row_indices = np.arange(n_seq).repeat(seq_l*2)
    col_indices = np.tile(np.repeat(np.arange(seq_l), 2), n_seq)
    data = np.hstack([paternal_sequences.raw().ravel().reshape(-1, 1), maternal_sequences.raw().ravel().reshape(-1, 1)]).ravel()
    return sparse.COO((data, (row_indices, col_indices)), shape=(n_seq, seq_l))


def write_one_hot(genome, intervals, all_variants):
    reference = genome.read_sequence()
    for i in range(len(intervals)):
        print(f'Running interval {i} of {len(intervals)}: {intervals[i].name}')
        interval = intervals[i:i+1]
        sequence = reference[interval][0]
        variants = all_variants[all_variants.chromosome == interval.name]
        maternal_mask = (variants.genotype == b'0|1') | (variants.genotype == b'1|1')
        paternal_mask = (variants.genotype == b'1|0') | (variants.genotype == b'1|1')
        maternal_sequences = []
        paternal_sequences = []
        for sample_id in range(10):#variants.genotype.shape[1]):
            maternal_sequences.append(apply_variants_to_sequence(sequence, variants[maternal_mask[:,i]]))
            paternal_sequences.append(apply_variants_to_sequence(sequence, variants[paternal_mask[:,i]]))
        one_hot = sparse_one_hot(bnp.as_encoded_array(paternal_sequences), bnp.as_encoded_array(maternal_sequences))
        sparse.save_npz(f'{interval.name}.npz', one_hot)


def write_new_sequences(genome, intervals, variants):
    sequences = SequenceEntry(intervals.name, genome.read_sequence()[intervals])
    for i, genotype in enumerate(variants.genotype.T):
        print(f'Running genotype {i} of {variants.genotype.shape[1]}')
        M_mask = (genotype == '1|0') | (genotype == '1|1')
        M_variants = variants[M_mask]
        bnp.open(f'tmp_M{i}.fasta', 'w', buffer_type=bnp.TwoLineFastaBuffer).write(
            apply_variants(sequences, M_variants))
        P_mask = (genotype == '0|1') | (genotype == '1|1')
        bnp.open(f'tmp_P{i}.fasta', 'w', buffer_type=bnp.TwoLineFastaBuffer).write(
            apply_variants(sequences, variants[P_mask]))


if __name__ == '__main__':
    get_one_hot('/home/knut/Data/hg38.fa', '/home/knut/Data/variants.vcf.gz', '/home/knut/Data/gencode.v43.annotation.gff3.gz')

