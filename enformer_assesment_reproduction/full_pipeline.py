import bionumpy as bnp
from bionumpy.arithmetics.intervals import clip
import pandas as pd
import sparse
from bionumpy import SequenceEntry
from bionumpy.variants import apply_variants
import numpy as np
from bionumpy.variants.consensus import apply_variants_to_sequence


def get_one_hot(genome_filename, variant_filename, annotation_filename, gene_list_name='gene_list.txt'):
    genome = bnp.Genome.from_file(genome_filename)
    intervals = write_gene_windows(annotation_filename, genome, gene_list_name)
    # intervals = bnp.open('tmp.bed', buffer_type=bnp.io.Bed6Buffer, lazy=False).read()
    map_vcf(genome.get_intervals(intervals), variant_filename)
    variants = bnp.open('tmp.vcf', buffer_type=bnp.io.vcf_buffers.VCFBuffer2).read_chunk()
    write_one_hot(genome, intervals, variants)


def get_genes(annotation_file_name, gene_list_name):
    gene_names = open(gene_list_name).read().strip().split()
    for chunk in bnp.open(annotation_file_name).read_chunks():
        genes = chunk.get_genes()
        mask = np.isin(genes.gene_id, gene_names)


def write_gene_windows(annotation_filename, genome, gene_list_name, out_filename='tmp.bed'):
    gene_names = open(gene_list_name).read().strip().split()
    genes_iter = (chunk.get_genes() for chunk in bnp.open(annotation_filename, lazy=False).read_chunks())
    genes_iter = (genes[np.isin(genes.gene_id, gene_names)] for genes in genes_iter)
    intervals = np.concatenate([bnp.datatypes.Bed6(genes.chromosome, genes.start, genes.stop, genes.gene_id,
                                                   np.zeros(len(genes), dtype=int), genes.strand) for genes in genes_iter])
    flank = 10000
    intervals.stop = intervals.start + flank + 1
    intervals.start = intervals.start - flank
    bnp.open(out_filename, 'w').write(intervals)
    genomic_intervals = genome.get_intervals(intervals)
    genomic_intervals = genomic_intervals.clip()
    intervals = genomic_intervals.data
    intervals = intervals[intervals.stop - intervals.start == 2 * flank + 1]
    bnp.open(out_filename, 'w').write(intervals)
    return intervals


def convert_bed_to_csv(bed_filename, csv_filename='tmp.csv'):
    intervals = bnp.read(bed_filename, buffer_type=bnp.io.Bed6Buffer)
    df = intervals.topandas()
    new_df = pd.DataFrame(
        {'ensg': df['name'],
         'chr': [s[3:] for s in df['chromosome']],
         'winS': df['start'],
         'winE': df['stop']})
    new_df.to_csv(csv_filename, index=False, header=False, sep='\t')


def map_vcf(genomic_intervals, variant_filename):
    with bnp.open('tmp.vcf', 'w') as out:
        for i, variants in enumerate(
            bnp.open(variant_filename, buffer_type=bnp.io.vcf_buffers.VCFBuffer2).read_chunks(min_chunk_size=1000000)):
            variants = variants[(variants.ref_seq.lengths == 1) & (variants.alt_seq.lengths == 1)]
            locations = genomic_intervals.map_locations(bnp.replace(variants, chromosome=bnp.as_encoded_array(
                ['chr' + c for c in variants.chromosome.tolist()])))
            out.write(locations)


def sparse_one_hot(paternal_sequences, maternal_sequences):
    if (np.any(paternal_sequences=='n') or np.any(maternal_sequences=='n')):
        print(np.sum(paternal_sequences=='n'), np.sum(maternal_sequences=='n'))
    n_seq = len(paternal_sequences)
    seq_l = len(paternal_sequences[0])
    subject_indices = np.arange(n_seq).repeat(seq_l*2)
    position_index = np.tile(np.repeat(np.arange(seq_l), 2), n_seq)
    one_hot_indices = np.hstack([paternal_sequences.raw().ravel().reshape(-1, 1),
                                 maternal_sequences.raw().ravel().reshape(-1, 1)+4]).ravel()
    return sparse.COO(data=np.full(subject_indices.size, 1), coords=(one_hot_indices, position_index, subject_indices))


def write_one_hot(genome, intervals, all_variants):
    reference = genome.read_sequence()
    for i in range(len(intervals)):
        interval = intervals[i:i+1]
        sequence = reference[interval][0]
        if np.any(sequence=='n'):
            print(f'Skipping interval {i} of {len(intervals)}: {intervals[i].name}')
            continue
        variants = all_variants[all_variants.chromosome == interval.name]
        if (not len(variants)) or np.all(variants.genotype[:, 3] == b'0|0'):
            paternal_sequences = [sequence for _ in range(3)]
            maternal_sequences = [sequence for _ in range(3)]
        else:
            maternal_mask = (variants.genotype == b'0|1') | (variants.genotype == b'1|1')
            paternal_mask = (variants.genotype == b'1|0') | (variants.genotype == b'1|1')
            maternal_sequences = []
            paternal_sequences = []
            for sample_id in range(3):
                maternal_sequences.append(apply_variants_to_sequence(sequence, variants[maternal_mask[:, i]]))
                paternal_sequences.append(apply_variants_to_sequence(sequence, variants[paternal_mask[:, i]]))
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
    get_one_hot('/home/knut/Data/hg38.fa',
                '/home/knut/Data/variants.vcf.gz',
                '/home/knut/Data/gencode.v43.annotation.gff3.gz',
                'gene_list.txt')

