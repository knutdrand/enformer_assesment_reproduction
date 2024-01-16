import os
from typing import Iterable, List
import bionumpy as bnp
import sparse
import numpy as np
import typer


def get_one_hot(genome_filename: str, variant_filename: str, annotation_filename: str, gene_list_name: str, out_folder: str, flank: int = 10000):
    genome = bnp.Genome.from_file(genome_filename)
    gene_names = open(gene_list_name).read().strip().split()
    intervals = get_tss_windows(bnp.open(annotation_filename).read_chunks(), gene_names, flank=flank)
    variants = map_vcf(genome.get_intervals(intervals),
                       bnp.open(variant_filename, buffer_type=bnp.io.vcf_buffers.VCFBuffer2).read_chunks(
                           min_chunk_size=1000000))
    write_one_hot(genome, intervals, variants, out_folder)


def get_tss_windows(gtf_entries: Iterable[bnp.datatypes.GTFEntry], gene_names: List[str], flank=10000) -> bnp.Bed6:
    all_transcripts = (entries.get_transcripts() for entries in gtf_entries)
    transcripts = np.concatenate([t[np.isin(t.gene_id, gene_names)] for t in all_transcripts])
    entry_tuples = []
    for gene_id, group in bnp.groupby(transcripts, 'gene_id'):
        tss = np.min(group.start) if group.strand[0] == '+' else np.max(group.stop) - 1
        entry_tuples.append((group.chromosome[0].to_string(), tss - flank, tss + flank + 1,
                             gene_id, 0, group.strand[0]))
    return bnp.datatypes.Bed6.from_entry_tuples(entry_tuples)


def map_vcf(genomic_intervals: bnp.GenomicIntervals, variants_iter: Iterable[bnp.datatypes.VCFEntry]) -> bnp.datatypes.VCFEntry:
    new_variants = []
    for i, variants in enumerate(variants_iter):
        variants = variants[(variants.ref_seq.lengths == 1) & (variants.alt_seq.lengths == 1)]
        locations = genomic_intervals.map_locations(bnp.replace(variants, chromosome=bnp.as_encoded_array(
            ['chr' + c for c in variants.chromosome.tolist()])))
        new_variants.append(locations)
    return np.concatenate(new_variants)


def sparse_one_hot(paternal_sequences: bnp.EncodedRaggedArray, maternal_sequences: bnp.EncodedRaggedArray) -> sparse.COO:
    n_seq = len(paternal_sequences)
    seq_l = len(paternal_sequences[0])
    subject_indices = np.arange(n_seq).repeat(seq_l * 2)
    position_index = np.tile(np.repeat(np.arange(seq_l), 2), n_seq)
    one_hot_indices = np.hstack([paternal_sequences.raw().ravel().reshape(-1, 1),
                                 maternal_sequences.raw().ravel().reshape(-1, 1) + 4]).ravel()
    return sparse.COO(data=np.full(subject_indices.size, 1), coords=(one_hot_indices, position_index, subject_indices))


def write_one_hot(genome: bnp.Genome, intervals: bnp.Bed6, all_variants: bnp.datatypes.VCFEntry, out_folder: str = ''):
    os.makedirs(out_folder, exist_ok=True)
    reference = genome.read_sequence()
    for gene_idx in range(len(intervals)):
        interval = intervals[gene_idx:gene_idx + 1]
        sequence = reference[interval][0]
        variants = all_variants[all_variants.chromosome == interval.name]
        masks = [np.isin(variants.genotype, genotypes)
                           for genotypes in (['0|1', '1|1'], ['1|0', '1|1'])]
        sequences = ([bnp.variants.apply_variants_to_sequence(sequence, variants[m]) for m in mask.T] for mask in masks)
        one_hot = sparse_one_hot(*(bnp.as_encoded_array(seq) for seq in sequences))
        sparse.save_npz(f'{out_folder}/{interval.name}.npz', one_hot)


def main():
    typer.run(get_one_hot)
