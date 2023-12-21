from typing import Iterable, List

import bionumpy as bnp
import numpy as np
from bionumpy import SequenceEntry
from bionumpy.variants import apply_variants

from enformer_assesment_reproduction.full_pipeline import get_tss_windows


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


def _map_vcf(genomic_intervals, variant_filename):
    with bnp.open('tmp.vcf', 'w') as out:
        for i, variants in enumerate(
            bnp.open(variant_filename, buffer_type=bnp.io.vcf_buffers.VCFBuffer2).read_chunks(min_chunk_size=1000000)):
            variants = variants[(variants.ref_seq.lengths == 1) & (variants.alt_seq.lengths == 1)]
            locations = genomic_intervals.map_locations(bnp.replace(variants, chromosome=bnp.as_encoded_array(
                ['chr' + c for c in variants.chromosome.tolist()])))
            out.write(locations)


def _get_tss_windows(chunks: Iterable[bnp.datatypes.GTFEntry], gene_names: List[str], genome: bnp.Genome) -> bnp.Bed6:
    genes_iter = (chunk.get_genes() for chunk in chunks)
    genes_iter = (genes[np.isin(genes.gene_id, gene_names)] for genes in genes_iter)
    intervals = np.concatenate([bnp.datatypes.Bed6(genes.chromosome, genes.start, genes.stop, genes.gene_id,
                                                   np.zeros(len(genes), dtype=int), genes.strand) for genes in
                                genes_iter])
    flank = 10000
    intervals.stop = intervals.start + flank + 1
    intervals.start = intervals.start - flank
    genomic_intervals = genome.get_intervals(intervals)
    genomic_intervals = genomic_intervals.clip()
    intervals = genomic_intervals.data
    intervals = intervals[intervals.stop - intervals.start == 2 * flank + 1]
    return intervals


def write_gene_windows(annotation_filename, genome, gene_list_name, out_filename='tmp.bed'):
    gene_names = open(gene_list_name).read().strip().split()
    chunks = bnp.open(annotation_filename, lazy=False).read_chunks()
    intervals = get_tss_windows(chunks, gene_names, genome)
    bnp.open(out_filename, 'w').write(intervals)
    return intervals
