import os
from pathlib import PurePath

import bionumpy as bnp
import pandas as pd
import pytest
from bionumpy.simulate.variants import simulate_variants
from bionumpy.io.vcf_buffers import VCFBuffer2
from enformer_assesment_reproduction.simulation import simulate_genome

gene_list = ['ENSG00000186092', 'ENSG00000237683', 'ENSG00000235249', 'ENSG00000185097', 'ENSG00000269831']
@pytest.fixture
def genome_folder():
    genome_folder = 'data/hg38/'
    bnp.MultiLineFastaBuffer.n_characters_per_line = 50
    if os.path.exists('data/hg38/'):
        return genome_folder
    else:
        os.makedirs('data/hg38/')
        genome = simulate_genome(23, [1000000 for i in range(23)])
        for i in range(len(genome)):
            with bnp.open(f'data/hg38/chr{i}.fa.gz', 'w', buffer_type=bnp.MultiLineFastaBuffer) as f:
                f.write(genome[i:i+1])
        with bnp.open(f'data/hg38/genome.fa', 'w', buffer_type=bnp.MultiLineFastaBuffer) as f:
            f.write(genome)
    return genome_folder

@pytest.fixture
def gene_list_path():
    path = 'data/to_process.txt'
    if os.path.exists(path):
        return path
    else:
        with open(path, 'w') as f:
            f.write('\t'.join(gene_list))
        return path

@pytest.fixture
def gene_table_path():
    path = 'data/geneWin100K.txt'
    if os.path.exists(path):
        return path
    else:
        table= pd.DataFrame({'gene_id': gene_list,
                             'chr': list(range(len(gene_list))),
                             'start': [10000*i for i in range(len(gene_list))],
                             'end': [10000*i+1000 for i in range(len(gene_list))]})
        table.to_csv(path, index=False, header=False, sep='\t')
        return path


@pytest.fixture
def variants_path(genome_folder):
    path = 'data/chrAll.phased.vcf.gz'
    if os.path.exists(path):
        return path
    else:
        genome = bnp.Genome.from_file(f'{genome_folder}/genome.fa').read_sequence()
        variants = simulate_variants(genome, 0.001, small_indel_prob=0, sv_prob=0, n_samples=3)
        with bnp.open(path, 'w', buffer_type=VCFBuffer2) as f:
            for v in variants:
                f.write(v)
        return path


