import os
from pathlib import PurePath

import bionumpy as bnp
import pytest

from enformer_assesment_reproduction.simulation import simulate_genome


@pytest.fixture
def genome_folder():
    genome_folder = 'data/hg38/'
    if os.path.exists('data/hg38/'):
        return PurePath(genome_folder)
    else:
        os.makedirs('data/hg38/')
        genome = simulate_genome(23, [1000000 for i in range(23)])
        for i in range(len(genome)):
            with bnp.open(f'data/hg38/chr{i}.fa', 'w', buffer_type=bnp.MultiLineFastaBuffer) as f:
                f.write(genome[i:i+1])
        with bnp.open(f'data/hg38/genome.fa', 'w', buffer_type=bnp.MultiLineFastaBuffer) as f:
            f.write(genome)
    return PurePath(genome_folder)


