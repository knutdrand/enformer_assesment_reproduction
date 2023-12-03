#!/usr/bin/env python
from pathlib import PurePath

import bionumpy as bnp
import pandas as pd

from .fixtures import *
from enformer_assesment_reproduction.extractSeq import extract_seq


def test_save_snp_info(variants_path, gene_list_path, gene_table_path):
    pass

def test_extract_sequences(genome_folder: PurePath, gene_table_path: str):
    data = pd.read_csv(gene_table_path, sep='\s+', header=None)
    for i in range(23):
        os.makedirs('data/sequences', exist_ok=True)
        extract_seq(data, i, genome_folder, 'data/sequences')



