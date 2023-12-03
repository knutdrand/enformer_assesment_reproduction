#!/usr/bin/env python
from pathlib import PurePath

import bionumpy as bnp
import pandas as pd

from enformer_assesment_reproduction.save_ref_seqs import save_ref_seqs, save_ref_seqs_new
from .fixtures import *
from enformer_assesment_reproduction.extractSeq import extract_seq
from enformer_assesment_reproduction.extract_seq_new import extract_seq as extract_seq_new


def test_save_snp_info(variants_path, gene_list_path, gene_table_path):
    pass

def test_extract_sequences(genome_folder: PurePath, gene_table_path: str):
    data = pd.read_csv(gene_table_path, sep='\s+', header=None)
    for i in range(23):
        os.makedirs('data/sequences', exist_ok=True)
        extract_seq(data, i, genome_folder, 'data/sequences')

def test_extract_sequences_new(genome_folder: PurePath, gene_table_path: str):
    data = pd.read_csv(gene_table_path, sep='\s+', header=None)
    for i in range(23):
        os.makedirs('data/sequences', exist_ok=True)
        extract_seq_new(data, i, genome_folder, 'data/sequences')

def test_save_ref_seqs(genome_folder: PurePath, gene_table_path: str):
    os.makedirs('data/ref_sequences', exist_ok=True)
    save_ref_seqs(gene_table_path, genome_folder, 'data/ref_sequences')

def test_save_ref_seqs_new(genome_folder: PurePath, gene_table_path: str):
    os.makedirs('data/ref_sequences', exist_ok=True)
    save_ref_seqs_new(gene_table_path, genome_folder, 'data/ref_sequences')


