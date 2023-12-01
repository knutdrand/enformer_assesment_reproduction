#!/usr/bin/env python
from pathlib import PurePath

import bionumpy as bnp
from .fixtures import genome_folder


def test_extract_sequences(genome_folder: PurePath):
    genome = bnp.open_indexed(str(genome_folder/'genome.fa'))


