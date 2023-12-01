import bionumpy as bnp
from bionumpy.simulate import simulate_sequence, simulate_sequences


def simulate_genome(n_chromosomes, lens):
    return simulate_sequences("ACGT", {f"chr{i}": length for i, length in enumerate(lens)})

