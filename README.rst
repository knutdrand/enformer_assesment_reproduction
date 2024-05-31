===============================
Enformer Assesment Reproduction
===============================




.. image:: https://readthedocs.org/projects/enformer-assesment-reproduction/badge/?version=latest
        :target: https://enformer-assesment-reproduction.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




Reproduce genomic processing using BioNumPy. For a set of genes and a set of variants this finds the maternal and paternal sequences in a window around the tss of the genes,
Writes the one-hot-encoded sequences in the same format as EnformerAsssesment.

.. code-block:: sh

    $ git clone git@github.com:knutdrand/enformer_assesment_reproduction.git
    $ cd enformer_assesment_reproduction
    $ pip install .
    $ cd Data
    $ chmod +x get_data.sh
    $ enformer_assesment_reproduction hg38.fa variants_cut.vcf.gz genocde.v45.gff.gz  gene_list.txt output/ --flank 10000

All code can be found here: `Code <enformer_assesment_reproduction/cli.py>`_
