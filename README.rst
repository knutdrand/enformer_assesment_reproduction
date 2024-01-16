===============================
Enformer Assesment Reproduction
===============================


.. image:: https://img.shields.io/pypi/v/enformer_assesment_reproduction.svg
        :target: https://pypi.python.org/pypi/enformer_assesment_reproduction

.. image:: https://img.shields.io/travis/knutdrand/enformer_assesment_reproduction.svg
        :target: https://travis-ci.com/knutdrand/enformer_assesment_reproduction

.. image:: https://readthedocs.org/projects/enformer-assesment-reproduction/badge/?version=latest
        :target: https://enformer-assesment-reproduction.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




Reproduce genomic processing using BioNumPy. For a set of genes and a set of variants this finds the maternal and paternal sequences in a window around the tss of the genes,
Writes the one-hot-encoded sequences in the same format as EnformerAsssesment.

.. code-block:: sh

    $ pip install enformer_assesment_reproduction
    $ enformer_assesment_reproduction ~/Data/hg38.fa ~/Data/variants_cut.vcf.gz ~/Data/gencode.v43.annotation.gff3.gz gene_list.txt output/ --flank 10000

All code can be found here: `Code <enformer_assesment_reproduction/cli.py>`_
