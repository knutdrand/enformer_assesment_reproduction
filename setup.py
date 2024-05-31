#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['bionumpy', 'sparse', 'numpy', 'typer']

test_requirements = ['pytest>=3', "hypothesis"]

setup(
    author="Knut Rand",
    author_email='knutdrand@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Reproduce genomic processing using bionumpy",
    entry_points={
        'console_scripts': [
            'enformer_assesment_reproduction=enformer_assesment_reproduction.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='enformer_assesment_reproduction',
    name='enformer_assesment_reproduction',
    packages=find_packages(include=['enformer_assesment_reproduction', 'enformer_assesment_reproduction.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/knutdrand/enformer_assesment_reproduction',
    version='0.0.1',
    zip_safe=False,
)
