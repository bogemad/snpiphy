# snpiphy
## An automated snp phylogeny pipeline

snpiphy is a workflow pipeline for the automated generation of SNP derived phylogenetic trees. This pipeline quickly generates a SNP-derived tree using a single reference sequence and established tools. Users need to provide a high-quality reference sequence and a directory containing interleaved fastq (or gzip-compressed fastq.gz) reads. Currently, split file, multi-file or single reads are not supported but these may be added in future versions. Assemblies in fasta format can also be provided but fastq reads are preferred.

###  Requirements

snpiphy requires a linux-based operating system and the following:

- python (>=3.5)
- biopython (>=1.71)
- pandas (>=0.23.0)
- numpy (>=1.14.3)
- [snippy](https://github.com/tseemann/snippy) (working with 4.0-dev)
- [gubbins](https://github.com/sanger-pathogens/gubbins) (>=2.3.1)

### Installation

To download and install snpiphy with all dependencies use conda:
```
conda install -c conda-forge -c bioconda snpiphy
```

or pip (requires manual dependency installation): 
```
pip install git+https://github.com/bogemad/snpiphy.git
```

or a manual install (again this requires you to manually install dependencies):
```
git clone https://github.com/bogemad/snpiphy.git
cd snpiphy
python setup.py install
```

### Usage

```
usage: snpiphy [-h] -q READS_DIR -o OUTDIR -r REFERENCE [-c CUTOFF]
               [-p PREFIX] [-t THREADS] [-f] [--version] [-v]

snpiphy - An automated snp phylogeny pipeline.

required arguments:
  -q READS_DIR, --fastq-dir READS_DIR
                        Path to a directory with your interleaved fastq
                        sequencing reads or fasta assemblies.
  -o OUTDIR, --output-directory OUTDIR
                        Path to the output directory. A directory will be
                        created if one does not exist.
  -r REFERENCE, --reference REFERENCE
                        Path to the reference sequence. Only fasta format is
                        accepted currently.

optional arguments:
  -h, --help            show this help message and exit
  -c CUTOFF, --cov-cutoff CUTOFF
                        Percent coverage of reference sequence (0-100%) used
                        to reject a sample. Samples lower than this threshold
                        will be excluded from phylogenetic pipeline steps.
                        Defaults to 85%.
  -p PREFIX, --prefix PREFIX
                        Prefix for output files
  -t THREADS, --threads THREADS
                        Number of threads to use for multiprocessing.
  -f, --force           Overwrite files in the output directories.
  --version             show program's version number and exit
  -v, --verbose         Increase verbosity on command line output (n.b.
                        verbose output is always saved to snpiphy.log in the
                        output directory).
```
