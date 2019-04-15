# snpiphy
## An automated snp phylogeny pipeline

snpiphy is a workflow pipeline for the automated generation of SNP derived phylogenetic trees. This pipeline quickly generates a SNP-derived tree using a single reference sequence and established tools. Users need to provide a high-quality reference sequence and a directory containing single-end (with the -s option) or interleaved fastq. Reads can be provided as uncompressed `.fastq` or gzip-compressed `.fastq.gz`. Currently, split file paired reads or reads from multiple runs are not supported. Assemblies in `fasta` format can also be provided but `fastq` reads are preferred.

###  Requirements

snpiphy requires a linux-based operating system and the following:

- python (>=3.5)
- biopython (>=1.72)
- pandas (>=0.23.4)
- numpy (>=1.16.0)
- [snippy](https://github.com/tseemann/snippy) (>=4.3.6)
- [gubbins](https://github.com/sanger-pathogens/gubbins) (>=2.3.4)

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
usage: snpiphy [-h] -o OUTDIR -r REFERENCE [-q READS_DIR] [-l READS_LIST]
               [-c CUTOFF] [-p PREFIX] [-t THREADS] [-j] [-s] [-m]
               [-b {raxml,fasttree}] [-f] [-n] [--version] [-v]

snpiphy - An automated snp phylogeny pipeline.

required arguments:
  -o OUTDIR, --output-directory OUTDIR
                        Path to the output directory. A directory will be
                        created if one does not exist.
  -r REFERENCE, --reference REFERENCE
                        Path to the reference sequence. Only fasta format is
                        accepted currently.

input arguments - required, provide only one:
  -q READS_DIR, --input_dir READS_DIR
                        Path to a directory with your interleaved fastq
                        sequencing reads, single-end sequencing reads (use
                        with the -s option) or fasta assemblies.
  -l READS_LIST, --input_list READS_LIST
                        Path to a tab-separated file containing read paths and
                        paired status. Best used when running a combination of
                        single-end and paired-end data or if your data is
                        spread across multiple directories. Run
                        'snpiphy_generate_input_list' to generate an example
                        list.

optional arguments:
  -h, --help            show this help message and exit
  -c CUTOFF, --cov_cutoff CUTOFF
                        Percent coverage of reference sequence (0-100%) used
                        to reject a sample. Samples lower than this threshold
                        will be excluded from phylogenetic pipeline steps.
                        Defaults to 85%.
  -p PREFIX, --prefix PREFIX
                        Prefix for output files
  -t THREADS, --threads THREADS
                        Number of threads to use for multiprocessing.
  -j, --parallel        Use GNU parallel to run multiple instances of snippy
                        (can speed things up if you have multiple cores
                        available)
  -s, --single_end      fastq reads are single end instead of paired-end. Use
                        for ion torrent or non-paired end illumina data
  -m, --gamma_model     Use GTRGAMMA model instead of GTRCAT during the
                        gubbins and RAxML tree building steps.
  -b {raxml,fasttree}, --tree_builder {raxml,fasttree}
                        Algorithm used for building the phylogenetic tree
                        (default: raxml)
  -f, --force           Overwrite files in the output directories.
  -n, --no_recombination_filter
                        Don't filter potential recombination events. Use for
                        organisms that are known undergo low amounts of
                        recombination.
  --version             show program's version number and exit
  -v, --verbose         Increase verbosity on command line output (n.b.
                        verbose output is always saved to snpiphy.log in the
                        output directory).
```
