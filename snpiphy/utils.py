#!/usr/bin/env python

import sys
import os
import gzip
import subprocess
import shutil
from contextlib import contextmanager
import logging
logger = logging.getLogger(__name__)


def configure_outdir(outdir, force):
	if os.path.isdir(outdir):
		if force == True:
			logger.debug("Force argument provided overwriting previous output directory...")
			shutil.rmtree(outdir)
			os.makedirs(outdir)
		else:
			logger.debug("Output directory present, continuing run with previous output directory...")
	else:
		logger.debug("Creating output directory...")
		os.makedirs(outdir)


def get_samplename(filename):
	basename = os.path.basename(filename)
	name, ext = os.path.splitext(basename)
	if ext == ".gz":
		name, ext = os.path.splitext(name)
	if ext == ".fastq" or ext == ".fq" or ext == ".fasta" or ext == ".fa" or ext == ".fna" or ext == ".fsa":
		name2, ext2 = os.path.splitext(name)
		if ext2 == ".lf":
			name = name2
		return name
	else:
		logger.error("Failed to determine samplename for file: {}. Check if this is a fastq|fastq.gz or fasta|fasta.gz file.".format(basename))
		sys.exit(1)


def id_read_format(file):
	logger.debug("Identifying format for file: {}".format(os.path.basename(file)))
	if os.path.splitext(file)[1] == ".gz":
		logger.debug("{} is gzipped".format(os.path.basename(file)))
		handle = gzip.open(file,'rt')
	else:
		logger.debug("{} is not gzipped".format(os.path.basename(file)))
		handle = open(file)
	for line in handle:
		if line.strip() == "":
			continue
		if line.strip().startswith('@'):
			logger.debug("{} is fastq".format(os.path.basename(file)))
			return "fastq"
		elif line.strip().startswith(">"):
			logger.debug("{} is fasta".format(os.path.basename(file)))
			return "fasta"
		else:
			logger.error("Can't properly determine read format for file {}. Please check your read files are in fastq for raw sequencing reads) or fasta (for assembled genomes) format".format(os.path.basename(file)))
			sys.exit(1)


def log_subprocess_output(pipe):
	for line in iter(pipe.readline, b''): # b'\n'-separated lines
		logger.debug(line.strip().decode('utf-8'))


def define_make_outdirs(dir):
	if os.path.isdir(dir) == False:
		os.makedirs(dir)
	return dir


def run_command(cmd):
	logger.info("Executing command: {}\n".format(" ".join(cmd)))
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	with process.stdout:
		log_subprocess_output(process.stdout)
	exitcode = process.wait()
	return exitcode


@contextmanager
def cd(newdir):
	prevdir = os.getcwd()
	os.chdir(os.path.expanduser(newdir))
	try:
		yield
	finally:
		os.chdir(prevdir)


def find_source_file(name, reads_dir):
	reads_files = os.listdir(reads_dir)
	for file in reads_files:
		searchObj = re.search(name,file)
		if searchObj:
			return os.path.join(reads_dir, file)


