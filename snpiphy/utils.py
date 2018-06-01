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



# def chunks(l, n):
	# """Yield successive n-sized chunks from l."""
	# for i in range(0, len(l), n):
		# yield l[:, i:i + n]


# def id_read_format(file):
	# if os.path.splitext(file)[1] == ".gz":
		# handle = gzip.open(file,'rt')
	# else:
		# handle = open(file)
	# for line in handle:
		# if line.strip() == "":
			# continue
		# if line.strip().startswith('@'):
			# return "fastq"
		# elif line.strip().startswith(">"):
			# return "fasta"
		# else:
			# print("Can't properly determine read format for file {}. Please check your read files are in fastq for raw sequencing reads) or fasta (for assembled genomes) format".format(os.path.basename(file)))
			# sys.exit(1)


# def find_source_file(name, reads_dir):
	# reads_files = os.listdir(reads_dir)
	# for file in reads_files:
		# searchObj = re.search(name,file)
		# if searchObj:
			# return os.path.join(reads_dir, file)







# def multijob(input):
	# query_id, subject_id, query, subject = input
	# snps = int(np.count_nonzero(query != subject))
	# print("{} SNPs found between isolates {} and {}".format(snps, query_id, subject_id))
	# return (query_id, subject_id, snps)


# def multirun(input_list, threads):
	# p = multiprocessing.Pool(threads)
	# result = p.map_async(multijob,input_list,chunksize=1)
	# prev_jobs_remaining = len(input_list)
	# while not result.ready():
		# jobs_remaining = result._number_left
		# if jobs_remaining != prev_jobs_remaining:
			# print("{} of {} sequence pairs remaining to be compaired".format(result._number_left, len(input_list)))
		# prev_jobs_remaining = jobs_remaining
	# results_list = result.get(999999999999999)
	# p.close()
	# p.join()
	# return results_list

# def gen_input_list(array, id_list):
	# input_list = []
	# for i, query in enumerate(array):
		# for j, subject in enumerate(array):
			# input_list.append((id_list[i], id_list[j], query, subject))
	# return input_list



# def generate_data_frame(snp_count_d,id_list):
	# df = pd.DataFrame(data=np.zeros((0,len(id_list))), columns=id_list)
	# for query_id in sorted(id_list):
		# print("Adding %s results to dataframe" % query_id)
		# df_1 = pd.DataFrame(snp_count_d[query_id], index=[query_id])
		# df = df.append(df_1)
	# return df


# def validate_multi_import(array, array2):
	# for i in range(len(array)):
		# if multijob((i ,i, array[i], array2[i]))[2] != 0:
			# print("Check multiprocessing array input is working")
			# sys.exit(1)
	# print("Multiprocessing array input is working fine")


# def import_array(aln, threads):
	# chunk_list = gen_chunks(aln, threads)
	# array_chunk_list = multi_import(chunk_list, threads)
	# return reassemble_array(array_chunk_list)







