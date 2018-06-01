import multiprocessing, logging, operator
import numpy as np
import pandas as pd
from Bio import AlignIO
from collections import defaultdict
logger = logging.getLogger(__name__)

def count_snps(alignment_file, outfile, statsfile, threads):
	logger.info("Importing alignment...")
	aln = AlignIO.read(alignment_file,"fasta")
	logger.info("Identifying core genome...")
	array = import_array(threads, aln)
	id_list = [x.id for x in aln]
	input_list = [ (id_list[i], id_list[j], query, subject) for i, query in enumerate(array) for j, subject in enumerate(array) ]
	results = multirun(threads, input_list)
	logger.info("Writing stats and results...")
	df = generate_data_frame(results, id_list)
	stats = df.describe()
	stats.to_csv(statsfile, sep="\t")
	df.to_csv(outfile, sep="\t")
	logger.info("Done.")

def import_array(threads, aln):
	chunk_num = threads
	chunk_size = int(len(aln[0])/chunk_num)+1
	logger.debug("Separating alignment into {}, {} bp sequence chunks for multiprocessing...".format(chunk_num, chunk_size))
	chunk_list = [(i, aln[:,i:i + chunk_size]) for i in range(0, len(aln[0]), chunk_size)]
	p = multiprocessing.Pool(threads)
	array_chunk_list = p.map_async(import_array_in_chunks,chunk_list,chunksize=1).get()
	p.close()
	p.join()
	logger.debug("Sorting processed sequence chunks into original order.")
	l = [array_chunk for i, array_chunk in sorted(array_chunk_list,key=operator.itemgetter(0))]
	logger.debug("Reassembling sequence chunks into alignment.")
	return np.concatenate(l , axis=1)

def import_array_in_chunks(indexed_chunk):
	index, aln = indexed_chunk
	chunk = np.array([list(rec) for rec in aln], np.character)
	indices = [i for i, v in enumerate(chunk[0]) if (b'N' in chunk[:,i]) or (b'-' in chunk[:,i])]
	chunk = np.delete(chunk, indices, axis=1)
	return index, chunk

def generate_data_frame(results, id_list):
	snp_count_d = defaultdict(dict)
	for query_id, subject_id, sum in results:
		snp_count_d[query_id][subject_id] = sum
	df = pd.DataFrame(data=np.zeros((0,len(id_list))), columns=id_list)
	for query_id in sorted(id_list):
		logger.debug("Adding %s results to dataframe" % query_id)
		df_1 = pd.DataFrame(snp_count_d[query_id], index=[query_id])
		df = df.append(df_1)
	return df

def multirun(threads, input_list):
	p = multiprocessing.Pool(threads)
	result = p.map_async(multijob,input_list,chunksize=1)
	prev_jobs_remaining = len(input_list)
	while not result.ready():
		jobs_remaining = result._number_left
		if prev_jobs_remaining - jobs_remaining == 100:
			logger.debug("{} of {} sequence pairs remaining to be compaired".format(result._number_left, len(input_list)))
			prev_jobs_remaining = jobs_remaining
	results_list = result.get(999999999999999)
	p.close()
	p.join()
	return results_list

def multijob(input):
	query_id, subject_id, query, subject = input
	snps = int(np.count_nonzero(query != subject))
	logger.debug("{} SNPs found between isolates {} and {}".format(snps, query_id, subject_id))
	return (query_id, subject_id, snps)

