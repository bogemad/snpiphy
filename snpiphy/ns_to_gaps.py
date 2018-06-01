import multiprocessing, logging, os, operator
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import numpy as np
logger = logging.getLogger(__name__)

def replace_Ns_with_gaps(outdir, threads):
	logger.info("Importing core alignment...")
	aln = AlignIO.read(os.path.join(outdir, "core.full.aln"),"fasta")
	logger.info("Filtering Ns...")
	array = filter_gaps_Ns(aln, int(threads))
	logger.info("Finished cleaning. Exporting to file...")
	AlignIO.write(reconstruct_alignment(aln, array), os.path.join(outdir, "core.full.trimmed.aln"), "fasta")
	logger.info("Done.")

def filter_gaps_Ns(aln, threads):
	chunk_num = threads
	chunk_size = int(len(aln[0])/chunk_num)+1
	logger.debug("Separating alignment into {}, {} bp sequence chunks for multiprocessing...".format(chunk_num, chunk_size))
	chunk_list = [(i, aln[:,i:i + chunk_size]) for i in range(0, len(aln[0]), chunk_size)]
	logger.debug("Removing Ns with {} threads...".format(threads))
	p = multiprocessing.Pool(threads)
	array_chunk_list = p.map_async(import_edit_array_in_chunks,chunk_list,chunksize=1).get()
	p.close()
	p.join()
	logger.debug("Ns have been converted to gaps. Sorting sequence chunks into original order.")
	l = [array_chunk for i, array_chunk in sorted(array_chunk_list,key=operator.itemgetter(0))]
	logger.debug("Reassembling sequence chunks into alignment.")
	return np.concatenate(l , axis=1)

def reconstruct_alignment(aln, array):
	for i, rec in enumerate(aln):
		sequence_bytes = b''.join(array[i])
		rec.seq = Seq(str(sequence_bytes,'utf-8'),IUPAC.unambiguous_dna)
	return aln

def import_edit_array_in_chunks(indexed_chunk):
	index, aln = indexed_chunk
	chunk = np.array([list(rec) for rec in aln], np.character)
	chunk[chunk == b'N'] = b'-'
	return index, chunk