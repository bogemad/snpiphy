import os, multiprocessing, logging
from Bio import SeqIO
logger = logging.getLogger(__name__)

def remove_degenerates(reference_file, outfile, threads):
	ref_type = check_ref(reference_file)
	outdir = os.path.dirname(outfile)
	outfile = outfile + ref_type
	logger.debug("Checking for degenerate bases and replacing with N...")
	if os.path.isfile(outfile) == False:
		sequences, annotations_present = parse_sequence_file(reference_file)
		cleaned_sequences = []
		for sequence in sequences:
			logger.debug("Processing reference contig {}, length = {} bp".format(sequence.id, len(sequence)))
			chunk_list = gen_ref_chunks(sequence, threads)
			logger.debug("Identifying degenerate positions...")
			p = multiprocessing.Pool(threads)
			result = p.map_async(filter_degenerates, chunk_list, chunksize=1)
			results_list = result.get()
			p.close()
			p.join()
			logger.debug("Filtering degenerate positions...")
			cleaned_sequences.append(clean_sequence(results_list, sequence))
			logger.debug("Done.")
		write_new_sequence_file(cleaned_sequences, annotations_present, outfile)
	else:
		logger.debug("Reference sequence has already been processed. Skipping this step...")
	return outfile


def check_ref(ref):
	logger.debug("Checking reference sequence type...")
	if os.path.splitext(ref)[1] == ".gbk" or os.path.splitext(ref)[1] == ".gb" or os.path.splitext(ref)[1] == ".gbff":
		logger.debug("Suffix = {}. Reference is a genbank file.".format(os.path.splitext(ref)[1]))
		return ".gbk"
	elif os.path.splitext(ref)[1] == ".fasta" or os.path.splitext(ref)[1] == ".fa" or os.path.splitext(ref)[1] == ".fna":
		logger.debug("Suffix = {}. Reference is a fasta file.".format(os.path.splitext(ref)[1]))
		return ".fasta"
	else:
		logger.error("Can't id reference sequence format. Ensure you have a genbank (.gb, .gbk, .gbff) or fasta (.fasta,.fa,.fna) formatted sequence.")
		sys.exit(1)


def parse_sequence_file(sequence_file):
	fasta_parse = [x for x in SeqIO.parse(sequence_file,'fasta')]
	gb_parse = [x for x in SeqIO.parse(sequence_file,'gb')]
	if len(fasta_parse) > 0:
		logger.debug("Reference is a fasta file, importing sequence only...")
		return fasta_parse, False
	elif len(gb_parse) > 0:
		logger.debug("Reference is a genbank file, importing sequence and annotations...")
		return gb_parse, True
	else:
		logger.error("Unable to parse reference sequence file, please check if it is a valid fasta or gb format.")
		sys.exit(1)


def gen_ref_chunks(sequence, threads):
	chunk_num = threads ** 2
	chunk_size = int(len(sequence)/chunk_num)
	if chunk_size < 1:
		chunk_size = 1
	logger.debug("Chopping sequence into {} chunks of size {} bp".format(chunk_num, chunk_size))
	return [(i, sequence[i:i + chunk_size]) for i in range(0, len(sequence), chunk_size)]


def check_degenerate(base):
	if base.upper() in ['A','T','G','C','N','-']:
		return False
	return True


def filter_degenerates(indexed_sequence_chunk):
	index, sequence = indexed_sequence_chunk
	degenerate_locations = [index+i for i, base in enumerate(sequence) if check_degenerate(base) == True]
	return index, degenerate_locations


def clean_sequence(results_list, sequence):
	result_dict = {i:j for i,j in results_list}
	keys = list(result_dict.keys())
	degenerate_locations = []
	for i in sorted(keys):
		degenerate_locations += result_dict[i]
	logger.info("Total %d degenerate bases removed from reference" % len(degenerate_locations))
	mutable_seq = sequence.seq.tomutable()
	for location in degenerate_locations:
		mutable_seq[location] = "N"
	new_seq = mutable_seq.toseq()
	sequence.seq = new_seq
	return sequence


def write_new_sequence_file(cleaned_sequences, annotations_present, outfile):
	if annotations_present == True:
		logger.info("Writing reference to %s" % outfile)
		null = SeqIO.write(cleaned_sequences, outfile, "gb")
	elif annotations_present == False:
		logger.info("Writing reference to %s" % outfile)
		null = SeqIO.write(cleaned_sequences, outfile, "fasta")

