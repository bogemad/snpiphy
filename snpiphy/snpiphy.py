#!/usr/bin/env python

import sys, os, re, shutil, subprocess, snpiphy, random, logging
from contextlib import contextmanager
from Bio import Phylo

class SnpiPhy:
	def __init__(self, threads, cutoff, reads_dir, outdir, reference, prefix):
		self.logger = logging.getLogger(__name__)
		self.outdir = os.path.abspath(outdir)
		self.reads_dir = os.path.abspath(reads_dir)
		self.reference = snpiphy.remove_degenerates(reference, os.path.join(outdir, ".".join([prefix, "reference"])), threads)
		self.threads = threads
		self.cutoff = cutoff
		self.temp_dir = os.path.join(outdir, '.temp')
		self.ref_aligns = snpiphy.define_make_outdirs(os.path.join(self.outdir, "reference_alignments"))
		self.sample_names = [ snpiphy.get_samplename(path) for path in os.listdir(self.reads_dir) if os.path.isfile(os.path.join(self.reads_dir, path)) ]
		self.snippy_not_done = [ os.path.join(self.reads_dir, path) for path in os.listdir(self.reads_dir) if (snpiphy.get_samplename(path) in os.listdir(self.ref_aligns)) == False ]
		self.core_align = snpiphy.define_make_outdirs(os.path.join(self.outdir, "core_genome_alignment"))
		self.excluded_seqs = snpiphy.define_make_outdirs(os.path.join(self.outdir, "excluded_sequences_and_alignments"))
		self.recomb_filter = snpiphy.define_make_outdirs(os.path.join(self.outdir, "recombination_filtration"))
		self.phylogenetic_trees = snpiphy.define_make_outdirs(os.path.join(self.outdir, "intermediate_phylogenetic_trees"))
		self.raw_core_aln = os.path.join(self.outdir, '.'.join((prefix,"raw_core.alignment.fasta")))
		self.N_trimmed_core_aln = os.path.join(self.outdir, '.'.join((prefix,"N_trimmed_core.alignment.fasta")))
		self.init_filtered_tree = os.path.join(self.outdir, '.'.join((prefix,"N_trimmed_core.tree.tre")))
		self.recomb_figure = os.path.join(self.outdir, '.'.join((prefix,"recombination_predictions.pdf")))
		self.bootstrap_tree = os.path.join(self.outdir, '.'.join((prefix,"bootstap_trees.nex")))
		self.final_unrooted_tree = os.path.join(self.outdir, '.'.join((prefix,"final_tree_with_bootstrap.nex")))
		self.core_snp_stats = os.path.join(self.outdir, '.'.join((prefix,"recombination_filtered_core.snp_counts.stats")))
		self.core_snp_matrix = os.path.join(self.outdir, '.'.join((prefix,"recombination_filtered_core.snp_counts.tsv")))
	
	def remove_results_without_reads(self):
		self.logger.debug("Removing results from previous runs that have been removed from reads directory...\n")
		for x in os.listdir(self.outdir):
			abspath = os.path.join(self.ref_aligns, x)
			if os.path.isdir(abspath) == True:
				matches = len([y for y in os.listdir(self.reads_dir) if re.match(x, y)])
				if matches == 1:
					self.logger.debug("Results directory {} is present, continuing...")
					continue
				elif matches > 1:
					self.logger.warn("Something strange going on here. More than one read file for results directory {}. You should check this.".format(x))
				else:
					self.logger.warn("Results directory {} has no associated reads file and has been moved to {}.".format(x, self.excluded_seqs))
					shutil.move(abspath, self.excluded_seqs)
	
	def run_snippy(self, read_path):
		name = snpiphy.get_samplename(read_path)
		self.logger.info("{}: Running snippy...\n".format(name))
		read_format = snpiphy.id_read_format(read_path)
		if read_format == 'fastq':
			ec = snpiphy.run_command([
									'snippy',
									'--cpus', str(self.threads),
									'--prefix', name,
									'--outdir', os.path.join(self.ref_aligns, name),
									'--ref', self.reference,
									'--peil', read_path
									])
			if ec != 0:
				self.logger.error("Error running snippy on fastq reads file: {}. Please check your files and error output.".format(os.path.basename(read_path)))
				sys.exit(1)
		else:
			ec = snpiphy.run_command([
									'snippy',
									'--cpus', str(self.threads),
									'--prefix', name,
									'--outdir', os.path.join(self.ref_aligns, name),
									'--ref', self.reference,
									'--ctgs', read_path
									])
			if ec != 0:
				self.logger.error("Error running snippy on fastq reads file: {}. Please check your files and error output.".format(os.path.basename(read_path)))
				sys.exit(1)
		self.logger.info("{}: snippy has completed successfully.\n".format(name))
	
	def run_snippy_core(self):
		self.logger.info("Building core genome alignment from individual reference alignments...")
		if os.path.isfile(self.raw_core_aln) == False:
			with snpiphy.cd(self.ref_aligns):
				self.logger.debug("Running snippy-core...")
				ec = snpiphy.run_command((
										[
										"snippy-core"
										]+self.sample_names
										))
				if ec != 0:
					self.logger.error("Error building core alignment before distant sequences removed. Please check your files and error output.")
					sys.exit(1)
				self.logger.debug("Examining snippy core output for poorly aligned sequences...")
				core_data_handle = open("core.txt", 'r')
				exclude_log = open(os.path.join(self.excluded_seqs, "removed_sequences.log"), 'w')
				bad_seqs = []
				for line in core_data_handle:
					if line.startswith("ID"):
						continue
					line_data = line.strip().split("\t")
					if float(line_data[3]) < self.cutoff:
						archive = "{}.alignment.tar.gz".format(line_data[0])
						reads_file = snpiphy.find_source_file(line_data[0], self.reads_dir)
						moved_reads_file = find_source_file(line_data[0], self.excluded_seqs)
						moved_archive = os.path.join(self.excluded_seqs, "{}.tar.gz".format(line_data[0]))
						self.logger.info("Sample {} coverage is too low ({}%), removing from core alignment. Reads and reference mapping data will be retained in archive: {}".format(line_data[0],line_data[3],archive))
						exclude_log.write("Sample {} coverage is too low ({}%), removing from core alignment. Reads and reference mapping data will be retained in archive: {}\n".format(line_data[0],line_data[3],archive))
						exitcode = snpiphy.run_command(["tar", "cvzf", archive, line_data[0]])
						if ec != 0:
							self.logger.error("Error compressing excluded alignment for sample: {}. Please check your files and error output.".format(line_data[0]))
							sys.exit(1)
						os.rename(reads_file, moved_reads_file)
						os.rename(archive, moved_archive)
						shutil.rmtree(line_data[0])
						bad_seqs.append(line_data[0])
					else:
						self.logger.debug("Sample {} coverage is ok ({}%)".format(line_data[0],line_data[3]))
				core_data_handle.close()
				exclude_log.close()
				self.logger.debug("Deleting first core alignment...")
				alignment_files = [os.path.join(self.ref_aligns, item) for item in os.listdir(self.ref_aligns) if os.path.isfile(item)]
				for file in alignment_files:
					if file.startswith("core"):
						os.remove(file)
				self.logger.debug("Running snippy-core again...")
				ec = snpiphy.run_command([
										"snippy-core"
										]+[
										x for x in self.sample_names if (x in bad_seqs) == False
										])
				if ec != 0:
					self.logger.error("Error building core alignment after distant sequences removed. Please check your files and error output.")
					sys.exit(1)
				self.logger.debug("Moving core alignments to {}".format(self.core_align))
				alignment_files = [os.path.join(self.ref_aligns, item) for item in os.listdir(self.ref_aligns) if os.path.isfile(item)]
				for file in alignment_files:
					if os.path.basename(file).startswith("core"):
						shutil.move(file, self.core_align)
				shutil.move(os.path.join(self.core_align, "core.full.aln"), self.raw_core_aln)
				os.symlink(self.raw_core_aln, os.path.join(self.core_align, "core.full.aln"))
		else:
			self.logger.info("Core alignment has already been generated. Skipping this step...")
	
	def Ns2gaps(self):
		if os.path.isfile(self.N_trimmed_core_aln) == False:
			self.logger.info("Converting Ns generated by snippy-core to gaps...")
			snpiphy.replace_Ns_with_gaps(self.core_align, self.threads)
			shutil.move(os.path.join(self.core_align, "core.full.trimmed.aln"), self.N_trimmed_core_aln)
			os.symlink(self.N_trimmed_core_aln, os.path.join(self.core_align, "core.full.trimmed.aln"))
		else:
			self.logger.info("snippy-core Ns have already been converted to gaps. Skipping this step...")
	
	def filter_recombinant_positions(self):
		if os.path.isfile(self.init_filtered_tree) == False:
			with snpiphy.cd(self.recomb_filter):
				self.logger.info("Scanning and filtering recombination positions with gubbins...")
				ec = snpiphy.run_command([
										"run_gubbins.py",
										"-v",
										'-i', '10',
										'-p', os.path.join(self.recomb_filter, 'filtered_core_aln'),
										'-c', str(self.threads),
										os.path.join(self.core_align, 'core.full.trimmed.aln')
										])
				if ec != 0:
					self.logger.warn("Recombination filtering using the RAxML only method has failed. Retrying with FastTree for first iteration.")
					for file in os.listdir(self.recomb_filter):
						if file.startswith('core.full.trimmed.aln.'):
							os.remove(file)
					ec = snpiphy.run_command([
											"run_gubbins.py",
											"-v", 
											'-i', '10', 
											'--tree_builder', 'hybrid', 
											'-p', os.path.join(self.recomb_filter, 'filtered_core_aln'), 
											'-c', str(self.threads), 
											os.path.join(self.core_align, 'core.full.trimmed.aln')
											])
					if ec != 0:
						self.logger.warn("Recombination filtering using hybrid RAxML/FastTree method has failed. Retrying with FastTree for all iterations.")
						for file in os.listdir(self.recomb_filter):
							if file.startswith('core.full.trimmed.aln.'):
								os.remove(file)
						ec = snpiphy.run_command([
												"run_gubbins.py",
												"-v", 
												'-i', '10', 
												'--tree_builder', 'fasttree', 
												'-p', os.path.join(self.recomb_filter, 'filtered_core_aln'), 
												'-c', str(self.threads), 
												os.path.join(self.core_align, 'core.full.trimmed.aln')
												])
						if ec != 0:
							self.logger.error("Running gubbins using all methods have failed. Please examine your alignment or consider removing highly divergent sequences. Looking at total.snp_counts.stats is a good place to start. Additionally consider using a different reference sequence.")
							sys.exit(1)
			Phylo.convert(os.path.join(self.recomb_filter, "filtered_core_aln.final_tree.tre"), 'newick', self.init_filtered_tree, 'nexus')
			self.logger.info("Recombination filtering by gubbins has completed successfully.")
		else:
			self.logger.info("Recombination filtering by gubbins has already been done. Skipping this step...")
	
	def visualise_recombination_predictions(self):
		if os.path.isfile(self.recomb_figure) == False:
			self.logger.info("Building a recombination prediction figure...")
			ec = snpiphy.run_command([
									'gubbins_drawer.py',
									'-o', self.recomb_figure,
									os.path.join(self.recomb_filter, "filtered_core_aln.final_tree.tre"),
									os.path.join(self.recomb_filter, "filtered_core_aln.recombination_predictions.embl")
									])
			if ec != 0:
				self.logger.warn("Recombination visualisation failed. This isn't critical so continuing to tree building...")
		else:
			self.logger.info("Recombination prediction visualisation has already been generated. Skipping this step...")
	
	def gen_bootstrap_tree(self):
		if os.path.isfile(self.bootstrap_tree) == False:
			self.logger.info("Running bootstrap analysis...")
			ec = snpiphy.run_command([
									"raxmlHPC-PTHREADS",
									"-T",str(self.threads),
									"-m","GTRCAT",
									"-V",
									"-p",str(random.randint(10000,99999)),
									"-b",str(random.randint(10000,99999)),
									"-#","100","-s",os.path.join(self.recomb_filter, "filtered_core_aln.filtered_polymorphic_sites.fasta"),
									"-n","filtered_core_aln.bootstrap",
									'-w', self.phylogenetic_trees
									])
			if ec != 0:
				self.logger.error("RAxML bootstrap has failed.")
				sys.exit(1)
			Phylo.convert(os.path.join(self.phylogenetic_trees, "RAxML_bootstrap.filtered_core_aln.bootstrap"), 'newick', self.bootstrap_tree, 'nexus')
		else:
			self.logger.info("Bootstrap RAxML trees have already been generated. Skipping this step...")
	
	def gen_final_tree(self):
		if os.path.isfile(self.final_unrooted_tree) == False:
			self.logger.info("Generating final tree...")
			ec = snpiphy.run_command([
									"raxmlHPC-PTHREADS",
									"-T", str(self.threads),
									"-m", "GTRCAT",
									"-V",
									"-p", str(random.randint(10000,99999)),
									"-f", "b",
									"-t", os.path.join(self.recomb_filter, "filtered_core_aln.final_tree.tre"),
									"-z", os.path.join(self.phylogenetic_trees, "RAxML_bootstrap.filtered_core_aln.bootstrap"),
									"-n", "filtered_core_aln.final",
									'-w', self.phylogenetic_trees
									])
			if ec != 0:
				self.logger.error("Final RAxML tree generation has failed.")
				sys.exit(1)
			Phylo.convert(os.path.join(self.phylogenetic_trees, 'RAxML_bipartitions.filtered_core_aln.final'), 'newick', self.final_unrooted_tree, 'nexus')
		else:
			self.logger.info("Final RAxML tree has already been generated. Skipping this step...")
	
	def count_core_snps(self):
		if os.path.isfile(os.path.join(self.outdir, "recombination_filtered_core.snp_counts.stats")) == False:
			self.logger.info("Calculating pairwise recombination-filtered_core SNPs...")
			snpiphy.count_snps(os.path.join(self.recomb_filter, 'filtered_core_aln.filtered_polymorphic_sites.fasta'), self.core_snp_matrix, self.core_snp_stats, self.threads)
		else:
			self.logger.info("Recombination filtered core snp count matrix has already been generated. Skipping this step...")



