#from Experiments.forontiers_jupyter.directory_structure_definer import DirectoryStructure
from Utility.aligner_wrapper import AlignerWrapper
import subprocess
from Utility.fastq_seq_size_filtering import filterFastqFileBySize
import os

bowtie_aligner = AlignerWrapper("bowtie {flags} --sam {reference_file} {lib} {positive_alignment} --un {negative_library} >> {log}.bowtie 2>&1", "bowtie-build {reference} {reference} {flags}",
		index_output_format_list=["{reference_file}"+suff for suff in [".1.ebwt",".2.ebwt",".3.ebwt",".4.ebwt",".rev.1.ebwt",".rev.2.ebwt"]])

def is_non_zero_file(fpath):  
	return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def preprocessing_mirna(mirna_filenames):
	'''
		filters out reads below 15 nts or above 25 nts
		'''
	for lib in mirna_filenames:
		filtered_filename = lib + ".filtered.fastq"
		print("filtered_filename")
		if is_non_zero_file(filtered_filename):
			print("skipping file: " + lib)
			continue
		else:
			print("filtering file: " + lib)
			filterFastqFileBySize(lib, 15,25, filtered_filename, lib+".garbage")

def align_mirna(reference_genome, sequence_directory, result_directory, wildtype_names, mutant_names, skip=False):
	aligned = []
	misaligned = []
	print("Starting WT alignments")
	for seq_lib in wildtype_names:
		fastq_lib_loc = sequence_directory + "demuxed/n2/"+ seq_lib + ".fastq" + ".filtered.fastq" 
		aligned.append(result_directory + seq_lib + "_aligned.sam")

		misaligned.append(result_directory + seq_lib + "_misaligned.fastq")
		if skip is False:
			print("WT alignment for: " + seq_lib)
			bowtie_aligner.align(reference=reference_genome, lib = fastq_lib_loc, flags = "", indexing_flags="",
					pos = aligned[-1], neg = misaligned[-1], log=res_directory+"wtlog",build_index=True)

	print("Starting mutant alignments: 0 and 3 mismatches")
	for seq_lib in mutant_names:
		fastq_lib_loc = sequence_directory + "demuxed/bb21/"+ seq_lib + ".fastq" + ".filtered.fastq" 
		aligned.append(result_directory + seq_lib + "_0mm_aligned.sam")
		misaligned.append(result_directory + seq_lib + "_0mm_misaligned.fastq")

		if skip is False:
			print("0mm alignment for: " + seq_lib)
			bowtie_aligner.align(reference=reference_genome, lib=fastq_lib_loc, flags = "-n 0",
					pos = aligned[-1], neg = misaligned[-1], log=res_directory+"0log")
		aligned.append(result_directory + seq_lib + "_3mm_aligned.sam")
		misaligned.append(result_directory + seq_lib + "_3mm_misaligned.fastq")
		if skip is False:
			print("3mm alignment for: " + seq_lib)
			bowtie_aligner.align(reference=reference_genome, lib=fastq_lib_loc, flags = "-n 3",
					pos = aligned[-1], neg = misaligned[-1], log=res_directory+"3log")

	return aligned, misaligned

def create_pileups(aligned, reference_genome):
	for pos in aligned:
		print("Creating pileup for: " + pos)
		sam_to_pileup(pos, reference_genome, pos+".bam", pos+".sorted.bam", pos+".pileup", pos+".sorted.pileup")

def create_genecount(aligned):
	for pos in aligned:
		print(f"generating genecount for {pos}")
		command = f"cat {pos}|grep -v \"^@\" | cut -f3 | sort | uniq -c | grep -v \"\*\" | awk  ' BEGIN {{  OFS=\",\" }}  {{ print $2,$1  }} ' | grep -v '^*'> {pos}.gene_count"
		try:
			subprocess.check_output(command, shell=True)
		except subprocess.CalledProcessError as e:
			print(f"genecount of {pos} failed")
			raise e
		#todo: check correctness

def sam_to_pileup(sam_name,fasta_name,bam_name,sorted_bam_name,pileup_name,sorted_pileup_name):
	command=(   f"samtools view -bS {sam_name} > {bam_name} && " +
			f"samtools sort {bam_name} -o {sorted_bam_name} && " +
			f"samtools mpileup -B -f {fasta_name} {sorted_bam_name} | tail -n +3 > {pileup_name} && "
			f"sed -r '/^[\t]*$/d' <{pileup_name} | LC_ALL=C sort -k1,1 -k2,2n -o {sorted_pileup_name} "
			)
	try:
		outout=subprocess.check_output(command,shell=True)
	
	except subprocess.CalledProcessError as e:
		print(f"pileup conversions of {sam_name} failed")
		raise e
	return

if __name__ == "__main__":
	#inputs:
		reference_genome = "/home/jupyter-mahmud/reference/mirbasev235.fasta"
		sequence_directory = "/home/jupyter-mahmud/libraries/"
		wt_seq_names = [	"11311_N2_Emb_Primary_siRNA_biorep1_1A","11321_N2_Emb_Primary_siRNA_biorep2_1B","11331_N2_Emb_Primary_siRNA_biorep3_1C_techrep_1"
				,"11332_N2_Emb_Primary_siRNA_biorep3_1C_techrep_2" ,"11333_N2_Emb_Primary_siRNA_biorep3_1C_techrep_3","11341_N2_Emb_Primary_siRNA_biorep4_1D"]
		mut_seq_names = ["21311_BB21_Emb_Primary_siRNA_biorep1_5A","21321_BB21_Emb_Primary_siRNA_biorep2_5B","21331_BB21_Emb_Primary_siRNA_biorep3_5C"]
		wt_seq_library = [sequence_directory + "demuxed/n2/"+ name + ".fastq" for name in wt_seq_names]
		mut_seq_library = [sequence_directory + "demuxed/bb21/" + name + ".fastq" for name in mut_seq_names]
		res_directory = "/home/jupyter-mahmud/mirna_results/"

		#preprocessing (adds .filtered.fastq before file extension)
		preprocessing_mirna(wt_seq_library+mut_seq_library)
		wt_seq_library = [s+".filtered.fastq" for s in wt_seq_library]
		mut_seq_library = [s+".filtered.fastq" for s in mut_seq_library]

		#do alignment
		aligned, misaligned = align_mirna(reference_genome, sequence_directory, res_directory, wt_seq_names, mut_seq_names, skip=True)	

		#pileup_creation:
		create_pileups(aligned, reference_genome)

		#generate genecount for differential expression:
		create_genecount(aligned)

		print("finished")




