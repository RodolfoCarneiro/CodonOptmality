import sys
import argparse
import random

parser = argparse.ArgumentParser(description = 'Creates a fasta-format random genome, relative to a model genome', usage = '[options]')
parser.add_argument('-F', '--fasta', type=str, required=True, help='model genome fasta archive')
parser.add_argument('-O', '--output', type=str, required=True, help='output file (fasta format)')
parser.add_argument('--keep_seq', action='store_const', const=True, default=False, help='Keep aminoacids sequence.')
parser.add_argument('-E', '--external', type=str, default=False, help='Use an external file for aminoacid frequence')
args = parser.parse_args()

codon_dictionary = {
	'TTT':'Phe', 'TTC':'Phe', 'TTA':'Leu', 'TTG':'Leu',
	'CTT':'Leu', 'CTC':'Leu', 'CTA':'Leu', 'CTG':'Leu',
	'ATT':'Ile', 'ATC':'Ile', 'ATA':'Ile', 'ATG':'Met',
	'GTT':'Val', 'GTC':'Val', 'GTA':'Val', 'GTG':'Val',

	'TCT':'Ser', 'TCC':'Ser', 'TCA':'Ser', 'TCG':'Ser',
	'CCT':'Pro', 'CCC':'Pro', 'CCA':'Pro', 'CCG':'Pro',
	'ACT':'Thr', 'ACC':'Thr', 'ACA':'Thr', 'ACG':'Thr',
	'GCT':'Ala', 'GCC':'Ala', 'GCA':'Ala', 'GCG':'Ala',

	'TAT':'Tyr', 'TAC':'Tyr',
	'CAT':'His', 'CAC':'His', 'CAA':'Gln', 'CAG':'Gln',
	'AAT':'Asn', 'AAC':'Asn', 'AAA':'Lys', 'AAG':'Lys',
	'GAT':'Asp', 'GAC':'Asp', 'GAA':'Glu', 'GAG':'Glu',

	'TGT':'Cys', 'TGC':'Cys',              'TGG':'Trp',
	'CGT':'Arg', 'CGC':'Arg', 'CGA':'Arg', 'CGG':'Arg',
	'AGT':'Ser', 'AGC':'Ser', 'AGA':'Arg', 'AGG':'Arg',
	'GGT':'Gly', 'GGC':'Gly', 'GGA':'Gly', 'GGG':'Gly'
	}

aminoacid_dictionary = {
	'Phe':0,  'Leu':1,  'Ile':2,  'Met':3,  'Val':4,  'Ser':5,  'Pro':6,  'Thr':7,  'Ala':8,  'Tyr':9, 
	'His':10, 'Gln':11, 'Asn':12, 'Lys':13, 'Asp':14, 'Glu':15, 'Cys':16, 'Trp':17, 'Arg':18, 'Gly':19
	}

#aminoacid_list = [
#	'Phe', 'Leu', 'Ile', 'Met', 'Val', 'Ser', 'Pro', 'Thr', 'Ala', 'Tyr', 
#	'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Gly'
#	]

def percent(index, genome_length): # Iterface updating the percentage of the ORFs read
	percentage = 100*index/genome_length
	interface = int(percentage)
	print str(interface)+'% complete\r',

def extract_sequence(genome_list, index, header=False):
		gene_code = ''
		genome_list[index] = genome_list[index].replace('\r', '')
		gene_by_line = str.split(genome_list[index], '\n') # Split each line
		if header == True:
			return str.split(gene_by_line[0], ' ')[0] # The first line is the header, the first word is its code
		else:
			for j in range(1, len(gene_by_line)):
				gene_code += gene_by_line[j] # Every line (except the first) is part of the code

			return list(gene_code) # List of all aminoacids in a gene

def get_codon(aminoacids_for_codon, index_for_codon, name, gene_length): # Writes the codon for the list
	if index_for_codon < gene_length: # Test if codon is not beyond the end of the gene 
		try: # Try to define the codon
			obtained_codon = aminoacids_for_codon[index_for_codon*3] + aminoacids_for_codon[index_for_codon*3+1] + aminoacids_for_codon[index_for_codon*3+2]
		except Exception: # If not possible, there is no enough aminoacids
			print "\r" + str.split(name, ' ')[0] + ": fragmented codon"
			return False # Errors return "False" value

		if obtained_codon == "TAG" or obtained_codon == "TGA" or obtained_codon == "TAA": # Test for nonsense stop codon
			if index_for_codon > gene_length - 1: # Nonsense only if not the last codon
				print "\r" + str.split(name, ' ')[0]+": nonsense stop codon"
				return False # Errors return "False" value
			else:
				return False # Errors return "False" value
		else:
			return obtained_codon # Return the codon code
	else: return False # Errors return "False" value

def make_external_list(freq_file, keep_sequence): # Makes a list containing external values.
	freq_archive = freq_file.replace('\r', '')
	freq_lines = str.split(freq_archive, '\n') # Split every line

	if keep_sequence == True: # If keeping aminoacid sequence
		freq_distribution = [[] for entry in range(20)] # Create a list containing 20 internal lists, for every aminoacid

		for entry in range(61): # Read every entry in external file
			current_codon_code = str.split(freq_lines[entry], ';')[0] # Define current codon
			current_codon_value = float(str.split(freq_lines[entry], ';')[1]) # Define current codon frequency value

			current_aminoacid = aminoacid_dictionary[codon_dictionary[current_codon_code]] # Define the aminoacid for that codon
			length = len(freq_distribution[current_aminoacid]) # Define the number of codons alredy written in current aminoacid list

			if length != 0: # If it is not the first codon to be written in the aminoacid list
				previous_value = freq_distribution[current_aminoacid][length - 1][0] # Value of the previous codon in the list
				current_codon_value += previous_value # The current value must be cumulative to the previous ones

			freq_distribution[current_aminoacid].append([current_codon_value, current_codon_code])

	else: # If not keeping aminoacid sequence
		freq_distribution = [[0,''] for entry in range(61)] # Create a list containing 61 internal lists, for every codon

		for entry in range(61): # read every entry in external file
			current_codon = str.split(freq_lines[entry], ';') # Split codon code from CSC value
			if entry == 0:
				freq_distribution[0][0] = float(current_codon[1]) # Write first value on list
			else:
				freq_distribution[entry][0] = float(current_codon[1]) + freq_distribution[entry-1][0] # Write comulative frequency
			freq_distribution[entry][1] = current_codon[0] # Write codon code on list
	
	return freq_distribution

def distribute(genome_list, keep_sequence): # Makes a list with the frequency of every codon found in the base fasta archive

	if keep_sequence == True: # If keeping aminoacid sequence
		freq_distribution = [[] for entry in range(20)] # Create a list containing 20 internal lists, for every aminoacid

		for i in range(1, len(genome_list)): # Repeat for each gene
			percent(i, 2*len(genome_list)) # Print percentage until 50% (first of 2 times the fasta file will be screened)
			bases_list = extract_sequence(genome_list, i) # Make a list with every base of the gene

			for j in range(len(bases_list)/3): # For every codon
				current_codon_code = get_codon(bases_list, j, extract_sequence(genome_list, i, header=True), len(bases_list))
				if current_codon_code == False:
					break

				current_aminoacid = aminoacid_dictionary[codon_dictionary[current_codon_code]] # Define the aminoacid for that codon
				length = len(freq_distribution[current_aminoacid]) # Define the number of codons alredy written in current aminoacid list

				match = False
				if length == 0: # If current aminoacid list is empty
					freq_distribution[current_aminoacid].append([0,current_codon_code]) # Add first value

				else:
					for k in range(length): # For each codon in list
						if current_codon_code == freq_distribution[current_aminoacid][k][1]: # If it is the current codon
							match = True
						if match == True: # Add 1 to current codon. Cumulative to the next codons
							freq_distribution[current_aminoacid][k][0] += 1
						if k == length-1 and match == False: # If the current codon is not in the list
							current_codon_value = freq_distribution[current_aminoacid][k][0] # Add the codon, current frequency +1
							freq_distribution[current_aminoacid].append([current_codon_value+1, current_codon_code])

	else: # If not keeping aminoacid sequence
		freq_distribution = [] # Create a list containing 61 internal lists, for every codon

		for i in range(1, len(genome_list)): # Repeat for each gene
			percent(i, 2*len(genome_list)) # Print percentage until 50% (first of 2 times the fasta file will be screened)
			bases_list = extract_sequence(genome_list, i) # Make a list with every base of the gene

			for j in range(len(bases_list)/3): # For every codon
				current_codon_code = get_codon(bases_list, j, extract_sequence(genome_list, i, header=True), len(bases_list))
				if current_codon_code == False:
					break

				match = False

				if len(freq_distribution) == 0: # If current aminoacid list is empty
					freq_distribution.append([0,current_codon_code]) # Add first value

				for k in range(len(freq_distribution)): # For each codon in list
					if current_codon_code == freq_distribution[k][1]: # If it is the current codon
						match = True
					if match == True: # Add 1 to current codon. Cumulative to the next codons
						freq_distribution[k][0] += 1
					if k == len(freq_distribution)-1 and match == False: # If the current codon is not in the list
						current_codon_value = freq_distribution[k][0] # Add the codon, current frequency +1
						freq_distribution.append([current_codon_value+1, current_codon_code])
	return freq_distribution

def main():
	file_in = open(args.fasta, 'r') # Input file
	in_file = file_in.read()
	fasta_list = str.split(in_file, '>') # Split genes

	out_file = open(args.output, 'w') # Output file

	keep_sequence = args.keep_seq

	if args.external == False:
		codons = distribute(fasta_list, keep_sequence) # Make codons list from base file

	else:
		file_freq = open(args.external, 'r') # External frequency file
		freq_file = file_freq.read()
		codons = make_external_list(freq_file, keep_sequence) # Make codons list from external file
		last = len(codons)-1 # Last value of codons list


	for i in range(len(fasta_list)): # Repeat for each gene

		if args.external == False: # If there is a external frequency file
			percent(i + len(fasta_list), 2*len(fasta_list)) # Print 50% to 100% (Last of 2 times the fasta file will be screened)
		else:
			percent(i, len(fasta_list)) # Print percentage (1% - 100%)

		out_file.write('>' + extract_sequence(fasta_list, i, header=True) + '\n')

		bases_list = extract_sequence(fasta_list, i) # Make a list with every base of the gene

		codons_in_line = 0

		for j in range(len(bases_list)/3): # For every codon

			if codons_in_line == 20: # Writes only 60 nitrogenous bases in each line
				out_file.write('\n')
				codons_in_line = 0

			current_codon_code = get_codon(bases_list, j, extract_sequence(fasta_list, i, header=True), len(bases_list))
			if current_codon_code == False:
				break

			if keep_sequence == True: # If keeping aminoacid sequence
				current_aminoacid = aminoacid_dictionary[codon_dictionary[current_codon_code]] # Define the aminoacid for that codon
				last = len(codons[current_aminoacid])-1 # Last value of codons list
			
				random_number = random.random()*codons[current_aminoacid][last][0] # Randomize codon to be written

				match = False
				for k in range(len(codons[current_aminoacid])):
					if random_number <= codons[current_aminoacid][k][0] and match == False: # Test if the current codon is the random one
						out_file.write(current_codon_code) # Write
						codons_in_line += 1 # Add 1 to count of codons in line
						match = True
					if match == True and args.external == False: # Deduct 1 from codon value for every codon after the match
						codons[current_aminoacid][k][0] -= 1

			else: # If not keeping aminoacid sequence
				last = len(codons)-1 # Last value of codons list
			
				random_number = random.random()*codons[last][0] # Randomize codon to be written

				match = False
				for k in range(len(codons)):
					if random_number <= codons[k][0] and match == False: # Test if the current codon is the random one
						out_file.write(current_codon_code) # Write
						codons_in_line += 1 # Add 1 to count of codons in line
						match = True
					if match == True and args.external == False: # Deduct 1 from codon value for every codon after the match
						codons[k][0] -= 1

		out_file.write('\n')

	print '100% complete'

main()