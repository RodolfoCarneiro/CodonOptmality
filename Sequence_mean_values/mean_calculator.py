import sys
import argparse

parser = argparse.ArgumentParser(description = 'Calculates CSCg (CSC mean value of a gene) for a list of genes.', usage = '[options]')
parser.add_argument('-F', '--fasta', type=str, required=True, help='fasta archive containing ORFs list')
parser.add_argument('-O', '--output', type=str, required=True, help='output file (csv format)')
parser.add_argument('-C', '--csc', type=str, required=True, help='.csv archive containing CSC codons and values')
args = parser.parse_args()

def make_csc(csc_archive): # Makes a dictionary containing input CSC values.
	csc_dictionary = {}
	csc_archive = csc_archive.replace('\r', '')
	csc_lines = str.split(csc_archive, '\n') # Splits every line
	for entry in range(61):
		current_codon = str.split(csc_lines[entry], ';') # Splits codon code from CSC value
		csc_dictionary[current_codon[0]] = float(current_codon[1]) # Writes on dictionary
	return csc_dictionary

def percent(index, genome_lenght): # Iterface updating the percentage of the ORFs read
	percentage = 100*index/genome_lenght
	interface = int(percentage)
	print str(interface)+'% complete\r',

def evaluate(aminoacids, gene_name, codons_dictionary, output_file): # Calculates CSCg values
	csc_count = 0.0
	for j in range(len(aminoacids)/3):
		codon = aminoacids[(j)*3] + aminoacids[(j)*3+1] + aminoacids[(j)*3+2] # Defines the current codon
		if codon == 'TGA' or codon == 'TAG' or codon == 'TAA': # Test for a stop codon
			if j != len(aminoacids)/3-1: # The last codon may be a stop codon
				print 'Stop codon in ' + str.split(gene_name, ' ')[0]
			break
		csc_count += codons_dictionary[codon] # Add codon value
	try:
		csc_count = csc_count*3/len(aminoacids) # Calculate the mean value of the codon
	except:
		return
	output_file.write(str.split(gene_name,' ')[0] + ',' + str(csc_count)) # Write gene and CSCg

def main():
	file_in = open(args.fasta, 'r') # Input file
	in_file = file_in.read()

	out_file = open(args.output, 'w') # Output file

	file_csc = open(args.csc, 'r') # Csc file
	csc_file = file_csc.read()
	codons = make_csc(csc_file) # Make a dictionary with csc values

	genome_list = str.split(in_file, '>') # Split genes by '>'

	for i in range(1, len(genome_list)): # Repeat for each gene
		percent(i, (len(genome_list)))

		gene_code = ''

		genome_list[i] = genome_list[i].replace('\r', '')
		gene_by_line = str.split(genome_list[i], '\n') # Split each line
		header = gene_by_line[0] # The first line is the header
		for j in range(1, len(gene_by_line)):
			gene_code += gene_by_line[j] # Every line (except the first) is part of the code

		gene_aminoacids = list(gene_code) # List of all aminoacids in a gene
		evaluate(gene_aminoacids, header, codons, out_file)
		if i < len(genome_list) - 1:
			out_file.write('\n')

	print '100% complete\nDone!'
	file_in.close()
	out_file.close()

main()