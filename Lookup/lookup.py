import sys
import argparse

parser = argparse.ArgumentParser(description='Screens a fasta archive containing ORFs, and looks for a stretch with a given CSC mean value.\nThe output file variables (gene lenght, start pos, end pos) are defined in codons.', usage='[options]')
parser.add_argument('-F', '--fasta', type=str, required=True, help='fasta archive containing ORFs list')
parser.add_argument('-O', '--output', type=str, required=True, help='output file (csv format)')
parser.add_argument('-C', '--csc', type=str, required=True, help='.csv archive containing CSC codons and values')
parser.add_argument('-L', '--length', type=int, required=True, help='length of the stretch')
parser.add_argument( '--min', type=float, default=-9999, help='minimum mean CSC value. (Optional)')
parser.add_argument( '--max', type=float, default=9999, help='maximum mean CSC value. (Optional)')
parser.add_argument( '--best', default=False, action='store_true', help='Writes only the best stretch (Optional)')
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
	print str(interface) + '% complete\r',

def get_codon(aminoacids_for_codon, index_for_codon, name, gene_length): # Writes the codon for the list
	if index_for_codon < gene_length: # Test if codon is not beyond the end of the gene 
		try: # Try to define the codon
			obtained_codon = aminoacids_for_codon[index_for_codon*3] + aminoacids_for_codon[index_for_codon*3+1] + aminoacids_for_codon[index_for_codon*3+2]
		except Exception: # If not possible, there is no enough aminoacids
			print '\r' + str.split(name, ' ')[0] + ': fragmented codon'
			return False # Errors return 'False' value

		if obtained_codon == 'TAG' or obtained_codon == 'TGA' or obtained_codon == 'TAA': # Test for nonsense stop codon
			if index_for_codon > gene_length - 1: # Nonsense only if not the last codon
				print '\r' + str.split(name, ' ')[0]+': nonsense stop codon'
				return False # Errors return 'False' value
			else:
				return False # Errors return 'False' value
		else:
			return obtained_codon # Return the codon code
	else: return False # Errors return 'False' value

def screen_gene(gene_aminoacids, codons, out_file, name): # Screen the given gene after the defined stretches

	number_of_codons = len(gene_aminoacids)/3

	j = 0
	CSC_sum = 0.0
	for k in range(args.length): # k value is the index of the codon inside the stretch
		codon = get_codon(gene_aminoacids, k, name, number_of_codons)
		if codon == False: # exit gene if error
			break
		CSC_sum += codons[codon] # Add the current codon value
	if codon == False: # exit gene if error
		return

	if args.best == False:
		no_end = False # True if a stretch has a start position but not a end position yet
		while j < number_of_codons - args.length: # j value is the index for the beginning of the stretch
			print str(args.max) + ' > ' + str(CSC_sum) + ' > ' + str(args.min) + '    ' + str(no_end) + ' = False'
			if CSC_sum >= args.min and CSC_sum <= args.max and no_end == False: # Test if the stretch corresponds to the wanted values
				out_file.write('\n' + str.split(name, ' ')[0] + ';' + str(number_of_codons) + ';' + str(j+1) + ';') # Write variables
				no_end = True
		
			j += 1 # Go to the next codon

			codon = get_codon(gene_aminoacids, j-1, name, number_of_codons) # Define the codon that will be excluded in the next stretch
			CSC_sum -= codons[codon] # Deduct its value
			codon = get_codon(gene_aminoacids, j-1+args.length, name, number_of_codons) # Define the codon that will me added in the next stretch
			if codon == False: # If there is no codon to be added
				if no_end == True:
					out_file.write(str(j + args.length)) # Write the last codon position
					no_end = False
				break
			CSC_sum += codons[codon] # Add the next codon value

			if (CSC_sum <= args.min or CSC_sum >= args.max) and no_end == True: # If it does not correspond to the wanted values anymore
				out_file.write(str(j + args.length)) # Write the end position
				no_end = False

	if args.best == True:
		current_gene_value = 0
		current_gene = ''
		while j < number_of_codons - args.length: # j value is the index for the beginning of the stretch

			# For positive values:
			#print str(args.max) + ' > ' + str(CSC_sum) + ' > ' + str(args.min) + ' > 0    !    ' + str(CSC_sum) + ' > ' + str(current_gene_value)
			if CSC_sum >= args.min and CSC_sum <= args.max and args.min > 0 and CSC_sum > current_gene_value: # Test if the stretch corresponds to the wanted values
				current_gene = '\n' + str.split(name, ' ')[0] + ';' + str(number_of_codons) + ';' + str(j+1) # Save output text structure
				current_gene_value = CSC_sum # Save value for future comparison

			# For negative values:
			if CSC_sum >= args.min and CSC_sum <= args.max and args.max < 0 and CSC_sum < current_gene_value: # Test if the stretch corresponds to the wanted values
				current_gene = '\n' + str.split(name, ' ')[0] + ';' + str(number_of_codons) + ';' + str(j+1) # Save output text structure
				current_gene_value = CSC_sum # Save value for future comparison
		
			j += 1 # Go to the next codon

			codon = get_codon(gene_aminoacids, j-1, name, number_of_codons) # Define the codon that will be excluded in the next stretch
			CSC_sum -= codons[codon] # Deduct its value
			codon = get_codon(gene_aminoacids, j-1+args.length, name, number_of_codons) # Define the codon that will me added in the next stretch
			if codon == False: # If there is no codon to be added
				break
			CSC_sum += codons[codon] # Add the next codon value

		out_file.write(current_gene) # Write line
		

def main():
	if args.min >= args.max:
		print 'Minimum value must be lower than maximum value'
		sys.exit()

	if args.best == True and args.min < 0 and args.max > 0:
		print 'Range must not contain zero when looking for best stretches'
		sys.exit()

	file_in = open(args.fasta, 'r') # Input file
	in_file = file_in.read()

	out_file = open(args.output, 'w') # Output file
	out_file.write('gene;gene length (in codons);start position (in codons)')
	if args.best == False:
		out_file.write(';end position (in codons)')

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
		name = str.split(header,' ')[0]
		for j in range(1, len(gene_by_line)):
			gene_code += gene_by_line[j] # Every line (except the first) is part of the code

		gene_aminoacids = list(gene_code) # List of all aminoacids in a gene
		screen_gene(gene_aminoacids, codons, out_file, name)
	print '\r100% complete'

main()
