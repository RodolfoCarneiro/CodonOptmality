import sys
import argparse



parser = argparse.ArgumentParser(description='Calculates the frequency of each codon by gene')
parser.add_argument('-F', '--fasta', type=str, required=True, help='fasta archive')
parser.add_argument('-O', '--output', type=str, required=True, help='output file (csv format)')
parser.add_argument('--global', action='store_const', dest='genomic', const=True, default=False, help='Use to calculate the global frequency of codons')
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()


codons_list = [                 # List of every codon
'ATG', 'TTT', 'TTC', 'TTA',
'TTG', 'CTT', 'CTC', 'CTA',
'CTG', 'ATT', 'ATC', 'ATA',
'GTT', 'GTC', 'GTA', 'GTG',

'TCT', 'TCC', 'TCA', 'TCG',
'AGT', 'AGC', 'CCT', 'CCC',
'CCA', 'CCG', 'ACT', 'ACC',
'ACA', 'ACG', 'GCT', 'GCC',

'GCA', 'GCG', 'TAT', 'TAC',
'CAT', 'CAC', 'CAA', 'CAG',
'AAT', 'AAC', 'AAA', 'AAG',
'GAT', 'GAC', 'GAA', 'GAG',

'TGT', 'TGC', 'TGG', 'CGT',
'CGC', 'CGA', 'CGG', 'AGA',
'AGG', 'GGT', 'GGC', 'GGA',
'GGG']

def percent(index, genome_lenght): # Iterface updating the percentage of the ORFs read
	percentage = 100*index/genome_lenght
	interface = int(percentage)
	print str(interface)+'% complete\r',

def evaluate(aminoacids, out_file, gene_name, codons): # Calculates the distribution of the codons
	if len(aminoacids) == 0:
		return 0

	if codons == 0:
		codons = {									# Initial values
		'ATG':0.0, 'TTT':0.0, 'TTC':0.0, 'TTA':0.0,
		'TTG':0.0, 'CTT':0.0, 'CTC':0.0, 'CTA':0.0,
		'CTG':0.0, 'ATT':0.0, 'ATC':0.0, 'ATA':0.0,
		'GTT':0.0, 'GTC':0.0, 'GTA':0.0, 'GTG':0.0,

		'TCT':0.0, 'TCC':0.0, 'TCA':0.0, 'TCG':0.0,
		'AGT':0.0, 'AGC':0.0, 'CCT':0.0, 'CCC':0.0,
		'CCA':0.0, 'CCG':0.0, 'ACT':0.0, 'ACC':0.0,
		'ACA':0.0, 'ACG':0.0, 'GCT':0.0, 'GCC':0.0,

		'GCA':0.0, 'GCG':0.0, 'TAT':0.0, 'TAC':0.0,
		'CAT':0.0, 'CAC':0.0, 'CAA':0.0, 'CAG':0.0,
		'AAT':0.0, 'AAC':0.0, 'AAA':0.0, 'AAG':0.0,
		'GAT':0.0, 'GAC':0.0, 'GAA':0.0, 'GAG':0.0,

		'TGT':0.0, 'TGC':0.0, 'TGG':0.0, 'CGT':0.0,
		'CGC':0.0, 'CGA':0.0, 'CGG':0.0, 'AGA':0.0,
		'AGG':0.0, 'GGT':0.0, 'GGC':0.0, 'GGA':0.0,
		'GGG':0.0
		}
		
	if args.genomic == False:
		out_file.write('\n' + gene_name) # Write header in output file

	for j in range(len(aminoacids)/3):
		codon = aminoacids[j*3] + aminoacids[j*3 + 1] + aminoacids[j*3 + 2] # Current codon being read
		try:
			codons[codon] += 1 # Add 1 to the counts of the current codon
		except:
			if j < len(aminoacids)/3 - 1: # If there is a stop codon, except in the end of the gene
				print 'Stop codon in ' + gene_name

	if args.genomic == False:
		for j in range(61): 
			current_codon = codons_list[j] 
			codon_relative_value = codons[current_codon]*3/(len(aminoacids)-3) # Calculate the codon distribution
			out_file.write(';' + str(codon_relative_value)) # Write the obtained relative values for every codon

	else:
		return codons

def main():
	file_in = open(args.fasta, 'r') # Input file
	in_file = file_in.read()

	out_file = open(args.output, 'w') # Output file

	list_of_genes = str.split(in_file, '>') # Split genes by '>'

	total_frequency = 0 # Define variable

	if args.genomic == False:
		for i in range(61):
			out_file.write(';' + codons_list[i]) # Write header with every codon

	for i in range(len(list_of_genes)):
		percent(i, len(list_of_genes))

		list_of_genes[i] = list_of_genes[i].replace("\r", "") 
		current_gene = str.split(list_of_genes[i], '\n') # Split the data by lines
		gene_name = str.split(current_gene[0], ' ')[0] # The first line is the header

		gene_code = ''
		for j in range(1, len(current_gene)):
			gene_code += current_gene[j] # Every line, except the first, are the code of the gene
		aminoacid_list = list(gene_code) # Make a list of every aminoacid in the gene

		total_frequency = evaluate(aminoacid_list, out_file, gene_name, total_frequency)

	if args.genomic == True:
		total_codons = 0.0

		for j in range(61): 
			current_codon = codons_list[j]
			total_codons += float(total_frequency[current_codon]) # Calculate the total number of codons

		for j in range(61): 
			current_codon = codons_list[j] 
			codon_relative_value = total_frequency[current_codon]/total_codons # Calculate the codon distribution
			out_file.write(current_codon + ';' + str(codon_relative_value) + '\n') # Write the obtained relative values for every codon


	print("100% complete\nDone!")
	sys.exit()

main()