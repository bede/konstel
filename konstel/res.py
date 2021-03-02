# uppercase by default
# remove terminal stop codons

from Bio import SeqIO


alphabets = {
	'dna': 'ACGT',
	'protein': '*ACDEFGHIKLMNPQRSTVWY',
	'dna-ambiguous': '-ABCDGHKMNRSTUVWYZ'
}

# transforms = {
# 	'lowercase': str.lowercase,
# 	'uppercase': str.uppercase
# }



def fasta(path):
	return str(SeqIO.read(path, 'fasta').seq).upper()


file_handlers = {
	'fasta': fasta
}

