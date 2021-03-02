# uppercase by default
# remove terminal stop codons
import io
import hashlib

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



def fasta(string):
	string_as_file = io.StringIO(string)
	record = SeqIO.read(string_as_file, 'fasta')
	assert len(record.seq) >= 1
	return str(record.seq).upper()


format_funcs = {
	'fasta': fasta
}

hash_funcs = {
	'sha1': hashlib.sha1,
	'sha256': hashlib.sha256
}