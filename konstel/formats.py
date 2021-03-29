import io

from Bio import SeqIO



alphabets = {
	'nucleotide': 'ACGTU',
	'nucleotide-ambiguous': '-ABCDGHKMNRSTUVWY',
	'protein': '*ACDEFGHIKLMNPQRSTVWY'
}



def fasta(string):
	'''
	Returns uppercased sequence from fasta string
	Prepends '>' if absent to enable naked sequence parsing
	'''
	if '>' not in string:
		string = f'>\n{string}'
	string_as_file = io.StringIO(string)
	record = SeqIO.read(string_as_file, 'fasta')
	assert len(record.seq) >= 1, 'Sequence is empty'
	return str(record.seq).upper()