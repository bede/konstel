import io

from Bio import SeqIO


def fasta(string):
	string_as_file = io.StringIO(string)
	record = SeqIO.read(string_as_file, 'fasta')
	assert len(record.seq) >= 1, 'Sequence is empty'
	return str(record.seq).upper()
