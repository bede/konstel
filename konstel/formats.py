import io

from Bio import SeqIO


def fasta(string):
	string_as_file = io.StringIO(string)
	record = SeqIO.read(string_as_file, 'fasta')
	assert len(record.seq) >= 1
	return str(record.seq).upper()

def phoneme_from_hash_b10(h_b10):
    '''Returns phoneme from base10 hash'''
    vowel_map = dict(zip(map(str, range(10)), 'aeiou'*2))
    consonant_map = dict(zip(map(str, range(10)), 'fhklmprstv'))
    phoneme = ''
    for char_i, char in enumerate(h_b10[:8]):
        if not char_i % 2:
            phoneme += consonant_map[char]
        else:
            phoneme += vowel_map[char]
    return phoneme