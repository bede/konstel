import os
import base64
import hashlib
import binascii

from pathlib import Path

from konstel.res import *

import fire
import strictyaml

from Bio import SeqIO
from Bio.Seq import Seq





def prepare_input(input, spec):
    if spec['remove_whitespace']:
        s = input.translate(str.maketrans('', '', ' \n\t\r'))




def test(scheme, string=None, file=None, format=None):
    library_yaml = '\n'.join([p.read_text() for p in Path('schemes').glob('*.yaml')])
    library = strictyaml.load(library_yaml)
    


    aliases_schemes = {**{k: k for k in library.data},
                       **{v['alias']: k for k, v in library.data.items()}}
    scheme = aliases_schemes[scheme]
    print(f"Using scheme {scheme}")


    for k, v in library[scheme]['input'].items():
        file_formats = library.data[scheme]['input'][k]['formats']
        print(type(file_formats))
        if file and not file_format:
            if len(file_formats == 1):
                file_format = file_formats[0]
            else:
                raise RunTimeError('Unspecified file format')
        if file:
            res.file_formats['fasta'](input)





        # prepare_input(input, library[scheme]['input'][k])

    # print(library[scheme]['input'])







def hash_b10(sequence):
    '''Returns base10 hash from string'''
    h = hashlib.sha1(sequence.encode())
    h_b10 = str(int(binascii.hexlify(h.digest()), 16))
    return h_b10

def hash_b32(sequence):
    '''Returns base32 hash from string'''
    h = hashlib.sha1(sequence.encode())
    h_b32 = base64.b32encode(h.digest()).decode().lower()
    return (h_b32)

def hash_prot_b10(sequence):
    '''Returns base10 hash of string comprising unambiguous IUPAC amino acids'''
    alphabet = set(list('ARNDCQEGHILKMFPSTWYV'))
    sequence_fmt = sequence.upper().strip('*').translate(str.maketrans('', '', ' \n\t\r'))
    assert(sequence_fmt != '' and set(sequence_fmt).issubset(alphabet))
    return hash_b10(sequence_fmt)

def hash_prot_b32(sequence):
    '''Returns base32 hash of string comprising unambiguous IUPAC amino acids'''
    alphabet = set(list('ARNDCQEGHILKMFPSTWYV'))
    sequence_fmt = sequence.upper().strip('*').translate(str.maketrans('', '', ' \n\t\r'))
    assert(sequence_fmt != '' and set(sequence_fmt).issubset(alphabet))
    return hash_b32(sequence_fmt)

def hash_nuc_b10(sequence):
    '''Returns base10 encoded hash of string comprising characters {A,C,G,T,U}'''
    alphabet = set(list('ACGTU'))
    sequence_fmt = sequence.upper().translate(str.maketrans('', '', ' \n\t\r'))
    assert(sequence_fmt != '' and set(sequence_fmt).issubset(alphabet))
    return hash_b10(sequence_fmt)

def hash_nuc_b32(sequence):
    '''Returns base32 encoded hash of string comprising characters {A,C,G,T,U}'''
    alphabet = set(list('ACGTU'))
    sequence_fmt = sequence.upper().translate(str.maketrans('', '', ' \n\t\r'))
    assert(sequence_fmt != '' and set(sequence_fmt).issubset(alphabet))
    return hash_b32(sequence_fmt)

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

def phoneme_generic(sequence):
    '''Returns 8 character phoneme from string'''
    return phoneme_from_hash_b10(hash_b10(sequence))

def phoneme_prot(sequence):
    '''Returns 8 character phoneme of string comprising unambiguous IUPAC amino acids'''
    return phoneme_from_hash_b10(hash_prot_b10(sequence))

def phoneme_nuc(sequence):
    '''Returns 8 character phoneme of string comprising characters {A,C,G,T,U}'''
    return phoneme_from_hash_b10(hash_nuc_b10(sequence))


# ---------- Main library functions ----------


def protein(sequence, hash_length=4):
    '''Returns dict of IDs given a string comprising unambiguous IUPAC amino acids'''
    h_b32 = hash_prot_b32(sequence)
    phoneme = phoneme_prot(sequence)
    return {'phoneme': phoneme,
            'hash': h_b32[:hash_length],
            'hash_full': h_b32,
            'sequence': sequence}

def nucleotide(sequence, hash_length=4):
    '''Returns dict of IDs given a string comprising characters {A,C,G,T,U}'''
    h_b32 = hash_nuc_b32(sequence)
    phoneme = phoneme_nuc(sequence)
    return {'phoneme': phoneme,
            'hash': h_b32[:hash_length],
            'hash_full': h_b32,
            'sequence': sequence}

def generic(sequence, hash_length=4):
    '''Returns dict of IDs given a string'''
    h_b32 = hash_b32(sequence)
    phoneme = phoneme_generic(sequence)
    return {'phoneme': phoneme,
            'hash': h_b32[:hash_length],
            'hash_full': h_b32,
            'sequence': sequence}


# ---------- SARS-CoV-2 functions ----------


def sars2_nuc_to_spike_prot(nuc_sequence):
    '''Returns translated SARS-CoV-2 spike contained in nucleotide string'''
    nuc_sequence_fmt = nuc_sequence.translate(str.maketrans('', '', ' \n\t\r'))
    seq = Seq(nuc_sequence_fmt)
    s_start_seq = 'atgtttgttttt'  # First 4 codons of Wuhan-Hu-1
    s_end_seq = 'ttacattacacataa'  # Last 5 codons of Wuhan-Hu-1
    s_start_pos = str(seq).lower().index(s_start_seq)
    s_end_pos = str(seq).lower().index(s_end_seq) + len(s_end_seq)
    spike_nucl = seq[s_start_pos:s_end_pos]
    spike_prot = str(spike_nucl.ungap().translate()).strip('*')
    assert(1200 < len(spike_prot) < 1300)  # SARS-CoV-2 is ~1270AAs long 
    return spike_prot

def sars2_nuc_fasta_to_spike_prot(fasta_path):
    '''Returns translated SARS-CoV-2 spike contained in nucleotide fasta'''
    record = SeqIO.read(fasta_path, 'fasta')
    spike_prot = sars2_nuc_to_spike_prot(str(record.seq))
    return spike_prot

def sars2_spike_from_nuc(nuc_sequence, hash_length=4):
    '''Returns hash of SARS-CoV-2 spike sequence contained in nucleotide string'''
    prot_sequence = sars2_nuc_to_spike_prot(nuc_sequence)
    h_b32 = hash_prot_b32(prot_sequence)
    phoneme = phoneme_prot(prot_sequence)
    return {
        'phoneme': phoneme,
        'hash': h_b32[:hash_length],
        'hash_full': h_b32,
        'sequence': prot_sequence}

def sars2_spike_from_nuc_fasta(fasta_path, hash_length=4):
    '''Returns dict of IDs of SARS-CoV-2 spike sequence contained in nucleotide fasta'''
    prot_sequence = sars2_nuc_fasta_to_spike_prot(fasta_path)
    h_b32 = hash_prot_b32(prot_sequence)
    phoneme = phoneme_prot(prot_sequence)
    return {
        'phoneme': phoneme,
        'hash': h_b32[:hash_length],
        'hash_full': h_b32,
        'sequence': prot_sequence}


def main():
    fire.Fire({
        'protein': protein,
        'nucleotide': nucleotide,
        'generic': generic,
        'sars2-spike-from-nuc': sars2_spike_from_nuc,
        'sars2-spike-from-nuc-fasta': sars2_spike_from_nuc_fasta,
        'test': test})