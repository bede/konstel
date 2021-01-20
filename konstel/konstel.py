import base64
import hashlib
import binascii

import fire

from Bio import SeqIO
from Bio.Seq import Seq


def hash_dna(sequence):
    '''Returns encoded hash of string comprising only characters {A,C,G,T,U}'''
    alphabet = set(list('ACGTU'))
    sequence_fmt = sequence.upper().translate(str.maketrans('', '', ' \n\t\r'))
    assert(set(sequence_fmt).issubset(alphabet))
    h = hashlib.sha1(sequence_fmt.encode())
    h_b32 = base64.b32encode(h.digest()).decode().lower()
    return (h_b32)

def hash_prot_b32(sequence):
    '''Returns encoded hash of string comprising only unambiguous IUPAC amino acids'''
    alphabet = set(list('ARNDCQEGHILKMFPSTWYV'))
    sequence_fmt = sequence.upper().strip('*').translate(str.maketrans('', '', ' \n\t\r'))
    assert(set(sequence_fmt).issubset(alphabet))
    h = hashlib.sha1(sequence_fmt.encode())
    h_b32 = base64.b32encode(h.digest()).decode().lower()
    return h_b32

def hash_prot_b10(sequence):
    '''Returns decimal hash of string comprising only unambiguous IUPAC amino acids'''
    alphabet = set(list('ARNDCQEGHILKMFPSTWYV'))
    sequence_fmt = sequence.upper().strip('*').translate(str.maketrans('', '', ' \n\t\r'))
    assert(set(sequence_fmt).issubset(alphabet))
    h = hashlib.sha1(sequence_fmt.encode())
    h_b10 = str(int(binascii.hexlify(h.digest()), 16))
    return h_b10

def prot_to_phoneme(sequence):
    h_b10 = hash_prot_b10(sequence)
    vowel_map = dict(zip(map(str, range(10)), 'aeiou'*2))
    consonant_map = dict(zip(map(str, range(10)), 'fhklmprstv'))
    phoneme = ''
    for char_i, char in enumerate(h_b10[:8]):
        if not char_i % 2:
            phoneme += consonant_map[char]
        else:
            phoneme += vowel_map[char]
    return phoneme


def sars2_nuc_to_spike_prot(nt_sequence):
    '''Returns translated SARS-CoV-2 spike contained in DNA string'''
    nt_sequence_fmt = nt_sequence.translate(str.maketrans('', '', ' \n\t\r'))
    seq = Seq(nt_sequence_fmt)
    s_start_seq = 'atgtttgtttttctt'  # First 5 codons of Wuhan-Hu-1
    s_end_seq = 'ttacattacacataa'  # Last 5 codons of Wuhan-Hu-1
    s_start_pos = str(seq).lower().index(s_start_seq)
    s_end_pos = str(seq).lower().index(s_end_seq) + len(s_end_seq)
    spike_nt = seq[s_start_pos:s_end_pos]
    spike_aa = str(spike_nt.ungap().translate()).strip('*')
    return spike_aa

def sars2_nuc_fasta_to_spike_prot(fasta_path):
    '''Returns translated SARS-CoV-2 spike contained in DNA fasta'''
    record = SeqIO.read(fasta_path, 'fasta')
    spike_aa = sars2_nuc_to_spike_prot(str(record.seq))
    return spike_aa

def sars2_nuc_to_spike_hash(nt_sequence):
    '''Returns qualified hash of SARS-CoV-2 spike sequence contained in DNA string'''
    return 'S:' + hash_prot_b32(sars2_nuc_to_spike_prot(nt_sequence))

def sars2_nuc_fasta_to_spike_hash(fasta_path):
    '''Returns qualified hash of SARS-CoV-2 spike sequence contained in DNA fasta'''
    return 'S:' + hash_prot_b32(sars2_nuc_fasta_to_spike_prot(fasta_path))


def main():
    fire.Fire({
        'hash-dna': hash_dna,
        'hash-prot': hash_prot_b32,
        'sars2-nuc-to-spike-prot': sars2_nuc_to_spike_prot,
        'sars2-nuc-fasta-to-spike-prot': sars2_nuc_fasta_to_spike_prot,
        'sars2-nuc-to-spike-hash': sars2_nuc_to_spike_hash,
        'sars2-nuc-fasta-to-spike-hash': sars2_nuc_fasta_to_spike_hash
    })