import base64
import hashlib

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
    return h_b32

def hash_prot(sequence):
    '''Returns encoded hash of string comprising only unambiguous IUPAC amino acids'''
    alphabet = set(list('ARNDCQEGHILKMFPSTWYV'))
    sequence_fmt = sequence.upper().strip('*').translate(str.maketrans('', '', ' \n\t\r'))
    assert(set(sequence_fmt).issubset(alphabet))
    h = hashlib.sha1(sequence_fmt.encode())
    h_b32 = base64.b32encode(h.digest()).decode().lower()
    return h_b32


def sars2_nuc_to_spike_prot(nt_sequence):
    '''Returns translated SARS-CoV-2 spike contained in DNA string'''
    seq = Seq(nt_sequence)
    s_start_seq = 'atgtttgtttttctt'
    s_end_seq = 'ttacattacacataa'
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
    '''Returns hash of SARS-CoV-2 spike sequence contained in DNA string'''
    return hash_prot(sars2_nuc_to_spike_prot(nt_sequence))


def sars2_nuc_fasta_to_spike_hash(fasta_path):
    '''Returns hash of SARS-CoV-2 spike sequence contained in DNA fasta'''
    return hash_prot(sars2_nuc_fasta_to_spike_prot(fasta_path))



def main():
    fire.Fire({
        'hash-dna': hash_dna,
        'hash-prot': hash_prot,
        'sars2-nuc-to-spike-prot': sars2_nuc_to_spike_prot,
        'sars2-nuc-fasta-to-spike-prot': sars2_nuc_fasta_to_spike_prot,
        'sars2-nuc-to-spike-hash': sars2_nuc_to_spike_hash,
        'sars2-nuc-fasta-to-spike-hash': sars2_nuc_fasta_to_spike_hash
    })