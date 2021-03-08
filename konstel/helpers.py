from Bio import SeqIO
from Bio.Seq import Seq

def sars_cov_2_s_b10_genome(nuc_sequence):
    '''Returns translated SARS-CoV-2 spike contained in nucleotide string'''
    nuc_sequence_fmt = nuc_sequence.translate(str.maketrans('', '', ' \n\t\r'))
    seq = Seq(nuc_sequence_fmt)
    s_start_seq = 'atgtttgttttt'  # First 4 codons of Wuhan-Hu-1
    s_end_seq = 'ttacattacacataa'  # Last 5 codons of Wuhan-Hu-1
    s_start_pos = str(seq).lower().index(s_start_seq)
    s_end_pos = str(seq).lower().index(s_end_seq) + len(s_end_seq)
    spike_nucl = seq[s_start_pos:s_end_pos]
    spike_prot = str(spike_nucl.ungap().translate()).strip('*')
    return spike_prot