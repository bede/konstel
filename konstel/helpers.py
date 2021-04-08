import os
import parasail

from Bio import SeqIO
from Bio.Seq import Seq


PACKAGE_PATH = os.path.dirname(__file__)


def sars_cov_2_s_genome(genome):
    '''Returns translated SARS-CoV-2 spike sequence contained in nucleotide string'''
    try:
        ref_gene = str(SeqIO.read(
            f'{PACKAGE_PATH}/schemes/sars-cov-2-s/assets/MN908947.3.S.fa', 'fasta').seq)
        result = parasail.sg_qx_trace(ref_gene, genome, 10, 1, parasail.dnafull)
        aln_start = result.traceback.query.rfind(ref_gene)
        aln_end = aln_start + len(ref_gene)
        gene = Seq(result.traceback.ref[aln_start:aln_end])
        protein = str(gene.ungap().translate()).strip('*')
    except Exception as e:
        raise RuntimeError('Problem extracting spike sequence').with_traceback(e.__traceback__)
    return protein


def sars_cov_2_s_legacy_genome(nuc_sequence):
    '''Returns translated SARS-CoV-2 spike sequence contained in nucleotide string'''
    nuc_sequence_fmt = nuc_sequence.translate(str.maketrans('', '', ' \n\t\r'))
    seq = Seq(nuc_sequence_fmt)
    s_start_seq = 'atgtttgttttt'  # First 4 codons of Wuhan-Hu-1
    s_end_seq = 'ttacattacacataa'  # Last 5 codons of Wuhan-Hu-1
    try:
        s_start_pos = str(seq).lower().index(s_start_seq)
        s_end_pos = str(seq).lower().index(s_end_seq) + len(s_end_seq)
        spike_nucl = seq[s_start_pos:s_end_pos]
        spike_prot = str(spike_nucl.ungap().translate()).strip('*')
    except Exception as e:
        raise RuntimeError('Problem extracting spike sequence').with_traceback(e.__traceback__)
    return spike_prot
