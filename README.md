# Konstel(lations)

CLI and Python library for creating short hash-based identifiers (e.g. `S:mfcq` or `S:papoheme`) for distinct SARS-CoV-2 spike proteins, or indeed any protein sequence, nucleotide sequence or generic string of interest. For further details and my SARS-CoV-2 naming proposal, please read my [blog post](https://log.bede.im/2021/01/19/covid-hashes).



## Install

```shell
# Python >= 3.6
pip install git+https://github.com/bede/konstel
```



## Usage 

### Command line

```
$ konstel protein PEACE
phoneme:   hosupuho
hash:      6gtq
hash_full: 6gtqivg444l3uvf6cy43pxyyjmvg3wpp
sequence:  PEACE

$ konstel sars2-spike-from-nuc-fasta genome.fa
phoneme:   papoheme
hash:      mfcq
hash_full: mfcqn6mh3bnp7vv6eirptvbqik5c65ip
sequence:  <truncated>
```



### Python

```python
from konstel import konstel

konstel.protein('PEACE')
{'phoneme': 'hosupuho', 'hash': '6gtq', 'hash_full': '6gtqivg444l3uvf6cy43pxyyjmvg3wpp', 'sequence': 'PEACE'}

konstel.sars2_spike_from_nuc_fasta('tests/data/test.fa')
{'phoneme': 'papoheme', 'hash': 'mfcq', 'hash_full': 'mfcqn6mh3bnp7vv6eirptvbqik5c65ip', 'sequence': '<truncated>'}
```