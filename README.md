# Konstel(lations)



**Currently under active development, proceed with caution**

CLI and Python library for creating short hash-based identifiers (e.g. `S:mfcq` or `S:papoheme`) for distinct SARS-CoV-2 spike proteins, or indeed any protein sequence, nucleotide sequence or generic string of interest. For further details and my SARS-CoV-2 naming proposal, please read my [blog post](https://log.bede.im/2021/01/19/covid-hashes).



## Install

```shell
# Python >= 3.6
pip install git+https://github.com/bede/konstel
```



## Usage 

### Command line

```bash
$ konstel gen sars-cov-2-s.genome tests/data/spike.genome.fa --output table
scheme: sars-cov-2-s (genome), input: tests/data/spike.genome.fa
hash                 S:mfcq         
hash-full            S:mfcqn6mh3bnp7vv6eirptvbqik5c65ip
id                   S:sapapag      

$ echo "ACGT" | konstel gen generic.nucl - --output table
scheme: generic (prot), input: stdin
hash                 4449jk         
hash-full            4449jkgqyv6akzs3aaptjav527dger1m
id                   mibofi    
```



### Python

```python
>>> from konstel import konstel
>>> konstel.generate('sars-cov-2-s.protein', 'tests/data/spike.prot.fa')
scheme: sars-cov-2-s (protein), input: tests/data/spike.prot.fa
{'hash': 'S:mfcq', 'hash-full': 'S:mfcqn6mh3bnp7vv6eirptvbqik5c65ip', 'id': 'S:sapapag'}
```