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
$ konstel scheme sars-cov-2-s.genome --file tests/data/spike.genome.fa --output table
Using scheme sars-cov-2-s (genome)
hash                 S:mfcq         
hash-full            S:mfcqn6mh3bnp7vv6eirptvbqik5c65ip
id                   S:sapapag           

```



### Python

```python
>>> from konstel import konstel
>>> konstel.generate_scheme('sars-cov-2-s.genome', file='tests/data/spike.genome.fa')
{'hash': 'S:mfcq', 'hash-full': 'S:mfcqn6mh3bnp7vv6eirptvbqik5c65ip', 'id': 'S:sapapag'}
```