# Konstel(lations)

[![Tests](https://github.com/bede/konstel/actions/workflows/test.yml/badge.svg)](https://github.com/bede/konstel/actions)
[![PyPI](https://badge.fury.io/py/konstel.svg)](https://badge.fury.io/py/konstel)

**Not yet stable, proceed with caution**

An extensible command line tool and library for generating memorable and pronounceable hash-based identifier schemes for sequences, biological or otherwise. For further details and my SARS-CoV-2 naming proposal, please read my [blog post](https://log.bede.im/2021/01/19/covid-hashes).



## Install

```shell
# Python >= 3.6
pip install konstel

# Latest
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