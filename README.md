# Konstel(lations)

[![Tests](https://img.shields.io/github/workflow/status/bede/konstel/tests)](https://github.com/bede/konstel/actions)
[![PyPI](https://img.shields.io/pypi/v/konstel.svg?color=brightgreen)](https://badge.fury.io/py/konstel)


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
$ konstel gen sars-cov-2-s.genome konstel/tests/data/spike2.genome.fa --output table
scheme               sars-cov-2-s   
hash                 S:w80qgz2k1fdds6x4mknxazm7psed5knd
hash-4               S:w80q         
id                   S:gofabil  

$ echo "ACGT" | konstel gen generic.nucl - --output table
scheme               generic        
hash                 4449jkgqyv6akzs3aaptjav527dger1m
id                   bodafanoja      
```


### Python

```python
>>> from konstel import konstel
>>> konstel.generate('sars-cov-2-s.protein', 'konstel/tests/data/spike.prot.fa')
{"scheme": "sars-cov-2-s", "hash": "S:c52gdyc7v1dfznny48hfkn1g8ax2yx8f", "hash-4": "S:c52g", "id": "S:dodidib"}
```