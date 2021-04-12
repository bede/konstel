# Konstel(lations)

[![Tests](https://img.shields.io/github/workflow/status/bede/konstel/tests)](https://github.com/bede/konstel/actions)
[![PyPI](https://img.shields.io/pypi/v/konstel.svg?color=brightgreen)](https://badge.fury.io/py/konstel)

**Not yet stable, proceed with caution**

An extensible command line tool and library for generating memorable and pronounceable hash-based identifier schemes for sequences, biological or otherwise. For further details and my SARS-CoV-2 naming proposal, please read my [blog post](https://log.bede.im/2021/01/19/covid-hashes). Requires Python 3.6+.

### SARS-CoV-2 naming

Phonemic and truncated cbase32 identifiers provide 36 and 40 bits of entropy respectively, producing no collisions within publicly deposited SARS-CoV-2 spike protein sequences as of 2021-04-12.

## Install

Ideally inside a new virtualenv or conda environment:

```shell
# Latest release
pip install konstel

# Development version
git clone https://github.com/bede/konstel
pip install --editable konstel
```


## Usage 

### Command line

```bash
$ konstel gen sars-cov-2-s.genome konstel/tests/data/spike.genome.fa --output table
scheme               sars-cov-2-s   
hash                 S:0k8n9hjh5xh5kbef1k6ye7e2d4brhpry5r985avrtf69v6amrbc0
hash-8               S:0k8n9hjh     
id                   S:huhiji-gakihi  

$ echo "ACGT" | konstel gen generic.nucl - --output table
scheme               generic        
hash                 3qzkx17yf1vy0ssvd6xxvkt02973jvhzk51xv28cj6va16pvkbr0
id                   bituzu-gupahu-zolodu-lumaki-suripi-rozitu-guhabi-figogo
```


### Python

```python
>>> from konstel import konstel
>>> konstel.generate('sars-cov-2-s.protein', 'konstel/tests/data/spike.prot.fa')
{'scheme': 'sars-cov-2-s', 'hash': 'S:0k8n9hjh5xh5kbef1k6ye7e2d4brhpry5r985avrtf69v6amrbc0', 'hash-8': 'S:0k8n9hjh', 'id': 'S:huhiji-gakihi'}
```

