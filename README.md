# Konstel(lations)

[![Tests](https://img.shields.io/github/workflow/status/bede/konstel/tests)](https://github.com/bede/konstel/actions)
[![PyPI](https://img.shields.io/pypi/v/konstel.svg?color=brightgreen)](https://badge.fury.io/py/konstel)

An extensible command line tool and library for generating memorable and pronounceable hash-based identifier schemes for sequences, biological or otherwise. Requires Python 3.6+.

### SARS-CoV-2 spike protein naming

Phonemic and truncated cbase32 identifiers provide 36 and 40 bits of entropy respectively, collision-free for publicly deposited SARS-CoV-2 spike protein sequences as of 2021-04-12. Phonemic identifiers include six consonant-vowel pairs with a separator after the fourth consonant (e.g. `dazator-isaki`). The first segment provides an empirically established useful compromise of identifier length and low collision rate, while inclusion of the second segment achieves collision resistance. Longer identifiers still may be minted by overriding the scheme's default length profile. For my original SARS-CoV-2 naming proposal, please refer to my [blog post](https://log.bede.im/2021/01/19/covid-hashes).

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

