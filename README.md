# Konstel(lations)

[![Tests](https://img.shields.io/github/workflow/status/bede/konstel/tests)](https://github.com/bede/konstel/actions)
[![PyPI](https://img.shields.io/pypi/v/konstel.svg?color=brightgreen)](https://badge.fury.io/py/konstel)

An extensible command line tool, library, and accompanying [web app](https://konstel.ew.r.appspot.com) for generating memorable and pronounceable hash-based identifiers for arbitrary sequences. Konstel normalises and hashes a given string, biological sequence or SARS-CoV-2 sequence before encoding the hash digest as a human-friendly phonemic word. This allows privacy preserving confirmation of input equality and can be thought of as the nomenclature equivalent of URL shortening. In the context of an emerging infectious disease outbreak like SARS-CoV-2, this approach can alleviate some of the challenges imposed restrictive data access agreements – I can tell that my sequence is the same as yours without you having to share it with me. Requires Python 3.6+. Presented at ABPHM '21.

### SARS-CoV-2 spike protein naming

Phonemic and truncated cbase32 identifiers provide 34 and 35 bits of entropy respectively, collision-free for publicly deposited SARS-CoV-2 spike protein sequences as of 2021-05-05. Phonemic identifiers include consonant-vowel pairs with a separator after the fourth consonant (e.g. `dazator-isak`). The first segment provides a useful compromise of identifier length and low collision rate, while inclusion of the second segment achieves collision resistance. Longer identifiers still may be minted by overriding the scheme's default length profile. For a discussion of SARS-CoV-2 naming, please refer to my [blog post](https://log.bede.im/2021/01/19/covid-hashes).

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

Arbitrary string

```bash
$ echo "test" | konstel gen string -  # Uses '-' to accept stdin
{"scheme": "string", "hash": "ky3d10c89hypb6hfxagcappg2phvykrv5c5r4b6hbnp1bc7g1840", "id": "fizodo-tafado-fahudu-tinino-mozupo-pagaji-kotabi"}

$ echo "test" | konstel gen string --length 6 --output table -  # Custom length, tabular output
scheme               string         
hash                 ky3d10c89hypb6hfxagcappg2phvykrv5c5r4b6hbnp1bc7g1840
id                   fizodo   
```

Arbitrary nucleotide sequence (alphabet `-ACGTU`):

```bash
$ konstel gen bio.nuc acgt.fa --output table  # Fasta containing ACGT
scheme               bio            
hash                 3qzkx17yf1vy0ssvd6xxvkt02973jvhzk51xv28cj6va16pvkbr0
id                   bituzu-gupahu-zolodu-lumaki-suripi-rozitu-guhabi

$ echo "ACGT" | konstel gen bio.nuc --output table -  # ACGT as stdin
scheme               bio        
hash                 3qzkx17yf1vy0ssvd6xxvkt02973jvhzk51xv28cj6va16pvkbr0
id                   bituzu-gupahu-zolodu-lumaki-suripi-rozitu-guhabi
```
Ambiguous arbitrary nucleotide sequence (alphabet `-ABCDGHKMNRSTUVWY`):
```bash
$ konstel gen bio.nuc-ambiguous acgtn.fa --output table  # Fasta containing ACGTN
scheme               bio            
hash                 t9a5abnf4nwtmbpb59b477218wqrwzf0hasz2qm9gw2ynpkpzgpg
id                   gifija-jihovo-rufiju-nopofu-rarapo-jinago-lahaja
```

Arbitrary protein sequence (alphabet `*-ACDEFGHIKLMNPQRSTVWY`):

```bash
$ konstel gen bio.pro taste.fa --output table
scheme               bio            
hash                 nr8npewt0bwamk8s3xwhgxnd47zn6rxsjtrjm4b3eqvvp40rp5g0
id                   fovahi-josuro-kobaru-mopohu-hinalu-lohimi-topuho
```

SARS-CoV-2 spike protein sequence:

```bash
$ konstel gen sars-cov-2-s.protein spike.prot.fa --output table
scheme               sars-cov-2-s   
hash                 S:0k8n9hjh5xh5kbef1k6ye7e2d4brhpry5r985avrtf69v6amrbc0
hash-7               S:0k8n9hj      
id                   S:huhijig-akih 
```

SARS-CoV-2 genome sequence (containing complete spike protein sequence)

```bash
$ konstel gen sars-cov-2-s.genome spike.genome.fa --output table
scheme               sars-cov-2-s   
hash                 S:0k8n9hjh5xh5kbef1k6ye7e2d4brhpry5r985avrtf69v6amrbc0
hash-7               S:0k8n9hj     
id                   S:huhijig-akih  
```

### Python

```python
>>> from konstel import konstel
```

Arbitrary string:

```python
>>> konstel.generate('string', '-', sequence='test')  # From string
{'scheme': 'string', 'hash': 'ky3d10c89hypb6hfxagcappg2phvykrv5c5r4b6hbnp1bc7g1840', 'id': 'fizodo-tafado-fahudu-tinino-mozupo-pagaji-kotabi'}
```

Arbitrary nucleotide sequence, length 6:

```python
>>> konstel.generate('bio.nuc', '-', sequence='ACGT', length=6)  # From string
{'scheme': 'bio', 'hash': '3qzkx17yf1vy0ssvd6xxvkt02973jvhzk51xv28cj6va16pvkbr0', 'id': 'bituzu'}
```

SARS-CoV-2 spike protein sequence:

```python
>>> konstel.generate('sars-cov-2-s.protein', 'spike.prot.fa')  # From fasta file
{'scheme': 'sars-cov-2-s', 'hash': 'S:0k8n9hjh5xh5kbef1k6ye7e2d4brhpry5r985avrtf69v6amrbc0', 'hash-7': 'S:0k8n9hj', 'id': 'S:huhijig-akih'}
```

## Issues

Please let me know if you run into problems by opening a GitHub issue, tweeting @beconstant or mailing me via b at bede dawt im.

## Contributing

If you would like to contribute to this project, please open an issue or contact the author directly using the details above.

Before issuing a pull request, please:

- Ensure tests pass by executing pytest inside the package directory (requires pytest package)
- Increment the version number inside `__init__.py` (Using SemVer as a guide)
- Update documentation and/or tests if possible

