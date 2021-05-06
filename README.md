# Konstel(lations)

[![Tests](https://img.shields.io/github/workflow/status/bede/konstel/tests)](https://github.com/bede/konstel/actions)
[![PyPI](https://img.shields.io/pypi/v/konstel.svg?color=brightgreen)](https://badge.fury.io/py/konstel)

An extensible command line tool and library for generating memorable and pronounceable hash-based identifiers for biological sequences and arbitrary strings. Requires Python 3.6+. A work-in-progress web app is hosted on App Engine: https://konstel.ew.r.appspot.com/

## Rationale

The challenges of biological taxonomy and nomenclature have in some respects become more acute in the era of abundantly available genome sequences, whose phylogenies continually reveal shortcomings in established taxonomy. A linked challenge lies in the sheer quantity of entities we wish to name, all while favouring names that are memorable, pronounceable and coherently organised. Nomenclature has received unusually widespread attention during the COVID-19 pandemic, with repeated independent emergence of variants of concern creating tension between lineage-focused and variant-focused naming for SARS-CoV-2. In this case, delays and uncertainty surrounding naming have led to dissemination of stigmatising variant names by media outlets and even public health agencies. A simple and arguably complementary approach to manually curated nomenclature in general is the use of random identifiers. If generated using an appropriate hash function, a distinct random identifier may be assigned to any given biological sequence. By subsequently transforming hash function output into pronouceable combinations of English letters, it is possible to generate information-dense identifiers that are pronouceable, memorable and unambiguous. This tool generalises the encoding of hashed sequences or arbitrary strings into pronounceable identifiers using customisable scheme definitions. In circumstances where sharing of the sequence data itself may be restricted—namely SARS-CoV-2—hash-based identifiers offer a convenient way to refer umabiguously to a sequence without breaching data access agreements, and konstel includes a scheme for generating identifiers linked to SARS-CoV-2 spike protein sequences given either a protein sequence or a genome sequence encoding the complete spike sequence.



### SARS-CoV-2 spike protein naming

Phonemic and truncated cbase32 identifiers provide 34 and 35 bits of entropy respectively, collision-free for publicly deposited SARS-CoV-2 spike protein sequences as of 2021-05-05. Phonemic identifiers include consonant-vowel pairs with a separator after the fourth consonant (e.g. `dazator-isak`). The first segment provides a useful compromise of identifier length and low collision rate, while inclusion of the second segment achieves collision resistance. Longer identifiers still may be minted by overriding the scheme's default length profile. For my original SARS-CoV-2 naming proposal, please refer to my [blog post](https://log.bede.im/2021/01/19/covid-hashes).

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

Arbitrary nucleotide sequence:

```bash
$ echo "ACGT" | konstel gen bio.nuc - --output table  # Use '-' to accept stdin 
scheme               bio        
hash                 3qzkx17yf1vy0ssvd6xxvkt02973jvhzk51xv28cj6va16pvkbr0
id                   bituzu-gupahu-zolodu-lumaki-suripi-rozitu-guhabi-figogo
```

Arbitrary protein sequence:

```bash
$ echo "TASTE" | konstel gen bio.pro - --output table
scheme               bio            
hash                 nr8npewt0bwamk8s3xwhgxnd47zn6rxsjtrjm4b3eqvvp40rp5g0
id                   fovahi-josuro-kobaru-mopohu-hinalu-lohimi-topuho-duzuja
```

SARS-CoV-2 spike protein sequence:

```bash
$ konstel gen sars-cov-2-s.protein konstel/tests/data/spike.prot.fa --output table
scheme               sars-cov-2-s   
hash                 S:0k8n9hjh5xh5kbef1k6ye7e2d4brhpry5r985avrtf69v6amrbc0
hash-7               S:0k8n9hj      
id                   S:huhijig-akih 
```

SARS-CoV-2 genome sequence (containing spike protein sequence)

```bash
$ konstel gen sars-cov-2-s.genome konstel/tests/data/spike.genome.fa --output table
scheme               sars-cov-2-s   
hash                 S:0k8n9hjh5xh5kbef1k6ye7e2d4brhpry5r985avrtf69v6amrbc0
hash-8               S:0k8n9hjh     
id                   S:huhijig-akih  
```



### Python

```python
>>> from konstel import konstel
>>> konstel.generate('sars-cov-2-s.protein', 'konstel/tests/data/spike.prot.fa')
{'scheme': 'sars-cov-2-s', 'hash': 'S:0k8n9hjh5xh5kbef1k6ye7e2d4brhpry5r985avrtf69v6amrbc0', 'hash-7': 'S:0k8n9hj', 'id': 'S:huhijig-akih'}
```

