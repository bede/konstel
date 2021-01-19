# Konstel(lations)

CLI and Python library underlying the [web service](https://konstel.ew.r.appspot.com/) for creating hash-based identifiers (e.g. `S:mfcq`) for distinct SARS-CoV-2 spike protein sequences, or indeed any protein or DNA sequence of interest. For my proposal and rationale, please read my [blog post](https://log.bede.im/2021/01/19/covid-hashes).



## Install

```shell
# Python >= 3.6
pip install git+https://github.com/bede/konstel
```



## Usage 

### Command line

```
$ konstel sars2-nuc-fasta-to-spike-hash genome.fasta
S:mfcqn6mh3bnp7vv6eirptvbqik5c65ip
```

### Python

```python
import konstel
print(konstel.sars2_nuc_fasta_to_spike_hash('genome.fasta'))
>> "S:mfcqn6mh3bnp7vv6eirptvbqik5c65ip"
```