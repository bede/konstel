import io
import sys
import json
import subprocess

import pytest

from konstel import konstel


data_dir = 'konstel/tests/data'


def run(cmd, cwd='./'):  # Helper for CLI testing
    return subprocess.run(cmd,
                          cwd=data_dir,
                          shell=True,
                          check=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def test_validation_alphabet():
    sys.stdin = io.StringIO('ACGTX')
    with pytest.raises(RuntimeError):
        konstel.generate('generic.nucl', file='-')

def test_sars2_generate_prot():
    result = konstel.generate('sars-cov-2-s.protein', f'{data_dir}/spike.prot.fa')
    assert result['id'] == 'S:huhijig-akihi'

def test_sars2_generate_nucl():
    result = konstel.generate('sars-cov-2-s.genome', f'{data_dir}/spike.genome.fa')
    assert result['id'] == 'S:huhijig-akihi'

def test_sars2_generate_length():
    result = konstel.generate('sars-cov-2-s.protein', f'{data_dir}/spike.prot.fa', length=21)
    assert result['id'] == 'S:huhijig-akihifu-tofikip'

def test_sars2_generate_stdin():
    with open(f'{data_dir}/spike.genome.fa') as sys.stdin:
        result = konstel.generate('sars-cov-2-s.genome', file='-')
    assert result['id'] == 'S:huhijig-akihi'

def test_sars2_regen():
    result = konstel.regenerate('sars-cov-2-s', '0k8n9hjh5xh5kbef1k6ye7e2d4brhpry5r985avrtf69v6amrbc0')
    assert result['id'] == 'S:huhijig-akihi'

def test_sars2_regen_directive():
    result = konstel.regenerate('sars-cov-2-s.protein', '0k8n9hjh5xh5kbef1k6ye7e2d4brhpry5r985avrtf69v6amrbc0')
    assert result['id'] == 'S:huhijig-akihi'

def test_sars2_regen_prefix():
    result = konstel.regenerate('sars-cov-2-s', 'S:0k8n9hjh5xh5kbef1k6ye7e2d4brhpry5r985avrtf69v6amrbc0')
    assert result['id'] == 'S:huhijig-akihi'

def test_sars2_legacy_regen():
    result = konstel.regenerate('sars-cov-2-s-legacy', 'mfcqn6mh3bnp7vv6eirptvbqik5c65ip')
    assert result['id'] == 'S:papoheme'

def test_sars2_legacy_regen_prefix():
    result = konstel.regenerate('sars-cov-2-s-legacy', 'S:mfcqn6mh3bnp7vv6eirptvbqik5c65ip')
    assert result['id'] == 'S:papoheme'

def test_sars2_regen_length():
    result = konstel.regenerate('sars-cov-2-s', 'S:0k8n9hjh5xh5kbef1k6ye7e2d4brhpry5r985avrtf69v6amrbc0', length=21)
    assert result['id'] == 'S:huhijig-akihifu-tofikip'

def test_sars2_legacy_protein():
    result = konstel.generate('sars-cov-2-s-legacy.protein', f'{data_dir}/spike.prot.fa')
    assert result['id'] == 'S:papoheme'

def test_sars2_legacy_genome():
    result = konstel.generate('sars-cov-2-s-legacy.genome', f'{data_dir}/spike.genome.fa')
    assert result['id'] == 'S:papoheme'

def test_stdin_generic_separator():
    sys.stdin = io.StringIO('ACGT')
    result = konstel.generate('generic.nucl', file='-')
    assert result['id'] == 'bituzu-gupahu-zolodu-lumaki-suripi-rozitu-guhabi-figogo'



# def test_cli_json():
#     cmd_run = run(f'konstel gen sars-cov-2-s.protein spike.prot.fa')
#     assert json.loads(cmd_run.stdout) == {"scheme": "sars-cov-2-s", "hash": "S:c52gdyc7v1dfznny48hfkn1g8ax2yx8f", "hash-4": "S:c52g", "id": "S:dodidib"}
