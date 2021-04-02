import io
import sys
import pytest
import subprocess

from konstel import konstel

data_dir = 'konstel/tests/data/'

def run(cmd, cwd='./'):  # Helper for CLI testing
    return subprocess.run(cmd,
                          cwd=cwd,
                          shell=True,
                          check=True,
                          universal_newlines=True,)


def test_sars2_generate_prot():
    result = konstel.generate('sars-cov-2-s.protein', f'{data_dir}/spike.prot.fa')
    print(result)
    assert result['id'] == 'S:sapapag'
    assert result['id-legacy'] == 'S:papoheme'


def test_sars2_generate_nucl():
    result = konstel.generate('sars-cov-2-s.genome', f'{data_dir}/spike.genome.fa')
    print(result)
    assert result['id'] == 'S:sapapag'
    assert result['id-legacy'] == 'S:papoheme'


def test_sars2_generate_stdin():
    with open(f'{data_dir}/spike.genome.fa') as sys.stdin:
        result = konstel.generate('sars-cov-2-s.genome', file='-')
    print(result)
    assert result['id'] == 'S:sapapag'


def test_validation_alphabet():
    sys.stdin = io.StringIO('ACGTX')
    with pytest.raises(RuntimeError):
        konstel.generate('generic.nucl', file='-')
