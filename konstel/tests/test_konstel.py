import io
import sys
import pytest
import subprocess

from konstel import konstel


data_dir = 'konstel/tests/data'


def run(cmd, cwd='./'):  # Helper for CLI testing
    return subprocess.run(cmd,
                          cwd=cwd,
                          shell=True,
                          check=True,
                          universal_newlines=True,)


def test_validation_alphabet():
    sys.stdin = io.StringIO('ACGTX')
    with pytest.raises(RuntimeError):
        konstel.generate('generic.nucl', file='-')


def test_sars2_generate_prot():
    result = konstel.generate('sars-cov-2-s.protein', f'{data_dir}/spike.prot.fa')
    print(result)
    assert result['id'] == 'S:dodidib'


def test_sars2_generate_nucl():
    result = konstel.generate('sars-cov-2-s.genome', f'{data_dir}/spike.genome.fa')
    print(result)
    assert result['id'] == 'S:dodidib'


def test_sars2_generate_stdin():
    with open(f'{data_dir}/spike.genome.fa') as sys.stdin:
        result = konstel.generate('sars-cov-2-s.genome', file='-')
    print(result)
    assert result['id'] == 'S:dodidib'


def test_sars2_legacy():
    result = konstel.generate('sars-cov-2-s-legacy.protein', f'{data_dir}/spike.prot.fa')
    assert result['id'] == 'S:papoheme'

    result = konstel.generate('sars-cov-2-s-legacy.genome', f'{data_dir}/spike.genome.fa')
    assert result['id'] == 'S:papoheme'