import io
import sys
import pytest
import subprocess

from konstel import konstel



def run(cmd, cwd='./'):  # Helper for CLI testing
    return subprocess.run(cmd,
                          cwd=cwd,
                          shell=True,
                          check=True,
                          universal_newlines=True,)


def test_sars2_generate_prot():
    result = konstel.generate('sars-cov-2-s.protein', 'tests/data/spike.prot.fa')
    print(result)
    assert result['id'] == 'S:sapapag'
    assert result['id-legacy'] == 'S:papoheme'


def test_sars2_generate_nucl():
    result = konstel.generate('sars-cov-2-s.genome', 'tests/data/spike.genome.fa')
    print(result)
    assert result['id'] == 'S:sapapag'
    assert result['id-legacy'] == 'S:papoheme'


def test_sars2_generate_stdin():
    with open('tests/data/spike.genome.fa') as sys.stdin:
        result = konstel.generate('sars-cov-2-s.genome', file='-')
    print(result)
    assert result['id'] == 'S:sapapag'

    run('echo CAT | konstel gen generic.nucl -')


def test_validation_alphabet():
    sys.stdin = io.StringIO('ACGTX')
    with pytest.raises(RuntimeError):
        konstel.generate('generic.nucl', file='-')
