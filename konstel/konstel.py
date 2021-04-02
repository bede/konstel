import os
import sys
import pathlib
import hashlib

import argh
from argh.decorators import named

from Bio import SeqIO
from Bio.Seq import Seq

import konstel.schema as schema
import konstel.formats as formats
import konstel.helpers as helpers
import konstel.encodings as encodings



def prepare(string, spec):
    ''''''
    if spec['remove_whitespace']:
        string = string.translate(str.maketrans('', '', ' \n\t\r'))
    if any(spec['strip_characters']):  # any() since schema default [''] is truthy
        string = string.strip(''.join(spec['strip_characters']))
    return string


def validate(string, spec):
    ''''''
    if 'alphabet' in spec:
        observed = set(string)
        permitted = set(formats.alphabets[spec['alphabet']])
        if any(char not in permitted for char in observed):
            illegal_fmt = ', '.join(observed.difference(permitted))
            raise RuntimeError(f"Validation failed; illegal characters: {illegal_fmt}")
    if 'min_length' in spec:
        if len(string) < spec['min_length']:
            raise RuntimeError(f'Validation failed: min_length')
    if 'max_length' in spec:
        if len(string) > spec['max_length']:
            raise RuntimeError(f'Validation failed: max_length')
    return True


def run_directive(string, scheme, directive, spec):
    ''''''
    if 'prepare' in spec[scheme]['directives'][directive]:
        string = prepare(string, spec[scheme]['directives'][directive]['prepare'])
    if 'validate' in spec[scheme]['directives'][directive]:
        validate(string, spec[scheme]['directives'][directive]['validate'])
    if spec[scheme]['directives'][directive]['helper']:
        helper_name = f"{scheme}_{directive}".replace('-','_')
        helper = getattr(helpers, helper_name)
        string = helper(string)
    return string


def generate_hash(string, algorithm):
    ''''''
    string_hash = getattr(hashlib, algorithm)(string.encode())
    return string_hash


def generate_output(string_hash, spec, hide_prefix):
    ''''''
    encodings_raw = {n: getattr(encodings, m['type'])(string_hash) for n, m in spec.items()}
    encodings_fmt = {}
    for name, encoding_raw in encodings_raw.items():
        prefix = spec[name]['prefix'] if not hide_prefix else ''
        length = spec[name]['length'] if 'length' in spec[name] else len(encoding_raw)
        encodings_fmt[name] = f"{prefix}{encoding_raw[:length]}"
        if spec[name].get('include_full'):
            encodings_fmt[f'{name}-full'] = f"{prefix}{encoding_raw}"
    return encodings_fmt


def format_output(outputs, output_type):
    '''Returns string format'''
    if output_type == 'dict':
        return outputs
    elif output_type == 'tab':
        outputs_fmt = ''
        for k, v in outputs.items():
            outputs_fmt += f'{k}\t{v}\n'
        return outputs_fmt
    elif output_type == 'table':
        outputs_fmt = ''
        for k, v in outputs.items():
            outputs_fmt += f'{k:<20} {v:<15}\n'
        return outputs_fmt


@named('gen')
def generate(
    scheme: 'scheme name; use {scheme}.{directive} if scheme has multiple directives',
    file: "input file path or '-' for stdin",
    format: 'input format; mandatory if scheme has multiple formats' = '',
    output: 'output format' = 'dict',
    hide_prefix: 'hide encoding prefix; overrides scheme' = False):
    '''Generate identifier(s) for input file path according to specified scheme'''
    PACKAGE_PATH = os.path.dirname(__file__)
    scheme, _, directive = scheme.partition('.')

    # Load scheme specification
    yaml_path = pathlib.Path(f'{PACKAGE_PATH}/schemes/{scheme}.yaml')
    if not os.path.exists(yaml_path):
        raise FileNotFoundError(f'Failed to open specification for scheme {scheme}')
    spec = schema.load_scheme(yaml_path.read_text()).data

    # Validate presence and unambiguity of scheme, directive and format
    if scheme not in spec:
        raise RuntimeError(f'Unrecognised scheme {scheme}.')
    if not directive:
        if len(spec[scheme]['directives']) == 1:  # One option
            directive = next(iter(spec[scheme]['directives']))
        else:
            raise RuntimeError(f"Ambiguous directive for scheme {scheme}. "
                               f"Options: {', '.join(spec[scheme]['directives'].keys())}")
    if directive not in spec[scheme]['directives']:
        raise RuntimeError(f'Unrecognised directive {directive} for scheme {scheme}')
    if not format:
        if len(spec[scheme]['directives'][directive]['formats']) == 1:
            format = next(iter(spec[scheme]['directives'][directive]['formats']))
        else:
            raise RuntimeError(f'Ambiguous format for directive {directive} of scheme {scheme}')
    
    print(f"scheme: {scheme} ({directive}),"
          f" input: {'stdin' if file == '-' else file}", file=sys.stderr)

    # Validate, read and format input file or stdin
    if file == '-':
        string = sys.stdin.read().strip()
    elif os.path.exists(file):
        string = pathlib.Path(file).read_text().strip()
    else:
        raise FileNotFoundError(f'File {file} not found ')
    string = getattr(formats, format)(string)

    # Validate output format
    if not output in schema.OUTPUT_TYPES:
        raise RuntimeError(f'Unrecognised output type {output}. Options: {schema.OUTPUT_TYPES}')

    # Handle chained directives
    dag = [directive]
    target = spec[scheme]['directives'][directive].get('target')
    if target:
        g = {s: spec[scheme]['directives'][s].get('target') for s in spec[scheme]['directives']}
        d = target
        while d:
            dag.append(d)
            d = g[d]
    for d in dag:
        string = run_directive(string, scheme, d, spec)
    string_hash = generate_hash(string, spec[scheme]['algorithm'])
    outputs = generate_output(string_hash, spec[scheme]['encodings'], hide_prefix)
    
    return format_output(outputs, output)



def main():
    argh.dispatch_commands([generate])


if __name__ == '__main__':
    main()