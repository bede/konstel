import os
import sys
import json
import typing
import hashlib
import pathlib

import defopt

from Bio import SeqIO
from Bio.Seq import Seq

import konstel.cli as cli
import konstel.schema as schema
import konstel.formats as formats
import konstel.helpers as helpers
import konstel.encodings as encodings


def load_scheme(scheme, validate_directive=True):
    ''''''
    PACKAGE_PATH = os.path.dirname(__file__)
    scheme, _, directive = scheme.partition('.')

    # Load scheme specification
    yaml_path = pathlib.Path(f'{PACKAGE_PATH}/schemes/{scheme}/scheme.yaml')
    if not os.path.exists(yaml_path):
        raise FileNotFoundError(f'Failed to open specification for scheme {scheme}')
    spec = schema.load_scheme(yaml_path.read_text()).data

    # Validate presence and unambiguity of scheme, directive and format
    if scheme not in spec:
        raise RuntimeError(f'Unrecognised scheme {scheme}.')
    
    if validate_directive:
        if not directive:
            if len(spec[scheme]['directives']) == 1:  # One option
                directive = next(iter(spec[scheme]['directives']))
            else:
                raise RuntimeError(f"Ambiguous directive for scheme {scheme}. "
                                   f"Options: {', '.join(spec[scheme]['directives'].keys())}")
        if directive not in spec[scheme]['directives']:
            raise RuntimeError(f'Unrecognised directive {directive} for scheme {scheme}')

    return scheme, directive, spec


def prepare(string, spec):
    ''''''
    if spec['remove_whitespace']:
        string = string.translate(str.maketrans('', '', ' \n\t\r'))
    if any(spec['strip_characters']):  # any() since schema default [''] is truthy
        string = string.strip(''.join(spec['strip_characters']))
    return string


def validate_string(string, spec):
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
        validate_string(string, spec[scheme]['directives'][directive]['validate'])
    if spec[scheme]['directives'][directive]['helper']:
        helper_name = f"{scheme}_{directive}".replace('-','_')
        helper = getattr(helpers, helper_name)
        string = helper(string)
    return string


def generate_hash(string, algorithm):
    ''''''
    hash_string = getattr(hashlib, algorithm)(string.encode()).hexdigest()
    return hash_string


def generate_encodings(hash_b16, spec, length, hide_prefix):
    ''''''
    encodings_raw = {n: getattr(encodings, m['type'])(hash_b16) for n, m in spec.items()}
    encodings_fmt = {}
    for name, encoding_raw in encodings_raw.items():
        prefix = spec[name]['prefix'] if not hide_prefix else ''
        
        # Determine length of each encoding
        if length:  # Manual override
            scheme_length = length
        elif 'length' in spec[name] and spec[name]['length']:  # In scheme and is truthy
            scheme_length = spec[name]['length']
        else:  # Fallback
            scheme_length = len(encoding_raw)

        # Handle 'hash' differently
        if name == 'hash':
            encodings_fmt[name] = f'{prefix}{encoding_raw}'
            if spec[name]['length']:
                encodings_fmt[f'{name}-{scheme_length}'] = f'{prefix}{encoding_raw[:scheme_length]}'
        else:
            encodings_fmt[name] = f'{prefix}{encoding_raw[:scheme_length]}'

        # Add separators
        if 'separator' in spec[name]:
            char = spec[name]['separator']['character']
            interval = spec[name]['separator']['interval']
            separated_chars = []
            for i in range(0, len(encodings_fmt[name])-2, interval):
                separated_chars.append(encodings_fmt[name][i+len(prefix):i+len(prefix)+interval])
            encodings_fmt[name] = prefix + char.join(separated_chars)

    return encodings_fmt


def format_encodings(outputs, output_type='dict'):
    '''Optionally prints encodings in tabular format'''
    if output_type == 'json':
        return json.dumps(outputs)  # Do nothing, return dict
    elif output_type == 'tsv':
        return '\t'.join(outputs.values())
    elif output_type == 'table':
        outputs_fmt = ''
        for k, v in outputs.items():
            outputs_fmt += f'{k:<20} {v:<15}\n'
        return outputs_fmt


def generate(scheme: str,
             file: str,
             format: typing.Union[str, None] = None,
             output: str = 'json',
             length: int = 0,
             hide_prefix: bool = False):
    '''
    Generate identifier(s) for input file path or stdin according to the specified scheme
    Returns dict, and prints format specified in OUTPUT

    :arg scheme: Scheme name; use {scheme}.{directive} if scheme specifies multiple directives
    :arg file: Input path or - for stdin
    :arg format: Input format; mandatory if scheme specifies multiple formats
    :arg output: Output format for CLI (json, tsv, table)
    :arg length: Encoding length; overrides scheme
    :arg hide_prefix: Hide encoding prefix; overrides scheme
    '''
    scheme, directive, spec = load_scheme(scheme)

    # Validate format
    if not format:
        if len(spec[scheme]['directives'][directive]['formats']) == 1:
            format = next(iter(spec[scheme]['directives'][directive]['formats']))
        else:
            raise RuntimeError(f'Ambiguous format for directive {directive} of scheme {scheme}')

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
    hash_b16 = generate_hash(string, spec[scheme]['algorithm'])
    encodings_ = generate_encodings(hash_b16, spec[scheme]['encodings'], length, hide_prefix)
    outputs = {**{'scheme': scheme}, **encodings_}

    return outputs


def regenerate(scheme: str,
               hash_string: str,
               output: str = 'json',
               length: int = None,
               hide_prefix: bool = False):
    '''
    Regenerate identifier(s) for an existing full length hash according to the specified scheme

    :arg scheme: Scheme name
    :arg hash_string: Full length hash digest
    :arg output: Output format (dict, tsv or table)
    :arg length: Encoding length (overrides scheme)
    :arg hide_prefix: Hide encoding prefix (overrides scheme)
    '''
    scheme, directive, spec = load_scheme(scheme, validate_directive=False)
    hash_type = spec[scheme]['encodings']['hash']['type']
    print(f"Scheme uses {hash_type}", file=sys.stderr)

    # Validate presence of 'hash' in scheme
    if not 'hash' in spec[scheme]['encodings']:
        raise RuntimeError(f'Scheme {scheme} must specify a hash encoding')

    # Validate output format
    if not output in schema.OUTPUT_TYPES:
        raise RuntimeError(f'Unrecognised output type {output}. Options: {schema.OUTPUT_TYPES}')

    # Normalise and decode supplied hash encoding
    prefix = spec[scheme]['encodings']['hash']['prefix']
    hash_string = hash_string[hash_string.startswith(prefix) and len(prefix):]  # Remove prefix
    if 'separator' in spec[scheme]['encodings']['hash']:
        sep_char = spec[scheme]['encodings']['hash']['separator']['character']
        hash_string = hash_string.replace(sep_char, '')
    hash_b16 = getattr(encodings, f'decode_{hash_type}')(hash_string)

    # Regenerate using scheme
    encodings_ = generate_encodings(hash_b16, spec[scheme]['encodings'], length, hide_prefix)
    outputs = {**{'scheme': scheme}, **encodings_}

    return outputs
