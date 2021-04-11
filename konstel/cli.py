import typing

import defopt

from konstel import konstel


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
    :arg output: Output format (json, tsv, table)
    :arg length: Encoding length; overrides scheme
    :arg hide_prefix: Hide encoding prefix; overrides scheme
    '''
    outputs = konstel.generate(scheme, file, format, output, length, hide_prefix)
    print(konstel.format_encodings(outputs, output))


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
    outputs = konstel.regenerate(scheme, hash_string, output, length, hide_prefix)
    print(konstel.format_encodings(outputs, output))


def main():
    defopt.run(
        {'gen': generate, 'regen': regenerate},
        strict_kwonly=False, no_negated_flags=True)