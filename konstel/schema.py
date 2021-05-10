import hashlib
import inspect

from strictyaml import load, Bool, Int, Str, Seq, Map, Enum, MapPattern, Optional, YAMLValidationError

from konstel import formats
from konstel import encodings


ALGORITHMS = hashlib.algorithms_available
ALPHABETS = formats.alphabets.keys()
BASE_ENCODINGS = {'base32', 'cbase32'}
ENCODINGS = {o[0] for o in inspect.getmembers(encodings, inspect.isfunction)}
FORMATS = {o[0] for o in inspect.getmembers(formats, inspect.isfunction)}
OUTPUT_TYPES = {'json', 'tsv', 'table'}


def load_scheme(yaml_text):
    '''
    Some optional keys have enforced default values, otherwise use dict.get()
    '''
    schema = MapPattern(
        Str(), Map({
            'description': Str(),
            Optional('alias'): Str(),
            'version': Str(),
            'directives': MapPattern(
                Str(), Map({
                    Optional('description'): Str(),
                    'formats': Seq(Enum(FORMATS)),
                    Optional('prepare'): Map({
                        Optional('remove_whitespace', default=False): Bool(),
                        Optional('remove_characters', default=['']): Seq(Str()),
                        Optional('strip_characters', default=['']): Seq(Str()),
                    }),
                    Optional('validate'): Map({
                        Optional('alphabet'): Enum(ALPHABETS),
                        Optional('min_length'): Int(),
                        Optional('max_length'): Int(),
                    }),
                    Optional('target'): Str(),
                    Optional('helper', default=False): Bool()
                }),
            ),
            'algorithm': Enum(ALGORITHMS),
            'encodings': MapPattern(
                Str(), Map({
                    'type': Enum(ENCODINGS),
                    Optional('length', default=0): Int(),
                    Optional('prefix', default=''): Str(),
                    Optional('separator'): Map({
                        'character': Str(),
                        'interval': Int()
                    })
                })
            )
        })
    )
    return load(yaml_text, schema)
