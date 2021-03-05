from strictyaml import load, Bool, Int, Str, Seq, Map, MapPattern, Optional, YAMLValidationError


def load_scheme(yaml_text):
    schema = MapPattern(
        Str(), Map({
            'description': Str(),
            Optional('alias'): Str(),
            'version': Str(),
            'directives': MapPattern(
                Str(), Map({
                    'description': Str(),
                    'class': Str(),
                    'formats': Seq(Str()),
                    Optional('prepare'): Map({
                        Optional('strip_characters'): Seq(Str()),
                        Optional('remove_whitespace'): Bool(),
                        Optional('custom_function'): Bool(),
                    }),
                    Optional('validate'): Map({
                        'min_length': Int(),
                        'max_length': Int(),
                        Optional('custom_function'): Bool(),
                    }),
                    Optional('target'): Str()
                }),
            ),
            'algorithm': Str(),
            'encodings': MapPattern(
                Str(), Map({
                    'type': Str(),
                    'length': Int(),
                    'include_full': Bool(),
                    Optional('prefix'): Str(),
                    Optional('custom_function', default=False): Bool(),
                })
            )
        })
    )
    return load(yaml_text, schema)
    
