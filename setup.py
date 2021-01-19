import re

from setuptools import setup

__version__ = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        open('konstel/__init__.py').read()).group(1)

setup(name = 'konstel',
    version = __version__,
    description = 'Hash-based identifiers for SARS-CoV-2 spike variants',
    url = 'http://github.com/bede/konstel',
    author = 'Bede Constantinides',
    author_email = 'bedeabc@gmail.com',
    license = 'license',
    packages=['konstel'],
    zip_safe=True,
    python_requires='>=3.6',
    install_requires=['biopython', 'fire'],
    entry_points = {'console_scripts':['konstel=konstel.konstel:main']},
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English'])
