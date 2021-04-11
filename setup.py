import konstel

from setuptools import setup


setup(
    name = 'konstel',
    version = konstel.__version__,
    description = 'Hash-based phonemic sequence identifiers',
    url = 'http://github.com/bede/konstel',
    author = 'Bede Constantinides',
    author_email = 'bedeabc@gmail.com',
    license = 'LICENSE',
    packages=['konstel'],
    include_package_data=True,  # MANIFEST.in
    python_requires='>=3.6',
    install_requires=[
        'defopt>=6.1.0',
        'strictyaml>=1.4.0',
        'biopython>=1.78',
        'parasail==1.2.4'],
    entry_points = {
        'console_scripts':['konstel=konstel.cli:main']},
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English'])
