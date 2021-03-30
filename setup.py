import konstel

from setuptools import setup


setup(name = 'konstel',
      version = konstel.__version__,
      description = 'Hash-based phonemic sequence identifiers',
      url = 'http://github.com/bede/konstel',
      author = 'Bede Constantinides',
      author_email = 'bedeabc@gmail.com',
      license = 'LICENSE',
      packages=['konstel'],
      python_requires='>=3.6',
      install_requires=['argh', 'biopython', 'strictyaml'],
      entry_points = {'console_scripts':['konstel=konstel.konstel:main']},
      classifiers=['Intended Audience :: Science/Research',
                   'License :: OSI Approved :: MIT License',
                   'Natural Language :: English'])
