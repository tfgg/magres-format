#!/usr/bin/env python
from distutils.core import setup

setup(name='Magres format',
      version='0.9',
      description='Ab-initio magnetic resonance format',
      author='Timothy Green',
      author_email='timothy.green@gmail.com',
      url='http://www.ccpnc.ac.uk/pmwiki.php/CCPNC/Fileformat',
      packages=['magres', 'magres.schema'],
      scripts=['scripts/convertoldmagres.py', 'scripts/magresjson.py', 'scripts/extract-all-jc-compare.py', 'scripts/compisc.py', 'scripts/magresmerge.py', 'scripts/extract-all-jc.py', 'scripts/absiscdev.py'],
      )
