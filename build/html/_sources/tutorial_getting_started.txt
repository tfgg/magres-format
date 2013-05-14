Getting started
===============

Installing
----------

You can either use git to clone the latest version of the repository or download the repository as a .zip.

>>> git clone git://github.com/tfgg/magres-format.git
>>> cd magres-format

or

>>> wget https://github.com/tfgg/magres-format/archive/master.zip
>>> unzip master.zip 
>>> cd magres-format-master

Now you run setup.py to install the module and associated scripts. If you want to install it system wide, run

>>> sudo python setup.py install

If you want to install it locally (e.g. if you don't have sudo), run

>>> python setup.py install --user

and it should now be available to you.

Python module
-------------

The package comes with a number of Python modules, documented separately.

.. toctree::
   :maxdepth: 2

   atoms
   format

Import the modules into your code with the standard import commands

>>> import magres.atoms
>>> import magres.format

or just import individual components

>>> from magres.atoms import MagresAtoms


Helper scripts
--------------

The package comes with a number of helper scripts that can be run from the command line.

Dump JSON: magresjson.py
~~~~~~~~~~~~~~~~~~~~~~~~

This will give you a JSON dump of a given .magres text file.

>>> magresjson.py sample.magres > sample.magres.json

Merge .magres files: magresmerge.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This will merge together several given .magres (e.g. different perturbing atoms in an indirect spin-spin coupling calculation) and output the result.

>>> magresmerge.py */*.magres > ethanol-jc.magres

Quickly analyse J-couplings: extract-all-jc-compare.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> extract-all-jc-compare.py [directory] [species] [index]

This will find all two-way indirect spin-spin coupling calculations below [directory] between atom [species] [index] and other atoms and print them out for analysis.

Samples and example scripts
---------------------------

There are some sample .magres files in the samples/ directory and some example Python processing scripts in the examples/ directory. Run the example scripts from the package root, like

>>> python examples/atoms_ms.py

or they will not find the correct .magres files.
