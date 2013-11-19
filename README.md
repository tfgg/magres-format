magres-format
=============

Code for parsing the CCP-NC ab-initio magnetic resonance file format as used in the latest version of [CASTEP](http://www.castep.org), coming soon to other codes such as [Quantum ESPRESSO](http://www.quantum-espresso.org). See [more on this page](http://www.ccpnc.ac.uk/pmwiki.php/CCPNC/Fileformat)

Documentation for the Python API is [available here](http://tfgg.github.io/magres-format/build/html/).

A few IPython notebooks have been written using the library as examples:

 * [Plotting bonding networks](http://nbviewer.ipython.org/7203658)
 * [Plotting glycine chemical shifts](http://nbviewer.ipython.org/6699984)

Installing
----------

Clone the repository or download and extract the .zip file somewhere. From the command line, run:

    sudo python setup.py install
    
to install it globally, or

    python setup.py install --user
    
to install it just for your user account ([read more here](http://docs.python.org/2/install/#alternate-installation)), and the Python module and associated scripts should now be installed.

If you have installed it locally with --user, you may have to add ~/.local/bin to your PATH. You can do this by adding

    export PATH=$HOME/.local/bin:$PATH
    
to your ~/.bashrc and restarting your session or running "source ~/.bashrc". If you use tcsh you do

    setenv PATH $HOME/.local/bin:$PATH
    
and then restart your session.

Conversion script usage
-----------------------

The convertoldmagres.py script installed by the above command will convert an old-style Castep magres file to
the new-style format for use with the new tools. You use it from the command line like:

    convertoldmagres.py sample.magres > sample.new.magres

and optionally with the associated job's .castep file, to capture the lattice information,

    convertoldmagres.py sample.magres sample.castep > sample.new.magres

Python module usage
-------------------

From inside Python you can import the magres.format module and use the MagresFile class.


    >>> from magres.format import MagresFile
    >>> magres_file = MagresFile(open('samples/T1Si0.magres'))
    >>> magres_file
    <magres.format.MagresFile object at 0x26ad690>

The 'data_dict' member variable contains all the information parsed from the file. You can also request the data
serialized as JSON by the .as_json() method and load JSON data with the .load_json() method.

There is also the magres.atoms.MagresAtoms structure, which adds a more user-friendly interface to manipulating magres data on
top of MagresFile. For example, loading and printing out the isotropic magnetic shieldings from an ethanol calculation:

    >>> from magres.atoms import MagresAtoms
    >>> atoms = MagresAtoms.load_magres('ethanol.magres')
    >>> for atom in atoms:
    >>>   print atom, atom.ms.iso
    1H1 29.5599376391
    1H2 30.2261866485
    1H3 30.0722544561
    1H4 26.9539908188
    1H5 27.3739467591
    1H6 31.8881193712
    13C1 156.12291535
    13C2 109.357530445
    17O1 267.012276599

More documentation is [available here](http://tfgg.github.io/magres-format/build/html/).

The module also include the useful magres.constants module, which gives the best-known gamma constants and quadrupole 
moments for all isotopes, the most common isotopes used in experiments.

JSON schema
-----------

We use the [JSONschema](http://json-schema.org/) definition to provide a specification for the internal datastructure used by the parser and the format of the JSON emitted and consumed by .as_json() and .load_json() on MagresFile.

To dump the JSON representation of a .magres simply do

    magresjson.py sample.magres > sample.magres.json

and sample.magres.json should now contain a schema-compliant JSON representation of sample.magres.
