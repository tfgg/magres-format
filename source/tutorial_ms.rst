Processing magnetic shieldings
========================

Suppose we have the following .magres file for an ethanol molecule in a box:

Example datafile: ethanol.magres
--------------------------------

.. code::

  #$magres-abinitio-v1.0
  <atoms>
  units lattice Angstrom
  lattice    6.0000000000000009E+00   0.0000000000000000E+00   0.0000000000000000E+00   0.0000000000000000E+00   6.0000000000000009E+00   0.0000000000000000E+00   0.0000000000000000E+00   0.0000000000000000E+00   6.0000000000000009E+00
  units atom Angstrom
  atom H       H                  1          3.9805990000000007E+00          4.1783420000000007E+00          3.2950789999999999E+00
  atom H       H                  2          5.0333940000000004E+00          3.4304300000000003E+00          4.5047590000000008E+00
  atom H       H                  3          5.7190700000000012E+00          4.5522570000000000E+00          3.3153530000000009E+00
  atom H       H                  4          3.7202350000000006E+00          5.3295050000000010E+00          5.5099090000000004E+00
  atom H       H                  5          4.4121710000000007E+00          6.4335720000000007E+00          4.3170010000000012E+00
  atom H       H                  6          5.9116109999999997E+00          5.0322839999999998E+00          6.2422020000000007E+00
  atom C       C                  1          4.8469400000000009E+00          4.3506309999999999E+00          3.9411360000000002E+00
  atom C       C                  2          4.6030250000000006E+00          5.5187379999999999E+00          4.8825320000000012E+00
  atom O       O                  1          5.7462540000000013E+00          5.8127050000000011E+00          5.6871000000000009E+00
  </atoms>
  <magres>
  units ms ppm
  ms H                  1          3.0275414382832704E+01          1.2343176147963255E+00          3.8055347702168261E+00          1.9319927686201195E+00          2.7549852274727169E+01          2.4936355990438019E+00          4.1246695670194242E+00          2.2382434391773800E+00          3.0854546259833729E+01
  ms H                  2          2.6965983864643054E+01         -3.8585624895769854E-01          7.4285824854937599E-01         -4.1196902839127542E-01          3.5292668226256616E+01         -1.8802552393259337E+00         -6.7122267950404291E-01         -1.3639578313254843E+00          2.8419907854724499E+01
  ms H                  3          2.9971223084896938E+01          8.3477502779773927E-01         -3.2689438749161925E+00         -5.7219252176672120E-01          2.7791184879829256E+01         -2.1461997294113141E-01         -3.6440239128383367E+00          3.7767506813290157E-02          3.2454355403489878E+01
  ms H                  4          3.2261642842180500E+01          6.2520948574546298E-01         -1.4918313319633616E+00          7.5053965528921041E-01          2.2658737518918223E+01          1.7257630753650515E+00         -1.6495892494737605E-01          2.0870476890698075E+00          2.5941592095402918E+01
  ms H                  5          2.5948305588592422E+01         -2.7316306883849486E+00          3.8807548928060029E+00         -1.8011484678760228E+00          2.9698619245848839E+01         -4.6789512034172942E-01          2.9480868200872661E+00         -1.3677707067034952E+00          2.6474915442944035E+01
  ms H                  6          2.8536732593448992E+01         -2.0217790318703996E+00          4.8738962044242289E+00         -3.4947137896321606E-01          3.2849570358829844E+01         -6.5479198378162398E+00          4.6010039021206746E+00         -5.8146999362447280E+00          3.4278055161276548E+01
  ms C                  1          1.4998628033046998E+02         -1.0199071049759386E+01         -8.4736870826874103E-02          2.3362956476235755E-01          1.6059638213098884E+02          2.3169348644807837E+01          7.5564002686506901E+00          1.5878380496240151E+01          1.5778608358849078E+02
  ms C                  2          1.3326536661448372E+02          3.2787517370578650E-01          3.0641881628898851E+01          3.6408749419944675E+00          8.6537759377111740E+01          1.3720991423773526E+01          3.3471610166046347E+01          1.2743115987611173E+01          1.0826946534305428E+02
  ms O                  1          2.4574869953866229E+02          4.8677155199684901E+00          2.3843699230434890E+01         -2.5684198667080988E+01          2.8842353056750665E+02         -2.1091258283138934E+01          1.8346324157049754E+01         -3.9130135989063421E+00          2.6686459969143857E+02
  </magres>

Playing with ethanol.magres
---------------------------

You can run all these Python commands from the interpreter, supposing that you have saved the above data into a file called *ethanol.magres* in the same directory.

We can load this into an :py:class:`magres.atoms.MagresAtoms` object using the :py:meth:`~magres.atoms.MagresAtoms.load_magres` method.

.. code:: python

  from magres.atoms import MagresAtoms

  atoms = MagresAtoms.load_magres('ethanol.magres')

If we run this code, *atoms* now contains a list of all the atoms in the system, plus a bunch of helpful functions as documented in :py:class:`magres.atoms.MagresAtomsView`. Each atom will have a *ms* attribute, corresponding to a :py:class:`magres.atoms.MagresAtomMs` object. This object contains lots of different ways to process the sigma tensor, including a number of different conventions.

For example, supposing that we want to iterate over each atom and print its name and isotropic shielding (:py:attr:`~magres.atoms.MagresAtomMs.iso`) we can simply do

.. code:: python

  for atom in atoms:
    print atom, atom.ms.iso

this will output

.. code::

  1H1 29.5599376391
  1H2 30.2261866485
  1H3 30.0722544561
  1H4 26.9539908188
  1H5 27.3739467591
  1H6 31.8881193712
  13C1 156.12291535
  13C2 109.357530445
  17O1 267.012276599

Maybe we want to print out the spans and skews for all the protons. To do this we select all atoms with the 'H' species with :py:meth:`~magres.atoms.MagresAtomsView.species` and then access the :py:attr:`~magres.atoms.MagresAtomMs.span` and :py:attr:`~magres.atoms.MagresAtomMs.skew` attributes.

.. code:: python

  for atom in atoms.species('H'):
    print atom, atom.ms.span, atom.ms.skew

this will output

.. code::

  1H1 9.38562488256 0.818487161614
  1H2 8.72896735943 0.7449106712
  1H3 7.36075836882 0.925944171311
  1H4 10.6945028349 0.0495013823109
  1H5 9.32371700876 -0.0606983834955
  1H6 16.2585008574 0.471974409294

Or perhaps we want to print out all isotropic and anisotropic shieldings for all atoms within 2 Angstrom of the C1 atom.

.. code:: python

  for atom in atoms.within(atoms.get_species('C', 1), 2.0):
    print atom, atom.ms.iso, atom.ms.aniso

this will output

.. code::

  1H1 29.5599376391 8.95972202944
  1H2 30.2261866485 8.17230075322
  1H3 30.0722544561 7.22448160362
  13C1 156.12291535 34.0294598396
  13C2 109.357530445 70.540995898

You can also directly access the tensor, eigenvectors and eigenvaleus via the :py:attr:`~magres.atoms.MagresAtomMs.sigma`, :py:attr:`~magres.atoms.MagresAtomMs.evecs` and :py:attr:`~magres.atoms.MagresAtomMs.evals`. The eigenvectors and eigenvalues are by default ordered according to the Haeberlen convention. Let's print out these for the C1 atom:

.. code:: python

  atom = atoms.get_species('C', 1)

  print atom
  print atom.ms.sigma
  print atom.ms.evecs
  print atom.ms.evals

outputting

.. code::

  13C1
  [[  1.49986280e+02  -1.01990710e+01  -8.47368708e-02]
   [  2.33629565e-01   1.60596382e+02   2.31693486e+01]
   [  7.55640027e+00   1.58783805e+01   1.57786084e+02]]
  [array([-0.42141732, -0.9060046 , -0.03953616]), array([-0.62585978,  0.25900804,  0.73567273]), array([ 0.65628269, -0.33476933,  0.67618232])]
  [136.76839422482036, 152.79112991538537, 178.809221909744]

Note that we're using numpy for array support.
