# Physical units lifted from castep's io.f90
from constants import *

# Lengths

bohr       = 1.0
metre      = electron_mass_si*speed_light_si*fine_structure_si/hbar_si
centimetre = metre * 1.0e-2
nanometre  = metre * 1.0e-9
angstrom   = metre * 1.0e-10

# Masses
 
electron_mass = 1.0
amu           = 1.0e-3/avogadro_si/electron_mass_si
kilogram      = 1.0/electron_mass_si
gram          = kilogram*1e-3

# Times

aut         = 1.0
second      = speed_light_si**2*fine_structure_si**2*electron_mass_si/hbar_si
millisecond = 1e-3*second
microsecond = 1e-6*second
nanosecond  = 1e-9*second
picosecond  = 1e-12*second
femtosecond = 1e-15*second

# Charges

elementary_charge = 1.0
coulomb           = 1.0/elementary_charge_si

# Electric Dipole moments

debye_si = 1.0e-21/speed_light_si # Defn is 1D = 10**-18 statcoulomb cm (cgs)
debye    = debye_si*coulomb*metre

# Energies

hartree            = 1.0
millihartree       = 1e-3
electron_volt      = elementary_charge_si/(fine_structure_si**2*electron_mass_si*speed_light_si**2)
millielectron_volt = electron_volt*1e-3
rydberg            = 0.5
millirydberg       = rydberg*1e-3
joule              = 1.0/(fine_structure_si**2*electron_mass_si*speed_light_si**2)
erg                = joule*1e-7
kilojoulepermole   = joule/avogadro_si*1e3
kilocalpermole     = kilojoulepermole*4.184 # Specific heat of water
hertz              = planck_si*joule
megahertz          = hertz*1e6
gigahertz          = hertz*1e9
terahertz          = hertz*1e12
wavenumber         = hertz*speed_light_si*1e2
kelvin             = boltzmann_si*joule

# Entropy
  
joulebymolebykelvin   = 1.0/molar_gas_si
caloriebymolebykelvin = 4.184*joulebymolebykelvin
  
# Forces

hartreebybohr = 1.0
eVbyang       = electron_volt/angstrom
newton        = joule/metre
dyne          = newton*1e-5

# Velocities

auv            = 1.0
angperps       = angstrom/picosecond
angperfs       = angstrom/femtosecond
bohrperps      = bohr/picosecond
bohrperfs      = bohr/femtosecond
metrepersecond = metre/second

# Pressures

hartreebybohr3 = 1.0
evbyang3       = electron_volt/angstrom**3
pascal         = newton/metre**2
megapascal     = pascal*1e6
gigapascal     = pascal*1e9
terapascal     = pascal*1e12
petapascal     = pascal*1e15
atmosphere     = pascal*101325.027 # Conversion to atmospheres
bar            = pascal*1e5
megabar        = bar*1e6

# Reciprocal length

invbohr      = 1.0
invmetre     = 1.0/metre
invnanometre = 1.0/nanometre
invangstrom  = 1.0/angstrom

# Force constant

hartreebybohr2   = 1.0
evbyang2         = electron_volt/angstrom**2
newtonbymetre    = newton/metre
dynebycentimetre = dyne/centimetre

# Volumes

bohr3       = 1.0
metre3      = metre**3
centimetre3 = (metre * 1e-2)**3
nanometre3  = (metre * 1e-9)**3
angstrom3   = (metre * 1e-10)**3

# Magres

acu         = 1.0
ampere      = (fine_structure_si**2*electron_mass_si*speed_light_si**2*planck_si)/elementary_charge_si
acd         = 1.0
amperemetre2 = ampere/(metre*metre)
amfd         = 1.0
tesla        = 1.0/planck_si
gauss        = (1.0/planck_si)*1.0e-4
agr          = 1.0
radsectesla  = electron_mass_si/elementary_charge_si
mhztesla     = 2.0*math.pi*radsectesla*1.0e-6
bohr2        = 1
fm2          = metre*metre*1.0e-30
barn         = metre*metre*1.0e-28
millibarn    = barn*1.0e-3

# IR Intensities

e2byamu     = 1.0/amu
d2byamuang2 = (debye/angstrom)**2/amu
kmbymol     = d2byamuang2/42.2561  # Strange unit favoured by spectroscopists.
                                      # Conversion factor is from Gaussian03 docs. 
# Efield
hartreebybohrbye = 1.0
eVbyangbye       = electron_volt/angstrom
newtonbycoulomb  = joule/metre/elementary_charge_si

