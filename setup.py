#!/usr/bin/env python
from distutils.core import setup
import os, site, sys

def bin():
  return os.path.join(site.USER_BASE, "bin")

def bin_on_path():
  return bin() in os.environ['PATH'].split(':')

def find_scripts():
  for f in os.listdir("scripts/"):
    if f.endswith(".py"):
      yield os.path.join("scripts", f)

setup(name='Magres format',
      version='0.9',
      description='Ab-initio magnetic resonance format',
      author='Timothy Green',
      author_email='timothy.green@gmail.com',
      url='http://www.ccpnc.ac.uk/pmwiki.php/CCPNC/Fileformat',
      packages=['magres', 'magres.schema'],
      scripts=list(find_scripts()),
      )

if not bin_on_path() and "--user" in sys.argv[1:]:
  print "\n\nWARNING: Your scripts directory (%s) is not on your PATH" % bin()

  shell = os.environ["SHELL"].split('/')[-1]

  if shell == "bash":
    print "You can fix this by adding the following to your .bashrc and reloading your session"
    print "export PATH=%s:$PATH" % bin()
  elif shell == "tcsh":
    print "You can fix this by adding the following to your .tcshrc and reloading your session"
    print "setenv PATH %s:$PATH" % bin()
  else:
    print "You can fix this by adding %s to your PATH environment variable" % bin()

