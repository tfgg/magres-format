#!/bin/bash
# The -e flag causes the script to exit as soon as one command returns a non-zero exit code.
# The -v flag makes the shell print all lines in the script before executing them.
set -ev

cd notebooks

FILES=*.ipynb
for f in $FILES
do
  python ipynbdoctest.py "$f"
done
